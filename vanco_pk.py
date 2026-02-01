# vanco_pk.py
import numpy as np
from dataclasses import dataclass
from datetime import datetime, timedelta

INFUSION_RATE = 1000.0  # mg/hour
LN2 = np.log(2)

@dataclass
class Dose:
    time: float
    amount: float
    ke: float = None

def pk_params_from_patient(
    age,
    sex,
    weight_kg,
    height_cm,
    cr_func,
    when,
):
    """
    Returns ke (1/h) and Vd (L) from patient characteristics.
    """
    # -------------------------
    # Creatinine at time point
    # -------------------------
    cr = cr_func(when)  # µmol/L

    # -------------------------
    # Weightless Cockcroft–Gault
    # -------------------------
    crcl = (140 - age) * 88.4 / cr
    if sex.lower().startswith("f"):
        crcl *= 0.85

    # Serum-to-plasma Cr correction
    crcl *= 0.9  # REQUIRED

    # -------------------------
    # Matzke elimination
    # -------------------------
    ke = 0.00083 * crcl + 0.0044

    KE_MAX = 0.693 / 5.5   # fastest allowed elimination
    KE_MIN = 0.693 / 120  # very slow elimination safety floor

    ke = min(ke, KE_MAX)
    ke = max(ke, KE_MIN)

    # -------------------------
    # Ideal & dosing body weight
    # -------------------------
    if sex.lower().startswith("m"):
        ibw = 50 + 0.9 * (height_cm - 152)
    else:
        ibw = 45.5 + 0.9 * (height_cm - 152)

    if weight_kg > 1.25 * ibw:
        # Obese → dosing body weight
        dbw = ibw + 0.4 * (weight_kg - ibw)
        weight_for_vd = dbw
    else:
        # Non-obese → actual body weight
        weight_for_vd = weight_kg

    # -------------------------
    # Volume of distribution
    # -------------------------
    vd = 0.75 * weight_for_vd  # L

    return ke, vd


class VancoPK:
    def __init__(self, ke, vd):
        self.ke = ke
        self.ke_prior = ke
        self.vd = vd
        self.ke_sd = None

    @property
    def clearance(self):
        """Clearance in L/h"""
        return self.ke * self.vd

    @property
    def half_life(self):
        """Half-life in hours"""
        return 0.693 / self.ke if self.ke > 0 else float("inf")

    def run(self, doses, duration_days=7, dt=0.1):
        """
        Simulates concentration over time.
        doses: list of (time_hours, dose_mg)
        returns dict with time, conc, auc24, ke, half_life, vd
        """
        t_end = duration_days * 24
        time = np.arange(0, t_end + dt, dt)
        conc = np.zeros_like(time)

        # Track active infusions
        infusions = []
        for t_dose, dose_mg in doses:
            duration = dose_mg / INFUSION_RATE  # hours
            infusions.append((t_dose, t_dose + duration, dose_mg))

        for i in range(1, len(time)):
            t = time[i]
            c_prev = conc[i - 1]

            # Elimination
            dc = -self.ke * c_prev * dt

            # Infusion input
            rate_in = 0.0
            for t_start, t_end_inf, dose_val in infusions:
                if t_start <= t < t_end_inf:
                    rate_in += INFUSION_RATE / self.vd  # mg/L/h

            dc += rate_in * dt
            conc[i] = c_prev + dc

        # AUC24 from last 24h
        mask24 = time >= (t_end - 24)
        if any(mask24):
            auc24 = np.trapezoid(conc[mask24], time[mask24])
        else:
            auc24 = 0.0

        return {
            "time": time,
            "conc": conc,
            "auc24": auc24,
            "ke": self.ke,
            "half_life": self.half_life,
            "vd": self.vd,
        }

    def simulate_regimen(self, dose_mg, interval_h, start_datetime, end_datetime, dt=0.1):
        """Simulate a standing regimen between two datetimes."""
        total_hours = (end_datetime - start_datetime).total_seconds() / 3600.0
        doses = []
        t = 0.0
        while t <= total_hours:
            doses.append((t, dose_mg))
            t += interval_h

        duration_days = total_hours / 24.0
        return self.run(doses, duration_days, dt)

    def predict_at_times(self, doses, times_h):
        """
        Predict concentrations at specific times using the SAME engine as run()
        """
        t_end = max(times_h) + 0.1
        time = np.arange(0, t_end + 0.1, 0.1)
        conc = np.zeros_like(time)

        # build infusions
        infusions = []
        for t_dose, dose_mg in doses:
            dur = dose_mg / INFUSION_RATE
            infusions.append((t_dose, t_dose + dur))

        for i in range(1, len(time)):
            c_prev = conc[i - 1]
            dc = -self.ke * c_prev * 0.1

            rate_in = 0.0
            for t0, t1 in infusions:
                if t0 <= time[i] < t1:
                    rate_in += INFUSION_RATE / self.vd

            conc[i] = c_prev + (dc + rate_in * 0.1)

        return np.interp(times_h, time, conc)


    def fit_ke_from_levels(self, doses, level_times, levels, sim_start, sigma_frac=0.15):
        """
        Bayesian MAP estimation of ke using measured levels
        """

        # Convert level times to hours since sim_start
        times_h = [
            (t - sim_start).total_seconds() / 3600
            for t in level_times
        ]

        ke_grid = np.linspace(self.ke_prior * 0.3, self.ke_prior * 3.0, 200)
        log_post = []

        obs = np.array(levels)
        sigma = sigma_frac * np.maximum(obs, 1.0)

        for k in ke_grid:
            self.ke = k
            preds = self.predict_at_times(doses, times_h)

            ll = -np.sum((obs - preds) ** 2 / (2 * sigma ** 2))
            prior = -((k - self.ke_prior) ** 2) / (2 * (0.3 * self.ke_prior) ** 2)

            log_post.append(ll + prior)

        log_post = np.array(log_post)
        idx = np.argmax(log_post)

        self.ke = ke_grid[idx]

        # curvature-based SD
        d2 = np.gradient(np.gradient(log_post, ke_grid), ke_grid)
        self.ke_sd = np.sqrt(-1 / d2[idx]) if d2[idx] < 0 else 0.25 * self.ke

        return self.ke


    def compute_ci(self, level=0.5):
        """Compute a confidence interval for ke."""
        z = 0.674 if level == 0.5 else 1.96
        sd = self.ke_sd if self.ke_sd is not None else 0.2 * self.ke
        lo = max(self.ke - z * sd, 0)
        hi = self.ke + z * sd
        return lo, hi


def dose_ci_from_ke(pk, target_auc, interval_h, n=1000, level=0.5):
    """
    Calculates dose confidence interval based on ke uncertainty.
    This is a standalone function called by streamlit_app.py.
    """
    sd = pk.ke_sd if pk.ke_sd else 0.2 * pk.ke
    ke_samples = np.random.normal(pk.ke, sd, n)
    ke_samples = ke_samples[ke_samples > 0]

    doses = []
    for k in ke_samples:
        cl = k * pk.vd
        auc_per_mg = interval_h / cl
        dose = target_auc / auc_per_mg
        doses.append(dose)

    return np.percentile(doses, 25), np.percentile(doses, 75)