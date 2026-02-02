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
        self.ke = ke  # Population base ke
        self.vd = vd
        self.ke_multiplier = 1.0  # NEW: The scaling factor from fitting
        self.ke_sd = 0.2 * ke     # Standard deviation

    @property
    def clearance(self):
        """Clearance in L/h"""
        return self.ke * self.vd

    @property
    def half_life(self):
        """Half-life in hours"""
        return 0.693 / self.ke if self.ke > 0 else float("inf")

    def run(self, doses, duration_days=7, sim_start=None, cr_func=None, patient_info=None):
        t_max = duration_days * 24
        # 0.1 hour steps for smooth curves and numerical stability
        t_grid = np.linspace(0, t_max, int(t_max * 10) + 1) 
        dt = t_grid[1] - t_grid[0]
        
        conc = np.zeros_like(t_grid)
        current_conc = 0.0
        
        # We track ke for the summary metrics at the end
        last_ke = self.ke 

        for i, t_h in enumerate(t_grid):
            # 1. DYNAMIC KE UPDATE
            if cr_func and patient_info and sim_start:
                current_dt = sim_start + timedelta(hours=t_h)
                # Recalculate based on patient characteristics at THIS specific time
                new_ke, _ = pk_params_from_patient(
                    patient_info['age'], patient_info['sex'], 
                    patient_info['weight'], patient_info['height'], 
                    cr_func, current_dt
                )
                # Apply the Bayesian multiplier (if any) to the new pop-ke
                active_ke = new_ke * getattr(self, 'ke_multiplier', 1.0)
                last_ke = active_ke
            else:
                active_ke = self.ke

            # 2. INFUSION LOGIC (Preserved)
            input_rate = 0.0
            for t_dose, amt in doses:
                infusion_duration = amt / INFUSION_RATE
                if t_dose <= t_h <= (t_dose + infusion_duration):
                    input_rate = INFUSION_RATE / self.vd
                    break

            # 3. DIFFERENTIAL EQUATION STEP (Euler Method)
            # dC/dt = (Infusion_In / Vd) - (ke * C)
            dC = (input_rate - active_ke * current_conc) * dt
            current_conc += dC
            conc[i] = max(current_conc, 0) # Prevent negative numbers

        return {
            "time": t_grid,
            "conc": conc,
            "ke": last_ke,
            "vd": self.vd,
            "half_life": np.log(2) / last_ke if last_ke > 0 else 0,
            "auc24": (conc[-240:].sum() * dt) if len(conc) >= 240 else 0
        }

    def simulate_regimen(self, dose_mg, interval_h, sim_start, sim_end, cr_func=None, patient_info=None):
        """
        Correctly simulates a regimen while preserving the Bayesian fit 
        and time-varying creatinine trends.
        """
        from dosing import build_ordered_doses
        
        # 1. Use the existing helper to correctly space doses relative to sim_start
        # We start the 'try' regimen at sim_start to show the full comparison
        doses = build_ordered_doses(
            dose_mg=dose_mg, 
            interval_h=interval_h, 
            start_dt=sim_start, 
            sim_start=sim_start, 
            sim_end=sim_end
        )

        # 2. Calculate duration
        duration_days = (sim_end - sim_start).total_seconds() / 86400.0

        # 3. Call run() with ALL the context required for dynamic ke
        return self.run(
            doses=doses, 
            duration_days=duration_days, 
            sim_start=sim_start, 
            cr_func=cr_func, 
            patient_info=patient_info
        )

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


    def fit_ke_from_levels(self, doses, times_dt, obs, sim_start, cr_func, patient_info):
        # Convert datetimes to hours since sim_start
        times_h = [(t - sim_start).total_seconds() / 3600 for t in times_dt]
        obs = np.array(obs)

        # We search for a multiplier between 0.3x and 3.0x of population ke
        mult_grid = np.linspace(0.3, 3.0, 100)
        log_post = []
        sigma = 2.0  # Level measurement error SD

        for m in mult_grid:
            self.ke_multiplier = m
            # Run a mini-sim to get predicted levels at the exact times
            # Note: You must ensure self.run() uses self.ke_multiplier inside its loop!
            res = self.run(doses, duration_days=7, sim_start=sim_start, 
                           cr_func=cr_func, patient_info=patient_info)
            
            # Interpolate predicted concentrations at the specific observation times
            preds = np.interp(times_h, res["time"], res["conc"])

            ll = -np.sum((obs - preds) ** 2 / (2 * sigma ** 2))
            # Prior: Expect the multiplier to be near 1.0 (population estimate)
            prior = -((m - 1.0) ** 2) / (2 * 0.3 ** 2)
            log_post.append(ll + prior)

        # Pick the best multiplier
        idx = np.argmax(log_post)
        self.ke_multiplier = mult_grid[idx]
        
        # Approximate the SD of the multiplier for the CI
        d2 = np.gradient(np.gradient(log_post, mult_grid), mult_grid)
        self.multiplier_sd = np.sqrt(-1 / d2[idx]) if d2[idx] < 0 else 0.2
        return self.ke_multiplier

    def compute_ci(self, level=0.5):
        z = 0.674 if level == 0.5 else 1.96
        sd = getattr(self, 'multiplier_sd', 0.2)
        lo_mult = max(self.ke_multiplier - z * sd, 0.1)
        hi_mult = self.ke_multiplier + z * sd
        return lo_mult, hi_mult


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