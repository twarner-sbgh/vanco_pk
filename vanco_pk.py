# vanco_pk.py
import numpy as np
from dataclasses import dataclass
from datetime import timedelta

LN2 = np.log(2)

KE_MAX = LN2 / 5.5        # 5.5 h half-life ceiling
KE_MIN = LN2 / 120        # safety floor (~120 h)

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
    Returns ke (1/h) and Vd (L) from patient characteristics
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
        # Dosing body weight with obesity factor 0.4
        dbw = ibw + 0.4 * (weight_kg - ibw)
    else:
        dbw = weight_kg

    # -------------------------
    # Volume of distribution
    # -------------------------
    vd = 0.75 * dbw  # L

    return ke, vd

class VancoPK:
    def __init__(self, ke, vd):
        self.ke = ke      # 1/h
        self.vd = vd      # L

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
        doses: list of (time_hours, dose_mg)
        returns dict with time, conc, auc24, ke, half_life, vd
        """

        t_end = duration_days * 24
        time = np.arange(0, t_end + dt, dt)
        conc = np.zeros_like(time)

        for t_dose, dose in doses:
            mask = time >= t_dose
            conc[mask] += (dose / self.vd) * np.exp(
                -self.ke * (time[mask] - t_dose)
            )

        mask24 = time >= (t_end - 24)
        auc24 = np.trapezoid(conc[mask24], time[mask24])

        return {
            "time": time,
            "conc": conc,
            "auc24": auc24,
            "ke": self.ke,
            "half_life": self.half_life,
            "vd": self.vd,
        }

    def simulate_regimen(
        self,
        dose_mg,
        interval_h,
        start_datetime,
        end_datetime,
        dt=0.1,
    ):
        """
        Simulate a standing regimen between two datetimes.
        """

        total_hours = (end_datetime - start_datetime).total_seconds() / 3600.0

        doses = []
        t = 0.0
        while t <= total_hours:
            doses.append((t, dose_mg))
            t += interval_h

        duration_days = total_hours / 24.0
        return self.run(doses, duration_days, dt)

def compute_ci(self, level=0.5):
    """
    Compute a confidence interval for ke based on fitted variance.
    Currently a simple symmetric CI around ke.
    """

    # Safety check
    if not hasattr(self, "ke"):
        raise ValueError("ke not calculated; cannot compute CI")

    # If no variance estimate exists, return a conservative interval
    if not hasattr(self, "ke_sd") or self.ke_sd is None:
        delta = 0.2 * self.ke
        return self.ke - delta, self.ke + delta

    z = 0.674 if level == 0.5 else 1.96
    lo = max(self.ke - z * self.ke_sd, 0)
    hi = self.ke + z * self.ke_sd
    return lo, hi
