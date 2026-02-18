import numpy as np
from dataclasses import dataclass
from datetime import datetime, timedelta # datetime used for type hinting or future extensions

INFUSION_RATE = 1000.0  # mg/hour - fixed infusion rate for all doses (1g/hr)
LN2 = np.log(2)

@dataclass
class Dose:
    time: float
    amount: float
    ke: float = None

def pk_params_from_patient(age, sex, weight_kg, height_cm, cr_func, when):
    """
    Returns ke (1/h) and Vd (L) from patient characteristics.
    """
    cr = cr_func(when)
    crcl = (140 - age) * 88.4 / cr
    if sex.lower().startswith("f"):
        crcl *= 0.85

    crcl *= 0.9  # Serum-to-plasma Cr correction
    ke = 0.00083 * crcl + 0.0044

    # Safety floors/ceilings
    ke = max(min(ke, 0.693 / 5.5), 0.693 / 120)

    if sex.lower().startswith("m"):
        ibw = 50 + 0.9 * (height_cm - 152)
    else:
        ibw = 45.5 + 0.9 * (height_cm - 152)

    weight_for_vd = ibw + 0.4 * (weight_kg - ibw) if weight_kg > 1.25 * ibw else weight_kg
    vd = 0.75 * weight_for_vd

    return ke, vd

class VancoPK:
    def __init__(self, ke, vd):
        self.ke = ke # Population base ke (will be updated to reflect current state)
        self.vd = vd
        self.ke_multiplier = 1.0
        self.multiplier_sd = 0.2

    def run(self, doses, duration_days=7, sim_start=None, cr_func=None, patient_info=None):
        t_max = duration_days * 24
        # CHANGE: Increase multiplier from 10 to 12 for 5-minute steps
        t_grid = np.linspace(0, t_max, int(t_max * 12) + 1)
        dt = t_grid[1] - t_grid[0] # This will now be 0.0833 (5 mins) instead of 0.1
        
        pop_ke_traj = np.full_like(t_grid, self.ke)
        if cr_func and patient_info and sim_start:
            for i, t_h in enumerate(t_grid):
                pop_ke_traj[i], _ = pk_params_from_patient(
                    patient_info['age'], patient_info['sex'], 
                    patient_info['weight'], patient_info['height'], 
                    cr_func, sim_start + timedelta(hours=t_h)
                )
        
        # UPDATED: Update self.ke to the LATEST population value so metrics/suggestions are current
        self.ke = pop_ke_traj[-1] 

        conc = np.zeros_like(t_grid)
        curr_c = 0.0
        for i, t_h in enumerate(t_grid):
            active_ke = pop_ke_traj[i] * self.ke_multiplier
            rate_in = 0.0
            for t_d, amt in doses:
                # IMPORTANT: Use a small epsilon (1e-9) to avoid floating point misses
                if t_d <= t_h <= (t_d + amt/INFUSION_RATE + 1e-9):
                    rate_in = INFUSION_RATE / self.vd
                    break
            curr_c += (rate_in - active_ke * curr_c) * dt
            conc[i] = max(curr_c, 0)

        return {
            "time": t_grid,
            "conc": conc,
            "ke": self.ke * self.ke_multiplier,
            "vd": self.vd,
            "half_life": LN2 / (self.ke * self.ke_multiplier),
            "auc24": (conc[-240:].sum() * dt) if len(conc) >= 240 else 0
        }

    def simulate_regimen(self, dose_mg, interval_h, start_dt, end_dt, cr_func, patient_info):
            """
            Calculates a sequence of doses based on an interval and runs the simulation.
            """
            # 1. Determine total simulation hours
            duration_h = (end_dt - start_dt).total_seconds() / 3600
            duration_days = duration_h / 24
            
            # 2. Build the list of doses [(time_in_hours, amount), ...]
            try_doses = []
            curr_t = 0.0
            while curr_t < duration_h:
                try_doses.append((curr_t, dose_mg))
                curr_t += interval_h
                
            # 3. Call the optimized run method with the necessary dynamic info
            return self.run(
                doses=try_doses, 
                duration_days=duration_days, 
                sim_start=start_dt, 
                cr_func=cr_func, 
                patient_info=patient_info
            )

    def _run_fast(self, doses, t_grid, ke_traj, multiplier):
        dt = t_grid[1] - t_grid[0]
        conc = np.zeros_like(t_grid)
        curr_c = 0.0
        for i, t_h in enumerate(t_grid):
            active_ke = ke_traj[i] * multiplier
            rate_in = 0.0
            for t_d, amt in doses:
                # Standardizing the check for 5-minute precision
                if t_d <= t_h <= (t_d + amt/INFUSION_RATE + 1e-9):
                    rate_in = INFUSION_RATE / self.vd
                    break
            curr_c += (rate_in - active_ke * curr_c) * dt
            conc[i] = max(curr_c, 0)
        return conc

    def fit_ke_from_levels(self, doses, times_dt, obs, sim_start, cr_func, patient_info, duration_days=7):
            """
            UPDATED: Accepts duration_days to ensure fitting window matches simulation window.
            """
            times_h = [(t - sim_start).total_seconds() / 3600 for t in times_dt]
            obs = np.array(obs)

            t_max = duration_days * 24 # Dynamic window
            t_grid = np.linspace(0, t_max, int(t_max * 12) + 1)
            
            # Pre-calculate once
            pop_ke_traj = np.array([
                pk_params_from_patient(
                    patient_info['age'], patient_info['sex'], 
                    patient_info['weight'], patient_info['height'], 
                    cr_func, sim_start + timedelta(hours=th)
                )[0] for th in t_grid
            ])

            mult_grid = np.linspace(0.3, 3.0, 100)
            log_post = []
            for m in mult_grid:
                res_conc = self._run_fast(doses, t_grid, pop_ke_traj, m)
                preds = np.interp(times_h, t_grid, res_conc)
                ll = -np.sum((obs - preds) ** 2 / (2 * 2.0 ** 2))
                prior = -((m - 1.0) ** 2) / (2 * 0.3 ** 2)
                log_post.append(ll + prior)

            idx = np.argmax(log_post)
            self.ke_multiplier = mult_grid[idx]
            d2 = np.gradient(np.gradient(log_post, mult_grid), mult_grid)
            self.multiplier_sd = np.sqrt(-1 / d2[idx]) if d2[idx] < 0 else 0.2
            return self.ke_multiplier

    # Confidence interval logic based on the fitted multiplier uncertainty
    def compute_ci(self, level=0.5):
            """
            Calculates the multiplier bounds for the confidence interval.
            level=0.5 corresponds to a 50% CI (Interquartile Range).
            """
            # z-score: 0.674 for 50% CI (IQR), 1.96 for 95% CI
            z = 0.674 if level == 0.5 else 1.96
            
            # Use the uncertainty (SD) calculated during the Bayesian fit
            # Default to 0.2 if no fit has been performed yet
            sd = getattr(self, 'multiplier_sd', 0.2)
            
            # Calculate bounds (ensuring they don't drop below a physical floor of 0.1)
            lo_mult = max(self.ke_multiplier - z * sd, 0.1)
            hi_mult = self.ke_multiplier + z * sd
            
            return lo_mult, hi_mult

# Align dose suggestions with fitted uncertainty
def dose_ci_from_ke(pk, target_auc, interval_h, n=1000):
    """
    Calculates dose confidence interval based on the Bayesian multiplier uncertainty.
    """
    mult_samples = np.random.normal(pk.ke_multiplier, pk.multiplier_sd, n)
    mult_samples = mult_samples[mult_samples > 0.1]

    doses = []
    for m in mult_samples:
        cl = (pk.ke * m) * pk.vd
        # Dose = (Target AUC * CL) * (interval / 24)
        dose_per_interval = (target_auc * cl) * (interval_h / 24)
        doses.append(dose_per_interval)

    return np.percentile(doses, 25), np.percentile(doses, 75)