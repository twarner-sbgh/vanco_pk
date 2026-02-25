from narwhals import when
import numpy as np
from dataclasses import dataclass
from datetime import datetime, timedelta # datetime used for type hinting or future extensions

INFUSION_RATE = 1000.0  # mg/hour - fixed infusion rate for all doses (1g/hr)
LN2 = np.log(2)

# Steady-state concentration calculations
def calculate_ss_conc(ke, vd, dose, interval, infusion_rate=INFUSION_RATE):
    """Calculates steady-state peak and trough concentrations based on PK parameters."""
    if ke <= 0 or vd <= 0 or interval <= 0:
        return 0.0, 0.0
    
    t_inf = dose / infusion_rate
    cl = ke * vd
    
    # SS Peak formula for intermittent infusion
    cpk_ss = (infusion_rate / cl) * (1 - np.exp(-ke * t_inf)) / (1 - np.exp(-ke * interval))
    
    # SS Trough formula
    ctr_ss = cpk_ss * np.exp(-ke * (interval - t_inf))
    
    return cpk_ss, ctr_ss

@dataclass
class Dose:
    time: float
    amount: float
    ke: float = None

def calculate_ss_conc(ke, vd, dose, interval, infusion_rate=INFUSION_RATE):
    """Calculates steady-state peak and trough concentrations (Intermittent Infusion)."""
    if ke <= 0 or vd <= 0 or interval <= 0:
        return 0.0, 0.0
    
    t_inf = dose / infusion_rate
    cl = ke * vd
    
    # SS Peak formula: Cpk = (R/CL) * (1 - e^-kTinf) / (1 - e^-kTau)
    cpk_ss = (infusion_rate / cl) * (1 - np.exp(-ke * t_inf)) / (1 - np.exp(-ke * interval))
    
    # SS Trough formula: Ctr = Cpk * e^-k(Tau - Tinf)
    ctr_ss = cpk_ss * np.exp(-ke * (interval - t_inf))
    
    return cpk_ss, ctr_ss

def pk_params_from_patient(age, sex, weight, height, cr_func, when):
    cr_data = cr_func(when)
    # Ensure we have a valid number for Creatinine to avoid division by zero
    cr = cr_data[0] if isinstance(cr_data, (tuple, list)) else cr_data
    # Avoid division by zero or extremely high values that crush CrCl to 0
    cr = max(cr, 10.0)
    if cr is None or cr <= 0:
        cr = 88.4 # Fallback to a standard 1.0 mg/dL (88.4 umol/L) if data are missing

    # CrCL estimation using Cockcroft-Gault with serum-to-plasma correction
    crcl = (140 - age) * 88.4 / cr * 0.9 # Serum-to-plasma Cr correction
    if sex.lower().startswith("f"):
        crcl *= 0.85

    ke = 0.00083 * crcl + 0.0044

    # Safety floors/ceilings
    ke = max(min(ke, 0.693 / 5.5), 0.693 / 120)

    if sex.lower().startswith("m"):
        ibw = 50 + 0.9 * (height - 152)
    else:
        ibw = 45.5 + 0.9 * (height - 152)

    weight_for_vd = ibw + 0.4 * (weight - ibw) if weight > 1.25 * ibw else weight
    vd = 0.8 * weight_for_vd

    return {"ke": ke, "vd": vd, "crcl": crcl}

class VancoPK:
    def __init__(self, ke, vd):
        self.ke = ke
        self.vd = vd
        self.ke_multiplier = 1.0
        self.multiplier_sd = 0.2

    def run(self, doses, duration_days=7, sim_start=None, cr_func=None, patient_info=None):
       
       # Extract parameters from the dictionary provided by the UI
        if patient_info is None:
                raise ValueError("patient_info dictionary must be provided to run simulation.")

        age = patient_info['age']
        sex = patient_info['sex']
        weight = patient_info['weight']
        height = patient_info['height']
                
        t_max = duration_days * 24
        # CHANGE: Increase multiplier from 10 to 12 for 5-minute steps
        t_grid = np.linspace(0, t_max, int(t_max * 12) + 1)
        dt = t_grid[1] - t_grid[0] # This will now be 0.0833 (5 mins) instead of 0.1
        # Calculate how many indices represent 24 hours
        steps_per_24h = int(24 / dt)    

        # PRE-CALCULATE infusion duration (Outside the loop)
        dose_intervals = []
        for t_d, amt in doses:
            duration = amt / INFUSION_RATE
            dose_intervals.append((t_d, t_d + duration))

        # Determine current baseline ke from renal function
        pop_ke_traj = np.full_like(t_grid, self.ke)
        if cr_func and patient_info and sim_start:
            for i, t_h in enumerate(t_grid):
                # FIX: Access ['ke'] from the dictionary
                p_dict = pk_params_from_patient(
                    patient_info['age'], patient_info['sex'], 
                    patient_info['weight'], patient_info['height'], 
                    cr_func, sim_start + timedelta(hours=t_h)
                )
                pop_ke_traj[i] = pk_params_from_patient(
                    patient_info['age'], patient_info['sex'], 
                    patient_info['weight'], patient_info['height'], 
                    cr_func, sim_start + timedelta(hours=t_h)
                )['ke']

        self.ke = pop_ke_traj[-1] 

        conc = np.zeros_like(t_grid)
        curr_c = 0.0
        active_ke = self.ke # Initialize

        for i, t_h in enumerate(t_grid):
            current_dt = sim_start + timedelta(hours=t_h)
            cr_val, kgfr = cr_func(current_dt)
            
            # Use vd_safe to prevent division by zero in dose entry
            vd_safe = self.vd if (self.vd and self.vd > 0) else 50.0
            
            # DETERMINING ELIMINATION (active_ke includes the Bayesian multiplier)
            if kgfr is not None:
                active_ke = ((kgfr * 0.06) / vd_safe) * self.ke_multiplier
            else:
                active_ke = pop_ke_traj[i] * self.ke_multiplier
            
            # Ensure ke stays within physiological limits (Half-life 4h to 120h)
            active_ke = max(min(active_ke, 0.17), 0.005)

#           FIX: Ultimate safety net against NaN poisoning
            if np.isnan(active_ke):
                active_ke = 0.05 # Fallback to a safe minimum if math fails

            # Store the unmultiplied base ke so it can be combined with 
            # the multiplier elsewhere without double-counting.
            self.ke = active_ke / max(self.ke_multiplier, 0.01)
            
            # ODE UPDATE
            # DOSE INFUSION
            rate_in = 0.0
            for t_d, amt in doses:
                duration = amt / INFUSION_RATE
                if t_d <= t_h <= (t_d + duration + 1e-9):
                    rate_in = INFUSION_RATE / vd_safe
                    break
            
            # Change in Conc = (Rate In) - (Elimination Rate * Current Conc)
            curr_c += (rate_in - active_ke * curr_c) * dt
            conc[i] = max(curr_c, 0)

        if len(conc) >= steps_per_24h:
            # Traversal sum of the last 24 hours of concentrations
            # Multiplying by dt converts the sum of concentrations into Area Under the Curve
            auc24 = conc[-steps_per_24h:].sum() * dt
        else:
            # Fallback if the simulation is shorter than 24 hours
            auc24 = conc.sum() * dt

        return {
            "time": t_grid, 
            "conc": conc, 
            "auc24": auc24, 
            "ke": active_ke, 
            "vd": self.vd,
            "half_life": (0.693 / active_ke) if active_ke > 0 else 0
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
        Fits the ke_multiplier based on measured serum levels using Bayesian estimation.
        """
        times_h = [(t - sim_start).total_seconds() / 3600 for t in times_dt]
        obs = np.array(obs)

        t_max = duration_days * 24 
        t_grid = np.linspace(0, t_max, int(t_max * 12) + 1)
        
        # NEW: Build trajectory using the exact same logic as the simulator (kGFR-aware)
        base_ke_traj = np.zeros_like(t_grid)
        for i, th in enumerate(t_grid):
            current_dt = sim_start + timedelta(hours=th)
            cr_val, kgfr = cr_func(current_dt)
            
            if kgfr is not None:
                # Match the simulator's logic for base k_e
                base_ke_traj[i] = (kgfr * 0.06) / self.vd
            else:
                p_dict = pk_params_from_patient(
                    patient_info['age'], patient_info['sex'], 
                    patient_info['weight'], patient_info['height'], 
                    cr_func, current_dt
                )
                base_ke_traj[i] = p_dict['ke']

        mult_grid = np.linspace(0.3, 3.0, 100)
        log_post = []
        for m in mult_grid:
            res_conc = self._run_fast(doses, t_grid, base_ke_traj, m)
            # MISSING LOGIC TO RESTORE:
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