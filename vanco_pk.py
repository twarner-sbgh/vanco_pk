from narwhals import when
import numpy as np
from dataclasses import dataclass
from datetime import datetime, timedelta

INFUSION_RATE = 1000.0
LN2 = np.log(2)

def calculate_ss_conc(ke, vd, dose, interval, infusion_rate=INFUSION_RATE):
    if ke <= 0 or vd <= 0 or interval <= 0:
        return 0.0, 0.0
    t_inf = dose / infusion_rate
    cl = ke * vd
    cpk_ss = (infusion_rate / cl) * (1 - np.exp(-ke * t_inf)) / (1 - np.exp(-ke * interval))
    ctr_ss = cpk_ss * np.exp(-ke * (interval - t_inf))
    return cpk_ss, ctr_ss

@dataclass
class Dose:
    time: float
    amount: float
    ke: float = None

def pk_params_from_patient(age, sex, weight, height, cr_func, when, muscle_factor=1.0):
    cr_data = cr_func(when)
    cr = cr_data[0] if isinstance(cr_data, (tuple, list)) else cr_data
    cr = max(cr, 10.0)
    if cr is None or cr <= 0:
        cr = 88.4
        
    # Apply muscle factor to the Cockcroft-Gault numerator
    k_sex = 0.85 if sex.lower().startswith("f") else 1.0
    crcl = ((140 - age) * 88.4 * k_sex * 0.9 * muscle_factor) / cr
    
    ke = 0.00083 * crcl + 0.0044
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

    def run(self, doses, duration_days=7, sim_start=None, cr_func=None, patient_info=None, mode="crcl"):
        if patient_info is None:
                raise ValueError("patient_info dictionary must be provided to run simulation.")

        age = patient_info['age']
        sex = patient_info['sex']
        weight = patient_info['weight']
        height = patient_info['height']
        muscle_factor = patient_info.get('muscle_factor', 1.0)
                
        t_max = duration_days * 24
        t_grid = np.linspace(0, t_max, int(t_max * 12) + 1)
        dt = t_grid[1] - t_grid[0]
        steps_per_24h = int(24 / dt)    

        dose_intervals = []
        for t_d, amt in doses:
            duration = amt / INFUSION_RATE
            dose_intervals.append((t_d, t_d + duration))

        pop_ke_traj = np.full_like(t_grid, self.ke)
        if cr_func and patient_info and sim_start:
            for i, t_h in enumerate(t_grid):
                pop_ke_traj[i] = pk_params_from_patient(
                    patient_info['age'], patient_info['sex'], 
                    patient_info['weight'], patient_info['height'], 
                    cr_func, sim_start + timedelta(hours=t_h),
                    muscle_factor=muscle_factor
                )['ke']

        self.ke = pop_ke_traj[-1] 

        conc = np.zeros_like(t_grid)
        curr_c = 0.0
        active_ke = self.ke

        ke_history = []

        for i, t_h in enumerate(t_grid):
            current_dt = sim_start + timedelta(hours=t_h)
            cr_val, kgfr = cr_func(current_dt)
            
            vd_safe = self.vd if (self.vd and self.vd > 0) else 50.0
            
            if mode == "kgfr" and kgfr is not None:
                active_ke = ((kgfr * 0.06) / vd_safe) * self.ke_multiplier
            else:
                active_ke = pop_ke_traj[i] * self.ke_multiplier
            
            active_ke = max(min(active_ke, 0.17), 0.005)
            
            if np.isnan(active_ke):
                active_ke = 0.05
                
            ke_history.append(active_ke) 

            self.ke = active_ke / max(self.ke_multiplier, 0.01)
            
            rate_in = 0.0
            for t_d, amt in doses:
                duration = amt / INFUSION_RATE
                if t_d <= t_h <= (t_d + duration + 1e-9):
                    rate_in = INFUSION_RATE / vd_safe
                    break
            
            curr_c += (rate_in - active_ke * curr_c) * dt
            conc[i] = max(curr_c, 0)

        last_24h_ke = np.mean(ke_history[-steps_per_24h:]) if len(ke_history) >= steps_per_24h else active_ke
        self.ke = last_24h_ke / max(self.ke_multiplier, 0.01)

        if len(conc) >= steps_per_24h:
            auc24 = conc[-steps_per_24h:].sum() * dt
        else:
            auc24 = conc.sum() * dt

        return {
            "time": t_grid, 
            "conc": conc, 
            "auc24": auc24, 
            "ke": active_ke, 
            "vd": vd_safe,
            "half_life": (np.log(2) / active_ke) if active_ke > 0 else 0
        }
                
    def simulate_regimen(self, dose_mg, interval_h, start_dt, end_dt, cr_func, patient_info, mode="crcl"):
            duration_h = (end_dt - start_dt).total_seconds() / 3600
            duration_days = duration_h / 24
            
            try_doses = []
            curr_t = 0.0
            while curr_t < duration_h:
                try_doses.append((curr_t, dose_mg))
                curr_t += interval_h
                
            return self.run(
                doses=try_doses, 
                duration_days=duration_days, 
                sim_start=start_dt, 
                cr_func=cr_func, 
                patient_info=patient_info,
                mode=mode
            )

    def _run_fast(self, doses, t_grid, ke_traj, multiplier):
        dt = t_grid[1] - t_grid[0]
        conc = np.zeros_like(t_grid)
        curr_c = 0.0
        for i, t_h in enumerate(t_grid):
            active_ke = ke_traj[i] * multiplier
            rate_in = 0.0
            for t_d, amt in doses:
                if t_d <= t_h <= (t_d + amt/INFUSION_RATE + 1e-9):
                    rate_in = INFUSION_RATE / self.vd
                    break
            curr_c += (rate_in - active_ke * curr_c) * dt
            conc[i] = max(curr_c, 0)
        return conc

    def fit_ke_from_levels(self, doses, times_dt, obs, sim_start, cr_func, patient_info, duration_days=7, mode="crcl"):
        times_h = [(t - sim_start).total_seconds() / 3600 for t in times_dt]
        obs = np.array(obs)

        t_max = duration_days * 24 
        t_grid = np.linspace(0, t_max, int(t_max * 12) + 1)
        
        base_ke_traj = np.zeros_like(t_grid)
        for i, th in enumerate(t_grid):
            current_dt = sim_start + timedelta(hours=th)
            cr_val, kgfr = cr_func(current_dt)
            muscle_factor = patient_info.get('muscle_factor', 1.0)
            
            if mode == "kgfr" and kgfr is not None:
                base_ke_traj[i] = (kgfr * 0.06) / self.vd
            else:
                p_dict = pk_params_from_patient(
                    patient_info['age'], patient_info['sex'], 
                    patient_info['weight'], patient_info['height'], 
                    cr_func, current_dt, muscle_factor=muscle_factor
                )
                base_ke_traj[i] = p_dict['ke']

        mult_grid = np.linspace(0.3, 3.0, 100)
        log_post = []
        for m in mult_grid:
            res_conc = self._run_fast(doses, t_grid, base_ke_traj, m)
            preds = np.interp(times_h, t_grid, res_conc)
            ll = -np.sum((obs - preds) ** 2 / (2 * 2.0 ** 2))
            prior = -((m - 1.0) ** 2) / (2 * 0.3 ** 2)
            log_post.append(ll + prior)

        idx = np.argmax(log_post)
        self.ke_multiplier = mult_grid[idx]
        d2 = np.gradient(np.gradient(log_post, mult_grid), mult_grid)
        self.multiplier_sd = np.sqrt(-1 / d2[idx]) if d2[idx] < 0 else 0.2
        return self.ke_multiplier

    def compute_ci(self, level=0.5):
            z = 0.674 if level == 0.5 else 1.96
            sd = getattr(self, 'multiplier_sd', 0.2)
            lo_mult = max(self.ke_multiplier - z * sd, 0.1)
            hi_mult = self.ke_multiplier + z * sd
            return lo_mult, hi_mult
