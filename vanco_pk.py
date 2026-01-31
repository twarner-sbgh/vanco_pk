import numpy as np
import datetime

UMOL_L_PER_MG_DL = 88.4

def cockcroft_gault_si(age, weight_kg, plasma_creatinine_umol_L, sex):
    plasma_creatinine_umol_L = max(plasma_creatinine_umol_L, 40.0) 
    serum_creatinine_umol_L = plasma_creatinine_umol_L * 0.9
    scr_mg_dl = serum_creatinine_umol_L / UMOL_L_PER_MG_DL
    crcl = ((140 - age) * weight_kg) / (72 * scr_mg_dl)
    if sex.lower() == "female":
        crcl *= 0.85
    return min(crcl, 150.0)

def calculate_single_dose_conc(t, dose_mg, dose_start_h, infusion_h, cl, v):
    k = cl / v
    rate = dose_mg / infusion_h
    t_dose = t - dose_start_h
    if t_dose < 0: return 0.0
    elif t_dose <= infusion_h:
        return (rate / cl) * (1 - np.exp(-k * t_dose))
    else:
        c_end = (rate / cl) * (1 - np.exp(-k * infusion_h))
        return c_end * np.exp(-k * (t_dose - infusion_h))

def simulate_vanco_timed(
    age, sex, weight_kg, 
    base_cr_val, base_cr_dt,
    new_cr_val, new_cr_dt,
    dose_history,
    standing_dose_mg, standing_interval_h,
    sim_start_dt, sim_end_dt,
    measured_level=None # {'val': float, 'dt': datetime}
):
    total_hours = (sim_end_dt - sim_start_dt).total_seconds() / 3600
    time_h = np.arange(0, total_hours, 0.1)
    
    # --- Creatinine Trajectory ---
    base_offset = (base_cr_dt - sim_start_dt).total_seconds() / 3600
    new_offset = (new_cr_dt - sim_start_dt).total_seconds() / 3600
    dur = max(0.1, new_offset - base_offset)
    delta = new_cr_val - base_cr_val
    
    creat = np.full_like(time_h, base_cr_val, dtype=float)
    mask_ramp = (time_h >= base_offset) & (time_h <= new_offset)
    mask_post = (time_h > new_offset)
    creat[mask_ramp] += delta * ((time_h[mask_ramp] - base_offset) / dur)
    creat[mask_post] += delta

    # --- PK Calculation ---
    v = 0.7 * weight_kg
    
    # Calculate population CL vector
    cl_vec_pop = np.array([4.5 * (cockcroft_gault_si(age, weight_kg, c, sex)/100)**0.75 for c in creat])
    
    # --- Best Fit Logic ---
    bias = 1.0
    if measured_level and measured_level['val'] > 0:
        # Find time index of the lab
        lab_offset = (measured_level['dt'] - sim_start_dt).total_seconds() / 3600
        idx = np.abs(time_h - lab_offset).argmin()
        
        # Iterative solver to find the bias (multiplier) for CL that hits the lab level
        best_bias = 1.0
        min_error = float('inf')
        for test_bias in np.linspace(0.3, 3.0, 100):
            test_cl = cl_vec_pop * test_bias
            test_conc = 0
            for d in dose_history:
                rel = (d['dt'] - sim_start_dt).total_seconds() / 3600
                test_conc += calculate_single_dose_conc(lab_offset, d['mg'], rel, 1.5, test_cl[idx], v)
            
            error = abs(test_conc - measured_level['val'])
            if error < min_error:
                min_error = error
                best_bias = test_bias
        bias = best_bias

    final_cl_vec = cl_vec_pop * bias
    conc = np.zeros_like(time_h)

    # Sum all doses (History + Standing)
    all_doses = []
    for d in dose_history:
        all_doses.append({'mg': d['mg'], 't': (d['dt'] - sim_start_dt).total_seconds() / 3600})
    
    if standing_dose_mg > 0:
        last_t = max([d['t'] for d in all_doses]) if all_doses else 0
        curr_t = last_t + standing_interval_h
        while curr_t < total_hours:
            all_doses.append({'mg': standing_dose_mg, 't': curr_t})
            curr_t += standing_interval_h

    for d in all_doses:
        for i, t in enumerate(time_h):
            conc[i] += calculate_single_dose_conc(t, d['mg'], d['t'], 1.5, final_cl_vec[i], v)

    # Derived PK constants (Final state)
    final_cl = final_cl_vec[-1]
    ke = final_cl / v
    thalf = 0.693 / ke
    mask_24 = time_h >= (time_h[-1] - 24)
    auc_24 = np.trapezoid(conc[mask_24], time_h[mask_24]) if any(mask_24) else 0
    
    return {
        "time_h": time_h, "conc": conc, "creat": creat, 
        "auc24": auc_24, "cl": final_cl, "v": v, 
        "ke": ke, "thalf": thalf
    }