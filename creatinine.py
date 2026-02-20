import numpy as np
from datetime import timedelta, datetime
from scipy.interpolate import interp1d

def calculate_kgfr(cr1, cr2, delta_t_hours, weight, age, sex):
    """Calculates Kinetic GFR (kGFR) for AKI settings."""
    k_sex = 0.85 if sex == "Female" else 1.0
    baseline_cl = ((140 - age) * 88.4) / (cr1) * k_sex * 0.9 
    
    v_dist = 0.6 # Fixed: Removed the * weight multiplier for dimensional consistency
    
    production_rate = baseline_cl * 0.06 * cr1
    delta_cr = cr2 - cr1
    
    if delta_t_hours <= 0 or production_rate <= 0:
        return baseline_cl
        
    kgfr = baseline_cl * (1 - (v_dist * delta_cr) / (production_rate * delta_t_hours))
    
    # Safety clamp to prevent negative or NaN kGFR
    if np.isnan(kgfr) or kgfr < 0.1:
        return 0.1
    return kgfr

def build_creatinine_function(cr_data, future_cr=None, modified_factor=1.0, patient_params=None):
    cr_data = sorted(list(set(cr_data)), key=lambda x: x[0])
    
    if modified_factor != 1.0:
        base_val = cr_data[0][1] * modified_factor
        return lambda t: (float(base_val), None)

    times = [d[0] for d in cr_data]
    values = [d[1] for d in cr_data]
    
    if future_cr is not None:
        times.append(times[-1] + timedelta(hours=48))
        values.append(future_cr)

    t_floats = [x.timestamp() for x in times]
    # If there is only 1 point, interp1d will crash, so we handle it gracefully
    if len(times) < 2:
        return lambda t: (float(values[0]), None)

    interp_func = interp1d(t_floats, values, kind='linear', fill_value="extrapolate")

    def cr_logic(t):
        if not isinstance(t, datetime):
            return float(values[0]), None
            
        current_val = float(interp_func(t.timestamp()))
        # Prevent extrapolated creatinine from dropping below physiological limits
        current_val = max(10.0, current_val)
        
        if patient_params is None:
            return current_val, None
            
        # FIX: Robust interval finding for kinetic GFR
        if t < times[0]:
            idx = 0
        elif t >= times[-1]:
            idx = len(times) - 2 # Use the most recent known interval for the future
        else:
            idx = 0
            for i in range(len(times)-1):
                if times[i] <= t < times[i+1]:
                    idx = i
                    break
        
        dt_hours = (times[idx+1] - times[idx]).total_seconds() / 3600
        kgfr = calculate_kgfr(values[idx], values[idx+1], dt_hours, 
                              patient_params['weight'], patient_params['age'], patient_params['sex'])
        return current_val, kgfr

    return cr_logic