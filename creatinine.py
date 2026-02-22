import numpy as np
from datetime import timedelta, datetime
from scipy.interpolate import interp1d

def calculate_kgfr(cr1, cr2, delta_t_hours, weight, height, age, sex):
    """Calculates Kinetic GFR (kGFR) for AKI settings using ABW for Volume of Distribution."""
    k_sex = 0.85 if sex == "Female" else 1.0
    # Baseline Cl (mL/min) - Cockcroft-Gault
    baseline_cl = ((140 - age) * 88.4) / (cr1) * k_sex * 0.9 # Serum-to-plasma correction
    
    # Calculate Ideal Body Weight (IBW) based on height (cm)
    if sex.lower().startswith("m") or sex == "Male":
        ibw = 50 + 0.9 * (height - 152)
    else:
        ibw = 45.5 + 0.9 * (height - 152)

    # Use adjusted body weight if total weight is > 1.25x IBW
    if weight > 1.25 * ibw:
        weight_for_vd = ibw + 0.4 * (weight - ibw)
    else:
        weight_for_vd = weight
        
    v_dist = 0.6 * weight_for_vd # Creatinine Vd in Liters
    
    # Production rate in umol/hr
    # Cl(mL/min) * 0.06 = L/hr. L/hr * umol/L = umol/hr
    production_rate = baseline_cl * 0.06 * cr1
    
    # Kinetic GFR formula
    delta_cr = cr2 - cr1
    if delta_t_hours <= 0 or production_rate <= 0:
        return baseline_cl
        
    kgfr = baseline_cl * (1 - (v_dist * delta_cr) / (production_rate * delta_t_hours))
    
    # Safety clamp to prevent negative or NaN kGFR
    if np.isnan(kgfr) or kgfr < 0.1:
        return 0.1
    return kgfr

def build_creatinine_function(cr_data, future_cr=None, modified_factor=1.0, patient_params=None):
    # Sort and remove duplicates to prevent interpolation errors
    cr_data = sorted(list(set(cr_data)), key=lambda x: x[0])
    
    # Modified factor (to allow for future projections or hypothetical Cr adjustments) left here in case UI element is added back in later
    if modified_factor != 1.0:
        base_val = cr_data[0][1] * modified_factor
        return lambda t: (float(base_val), None)

    times = [d[0] for d in cr_data]
    values = [d[1] for d in cr_data]
    
    t_floats = [x.timestamp() for x in times]
    
    # If there is only 1 point, interp1d will crash, so we handle it gracefully
    if len(times) < 2:
        return lambda t: (float(values[0]), None)
        
    interp_func = interp1d(t_floats, values, kind='linear', fill_value="extrapolate")

    def cr_logic(t):
        # Robustness check: if t is not a datetime, return first value
        if not isinstance(t, datetime):
            return float(values[0]), None
            
        current_val = float(interp_func(t.timestamp()))
        # Prevent extrapolated creatinine from dropping below physiological limits
        current_val = max(10.0, current_val)
        
        if patient_params is None:
            return current_val, None
            
        # Find the interval to calculate kGFR slope
        if t < times[0]:
            idx = 0
        elif t >= times[-1]:
            idx = len(times) - 2 # Use the most recent interval for the future
        else:
            idx = 0
            for i in range(len(times)-1):
                if times[i] <= t < times[i+1]:
                    idx = i
                    break
        
        dt_hours = (times[idx+1] - times[idx]).total_seconds() / 3600
        
        # UPDATED: Pass height down to calculate_kgfr for the IBW calculation
        kgfr = calculate_kgfr(
            values[idx], values[idx+1], dt_hours, 
            patient_params['weight'], patient_params['height'], 
            patient_params['age'], patient_params['sex']
        )
        return current_val, kgfr

    return cr_logic