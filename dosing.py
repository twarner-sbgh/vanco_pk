from datetime import timedelta
import numpy as np

# --- Existing manual/ordered dose functions remain the same ---
def build_manual_doses(dose_list, time_list, sim_start):
    doses = []
    for dose, dt in zip(dose_list, time_list):
        t_hours = (dt - sim_start).total_seconds() / 3600
        doses.append((t_hours, dose))
    return doses

def build_ordered_doses(dose_mg, interval_h, start_dt, sim_start, sim_end):
    doses = []
    t = start_dt
    while t <= sim_end:
        hours_since_start = (t - sim_start).total_seconds() / 3600
        doses.append((hours_since_start, dose_mg))
        t += timedelta(hours=interval_h)
    return doses

def suggest_regimen(pk, target_auc=500, patient_info=None):
    """
    Suggests a dose + interval that targets AUC24 ≈ target_auc.
    """
    allowed_intervals = [6, 8, 12, 18, 24, 36, 48, 72]
    allowed_doses = np.array([250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500])

    # FIX: Use 'pk' (the variable passed in) instead of 'pk_model'
    # Use the current multiplier and Vd from the object
    cl = (pk.ke * pk.ke_multiplier) * pk.vd  
    vd = pk.vd                                          

    # Calculate half-life for clinical interval selection
    half_life = (np.log(2) * vd / cl) if cl > 0 else 24

    best = (1000, 12, 500) # Default fallback to avoid NoneType unpacking error
    best_err = float("inf")

    for interval in allowed_intervals:
        # AUC24 formula: (Dose * (24/Interval)) / Clearance
        target_daily_dose = target_auc * cl
        doses_per_day = 24 / interval
        dose_est = target_daily_dose / doses_per_day

        # Round to nearest allowed increment
        dose_rounded = int(allowed_doses[np.argmin(np.abs(allowed_doses - dose_est))])
        
        # Calculate predicted AUC for the rounded dose
        predicted_auc = (dose_rounded * doses_per_day) / cl

        # Penalty logic: 1. Error from target AUC, 2. Clinically odd intervals
        interval_penalty = abs(interval - half_life) * 2
        err = abs(predicted_auc - target_auc) + interval_penalty

        if err < best_err:
            best_err = err
            best = (dose_rounded, interval, predicted_auc)

    return best