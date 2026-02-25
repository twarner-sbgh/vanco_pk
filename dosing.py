from datetime import timedelta
import numpy as np

# --- Existing manual/ordered dose functions ---
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
    Suggests a dose + interval based on half-life and target AUC.
    """
    allowed_intervals = [6, 8, 12, 18, 24, 36, 48, 72]
    allowed_doses = np.array([500, 750, 1000, 1250, 1500, 1750, 2000, 2500])

    # Calculate current clearance (Cl = Ke * Vd)
    # Using pk.ke (baseline) * pk.ke_multiplier (Bayesian fit)
    cl = (pk.ke * pk.ke_multiplier) * pk.vd  
    vd = pk.vd                                          

    # 1. Calculate half-life
    half_life = (np.log(2) * vd / cl) if cl > 0 else 24
    
    # 2. Find interval closest to half-life
    idx = np.argmin([abs(i - half_life) for i in allowed_intervals])
    best_interval = allowed_intervals[idx]
    
    # 3. Rule: Extend to next interval if half-life > 1.2x the tested interval
    if half_life > 1.2 * best_interval and idx < len(allowed_intervals) - 1:
        best_interval = allowed_intervals[idx + 1]

    # 4. Find the standard dose closest to target AUC (500) for this interval
    # AUC24 = (Daily Dose) / Clearance
    doses_per_day = 24 / best_interval
    target_daily_dose = target_auc * cl
    dose_est = target_daily_dose / doses_per_day
    
    # Round to nearest allowed increment
    dose_rounded = int(allowed_doses[np.argmin(np.abs(allowed_doses - dose_est))])
    
    # Recalculate actual predicted AUC for the rounded dose
    predicted_auc = (dose_rounded * doses_per_day) / cl

    return (dose_rounded, best_interval, predicted_auc)