from datetime import timedelta
import numpy as np

# ----------------------------
# Manual doses
# ----------------------------

def build_manual_doses(dose_list, time_list, sim_start):
    """
    dose_list: list of dose amounts (mg)
    time_list: list of datetime objects
    sim_start: datetime of simulation start

    returns list of (time_hours, dose_mg)
    """
    doses = []

    for dose, dt in zip(dose_list, time_list):
        t_hours = (dt - sim_start).total_seconds() / 3600
        doses.append((t_hours, dose))

    return doses

# ----------------------------
# Ordered doses
# ----------------------------

def build_ordered_doses(dose_mg, interval_h, start_dt, sim_start, sim_end):
    doses = []

    t = start_dt
    while t <= sim_end:
        hours_since_start = (t - sim_start).total_seconds() / 3600
        doses.append((hours_since_start, dose_mg))
        t += timedelta(hours=interval_h)

    return doses

# ----------------------------
# Try regimen
# ----------------------------



# ----------------------------
# Suggested regimen
# ----------------------------

def suggest_regimen(pk_model, target_auc=500):
    """
    Suggests a dose + interval that targets AUC24 ≈ target_auc.
    Uses the individual fitted clearance derived from population ke 
    and the Bayesian multiplier.

    Returns:
        dose_mg, interval_h, predicted_auc24
    """

    allowed_intervals = [6, 8, 12, 18, 24, 36, 48, 72]
    allowed_doses = np.array([250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500])

    # FIX: Calculate current individual clearance
    # cl = (Base Population ke * Multiplier) * Vd
    cl = (pk_model.ke * pk_model.ke_multiplier) * pk_model.vd  # L/h
    vd = pk_model.vd                                          # L

    # Calculate half-life for interval selection logic
    half_life = np.log(2) * vd / cl

    best = None
    best_err = float("inf")

    for interval in allowed_intervals:
        # Favor intervals close to half-life to minimize peak/trough swings
        interval_penalty = abs(interval - half_life)

        # AUC24 ≈ daily dose / clearance
        target_daily_dose = target_auc * cl

        # Convert to per-dose amount based on the interval
        doses_per_day = 24 / interval
        dose_est = target_daily_dose / doses_per_day

        # Round to nearest allowed dose (250 mg increments)
        dose_rounded = int(allowed_doses[np.argmin(np.abs(allowed_doses - dose_est))])

        # Predict AUC24 for this candidate regimen
        predicted_auc = (dose_rounded * doses_per_day) / cl

        # Scoring function: Balance AUC accuracy with clinical appropriateness (interval penalty)
        err = abs(predicted_auc - target_auc) + interval_penalty * 5

        if err < best_err:
            best_err = err
            best = (dose_rounded, interval, predicted_auc)

    return best
