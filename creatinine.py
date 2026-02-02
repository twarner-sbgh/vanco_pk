import numpy as np
from datetime import timedelta


def build_creatinine_function(
    cr1,
    t1,
    cr2=None,
    t2=None,
    use_multiplier=False,
    multiplier=1.0,
):
    """
    Returns a function cr_func(t: datetime) -> creatinine (µmol/L)

    Rules:
    - If only cr1 is provided: constant creatinine
    - If cr1 + cr2 provided: linear trend between them, then slope-preserving
      asymptotic leveling after 4 days from t2
    - If use_multiplier=True (and no cr2): crude ~4-day trajectory using multiplier
    """

    # -------------------------
    # Case 1: single creatinine
    # -------------------------
    if cr2 is None and not use_multiplier:
        def cr_func(t):
            return float(cr1)
        return cr_func

    # -------------------------
    # Case 2: Crude multiplier trajectory (Asymptotic)
    # -------------------------
    if cr2 is None and use_multiplier:
        cr_target = cr1 * multiplier
        # k=0.03 provides a curve that is ~95% complete at 4 days (96 hours)
        k_decay = 0.03 

        def cr_func(t):
            dt_hours = (t - t1).total_seconds() / 3600
            if dt_hours <= 0:
                return float(cr1)
            # Asymptotic approach formula
            return cr_target + (cr1 - cr_target) * np.exp(-k_decay * dt_hours)

        return cr_func

    # -------------------------
    # Case 3: two creatinine values
    # -------------------------

    # Convert times to days
    total_days = (t2 - t1).total_seconds() / 86400
    if total_days <= 0:
        raise ValueError("Second creatinine must be after first")

    slope = (cr2 - cr1) / total_days  # µmol/L per day

    # Choose k so slope ~0 after ~4 days
    k = 1.5 / 4.0  # per day

    def cr_func(t):
        dt_days = (t - t1).total_seconds() / 86400

        # Before first value
        if dt_days <= 0:
            return cr1

        # Between values → linear
        if t <= t2:
            return cr1 + slope * dt_days

        # After second value → slope-preserving taper
        dt2 = (t - t2).total_seconds() / 86400
        return cr2 + (slope / k) * (1 - np.exp(-k * dt2))

    return cr_func
