import numpy as np
from datetime import timedelta
from scipy.interpolate import interp1d

def build_creatinine_function(cr_data, future_cr=None, modified_factor=1.0):
    """
    cr_data: List of tuples [(datetime, value), ...]
    future_cr: Optional float for a future estimate
    modified_factor: Multiplier for the first measured value (if requested)
    """
    # Sort data by time
    cr_data = sorted(cr_data, key=lambda x: x[0])
    
    # Case: User requested a "Modified" creatinine using the slider
    if modified_factor != 1.0:
        base_cr = cr_data[0][1]
        mod_cr = base_cr * modified_factor
        return lambda t: float(mod_cr)

    # Case: Multiple values or Future Estimate
    times = [d[0] for d in cr_data]
    values = [d[1] for d in cr_data]
    
    # If a future estimate is provided, anchor it at +48h from last measured
    if future_cr is not None:
        times.append(times[-1] + timedelta(hours=48))
        values.append(future_cr)

    def cr_func(t):
        # Before first measurement
        if t <= times[0]:
            return float(values[0])
        # After last measurement (carry forward)
        if t >= times[-1]:
            return float(values[-1])
        
        # Interpolate between values
        # Convert times to floats for interpolation
        t_floats = [x.timestamp() for x in times]
        f = interp1d(t_floats, values, kind='linear')
        return float(f(t.timestamp()))

    return cr_func