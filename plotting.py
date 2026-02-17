import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import timedelta

def plot_vanco_simulation(
    sim_start, 
    results, 
    cr_func, 
    levels=None, 
    level_times=None, 
    try_results=None,
    ci_bounds=None
):
    """
    Handles all visual logic for the Vancomycin PK simulation.
    """
    fig, ax1 = plt.subplots(figsize=(12, 6))
    t_dates = [sim_start + timedelta(hours=h) for h in results["time"]]

    # 1. Plot Main Vancomycin Curve
    ax1.plot(t_dates, results["conc"], color="blue", linewidth=2, label="Vanco (Entered)")

    # 2. Plot Confidence Intervals if provided
    if ci_bounds:
        res_lo, res_hi = ci_bounds
        ax1.fill_between(t_dates, res_lo["conc"], res_hi["conc"], 
                         color="blue", alpha=0.1, label="50% CI (Fitted)")

    # 3. Plot "Try" Regimen if provided
    if try_results:
        try_dates = [sim_start + timedelta(hours=h) for h in try_results["time"]]
        ax1.plot(try_dates, try_results["conc"], color="green", 
                 linestyle="--", alpha=0.6, label="Try Regimen")

    # 4. Plot Measured Levels
    if levels and level_times:
        for lvl, t in zip(levels, level_times):
            ax1.scatter(t, lvl, color="#ff5733", s=60, zorder=5, label="Measured Level" if lvl == levels[0] else "")

    ax1.set_ylabel("Vancomycin (mg/L)", color="blue", fontweight='bold')
    ax1.tick_params(axis='y', labelcolor="blue")

    # 5. Secondary Axis for Creatinine
    ax2 = ax1.twinx()
    cr_vals = [cr_func(sim_start + timedelta(hours=h)) for h in results["time"]]
    ax2.plot(t_dates, cr_vals, color="purple", linestyle=":", linewidth=1.5, alpha=0.6, label="Creatinine")
    
    ax2.set_ylabel("Creatinine (µmol/L)", color="purple", fontweight='bold')
    ax2.tick_params(axis='y', labelcolor="purple")
    
    # Adjust y-limit for Cr to keep the line visible but not overwhelming
    if cr_vals:
        ax2.set_ylim(min(cr_vals)*0.8, max(cr_vals)*1.2)

    # 6. Formatting
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.grid(alpha=0.2)
    
    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', frameon=True, framealpha=0.5)

    return fig