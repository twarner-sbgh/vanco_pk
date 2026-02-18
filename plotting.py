import plotly.graph_objects as go
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
    Plotly-based Vancomycin PK simulation plot.
    Enhanced for interactive teaching with fixed white background, 
    visible legend in dark mode, and labeled measured levels.
    """

    # Convert simulation time to datetimes
    t_dates = [sim_start + timedelta(hours=h) for h in results["time"]]
    cr_vals = [cr_func(sim_start + timedelta(hours=h)) for h in results["time"]]

    fig = go.Figure()

    # 1. Confidence Interval Shading (Plotted first so it's in the background)
    if ci_bounds:
        res_lo, res_hi = ci_bounds
        # Upper anchor (invisible)
        fig.add_trace(go.Scatter(
            x=t_dates, y=res_hi["conc"], mode="lines", line=dict(width=0),
            showlegend=False, legendgroup="ci", hoverinfo="skip"
        ))
        # Lower bound + Fill
        fig.add_trace(go.Scatter(
            x=t_dates, y=res_lo["conc"], mode="lines", fill="tonexty",
            fillcolor="rgba(0, 0, 255, 0.15)", line=dict(width=0),
            name="50% CI (Fitted)", legendgroup="ci",
            text=res_hi["conc"],
            hovertemplate="<b>CI Range</b>: %{y:.1f} - %{text:.1f} mg/L<extra></extra>"
        ))

    # 2. Main Vancomycin Curve
    fig.add_trace(go.Scatter(
        x=t_dates, y=results["conc"], mode="lines",
        name="Vanco (Entered)", line=dict(color="blue", width=3),
        hovertemplate="Time: %{x}<br>Vanco: %{y:.2f} mg/L<extra></extra>"
    ))

    # 3. "Try" Regimen Curve
    if try_results:
        try_dates = [sim_start + timedelta(hours=h) for h in try_results["time"]]
        fig.add_trace(go.Scatter(
            x=try_dates, y=try_results["conc"], mode="lines",
            name="Try Regimen", line=dict(color="green", width=2, dash="dash"),
            hovertemplate="Time: %{x}<br>Try Conc: %{y:.2f} mg/L<extra></extra>"
        ))

    # 4. Measured Levels (WITH LABELS)
    if levels and level_times:
        # Create labels: "Value @ Time"
        level_labels = [f"<b>{val}</b><br>{t.strftime('%H:%M')}" for val, t in zip(levels, level_times)]
        
        fig.add_trace(go.Scatter(
            x=level_times,
            y=levels,
            mode="markers+text",  # Added +text
            name="Measured Level",
            marker=dict(color="#ff5733", size=12, line=dict(color="white", width=2)),
            text=level_labels,    # Data labels
            textposition="top center",
            textfont=dict(color="black", size=11),
            hovertemplate="Measured Level<br>Time: %{x}<br>Value: %{y} mg/L<extra></extra>"
        ))

    # 5. Creatinine Trend (Secondary Axis)
    fig.add_trace(go.Scatter(
        x=t_dates, y=cr_vals, mode="lines",
        name="Creatinine", line=dict(color="purple", width=1.5, dash="dot"),
        yaxis="y2", hovertemplate="Time: %{x}<br>Cr: %{y:.1f} µmol/L<extra></extra>"
    ))

    # 6. Layout Fixes (White BG, Visible Legend, Properly Scaled Axes)
    fig.update_layout(
        height=500,
        paper_bgcolor="white",
        plot_bgcolor="white",
        hovermode="x unified",
        # FIX: Explicitly set legend font color for Dark Mode compatibility
        legend=dict(
            orientation="h", yanchor="bottom", y=1.05, xanchor="left", x=0,
            font=dict(color="black"),
            bgcolor="rgba(255, 255, 255, 0.6)"
        ),
        xaxis=dict(
            title=dict(text="Date", font=dict(color="black")),
            tickfont=dict(color="black"),
            gridcolor="lightgrey",
            linecolor="black",
            range=[t_dates[0], t_dates[-1]], # Limit start to T=0
            rangemode="nonnegative"
        ),
        yaxis=dict(
            title=dict(text="Vancomycin (mg/L)", font=dict(color="blue")),
            tickfont=dict(color="blue"),
            gridcolor="lightgrey",
            linecolor="blue",
            rangemode="nonnegative",
            range=[0, max(results["conc"] + (levels or [0])) * 1.1] # Buffer for labels
        ),
        yaxis2=dict(
            title=dict(text="Creatinine (µmol/L)", font=dict(color="purple")),
            tickfont=dict(color="purple"),
            overlaying="y", side="right",
            showgrid=False,
            linecolor="purple",
            range=[min(cr_vals)*0.8, max(cr_vals)*1.2] if cr_vals else None
        ),
        margin=dict(l=50, r=50, t=80, b=50)
    )

# 7. Add subtle vertical lines to divide the days
    # Find the start of the simulation day
    sim_start_day = sim_start.replace(hour=0, minute=0, second=0, microsecond=0)
    
    # Calculate midnight for each day in the simulation range
    current_divider = sim_start_day + timedelta(days=1)
    last_date = t_dates[-1]
    
    while current_divider <= last_date:
        fig.add_vline(
            # FIX: Pass as string (ISO format) to avoid UTC offset issues in Plotly
            x=current_divider.strftime("%Y-%m-%d %H:%M:%S"),
            line_width=1,
            line_dash="dash",
            line_color="rgba(180, 180, 180, 0.3)", # Very subtle
            layer="below"
        )
        current_divider += timedelta(days=1)

    return fig