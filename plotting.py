import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import timedelta
import plotly.graph_objects as go

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
    Signature preserved for compatibility with existing codebase.
    """

    # Convert simulation time to datetimes
    t_dates = [sim_start + timedelta(hours=h) for h in results["time"]]

    fig = go.Figure()

    # 1. Main Vancomycin Curve
    fig.add_trace(
        go.Scatter(
            x=t_dates,
            y=results["conc"],
            mode="lines",
            name="Vanco (Entered)",
            line=dict(color="blue", width=3),
            hovertemplate="Time: %{x}<br>Vanco: %{y:.2f} mg/L<extra></extra>"
        )
    )

    # 2. Confidence Interval (if provided)
    if ci_bounds:
        res_lo, res_hi = ci_bounds

        # Upper bound
        fig.add_trace(
            go.Scatter(
                x=t_dates,
                y=res_hi["conc"],
                mode="lines",
                line=dict(width=0),
                showlegend=False,
                hoverinfo="skip"
            )
        )

        # Lower bound + fill
        fig.add_trace(
            go.Scatter(
                x=t_dates,
                y=res_lo["conc"],
                mode="lines",
                fill="tonexty",
                fillcolor="rgba(0, 0, 255, 0.15)",
                line=dict(width=0),
                name="50% CI (Fitted)",
                hovertemplate=(
                    "Time: %{x}<br>"
                    "CI Low: %{y:.2f} mg/L<extra></extra>"
                )
            )
        )

    # 3. Try Regimen
    if try_results:
        try_dates = [sim_start + timedelta(hours=h) for h in try_results["time"]]
        fig.add_trace(
            go.Scatter(
                x=try_dates,
                y=try_results["conc"],
                mode="lines",
                name="Try Regimen",
                line=dict(color="green", dash="dash"),
                hovertemplate="Time: %{x}<br>Vanco: %{y:.2f} mg/L<extra></extra>"
            )
        )

    # 4. Measured Levels
    if levels and level_times:
        fig.add_trace(
            go.Scatter(
                x=level_times,
                y=levels,
                mode="markers",
                name="Measured Level",
                marker=dict(color="#ff5733", size=10),
                hovertemplate=(
                    "Measured Level<br>"
                    "Time: %{x}<br>"
                    "Value: %{y:.2f} mg/L<extra></extra>"
                )
            )
        )

    # 5. Creatinine (Secondary Axis)
    cr_vals = [cr_func(sim_start + timedelta(hours=h)) for h in results["time"]]

    fig.add_trace(
        go.Scatter(
            x=t_dates,
            y=cr_vals,
            mode="lines",
            name="Creatinine",
            line=dict(color="purple", dash="dot"),
            yaxis="y2",
            hovertemplate="Time: %{x}<br>Cr: %{y:.1f} µmol/L<extra></extra>"
        )
    )

    # 6. Layout
    fig.update_layout(
        height=500,
        hovermode="x unified",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0
        ),
        xaxis=dict(
            title="Date",
            showgrid=True
        ),
        yaxis=dict(
            title="Vancomycin (mg/L)",
            showgrid=True
        ),
        yaxis2=dict(
            title="Creatinine (µmol/L)",
            overlaying="y",
            side="right",
            showgrid=False
        ),
        margin=dict(l=40, r=40, t=40, b=40)
    )

    return fig