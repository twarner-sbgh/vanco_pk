from datetime import timedelta
import plotly.graph_objects as go

def plot_vanco_simulation(sim_start, results, cr_func, levels=None, level_times=None, try_results=None, ci_bounds=None):
    """
    Plots Vancomycin simulation with separate colors for Cr and kGFR.
    Legend is placed inside the top-left corner of the plot area.
    """
    # Primary time grid for the main simulation and Cr/kGFR lines
    t_dates_main = [sim_start + timedelta(hours=h) for h in results["time"]]
    
    # Get dual data from cr_func (Expected return: (Cr, kGFR))
    cr_plot_data = [cr_func(d) for d in t_dates_main]
    cr_vals = [d[0] for d in cr_plot_data]
    kgfr_vals = [d[1] for d in cr_plot_data]

    fig = go.Figure()

    # 1. Confidence Interval (IQR Shadow)
    if ci_bounds:
        res_lo, res_hi = ci_bounds
        t_dates_ci = [sim_start + timedelta(hours=h) for h in res_lo["time"]]
        fig.add_trace(go.Scatter(
            x=t_dates_ci + t_dates_ci[::-1],
            y=list(res_hi["conc"]) + list(res_lo["conc"])[::-1],
            fill='toself',
            fillcolor='rgba(0, 0, 255, 0.1)',
            line=dict(color='rgba(255,255,255,0)'),
            name='Confidence Interval',
            showlegend=True
        ))

    # 2. Vanco Concentration Line (Main)
    fig.add_trace(go.Scatter(
        x=t_dates_main, 
        y=results["conc"], 
        name="Vanco (Current)", 
        line=dict(color="blue", width=3)
    ))

    # 3. Measured Levels (High-contrast markers)
    if levels is not None and len(levels) > 0:
        fig.add_trace(go.Scatter(
            x=level_times, 
            y=levels, 
            mode="markers", 
            name="Measured Levels", 
            marker=dict(
                color="red", 
                size=12, 
                symbol="diamond",
                line=dict(width=2, color="black")
            )
        ))

    # 4. Try Regimen Comparison
    if try_results is not None:
        t_dates_try = [sim_start + timedelta(hours=h) for h in try_results["time"]]
        fig.add_trace(go.Scatter(
            x=t_dates_try, 
            y=try_results["conc"], 
            name="Vanco ('Try' Regimen))", 
            line=dict(color="green", width=2, dash="dash")
        ))

    # 5. Creatinine (Right Axis - Orchid)
    fig.add_trace(go.Scatter(
        x=t_dates_main, 
        y=cr_vals, 
        name="Creatinine", 
        line=dict(color="orchid", width=1.5, dash="dot"), 
        yaxis="y2"
    ))

    # 6. Kinetic GFR (Right Axis - Indigo)
    if any(k is not None for k in kgfr_vals):
        fig.add_trace(go.Scatter(
            x=t_dates_main, 
            y=kgfr_vals, 
            name="Kinetic GFR", 
            line=dict(color="indigo", width=2, dash="dash"), 
            yaxis="y2"
        ))

    # --- STYLE & LEGEND ---
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(
            title="Date/Time",
            showline=True, 
            linecolor='black', 
            gridcolor='lightgrey',
            showgrid=True
        ),
        yaxis=dict(
            title=dict(text="Vanco (mg/L)", font=dict(color="blue")),
            tickfont=dict(color="blue"),
            showline=True, 
            linecolor='black',
            gridcolor='lightgrey',
            showgrid=True,
            zeroline=True,
            zerolinecolor='lightgrey'
        ),
        yaxis2=dict(
            title=dict(text="Cr (µmol/L) / kGFR (mL/min)", font=dict(color="indigo")),
            tickfont=dict(color="indigo"),
            overlaying="y",
            side="right",
            showline=True, 
            linecolor='black',
            showgrid=False, # Keeps grid tied only to Vanco scale
            rangemode='tozero'
        ),
        # Legend positioned inside the plot area at the top-left
        legend=dict(
            orientation="v", 
            yanchor="top",
            y=0.99, 
            xanchor="left",
            x=0.01,
            font=dict(color="black", size=11),
            bgcolor="rgba(255, 255, 255, 0.8)", # Slightly more opaque for readability
            bordercolor="black",
            borderwidth=1
        ),
        margin=dict(l=50, r=50, t=50, b=50) # Standard margins
    )
    
    return fig