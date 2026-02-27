from datetime import timedelta
import plotly.graph_objects as go

def plot_vanco_simulation(sim_start, results, cr_func, levels=None, level_times=None, try_results=None, ci_bounds=None, static_crcl=None):
    """
    Plots Vancomycin simulation with separate colors for CrCl and kGFR on the Y2 axis.
    """
    # Primary time grid for the main simulation and Cr/kGFR lines
    t_dates_main = [sim_start + timedelta(hours=h) for h in results["time"]]
    
    # Get dual data from cr_func (Expected return: (Cr, kGFR))
    cr_plot_data = [cr_func(d) for d in t_dates_main]
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
            name='CI (50%)',
            showlegend=True
        ))

    # 2. Vanco Concentration Line (Main)
    fig.add_trace(go.Scatter(
        x=t_dates_main, 
        y=results["conc"], 
        name="Vanco (Current)", 
        line=dict(color="blue", width=3)
    ))

    # 3. Measured Levels (High-contrast markers with text labels)
    if levels is not None and len(levels) > 0:
        # Create a formatted list of strings for the label showing value and time
        # Format: "15.0 mg/L<br>Feb 21, 14:30"
        level_texts = [f"{lvl:.1f} mg/L<br>{t.strftime('%b %d, %H:%M')}" for lvl, t in zip(levels, level_times)]
        
        fig.add_trace(go.Scatter(
            x=level_times, 
            y=levels, 
            mode="markers+text",       # ENABLE TEXT
            name="Measured Levels", 
            text=level_texts,          # APPLY THE FORMATTED STRINGS
            textposition="bottom center", # PLACE IT DIRECTLY UNDER THE MARKER
            textfont=dict(color="black", size=10), # KEEP IT READABLE BUT SMALL
            marker=dict(
                color="red", 
                size=8,                # CHANGED FROM 12 TO 8
                symbol="diamond",
                line=dict(width=1.5, color="black") # Slightly thinner outline to match the smaller size
            )
        ))

    # 4. Try Regimen Comparison
    if try_results is not None:
        t_dates_try = [sim_start + timedelta(hours=h) for h in try_results["time"]]
        fig.add_trace(go.Scatter(
            x=t_dates_try, 
            y=try_results["conc"], 
            name="Vanco ('Try' Regimen)", 
            line=dict(color="teal", width=2, dash="dot")
        ))

    # 5. Static Estimated CrCl 
    # This plots a horizontal reference line based on the last Cr lab value entered
    if static_crcl is not None:
        fig.add_trace(go.Scatter(
            x=[t_dates_main[0], t_dates_main[-1]],
            y=[static_crcl, static_crcl],
            name=f"Est CrCl (CG): {static_crcl:.0f} mL/min",
            mode="lines",  
            line=dict(color='indigo', width=1, dash='dot'),
            yaxis="y2"
        ))

    # 6. Kinetic GFR (Right Axis - darkorchid)
    if any(k is not None for k in kgfr_vals):
        fig.add_trace(go.Scatter(
            x=t_dates_main, 
            y=kgfr_vals, 
            name="Kinetic GFR", 
            line=dict(color="darkorchid", width=2, dash="dash"), 
            yaxis="y2"
        ))

    # --- STYLE & LEGEND ---
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(
            title=dict(text="Date/Time", font=dict(color="black")),
            tickfont=dict(color="black"),            
            type="date",     # Force Plotly to treat this as a timeline
            tickformat="%d %b", # Can return to date & time if needed with: "%d %b %H:%M"
            hoverformat="%b %d, %H:%M", # to display time on hover
            dtick=86400000, # ticks at 24 hr intervals
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
            zerolinecolor='lightgrey',
            rangemode='tozero'         # ALIGNS THE BOTTOM TO ZERO
        ),
        yaxis2=dict(
            title=dict(text="Estimated CrCl / Kinetic GFR (mL/min)", font=dict(color="indigo")),
            tickfont=dict(color="indigo"),
            overlaying="y",
            side="right",
            showline=True, 
            linecolor='black',
            showgrid=False,            # Keeps grid tied only to Vanco scale
            rangemode='tozero'         # MATCHES PRIMARY AXIS
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
        margin=dict(l=50, r=50, t=50, b=80) # Standard margins
    )
    
    return fig