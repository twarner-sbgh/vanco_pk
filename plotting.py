from datetime import timedelta
import plotly.graph_objects as go

def plot_vanco_simulation(sim_start, results, cr_func, levels=None, level_times=None, try_results=None, ci_bounds=None, static_crcl=None, results_kgfr=None):
    """
    Plots Vancomycin simulation with separate colors for CrCl and kGFR on the Y2 axis.
    """
    t_dates_main = [sim_start + timedelta(hours=h) for h in results["time"]]
    
    cr_plot_data = [cr_func(d) for d in t_dates_main]
    kgfr_vals = [d[1] for d in cr_plot_data]
    valid_kgfrs = [k for k in kgfr_vals if k is not None]

    fig = go.Figure()

    # 1. Plot the Shadow Simulation (kGFR)
    if results_kgfr is not None:
        t_dates_kgfr = [sim_start + timedelta(hours=h) for h in results_kgfr["time"]]
        fig.add_trace(go.Scatter(
            x=t_dates_kgfr, y=results_kgfr["conc"],
            mode='lines', name='Vanco (Predicted with kGFR)',
            line=dict(color='rgba(173, 216, 230, 0.6)', width=4), # Pale Blue (LightBlue)
            hoverinfo='none' # Keep tooltip clean
        ))

    # 2. Confidence Interval (IQR Shadow)
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

    # 3. Vanco Concentration Line (Main/CrCl)
    fig.add_trace(go.Scatter(
        x=t_dates_main, 
        y=results["conc"], 
        name="Vanco (Current)", 
        mode='lines',
        line=dict(color="blue", width=3)
    ))

    # 4. Measured Levels (High-contrast markers with text labels)
    if levels is not None and len(levels) > 0:
        level_texts = [f"{lvl:.1f} mg/L<br>{t.strftime('%b %d, %H:%M')}" for lvl, t in zip(levels, level_times)]
        
        fig.add_trace(go.Scatter(
            x=level_times, 
            y=levels, 
            mode="markers+text",       
            name="Measured Levels", 
            text=level_texts,          
            textposition="bottom center", 
            textfont=dict(color="black", size=10), 
            marker=dict(
                color="red", 
                size=8,                
                symbol="diamond",
                line=dict(width=1.5, color="black") 
            )
        ))

    # 5. Try Regimen Comparison
    if try_results is not None:
        t_dates_try = [sim_start + timedelta(hours=h) for h in try_results["time"]]
        fig.add_trace(go.Scatter(
            x=t_dates_try, 
            y=try_results["conc"], 
            name="Vanco ('Try' Regimen)", 
            mode='lines',
            line=dict(color="#29b09d", width=1.5)
        ))

    # 6. Static Estimated CrCl 
    if static_crcl is not None:
        fig.add_trace(go.Scatter(
            x=[t_dates_main[0], t_dates_main[-1]],
            y=[static_crcl, static_crcl],
            name=f"Est CrCl (CG): {static_crcl:.0f} mL/min",
            mode="lines",  
            line=dict(color='indigo', width=1, dash='dash'),
            yaxis="y2"
        ))

    # 7. Kinetic GFR Line 
    valid_kgfrs = [k for k in kgfr_vals if k is not None]
    if valid_kgfrs:
        latest_kgfr = valid_kgfrs[-1]
        fig.add_trace(go.Scatter(
            x=t_dates_main, 
            y=kgfr_vals, 
            name=f"Kinetic GFR: {latest_kgfr:.0f} mL/min", 
            line=dict(color="darkorchid", width=2, dash="dot"), 
            yaxis="y2"
        ))

    # --- STYLE & LEGEND ---

    # Dynamic Y2 Axis Title
    y2_title = "Kinetic GFR (mL/min)" if len(valid_kgfrs) > 0 else "Est CrCl (mL/min)"

    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(
            title=dict(text="Date/Time", font=dict(color="black")),
            tickfont=dict(color="black"),            
            type="date",     
            tickformat="%d %b", 
            hoverformat="%b %d, %H:%M", 
            dtick=86400000, 
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
            rangemode='tozero'         
        ),
        yaxis2=dict(
            title=dict(text="Estimated CrCl / Kinetic GFR (mL/min)", font=dict(color="indigo")),
            tickfont=dict(color="indigo"),
            overlaying="y",
            side="right",
            showline=True, 
            linecolor='black',
            showgrid=False,            
            rangemode='tozero'         
        ),
        legend=dict(
            orientation="v", 
            yanchor="top",
            y=0.99, 
            xanchor="left",
            x=0.01,
            font=dict(color="black", size=11),
            bgcolor="rgba(255, 255, 255, 0.7)", 
            bordercolor="black",
            borderwidth=1
        ),
        margin=dict(l=40, r=40, t=40, b=40) 
    )
    
    return fig