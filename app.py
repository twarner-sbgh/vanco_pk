import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
from vanco_pk import simulate_vanco_timed, cockcroft_gault_si

st.set_page_config(layout="wide", page_title="Vanco Precision Planner")

with st.sidebar:
    st.header("1. Patient & Standing Dose")
    age = st.number_input("Age", 18, 100, 65)
    sex = st.selectbox("Sex", ["Male", "Female"])
    weight = st.number_input("Weight (kg)", 40, 200, 80)
    std_mg = st.number_input("Standing Dose (mg)", 0, 3000, 1000, step=250)
    std_int = st.selectbox("Frequency (h)", [8, 12, 24, 48], index=1)

    st.header("2. Renal Trajectory")
    st.caption("Base Measurement")
    bc_val = st.number_input("Base Cr", 40, 600, 100)
    bc_dt = datetime.datetime.combine(st.date_input("Base Date", datetime.date.today() - datetime.timedelta(days=1)), 
                                     st.time_input("Base Time", datetime.time(8, 0)))
    
    st.caption("New Measurement")
    nc_val = st.number_input("New Cr", 40, 600, 100)
    nc_dt = datetime.datetime.combine(st.date_input("New Date", datetime.date.today()), 
                                     st.time_input("New Time", datetime.time(8, 0)))

col_inputs, col_graph = st.columns([1, 2])

with col_inputs:
    st.subheader("Historical Doses")
    hist_doses = []
    for i in range(4):
        with st.expander(f"Dose {i+1}"):
            d_dt = datetime.datetime.combine(st.date_input("Date", datetime.date.today(), key=f"d{i}"),
                                           st.time_input("Time", datetime.time(8,0), key=f"t{i}"))
            d_mg = st.number_input("mg", 0, 3000, 0, step=250, key=f"m{i}")
            if d_mg > 0: hist_doses.append({'dt': d_dt, 'mg': d_mg})

    st.subheader("Measured Serum Level")
    lab_val = st.number_input("Level (mg/L)", 0.0, 100.0, 0.0)
    lab_dt = datetime.datetime.combine(st.date_input("Lab Date", datetime.date.today()),
                                      st.time_input("Lab Time", datetime.time(12, 0)))

# Simulation Window
sim_start = min([bc_dt, nc_dt] + ([d['dt'] for d in hist_doses] if hist_doses else [datetime.datetime.now()]))
sim_end = datetime.datetime.now() + datetime.timedelta(days=2)

res = simulate_vanco_timed(
    age, sex, weight, bc_val, bc_dt, nc_val, nc_dt,
    hist_doses, std_mg, std_int, sim_start, sim_end,
    measured_level={'val': lab_val, 'dt': lab_dt} if lab_val > 0 else None
)

with col_graph:
    # Graph
    fig, ax1 = plt.subplots(figsize=(10, 5))
    plot_times = [sim_start + datetime.timedelta(hours=h) for h in res["time_h"]]
    ax1.plot(plot_times, res["conc"], color="#1f77b4", label="Predicted Level")
    if lab_val > 0:
        ax1.scatter([lab_dt], [lab_val], color="red", label="Measured Level (Fit Target)", zorder=5)
    
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d %H:%M'))
    ax1.set_ylabel("Vanco (mg/L)")
    ax1.axhspan(10, 20, color='green', alpha=0.05)
    plt.xticks(rotation=45)
    st.pyplot(fig)

    # Metrics
    st.subheader("Predicted Outcomes")
    m1, m2, m3 = st.columns(3)
    
    # Color coding for AUC
    auc = res['auc24']
    auc_color = "normal" if 400 <= auc <= 600 else "inverse"
    auc_label = "✅ In Range" if 400 <= auc <= 600 else "⚠️ Out of Range"
    
    # Use Markdown for colored metrics because st.metric color is limited
    m1.metric("Predicted AUC24", f"{auc:.0f}", delta=auc_label, delta_color=auc_color)
    m2.metric("Est. Vd (L)", f"{res['v']:.1f} L")
    m3.metric("Elimination (Ke)", f"{res['ke']:.3f} h⁻¹")

    m4, m5, m6 = st.columns(3)
    m4.metric("Half-life (t½)", f"{res['thalf']:.1f} h")
    m5.metric("Final CL", f"{res['cl']:.2f} L/h")
    m6.metric("Final CrCl", f"{cockcroft_gault_si(age, weight, nc_val, sex):.1f} mL/min")