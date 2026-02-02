import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from vanco_pk import VancoPK, pk_params_from_patient
from creatinine import build_creatinine_function
from dosing import (
    build_manual_doses,
    build_ordered_doses,
    suggest_regimen
)

# ---------------------------
# Streamlit setup
# ---------------------------
st.set_page_config(layout="centered")
st.title("Vancomycin PK Simulator")

# ---------------------------
# Patient inputs
# ---------------------------
st.header("Patient")

age = st.number_input("Age (years)", 18, 100, 65)
sex = st.radio("Sex", ["Male", "Female"], horizontal=True)
weight = st.slider("Weight (kg)", 30.0, 200.0, 75.0, 0.5)
height = st.slider("Height (cm)", 140.0, 230.0, 175.0, 0.5)

# ---------------------------
# Creatinine
# ---------------------------
st.header("Creatinine")

cr1 = st.slider("Initial Creatinine (µmol/L)", 35, 500, 100)
cr1_time = st.datetime_input("Initial Creatinine time", datetime.now() - timedelta(days=1))

# Initialize variables to prevent NameErrors
cr2 = None
cr2_time = None
multiplier = 1.0
use_cr_traj = False

cr_mode = st.radio(
    "Creatinine Simulation Type",
    ["Constant (Single Value)", "Two-Point Linear Trend", "Predicted Multiplier Trajectory"]
)

if cr_mode == "Two-Point Linear Trend":
    cr2 = st.slider("Second Creatinine (µmol/L)", 35, 500, 100)
    cr2_time = st.datetime_input("Second Creatinine time", datetime.now())
elif cr_mode == "Predicted Multiplier Trajectory":
    use_cr_traj = True
    multiplier = st.selectbox("Predicted change", [3, 2, 1.5, 1, 0.75, 0.5, 0.25], index=3)

cr_func = build_creatinine_function(
    cr1, cr1_time, cr2, cr2_time, use_cr_traj, multiplier
)

# ---------------------------
# Simulation settings
# ---------------------------
st.header("Simulation Start & Duration")
sim_start_date = st.date_input(
    "Simulation Start Date",
    datetime.now().date() - timedelta(days=1)
)
sim_start = datetime.combine(sim_start_date, datetime.min.time())

sim_duration_days = st.number_input(
    "Simulation duration (days)",
    min_value=1,
    max_value=30,
    value=7,
    step=1,
)
sim_end = sim_start + timedelta(days=sim_duration_days)

# ---------------------------
# Dose List Construction
# ---------------------------
manual_dose_inputs = []
manual_time_inputs = []

st.header("Manual Doses")
for i in range(5):
    if st.checkbox(f"Enable manual dose {i+1}", key=f"md_on_{i}"):
        dose = st.selectbox("Dose (mg)", [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500], index=3, key=f"md_dose_{i}")
        time = st.datetime_input("Time", sim_start + timedelta(hours=9, minutes=30), key=f"md_time_{i}")
        manual_dose_inputs.append(dose)
        manual_time_inputs.append(time)

st.header("Ordered Regimen")
show_ordered_dose = st.checkbox("Display ordered regimen", value=False)
ordered_dose, ordered_interval, ordered_start = None, None, None

if show_ordered_dose:
    ordered_dose = st.selectbox("Ordered dose (mg)", [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500], index=3)
    ordered_interval = st.selectbox("Interval (h)", [6, 8, 12, 18, 24, 36, 48, 72], index=2)
    ordered_start = st.datetime_input("Ordered dose start time", sim_start + timedelta(days=1, hours=9, minutes=30))
    if ordered_start < sim_start:
        st.error("Ordered dose start must be on or after simulation start.")
        st.stop()

doses = build_manual_doses(manual_dose_inputs, manual_time_inputs, sim_start)
if show_ordered_dose and ordered_dose:
    doses += build_ordered_doses(ordered_dose, ordered_interval, ordered_start, sim_start, sim_end)

# ---------------------------
# Measured levels
# ---------------------------
st.header("Measured Levels")
levels, level_times = [], []
for i in range(5):
    if st.checkbox(f"Enable level {i+1}", key=f"lvl_on_{i}"):
        lvl = st.number_input("Level (mg/L)", 0.0, 100.0, 15.0, step=0.1, key=f"lvl_val_{i}")
        t = st.datetime_input("Time", sim_start + timedelta(hours=9), key=f"lvl_time_{i}")
        levels.append(lvl)
        level_times.append(t)

# ---------------------------
# Core Simulation Logic
# ---------------------------
ke_pop, vd = pk_params_from_patient(age, sex, weight, height, cr_func, sim_start)
pk = VancoPK(ke_pop, vd)

if len(levels) >= 1:
    pk.fit_ke_from_levels(doses, level_times, levels, sim_start)
    st.info(f"Model fitted to {len(levels)} level(s).")
else:
    st.warning("Using population PK estimates.")

# Pack patient data for the simulation loop
p_info = {
    'age': age, 'sex': sex, 
    'weight': weight, 'height': height
}

# Run simulation
results = pk.run(
    doses, 
    duration_days=sim_duration_days, 
    sim_start=sim_start, 
    cr_func=cr_func, 
    patient_info=p_info
)

# ---------------------------
# Suggested Regimen
# ---------------------------
st.header("Try Regimen / Suggested Regimen")
suggested_dose, suggested_interval, predicted_auc = suggest_regimen(pk, target_auc=500)
st.markdown(f"**Suggested: {suggested_dose} mg q{suggested_interval}h** (AUC24 ≈ {predicted_auc:.0f})")

show_try_regimen = st.checkbox("Show try/suggested regimen on graph", value=False)
try_dose = st.selectbox("Try dose (mg)", [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500], 
                        index=[250,500,750,1000,1250,1500,1750,2000,2500].index(suggested_dose))
try_interval = st.selectbox("Try interval (h)", [6, 8, 12, 18, 24, 36, 48, 72], 
                            index=[6,8,12,18,24,36,48,72].index(suggested_interval))

# ---------------------------
# Plotting (Dual Axis)
# ---------------------------
fig, ax1 = plt.subplots(figsize=(12, 6))

t_dates = [sim_start + timedelta(hours=h) for h in results["time"]]
ax1.plot(t_dates, results["conc"], color="blue", linewidth=2, label="Vanco (Entered)")

# CI shading logic
# 1. Get the CI Multipliers (relative to the patient's curve)
mult_lo, mult_hi = pk.compute_ci(level=0.5)
current_fitted_mult = pk.ke_multiplier # Save the best fit to restore later

# 2. Simulate the LOWER bound (Higher clearance = Lower levels)
pk.ke_multiplier = mult_hi 
res_lo_bound = pk.run(doses, duration_days=sim_duration_days, sim_start=sim_start, 
                      cr_func=cr_func, patient_info=p_info)

# 3. Simulate the UPPER bound (Lower clearance = Higher levels)
pk.ke_multiplier = mult_lo
res_hi_bound = pk.run(doses, duration_days=sim_duration_days, sim_start=sim_start, 
                      cr_func=cr_func, patient_info=p_info)

# 4. RESTORE the best fit multiplier for the rest of the app logic
pk.ke_multiplier = current_fitted_mult

# 5. Plot the dynamic shaded area
ax1.fill_between(t_dates, res_lo_bound["conc"], res_hi_bound["conc"], 
                 color="blue", alpha=0.1, label="50% CI (Fitted Multiplier)")
# -------------------------------

try_results = None
if show_try_regimen:
    try_results = pk.simulate_regimen(
        try_dose, 
        try_interval, 
        sim_start, 
        sim_end, 
        cr_func=cr_func, 
        patient_info=p_info
    )
    
    ax1.plot([sim_start + timedelta(hours=h) for h in try_results["time"]], 
             try_results["conc"], color="green", linestyle="--", alpha=0.6, label="Try Regimen")

for lvl, t in zip(levels, level_times):
    ax1.scatter(t, lvl, color="#ff5733", s=60, zorder=5)

ax1.set_ylabel("Vancomycin (mg/L)", color="blue", fontweight='bold')
ax1.tick_params(axis='y', labelcolor="blue")

ax2 = ax1.twinx()
cr_vals = [cr_func(sim_start + timedelta(hours=h)) for h in results["time"]]
ax2.plot(t_dates, cr_vals, color="purple", linestyle=":", linewidth=1.5, alpha=0.6, label="Creatinine")

# Always plot the first creatinine point
ax2.scatter(cr1_time, cr1, color="purple", marker="x", s=60, label="Cr Point 1")

# Only plot the second point if the user is actually in "Two-Point" mode
if cr_mode == "Two-Point Linear Trend" and cr2_time is not None:
    ax2.scatter(cr2_time, cr2, color="purple", marker="x", s=60, label="Cr Point 2")

ax2.set_ylabel("Creatinine (µmol/L)", color="purple", fontweight='bold')
ax2.tick_params(axis='y', labelcolor="purple")
ax2.set_ylim(min(cr_vals)*0.8, max(cr_vals)*1.2)

ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax1.grid(alpha=0.2)
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1+h2, l1+l2, loc='upper left', frameon=True, framealpha=0.5)

st.pyplot(fig)

# ---------------------------
# Metrics with Color-Coding
# ---------------------------
def show_metrics(label, res):
    st.subheader(label)
    cols = st.columns(4)
    cols[0].metric("ke (1/h)", f"{res['ke']:.3f}")
    cols[1].metric("Half-life (h)", f"{res['half_life']:.1f}")
    cols[2].metric("Vd (L)", f"{res['vd']:.1f}")
    cols[3].metric("AUC24", f"{res['auc24']:.0f}")

    # Color coding logic for AUC24
    auc = res['auc24']
    if 400 <= auc <= 600:
        st.success(f"AUC24 of {auc:.0f} is within target range (400-600).")
    elif auc < 400:
        st.error(f"AUC24 of {auc:.0f} is below target range (< 400).")
    else:
        st.error(f"AUC24 of {auc:.0f} is above target range (> 600).")

show_metrics("Summary: Entered Regimen", results)
if try_results:
    show_metrics("Summary: Try Regimen", try_results)