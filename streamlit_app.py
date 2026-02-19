import streamlit as st
from datetime import datetime, timedelta
from vanco_pk import VancoPK, pk_params_from_patient
from creatinine import build_creatinine_function
from dosing import (
    build_manual_doses,
    build_ordered_doses,
    suggest_regimen
)
from plotting import plot_vanco_simulation

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
# Creatinine History & Modeling
# ---------------------------
st.header("Creatinine")

# 1. Initialize session state using the list-based approach
if 'cr_entries' not in st.session_state:
    # Default to one entry from 24 hours ago
    st.session_state.cr_entries = [{'val': 100.0, 'time': datetime.now() - timedelta(days=1)}]

# 2. UI to display and manage entries
st.subheader("Measured Levels")
for i, entry in enumerate(st.session_state.cr_entries):
    cols = st.columns([2, 2, 0.5])
    with cols[0]:
        # Update values directly in session state
        entry['val'] = st.number_input(f"Cr {i+1} (µmol/L)", value=entry['val'], key=f"cr_val_input_{i}")
    with cols[1]:
        # datetime_input is essential for PK accuracy
        entry['time'] = st.datetime_input(f"Time {i+1}", value=entry['time'], key=f"cr_time_input_{i}")
    with cols[2]:
        # Allow deletion of specific rows
        st.write("##") # Visual alignment for the button
        if st.button("🗑️", key=f"del_btn_{i}"):
            if len(st.session_state.cr_entries) > 1:
                st.session_state.cr_entries.pop(i)
                st.rerun()

if st.button("➕ Add Measured Creatinine"):
    # Append a new default entry based on the last entry's value
    last_val = st.session_state.cr_entries[-1]['val']
    st.session_state.cr_entries.append({'val': last_val, 'time': datetime.now()})
    st.rerun()

# 3. Overrides and Projections
st.divider()
st.subheader("Scenario Modeling")
col_a, col_b = st.columns(2)

with col_a:
    use_mod = st.toggle("Override with 'Modified' Creatinine", help="Ignores history and uses a single adjusted value.")
    # Uses 0.1x increments from 0.5x to 4x as requested
    mod_factor = st.slider("Multiplier Factor", 0.5, 4.0, 1.0, step=0.1, disabled=not use_mod)

with col_b:
    use_future = st.toggle("Project Future Trend", help="Extrapolates from the last measurement to a predicted future value.")
    future_val = st.number_input("Future Estimated Cr", value=100.0, disabled=not use_future)

# 4. Final Logic Assembly
# Format data for the build_creatinine_function
cr_data = [(e['time'], e['val']) for e in st.session_state.cr_entries]

# Call your updated backend function
cr_func = build_creatinine_function(
    cr_data=cr_data, 
    future_cr=future_val if use_future else None, 
    modified_factor=mod_factor if use_mod else 1.0
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
# Pack patient data for the simulation loop
# Ensure p_info is defined for the fitting logic
p_info = {
    'age': age,
    'sex': sex,
    'weight': weight,
    'height': height
}

ke_pop, vd = pk_params_from_patient(age, sex, weight, height, cr_func, sim_start)
pk = VancoPK(ke_pop, vd)

if len(levels) >= 1:
    pk.fit_ke_from_levels(doses, level_times, levels, sim_start, cr_func=cr_func, patient_info=p_info)
    st.info(f"Model fitted to {len(levels)} level(s).")
else:
    st.warning("Using population PK estimates.")

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

# Initialize try_results to prevent NameError
try_results = None 

if show_try_regimen:
    # end_dt is usually the start + the simulation duration
    sim_end = sim_start + timedelta(days=sim_duration_days)
    
    try_results = pk.simulate_regimen(
        dose_mg=try_dose, 
        interval_h=try_interval, 
        start_dt=sim_start, 
        end_dt=sim_end,
        cr_func=cr_func, 
        patient_info=p_info
    )

# ---------------------------
# Plotting (Dual Axis)
# ---------------------------
# 1. Prepare CI bounds if necessary
ci_bounds = None

if len(levels) >= 1:
    mult_lo, mult_hi = pk.compute_ci(level=0.5)
    
    current_fitted_mult = pk.ke_multiplier
    
    # Simulate bounds for shading
    pk.ke_multiplier = mult_hi 
    res_lo = pk.run(doses, duration_days=sim_duration_days, sim_start=sim_start, cr_func=cr_func, patient_info=p_info)
    
    pk.ke_multiplier = mult_lo
    res_hi = pk.run(doses, duration_days=sim_duration_days, sim_start=sim_start, cr_func=cr_func, patient_info=p_info)
    
    pk.ke_multiplier = current_fitted_mult
    ci_bounds = (res_lo, res_hi)

# 2. Call the modularized plotting function
fig = plot_vanco_simulation(
    sim_start=sim_start,
    results=results,
    cr_func=cr_func,
    levels=levels,
    level_times=level_times,
    try_results=try_results, # No longer needs 'if show_try_regimen else None'
    ci_bounds=ci_bounds
)

st.plotly_chart(fig, use_container_width=True)

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