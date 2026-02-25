import streamlit as st
from datetime import datetime, timedelta
from vanco_pk import VancoPK, pk_params_from_patient, calculate_ss_conc
from creatinine import build_creatinine_function
from dosing import (
    build_manual_doses,
    build_ordered_doses,
    suggest_regimen
)
from plotting import plot_vanco_simulation
import uuid

# ---------------------------
# Streamlit setup
# ---------------------------
st.set_page_config(layout="centered")
st.title("Vancomycin PK Simulator (under construction)")

# ---------------------------
# Simulation settings
# ---------------------------
sim_start_date = st.date_input(
    "Simulation Start Date",
    datetime.now().date() - timedelta(days=1)
)
sim_start = datetime.combine(sim_start_date, datetime.min.time())

duration_days = st.slider("Simulation Duration (Days)", 1, 14, 7)
sim_end = sim_start + timedelta(days=duration_days)

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
st.header("Plasma Creatinine")

# 1. Initialize session state
if 'cr_entries' not in st.session_state:
    st.session_state.cr_entries = [{
        'id': str(uuid.uuid4()),
        'val': 100,
        'time': datetime.now() - timedelta(days=1)
    }]
else:
    # Backfill IDs for old entries
    for entry in st.session_state.cr_entries:
        if 'id' not in entry:
            entry['id'] = str(uuid.uuid4())

# 2. UI with Sliders for multi-point entry
st.subheader("Measured PCr")
for i, entry in enumerate(st.session_state.cr_entries):
    cols = st.columns([3, 2, 2, 0.6])

    with cols[0]:
        st.markdown(f"**PCr {i+1} (µmol/L)**")
        entry['val'] = st.slider(
            "",
            min_value=35,
            max_value=500,
            value=int(entry['val']),
            step=1,
            key=f"cr_slider_input_{entry['id']}",
            label_visibility="collapsed"
        )

    with cols[1]:
        d = st.date_input(
            f"Date {i+1}",
            value=entry['time'].date(),
            key=f"cr_date_{entry['id']}"
        )

    with cols[2]:
        t = st.time_input(
            f"Time {i+1}",
            value=entry['time'].time(),
            key=f"cr_time_{entry['id']}"
        )
        entry['time'] = datetime.combine(d, t)

    with cols[3]:
        if i > 0:
            st.markdown("<div style='height:28px'></div>", unsafe_allow_html=True)
            if st.button("✖", key=f"del_{entry['id']}", help="Remove this measurement"):
                st.session_state.cr_entries.pop(i)
                st.rerun()

# Using a heavy plus sign (✚) and text. 
# It will use your Blue theme color if you click it, and look clean and native.
if st.button("✚ Add Measured PCr", key="add_cr_btn"):
    last_val = st.session_state.cr_entries[-1]['val']
    st.session_state.cr_entries.append({
        'id': str(uuid.uuid4()),
        'val': last_val,
        'time': datetime.now()
    })
    st.rerun()

# 3. Final Logic Assembly
cr_data = [(e['time'], e['val']) for e in st.session_state.cr_entries]

# We set future_cr to None and modified_factor to 1.0 (no change)
cr_func = build_creatinine_function(
    cr_data=cr_data, 
    future_cr=None, 
    modified_factor=1.0
)

# ---------------------------
# Dose List Construction
# ---------------------------
manual_dose_inputs = []
manual_time_inputs = []

st.header("Manual Vancomycin Doses")
for i in range(5):
    if st.checkbox(f"Enable manual dose {i+1}", key=f"md_on_{i}"):
        dose = st.selectbox("Dose (mg)", [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500], index=3, key=f"md_dose_{i}")
        
        # Fix for datetime input
        col1, col2 = st.columns(2)
        d_md = col1.date_input("Date", sim_start.date(), key=f"md_date_{i}")
        t_md = col2.time_input("Time", (sim_start + timedelta(hours=9, minutes=30)).time(), key=f"md_time_{i}")
        time = datetime.combine(d_md, t_md)
        
        manual_dose_inputs.append(dose)
        manual_time_inputs.append(time)

st.header("Ordered Vancomycin Regimen")
show_ordered_dose = st.checkbox("Display ordered regimen", value=False)
ordered_dose, ordered_interval, ordered_start = None, None, None

if show_ordered_dose:
    ordered_dose = st.selectbox("Ordered dose (mg)", [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500], index=3)
    ordered_interval = st.selectbox("Interval (h)", [6, 8, 12, 18, 24, 36, 48, 72], index=2)
    
    # Fix for datetime input
    col1, col2 = st.columns(2)
    d_ord = col1.date_input("Ordered Start Date", (sim_start + timedelta(days=1)).date())
    t_ord = col2.time_input("Ordered Start Time", (sim_start + timedelta(hours=9, minutes=30)).time())
    ordered_start = datetime.combine(d_ord, t_ord)

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
        
        # Fix for datetime input
        col1, col2 = st.columns(2)
        d_lvl = col1.date_input("Date", sim_start.date(), key=f"lvl_date_{i}")
        t_lvl = col2.time_input("Time", (sim_start + timedelta(hours=9)).time(), key=f"lvl_time_{i}")
        t = datetime.combine(d_lvl, t_lvl)
        
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

# Parameters for PK model are now returned as a dictionary
params = pk_params_from_patient(age, sex, weight, height, cr_func, sim_start)
ke_pop = params['ke']
vd_pop = params['vd']

cr_func = build_creatinine_function(
    cr_data=cr_data, 
    future_cr=None,           # Removes the "Future" projection
    modified_factor=1.0,      # Keeps measurements at 1:1 (no multiplier)
    patient_params=p_info     # Keeps the patient demographics/context
)

# Initialize the PK engine
pk = VancoPK(ke_pop, vd_pop)

if len(levels) >= 1:
    pk.fit_ke_from_levels(doses, level_times, levels, sim_start, cr_func=cr_func, patient_info=p_info)
    st.info(f"Model fitted to {len(levels)} level(s).")
else:
    st.warning("Using population PK estimates.")

# Run the simulation
results = pk.run(
    doses=doses,
    duration_days=duration_days,
    sim_start=sim_start,
    cr_func=cr_func,
    patient_info=p_info 
)

# ---------------------------
# Suggested Regimen
# ---------------------------
st.header("Try Regimen / Suggested Regimen")
# Pass the p_info dictionary we created earlier to the suggester
suggested_dose, suggested_interval, predicted_auc = suggest_regimen(pk, target_auc=500, patient_info=p_info)
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
    sim_end = sim_start + timedelta(days=duration_days)
    
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
    res_hi = pk.run(doses, duration_days=duration_days, sim_start=sim_start, cr_func=cr_func, patient_info=p_info)
    
    pk.ke_multiplier = mult_lo
    res_lo = pk.run(doses, duration_days=duration_days, sim_start=sim_start, cr_func=cr_func, patient_info=p_info)
    
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
def show_metrics(label, res, dose=None, interval=None):
    st.subheader(label)
    cols = st.columns(6) # Increased to 6 columns
    cols[0].metric("ke (1/h)", f"{res['ke']:.3f}")
    cols[1].metric("Half-life (h)", f"{res['half_life']:.1f}")
    cols[2].metric("Vd (L)", f"{res['vd']:.1f}")
    cols[3].metric("AUC24", f"{res['auc24']:.0f}")
    
    if dose and interval:
        # Use the finalized ke from the simulation results
        cpk, ctr = calculate_ss_conc(res['ke'], res['vd'], dose, interval)
        cols[4].metric("Cpkss (mg/L)", f"{cpk:.1f}")
        cols[5].metric("Ctrss (mg/L)", f"{ctr:.1f}")
    else:
        cols[4].metric("Cpkss", "N/A")
        cols[5].metric("Ctrss", "N/A")

    # Color coding logic for AUC24
    auc = res['auc24']
    if 400 <= auc <= 600:
        st.success(f"AUC24 of {auc:.0f} is within target range (400-600).")
    elif auc < 400:
        st.error(f"AUC24 of {auc:.0f} is below target range (< 400).")
    else:
        st.error(f"AUC24 of {auc:.0f} is above target range (> 600).")

show_metrics("Summary: Entered Regimen", results, 
             dose=ordered_dose if show_ordered_dose else None, 
             interval=ordered_interval if show_ordered_dose else None)

if try_results:
    show_metrics("Summary: Try Regimen", try_results, 
                 dose=try_dose, interval=try_interval)