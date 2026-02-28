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

# DISCLAIMER
with st.expander("⚖️ Legal Disclaimer & Terms of Use"):
    st.caption("""
    By using this application, you acknowledge that:
    1. This tool is for **educational and informational purposes only**. It is not a substitute for professional clinical judgment, institutional protocols, or the advice of a qualified healthcare provider.
    2. This software is provided "as is" without warranties of any kind.
    3. Final dosing decisions are the sole responsibility of the prescribing clinician. The developer assumes no liability for errors, omissions, or any clinical outcomes resulting from the use of this software.
    4. Pharmacokinetic models are mathematical approximations and may not account for all patient-specific variables. **Always verify dosing calculations and monitor serum levels according to your local guidelines and directives.**
    """)

# ---------------------------
# Streamlit setup
# ---------------------------
st.set_page_config(layout="centered")
st.title("Vancomycin PK Simulator")

# Formatting of tabs to make them larger and more distinct
st.markdown("""
    <style>
    /* Force tabs to take up equal width (50% each) */
    div[data-testid="stTabs"] button {
        flex: 1;
        width: 100%;
    }

    /* Style for ALL tabs (Inactive state) */
    button[data-baseweb="tab"] {
        font-size: 22px !important;    /* Larger text */
        font-weight: 800 !important;   /* Extra bold */
        background-color: #707070 !important; /* Darker medium-grey */
        color: #FFFFFF !important;     /* White text for contrast */
        border-radius: 8px 8px 0 0 !important;
        margin: 4px !important;
        transition: background-color 0.3s ease;
    }

    /* Style for the ACTIVE tab */
    button[aria-selected="true"] {
        background-color: #e1f5fe !important; /* Soft blue highlight */
        color: #007bff !important;           /* Bright blue text */
        border-bottom: 5px solid #007bff !important; /* Thicker accent line */
    }

    /* Hover effect for the inactive tabs */
    button[data-baseweb="tab"]:hover {
        background-color: #505050 !important; /* Even darker on hover */
        color: #007bff !important;
    }
    
    /* Ensure the text stays visible and large inside the button container */
    div[data-testid="stTabs"] p {
        font-size: 20px !important;
        font-weight: 800 !important;
    }
    </style>
    """, unsafe_allow_html=True)

tab1, tab2 = st.tabs(["Patient Data & Dosing", "Results & Simulation"])

with tab1:
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
    with st.container(border=True):
        st.header("Patient")

        age = st.number_input("Age (years)", 18, 100, 65)
        sex = st.radio("Sex", ["Male", "Female"], horizontal=True)
        weight = st.slider("Weight (kg)", 30.0, 200.0, 75.0, 0.5)
        height = st.slider("Height (cm)", 140.0, 230.0, 175.0, 0.5)

    # ---------------------------
    # Plasma Creatinine
    # ---------------------------
    with st.container(border=True):
        st.header("Plasma Creatinine")

        # Muscle Mass Factor - This is a simple categorical choice that adjusts the creatinine production rate in the model.
        MUSCLE_FACTORS = {
            "High (Athletic / High Muscle), 1.25x": 1.25,
            "Average": 1.0,
            "Low (Frail / Elderly / Mildly Cachectic), 0.75x": 0.75,
            "Very Low (Severe Sarcopenia / Paralysis / Bed-Bound), 0.5x": 0.5
        }

        muscle_mass_choice = st.selectbox(
            "Presumed Muscle Mass", 
            options=list(MUSCLE_FACTORS.keys()), 
            index=1,
            help="Adjusts presumed creatinine production used to generate kinetic GFR from two creatinine values."
        )
        selected_factor = MUSCLE_FACTORS[muscle_mass_choice]

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

        # Updated "add Cr" button with tooltip
        if st.button(
            "✚ Add Measured PCr", 
            key="add_cr_btn", 
            help="For best results, enter two measurements at least a day apart"
        ):
            last_val = st.session_state.cr_entries[-1]['val']
            st.session_state.cr_entries.append({
                'id': str(uuid.uuid4()),
                'val': last_val,
                'time': datetime.now(),         
            })
            st.rerun()

        # 3. Final Logic Assembly
        cr_data = [(e['time'], e['val']) for e in st.session_state.cr_entries]

        # --- Core Patient Data ---
        p_info = {
            'age': age,
            'sex': sex,
            'weight': weight,
            'height': height,
            'muscle_factor': selected_factor
        }

        cr_func = build_creatinine_function(
            cr_data=cr_data, 
            future_cr=None, 
            modified_factor=1.0,
            patient_params=p_info
        )

    # ---------------------------
    # Dose List Construction
    # ---------------------------
    with st.container(border=True):
        manual_dose_inputs = []
        manual_time_inputs = []

        st.header("Individual Vancomycin Doses")

        # Initialize session state if not present
        if 'manual_doses' not in st.session_state:
            st.session_state.manual_doses = [{'id': str(uuid.uuid4())}]

        manual_dose_inputs, manual_time_inputs = [], []

        for i, entry in enumerate(st.session_state.manual_doses):
            # Render the checkbox
            is_enabled = st.checkbox(f"Enable individual dose {i+1}", key=f"md_on_{entry['id']}")
            
            if is_enabled:
                c1, c2, c3 = st.columns([1, 1, 1])
                dose = c1.selectbox("Dose (mg)", [500, 750, 1000, 1250, 1500, 1750, 2000, 2500], index=2, key=f"md_val_{entry['id']}")
                d_md = c2.date_input("Date", sim_start.date(), key=f"md_d_{entry['id']}")
                t_md = c3.time_input("Time", (sim_start + timedelta(hours=9, minutes=30)).time(), key=f"md_t_{entry['id']}")
                
                # Collect the data for simulation
                manual_dose_inputs.append(dose)
                manual_time_inputs.append(datetime.combine(d_md, t_md))

                # If this is the last entry and it's checked, add a new empty one
                if i == len(st.session_state.manual_doses) - 1:
                    st.session_state.manual_doses.append({'id': str(uuid.uuid4())})
                    st.rerun()
                    
            else:
                # If a checkbox is unchecked AND it's not the only one in the list,
                # AND there are extra "empty" slots following it, prune the list.
                if len(st.session_state.manual_doses) > 1 and i < len(st.session_state.manual_doses) - 1:
                    # Remove all entries after this one to "collapse" the list
                    st.session_state.manual_doses = st.session_state.manual_doses[:i+1]
                    st.rerun()

        st.header("Ordered Vancomycin Regimen")
        show_ordered_dose = st.checkbox("Display ordered regimen", value=False)
        ordered_dose, ordered_interval, ordered_start = None, None, None

        if show_ordered_dose:
            ordered_dose = st.selectbox("Ordered dose (mg)", [500, 750, 1000, 1250, 1500, 1750, 2000, 2500], index=2)
            ordered_interval = st.selectbox("Interval (h)", [6, 8, 12, 18, 24, 36, 48, 72], index=2)
            
            # Datetime input
            col1, col2 = st.columns(2)
            d_ord = col1.date_input(
                "Ordered Start Date", 
                (sim_start + timedelta(days=1)).date(),
            )
            t_ord = col2.time_input(
                "Ordered Start Time", 
                (sim_start + timedelta(hours=9, minutes=30)).time(),
                help="Doses will be scheduled according to the entered interval starting from this time"
            )
            ordered_start = datetime.combine(d_ord, t_ord)

        doses = build_manual_doses(manual_dose_inputs, manual_time_inputs, sim_start)
        if show_ordered_dose and ordered_dose:
            doses += build_ordered_doses(ordered_dose, ordered_interval, ordered_start, sim_start, sim_end)

    # ---------------------------
    # Measured levels
    # ---------------------------
    with st.container(border=True):
        st.header("Measured Vancomycin Levels")

        # 1. Initialize session state for levels if not present
        if 'level_entries' not in st.session_state:
            st.session_state.level_entries = [{'id': str(uuid.uuid4())}]

        levels, level_times = [], []

        # 2. Iterate through the entries in session state
        for i, entry in enumerate(st.session_state.level_entries):
            # Render the checkbox
            is_enabled = st.checkbox(f"Enable level {i+1}", key=f"lvl_on_{entry['id']}")
            
            if is_enabled:
                lvl = st.number_input("Level (mg/L)", 0.0, 100.0, 15.0, step=0.1, key=f"lvl_val_{entry['id']}")
                
                col1, col2 = st.columns(2)
                d_lvl = col1.date_input("Date", sim_start.date(), key=f"lvl_date_{entry['id']}")
                t_lvl = col2.time_input("Time", (sim_start + timedelta(hours=9)).time(), key=f"lvl_time_{entry['id']}")
                t = datetime.combine(d_lvl, t_lvl)
                
                # Collect data for simulation logic
                levels.append(lvl)
                level_times.append(t)

                # EXPANSION: If this is the last entry and it was just checked, add a new empty slot
                if i == len(st.session_state.level_entries) - 1:
                    st.session_state.level_entries.append({'id': str(uuid.uuid4())})
                    st.rerun()
                    
            else:
                # CONTRACTION: If unchecked and it's not the only box, remove trailing empty boxes
                if len(st.session_state.level_entries) > 1 and i < len(st.session_state.level_entries) - 1:
                    st.session_state.level_entries = st.session_state.level_entries[:i+1]
                    st.rerun()

    # ---------------------------
    # Core Simulation Logic
    # ---------------------------
    # Parameters for PK model are now returned as a dictionary
    params = pk_params_from_patient(age, sex, weight, height, cr_func, sim_start, muscle_factor=selected_factor)
    ke_pop = params['ke']
    vd_pop = params['vd']

    # Initialize the PK engine
    pk = VancoPK(ke_pop, vd_pop)

    # 1. Primary Simulation (The Blue Plot - CrCl based)
    if len(levels) >= 1:
        pk.fit_ke_from_levels(doses, level_times, levels, sim_start, cr_func=cr_func, patient_info=p_info, mode="crcl")
        # Store message for later display
        fit_status_msg = f"Model fitted to {len(levels)} level(s)."
        is_fitted = True
    else:
        # Store message for later display
        fit_status_msg = "Using population PK estimates (no levels entered)."
        is_fitted = False

    # 2. Kinetic Simulation (The Shadow - kGFR based)
    results_kgfr = None
    if len(st.session_state.cr_entries) >= 2:
        results_kgfr = pk.run(doses=doses, duration_days=duration_days, sim_start=sim_start, cr_func=cr_func, patient_info=p_info, mode="kgfr")

    # 3. Standard Simulation (The Results - CrCl based)
    results = pk.run(
        doses=doses,
        duration_days=duration_days,
        sim_start=sim_start,
        cr_func=cr_func,
        patient_info=p_info,
        mode="crcl" 
    )


with tab2:
    # ---------------------------
    # Suggestion & Alignment
    # ---------------------------
    with st.container(border=True):
        st.header("Try Regimen / Suggested Regimen")

        # 1. Get the raw suggestion (Dose and Interval)
        suggested_dose, suggested_interval, _ = suggest_regimen(pk, target_auc=500, patient_info=p_info)

        # 2. RUN A SIMULATION for that specific suggestion to get the "Summary-style" AUC
        suggestion_sim = pk.simulate_regimen(
            suggested_dose, 
            suggested_interval, 
            sim_start, 
            sim_end, 
            cr_func, 
            p_info
        )
        simulated_suggested_auc = suggestion_sim['auc24']

        # 3. Display the label using the simulated value
        st.markdown(f"**Suggested: {suggested_dose} mg q{suggested_interval}h** (Simulated AUC24 ≈ {simulated_suggested_auc:.0f})")

        # ---------------------------
        # Try Regimen UI
        # ---------------------------
        show_try_regimen = st.checkbox("Show try/suggested regimen on graph", value=False)

        # Sync the 'Try' boxes to the suggestion by default
        try_dose = st.selectbox("Try dose (mg)", [500, 750, 1000, 1250, 1500, 1750, 2000, 2500], 
                                index=[500,750,1000,1250,1500,1750,2000,2500].index(suggested_dose))
        try_interval = st.selectbox("Try interval (h)", [6, 8, 12, 18, 24, 36, 48, 72], 
                                    index=[6,8,12,18,24,36,48,72].index(suggested_interval))

        if show_try_regimen:
            try_results = pk.simulate_regimen(try_dose, try_interval, sim_start, sim_end, cr_func, p_info)
        else:
            try_results = None

    # ---------------------------
    # Plotting (Dual Axis)
    # ---------------------------
    ci_bounds = None

    if len(levels) >= 1:
        mult_lo, mult_hi = pk.compute_ci(level=0.5)
        current_fitted_mult = pk.ke_multiplier
        
        pk.ke_multiplier = mult_hi 
        res_hi = pk.run(doses, duration_days=duration_days, sim_start=sim_start, cr_func=cr_func, patient_info=p_info, mode="crcl")
        
        pk.ke_multiplier = mult_lo
        res_lo = pk.run(doses, duration_days=duration_days, sim_start=sim_start, cr_func=cr_func, patient_info=p_info, mode="crcl")
        
        pk.ke_multiplier = current_fitted_mult
        ci_bounds = (res_lo, res_hi)

    # --- Calculate Static CrCl for inclusion on the Y2 Axis ---
    if st.session_state.cr_entries:
        last_entry = st.session_state.cr_entries[-1]
        static_params = pk_params_from_patient(
            age=age, sex=sex, weight=weight, height=height, 
            cr_func=cr_func, when=last_entry['time'], muscle_factor=selected_factor
        )
        current_static_crcl = static_params['crcl']
    else:
        current_static_crcl = None

    # Pass everything to the chart
    fig = plot_vanco_simulation(
        sim_start, 
        results, 
        cr_func, 
        levels, 
        level_times, 
        try_results, 
        ci_bounds, 
        static_crcl=current_static_crcl, 
        results_kgfr=results_kgfr  
    )

    st.plotly_chart(fig, use_container_width=True)

    # Disclaimer re pop PK vs fitted model after the plot
    if is_fitted:
        st.info(fit_status_msg)
    else:
        st.warning(fit_status_msg)

    # ---------------------------
    # Metrics with Color-Coding
    # ---------------------------
    def show_metrics(label, res, dose=None, interval=None):
        if dose and interval:
            st.subheader(f"{label} ({dose:.0f} mg q{interval:.0f}h)")
        else:
            st.subheader(label)
            
        cols = st.columns(6) 
        cols[0].metric("ke (1/h)", f"{res['ke']:.3f}")
        cols[1].metric("Half-life (h)", f"{res['half_life']:.1f}")
        cols[2].metric("Vd (L)", f"{res['vd']:.1f}")
        cols[3].metric("AUC24", f"{res['auc24']:.0f}")
        
        if dose and interval:
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

    show_metrics(
        "Summary: Ordered Regimen", 
        results, 
        dose=ordered_dose if show_ordered_dose else None, 
        interval=ordered_interval if show_ordered_dose else None
    )

    if try_results:
        show_metrics(
            "Summary: Try Regimen", 
            try_results, 
            dose=try_dose, 
            interval=try_interval
        )
