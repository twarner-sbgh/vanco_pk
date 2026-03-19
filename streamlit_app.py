import streamlit as st
import streamlit.components.v1 as components
import numpy as np
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
st.markdown("""
    <div style='text-align: center;'>
        <h1>Vancomycin PK Simulator</h1>
        <p style='color: grey;'>
    </div>
    """, unsafe_allow_html=True)

# DISCLAIMER
with st.expander("⚖️ Legal Disclaimer & Terms of Use"):
    st.caption("""
    By using this application, you acknowledge that:
    1. This tool is for educational and informational purposes only.
    2. This software is provided "as is" without warranties of any kind.
    3. Final dosing decisions are the sole responsibility of the prescribing clinician.
    4. Pharmacokinetic models are mathematical approximations. Always verify dosing calculations.
    """)

# Formatting of tabs
st.markdown("""
    <style>
    div[data-testid="stTabs"] button { flex: 1; width: 100%; }
    button[data-baseweb="tab"] {
        font-size: 20px !important; font-weight: 800 !important;
        background-color: #707070 !important; color: #FFFFFF !important;
        border-radius: 8px 8px 0 0 !important; margin: 4px !important;
        transition: background-color 0.3s ease;
    }
    button[aria-selected="true"] {
        background-color: #e1f5fe !important; color: #007bff !important;
        border-bottom: 5px solid #007bff !important;
    }
    button[data-baseweb="tab"]:hover {
        background-color: #505050 !important; color: #007bff !important;
    }
    div[data-testid="stTabs"] p { font-size: 19px !important; font-weight: 800 !important; }
    </style>
    """, unsafe_allow_html=True)

tab1, tab2 = st.tabs(["Patient Data & Dosing", "Results & Simulation"])

with tab1:
    # ---------------------------
    # Simulation settings
    # ---------------------------
    sim_start_date = st.date_input("Simulation Start Date", datetime.now().date() - timedelta(days=1))
    sim_start = datetime.combine(sim_start_date, datetime.min.time())

    # ---------------------------
    # Patient inputs
    # ---------------------------
    with st.container(border=True):
        st.header("Patient")
        age = st.slider("Age (years)", 17, 100, 65)
        sex = st.radio("Sex", ["Male", "Female"], horizontal=True)
        weight = st.slider("Weight (kg)", 30.0, 200.0, 75.0, 0.5)
        height = st.slider("Height (cm)", 140.0, 230.0, 175.0, 0.5)

    # ---------------------------
    # Plasma Creatinine
    # ---------------------------
    with st.container(border=True):
        st.header("Plasma Creatinine")
        MUSCLE_FACTORS = {
            "High (Athletic / High Muscle), 1.25x": 1.25,
            "Average": 1.0,
            "Low (Frail / Elderly / Mildly Cachectic), 0.75x": 0.75,
            "Very Low (Severe Sarcopenia / Paralysis / Bed-Bound), 0.5x": 0.5
        }
        muscle_mass_choice = st.selectbox("Presumed Muscle Mass", options=list(MUSCLE_FACTORS.keys()), index=1)
        selected_factor = MUSCLE_FACTORS[muscle_mass_choice]

        if 'cr_entries' not in st.session_state:
            st.session_state.cr_entries = [{'id': str(uuid.uuid4()), 'val': 100, 'time': datetime.now() - timedelta(days=1)}]

        st.subheader("Measured PCr")
        for i, entry in enumerate(st.session_state.cr_entries):
            st.markdown(f"**PCr {i+1} (µmol/L)**")
            entry['val'] = st.slider("", 35, 500, int(entry['val']), 1, key=f"cr_slider_input_{entry['id']}", label_visibility="collapsed")
            cols = st.columns([2, 2, 0.5])
            with cols[0]: d = st.date_input(f"Date {i+1}", value=entry['time'].date(), key=f"cr_date_{entry['id']}")
            with cols[1]: 
                t = st.time_input(f"Time {i+1}", value=entry['time'].time(), key=f"cr_time_{entry['id']}")
                entry['time'] = datetime.combine(d, t)
            with cols[2]:
                if i > 0:
                    st.markdown("<div style='height:28px'></div>", unsafe_allow_html=True)
                    if st.button("✖", key=f"del_{entry['id']}", help="Remove this measurement"):
                        st.session_state.cr_entries.pop(i)
                        st.rerun()
            st.divider()

        if st.button("✚ Add Measured PCr", key="add_cr_btn", help="Allows for estimation of kinetic GFR with changing renal function. For best results, add one additional PCr measurement at least 24 hours after the first"):
            last_val = st.session_state.cr_entries[-1]['val']
            st.session_state.cr_entries.append({'id': str(uuid.uuid4()), 'val': last_val, 'time': datetime.now()})
            st.rerun()

        cr_data = [(e['time'], e['val']) for e in st.session_state.cr_entries]
        p_info = {'age': age, 'sex': sex, 'weight': weight, 'height': height, 'muscle_factor': selected_factor}
        cr_func = build_creatinine_function(cr_data=cr_data, future_cr=None, modified_factor=1.0, patient_params=p_info)

    # ---------------------------
    # Dose List Construction
    # ---------------------------
    with st.container(border=True):
        st.header("Individual Vancomycin Doses")
        if 'manual_doses' not in st.session_state:
            st.session_state.manual_doses = [{'id': str(uuid.uuid4())}]

        manual_dose_inputs, manual_time_inputs = [], []
        for i, entry in enumerate(st.session_state.manual_doses):
            is_enabled = st.checkbox(f"Enable individual dose {i+1}", key=f"md_on_{entry['id']}")
            if is_enabled:
                c1, c2, c3 = st.columns([1, 1, 1])
                dose = c1.selectbox("Dose (mg)", [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500], index=3, key=f"md_val_{entry['id']}")
                d_md = c2.date_input("Date", sim_start.date(), key=f"md_d_{entry['id']}")
                t_md = c3.time_input("Time", (sim_start + timedelta(hours=9, minutes=30)).time(), key=f"md_t_{entry['id']}")
                manual_dose_inputs.append(dose)
                manual_time_inputs.append(datetime.combine(d_md, t_md))
                if i == len(st.session_state.manual_doses) - 1:
                    st.session_state.manual_doses.append({'id': str(uuid.uuid4())})
                    st.rerun()
            else:
                if len(st.session_state.manual_doses) > 1 and i < len(st.session_state.manual_doses) - 1:
                    st.session_state.manual_doses = st.session_state.manual_doses[:i+1]
                    st.rerun()

        st.header("Ordered Vancomycin Regimen")
        show_ordered_dose = st.checkbox("Display ordered regimen", value=False)
        ordered_dose, ordered_interval, ordered_start = None, None, None

        if show_ordered_dose:
            ordered_dose = st.selectbox("Ordered dose (mg)", [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500], index=3)
            ordered_interval = st.selectbox("Interval (h)", [6, 8, 12, 18, 24, 36, 48, 72], index=2)
            col1, col2 = st.columns(2)
            d_ord = col1.date_input("Ordered Start Date", (sim_start + timedelta(days=1)).date())
            t_ord = col2.time_input("Ordered Start Time", (sim_start + timedelta(hours=9, minutes=30)).time())
            ordered_start = datetime.combine(d_ord, t_ord)

        # Build doses using a 30-day max window so it supports the auto-extend
        max_sim_end = sim_start + timedelta(days=30)
        doses = build_manual_doses(manual_dose_inputs, manual_time_inputs, sim_start)
        if show_ordered_dose and ordered_dose:
            doses += build_ordered_doses(ordered_dose, ordered_interval, ordered_start, sim_start, max_sim_end)

    # ---------------------------
    # Measured levels
    # ---------------------------
    with st.container(border=True):
        st.header("Measured Vancomycin Levels")
        if 'level_entries' not in st.session_state:
            st.session_state.level_entries = [{'id': str(uuid.uuid4())}]

        levels, level_times = [], []
        for i, entry in enumerate(st.session_state.level_entries):
            is_enabled = st.checkbox(f"Enable level {i+1}", key=f"lvl_on_{entry['id']}")
            if is_enabled:
                lvl = st.number_input("Level (mg/L)", 0.0, 100.0, 15.0, step=0.1, key=f"lvl_val_{entry['id']}")
                col1, col2 = st.columns(2)
                d_lvl = col1.date_input("Date", sim_start.date(), key=f"lvl_date_{entry['id']}")
                t_lvl = col2.time_input("Time", (sim_start + timedelta(hours=9)).time(), key=f"lvl_time_{entry['id']}")
                levels.append(lvl)
                level_times.append(datetime.combine(d_lvl, t_lvl))
                if i == len(st.session_state.level_entries) - 1:
                    st.session_state.level_entries.append({'id': str(uuid.uuid4())})
                    st.rerun()
            else:
                if len(st.session_state.level_entries) > 1 and i < len(st.session_state.level_entries) - 1:
                    st.session_state.level_entries = st.session_state.level_entries[:i+1]
                    st.rerun()

    # --- PROCEED BUTTON (Using JS to switch tabs without breaking CSS) ---
    st.markdown("<br>", unsafe_allow_html=True)
    if st.button("Proceed to Results & Simulation ➡️", use_container_width=True, type="primary"):
        js = '''
        <script>
            var tabs = window.parent.document.querySelectorAll('button[data-baseweb="tab"]');
            if (tabs.length > 1) {
                tabs[1].click();
            }
        </script>
        '''
        components.html(js, height=0, width=0)

# ---------------------------
# Core Simulation Logic (Runs before Tab 2 UI)
# ---------------------------
params = pk_params_from_patient(age, sex, weight, height, cr_func, sim_start, muscle_factor=selected_factor)
pk = VancoPK(params['ke'], params['vd'])

if len(levels) >= 1:
    pk.fit_ke_from_levels(doses, level_times, levels, sim_start, cr_func=cr_func, patient_info=p_info, mode="crcl")
    fit_status_msg = f"Model fitted to {len(levels)} level(s)."
    is_fitted = True
else:
    fit_status_msg = "Using population PK estimates (no levels entered)."
    is_fitted = False

# Calculate max half-life for auto-duration calculation
effective_ke_crcl = pk.ke * pk.ke_multiplier
hl_crcl = (np.log(2) / effective_ke_crcl) if effective_ke_crcl > 0 else 24

hl_kgfr = 0
if len(st.session_state.cr_entries) >= 2:
    last_cr_time = st.session_state.cr_entries[-1]['time']
    _, latest_kgfr = cr_func(last_cr_time)
    if latest_kgfr is not None:
        vd_safe = pk.vd if (pk.vd and pk.vd > 0) else 50.0
        effective_ke_kgfr = ((latest_kgfr * 0.06) / vd_safe) * pk.ke_multiplier
        hl_kgfr = (np.log(2) / effective_ke_kgfr) if effective_ke_kgfr > 0 else 24

max_hl = max(hl_crcl, hl_kgfr)
auto_duration_days = int(np.ceil((5 * max_hl) / 24.0))
auto_duration_days = max(7, min(auto_duration_days, 30)) # Restrict between 7 and 30 days

# ---------------------------
# Results Tab
# ---------------------------
with tab2:
    # Full-width slider
    duration_days = st.slider(
        "Simulation Duration (Days)", 
        min_value=1, max_value=30, 
        value=auto_duration_days,
        help="Defaults to capturing at least 5 half-lives to show steady state (max 30 days)."
    )
    
    # Notification appearing below the slider
    if auto_duration_days > 7 and duration_days == auto_duration_days:
        st.info(f"⏳ Auto-extended to **{auto_duration_days} days** (5 × t½ of ~{max_hl:.1f}h).")

    sim_end = sim_start + timedelta(days=duration_days)

    # 3. Final Simulation Runs
    results_kgfr = None
    if len(st.session_state.cr_entries) >= 2:
        results_kgfr = pk.run(doses=doses, duration_days=duration_days, sim_start=sim_start, cr_func=cr_func, patient_info=p_info, mode="kgfr")

    results = pk.run(doses=doses, duration_days=duration_days, sim_start=sim_start, cr_func=cr_func, patient_info=p_info, mode="crcl")

    # ---------------------------
    # Suggestion & Alignment
    # ---------------------------
    with st.container(border=True):
        st.header("Try Regimen / Suggested Regimen")

        show_try_regimen = st.checkbox("Show try/suggested regimen on graph", value=False)
        
        use_kgfr_suggestion = False
        if results_kgfr is not None:
            use_kgfr_suggestion = st.checkbox("Use estimated PK parameters from kGFR", value=False)

        if use_kgfr_suggestion:
            base_ke_kgfr = results_kgfr['ke'] / max(pk.ke_multiplier, 0.01)
            pk_sugg = VancoPK(base_ke_kgfr, results_kgfr['vd'])
            pk_sugg.ke_multiplier = pk.ke_multiplier
            suggestion_mode = "kgfr"
        else:
            pk_sugg = pk
            suggestion_mode = "crcl"

        suggested_dose, suggested_interval, _ = suggest_regimen(pk_sugg, target_auc=500, patient_info=p_info)
        
        suggestion_sim = pk_sugg.simulate_regimen(
            suggested_dose, suggested_interval, sim_start, sim_end, cr_func, p_info, mode=suggestion_mode
        )
        simulated_suggested_auc = suggestion_sim['auc24']

        st.markdown(f"**Suggested: {suggested_dose} mg q{suggested_interval}h** (Simulated AUC24 ≈ {simulated_suggested_auc:.0f})")

        col1, col2 = st.columns(2)
        try_dose = col1.selectbox("Try dose (mg)", [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500], 
                                index=[250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500].index(suggested_dose))
        try_interval = col2.selectbox("Try interval (h)", [6, 8, 12, 18, 24, 36, 48, 72], 
                                    index=[6, 8, 12, 18, 24, 36, 48, 72].index(suggested_interval))

        if show_try_regimen:
            try_results = pk_sugg.simulate_regimen(try_dose, try_interval, sim_start, sim_end, cr_func, p_info, mode=suggestion_mode)
        else:
            try_results = None

    # ---------------------------
    # Plotting
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

    if st.session_state.cr_entries:
        last_entry = st.session_state.cr_entries[-1]
        static_params = pk_params_from_patient(age, sex, weight, height, cr_func, when=last_entry['time'], muscle_factor=selected_factor)
        current_static_crcl = static_params['crcl']
    else:
        current_static_crcl = None

    fig = plot_vanco_simulation(sim_start, results, cr_func, levels, level_times, try_results, ci_bounds, static_crcl=current_static_crcl, results_kgfr=results_kgfr)
    st.plotly_chart(fig, use_container_width=True)

    if is_fitted:
        st.info(fit_status_msg)
    else:
        st.warning(fit_status_msg)

    # ---------------------------
    # Metrics
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

        auc = res['auc24']
        if 400 <= auc <= 600:
            st.success(f"AUC24 of {auc:.0f} is within target range (400-600).")
        elif auc < 400:
            st.error(f"AUC24 of {auc:.0f} is below target range (< 400).")
        else:
            st.error(f"AUC24 of {auc:.0f} is above target range (> 600).")
    
    with st.container(border=True):
        show_metrics("Summary: Ordered Regimen", results, dose=ordered_dose if show_ordered_dose else None, interval=ordered_interval if show_ordered_dose else None)

    if try_results:
        with st.container(border=True):
            show_metrics("Summary: Try Regimen", try_results, dose=try_dose, interval=try_interval)

    if results_kgfr is not None:
        with st.container(border=True):
            show_metrics("Summary: Kinetic GFR", results_kgfr, dose=ordered_dose if show_ordered_dose else None, interval=ordered_interval if show_ordered_dose else None)
    
    st.markdown("<br>", unsafe_allow_html=True)
    if st.button("⬅️ Back to Patient Data & Dosing", use_container_width=True):
        js_back = '''
        <script>
            var tabs = window.parent.document.querySelectorAll('button[data-baseweb="tab"]');
            if (tabs.length > 0) {
                tabs[0].click();
            }
        </script>
        '''
        components.html(js_back, height=0, width=0)