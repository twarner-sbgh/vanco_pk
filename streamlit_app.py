"""
Vancomycin PK Simulator — lightweight clinical practice tool.

Architecture note:
    The app is a two-step wizard (Patient Data -> Results) driven entirely by
    st.session_state["view"]. Only the active view's widgets are rendered on any
    given run, so the (relatively expensive) PK simulations execute *only* when
    the Results view is shown — not on every input keystroke.

    Persistence across the view switch is handled by keeping the authoritative
    data in our own session_state structures (ss.patient, ss.cr_entries,
    ss.dose_entries, ss.level_entries, ss.ordered) rather than relying on widget
    state, which Streamlit clears for widgets that aren't rendered.
"""

import uuid
from datetime import datetime, timedelta

import numpy as np
import streamlit as st

from vanco_pk import VancoPK, pk_params_from_patient, calculate_ss_conc
from creatinine import build_creatinine_function
from dosing import build_manual_doses, build_ordered_doses, suggest_regimen
from plotting import plot_vanco_simulation

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DOSE_OPTIONS = [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500]
INTERVAL_OPTIONS = [6, 8, 12, 18, 24, 36, 48, 72]
MUSCLE_FACTORS = {
    "High (Athletic / High Muscle), 1.25x": 1.25,
    "Average": 1.0,
    "Low (Frail / Elderly / Mildly Cachectic), 0.75x": 0.75,
    "Very Low (Severe Sarcopenia / Paralysis / Bed-Bound), 0.5x": 0.5,
}
AUC_LOW, AUC_HIGH = 400, 600      # target AUC24 window
AUC_TARGET = 500                  # midpoint used for regimen suggestion
DEFAULT_DOSE_HOUR = timedelta(hours=9, minutes=30)
DEFAULT_LEVEL_HOUR = timedelta(hours=9)

st.set_page_config(layout="centered")
ss = st.session_state

# ---------------------------------------------------------------------------
# Session-state defaults (authoritative data lives here, not in widget state)
# ---------------------------------------------------------------------------
# Streamlit preserves session_state across code reloads/redeploys. If an older
# version of the app populated these collections with a different entry schema,
# reusing them would crash (missing keys) or, worse, silently feed phantom
# measured levels into the Bayesian fit. Bumping SCHEMA_VERSION clears the
# affected collections once, so upgrades always start from a clean, valid state.
SCHEMA_VERSION = 2
if ss.get("_schema") != SCHEMA_VERSION:
    for _stale in ("patient", "cr_entries", "dose_entries", "level_entries", "ordered"):
        ss.pop(_stale, None)
    ss["_schema"] = SCHEMA_VERSION

ss.setdefault("view", "input")
ss.setdefault("sim_start_date", datetime.now().date() - timedelta(days=1))
ss.setdefault("patient", {
    "age": 65, "sex": "Male", "weight": 75.0, "height": 175.0,
    "muscle_choice": "Average",
})
ss.setdefault("cr_entries", [
    {"id": str(uuid.uuid4()), "val": 100, "time": datetime.now() - timedelta(days=1)},
])
ss.setdefault("dose_entries", [{"id": str(uuid.uuid4()), "dose": 1000, "time": None}])
ss.setdefault("level_entries", [])          # measured levels are optional
ss.setdefault("ordered", {
    "show": False, "dose": 1000, "interval": 12, "start": None,
})


def go_to(view):
    ss.view = view
    st.rerun()


def _normalize_entries():
    """Guarantee every persisted entry dict has its required keys.

    Safety net against partially-formed entries (e.g. state left by an older
    build). Time defaults are left as None here and filled lazily where
    sim_start is known.
    """
    for e in ss.cr_entries:
        e.setdefault("id", str(uuid.uuid4()))
        e.setdefault("val", 100)
        e.setdefault("time", datetime.now() - timedelta(days=1))
    for e in ss.dose_entries:
        e.setdefault("id", str(uuid.uuid4()))
        e.setdefault("dose", 1000)
        e.setdefault("time", None)
    for e in ss.level_entries:
        e.setdefault("id", str(uuid.uuid4()))
        e.setdefault("lvl", 15.0)
        e.setdefault("time", None)


# ===========================================================================
# HEADER (shown on both views)
# ===========================================================================
st.markdown(
    "<div style='text-align:center'><h1>Vancomycin PK Simulator</h1></div>",
    unsafe_allow_html=True,
)

with st.expander("⚖️ Legal Disclaimer & Terms of Use"):
    st.caption(
        "By using this application, you acknowledge that:\n"
        "1. This tool is for educational and informational purposes only.\n"
        "2. This software is provided \"as is\" without warranties of any kind.\n"
        "3. Final dosing decisions are the sole responsibility of the prescribing clinician.\n"
        "4. Pharmacokinetic models are mathematical approximations. Always verify dosing calculations."
    )


# ===========================================================================
# VIEW 1 — PATIENT DATA & DOSING
# ===========================================================================
def render_input_view():
    st.caption("Step 1 of 2 — Patient Data & Dosing")

    # --- Simulation start date -------------------------------------------
    sim_start_date = st.date_input("Simulation Start Date", value=ss.sim_start_date)
    ss.sim_start_date = sim_start_date
    sim_start = datetime.combine(sim_start_date, datetime.min.time())

    # --- Patient ----------------------------------------------------------
    pt = ss.patient
    with st.container(border=True):
        st.header("Patient")
        pt["age"] = st.slider("Age (years)", 17, 100, pt["age"])
        pt["sex"] = st.radio("Sex", ["Male", "Female"],
                             index=["Male", "Female"].index(pt["sex"]), horizontal=True)
        pt["weight"] = st.slider("Weight (kg)", 30.0, 200.0, pt["weight"], 0.5)
        pt["height"] = st.slider("Height (cm)", 140.0, 230.0, pt["height"], 0.5)

    # --- Plasma creatinine ------------------------------------------------
    with st.container(border=True):
        st.header("Plasma Creatinine")
        pt["muscle_choice"] = st.selectbox(
            "Presumed Muscle Mass", options=list(MUSCLE_FACTORS),
            index=list(MUSCLE_FACTORS).index(pt["muscle_choice"]),
        )
        muscle_factor = MUSCLE_FACTORS[pt["muscle_choice"]]

        st.subheader("Measured PCr")
        for i, e in enumerate(ss.cr_entries):
            st.markdown(f"**PCr {i + 1} (µmol/L)**")
            e["val"] = st.slider(f"PCr {i + 1} (µmol/L)", 35, 500, int(e["val"]), 1,
                                 key=f"cr_val_{e['id']}", label_visibility="collapsed")
            c1, c2, c3 = st.columns([2, 2, 0.5])
            d = c1.date_input(f"Date {i + 1}", value=e["time"].date(), key=f"cr_d_{e['id']}")
            t = c2.time_input(f"Time {i + 1}", value=e["time"].time(), key=f"cr_t_{e['id']}")
            e["time"] = datetime.combine(d, t)
            if i > 0:
                c3.markdown("<div style='height:28px'></div>", unsafe_allow_html=True)
                if c3.button("✖", key=f"cr_del_{e['id']}", help="Remove this measurement"):
                    ss.cr_entries.pop(i)
                    st.rerun()
            st.divider()

        if st.button("✚ Add Measured PCr", key="cr_add",
                     help="Enables kinetic GFR estimation with changing renal function. "
                          "For best results add a second PCr at least 24 h after the first."):
            ss.cr_entries.append({"id": str(uuid.uuid4()),
                                  "val": ss.cr_entries[-1]["val"], "time": datetime.now()})
            st.rerun()

    # --- Individual doses -------------------------------------------------
    with st.container(border=True):
        st.header("Individual Vancomycin Doses")
        for i, e in enumerate(ss.dose_entries):
            if e["time"] is None:
                e["time"] = sim_start + DEFAULT_DOSE_HOUR
            c1, c2, c3, c4 = st.columns([1.2, 1.2, 1.2, 0.4])
            e["dose"] = c1.selectbox("Dose (mg)", DOSE_OPTIONS,
                                     index=DOSE_OPTIONS.index(e["dose"]), key=f"md_v_{e['id']}")
            d = c2.date_input("Date", e["time"].date(), key=f"md_d_{e['id']}")
            t = c3.time_input("Time", e["time"].time(), key=f"md_t_{e['id']}")
            e["time"] = datetime.combine(d, t)
            c4.markdown("<div style='height:28px'></div>", unsafe_allow_html=True)
            if c4.button("✖", key=f"md_del_{e['id']}", help="Remove this dose"):
                ss.dose_entries.pop(i)
                st.rerun()

        if st.button("✚ Add Dose", key="md_add"):
            last = ss.dose_entries[-1]["time"] if ss.dose_entries else sim_start + DEFAULT_DOSE_HOUR
            ss.dose_entries.append({"id": str(uuid.uuid4()), "dose": 1000,
                                    "time": last + timedelta(hours=12)})
            st.rerun()

    # --- Ordered regimen --------------------------------------------------
    with st.container(border=True):
        st.header("Ordered Vancomycin Regimen")
        od = ss.ordered
        od["show"] = st.checkbox("Display ordered regimen", value=od["show"])
        if od["show"]:
            od["dose"] = st.selectbox("Ordered dose (mg)", DOSE_OPTIONS,
                                      index=DOSE_OPTIONS.index(od["dose"]))
            od["interval"] = st.selectbox("Interval (h)", INTERVAL_OPTIONS,
                                          index=INTERVAL_OPTIONS.index(od["interval"]))
            if od["start"] is None:
                od["start"] = sim_start + timedelta(days=1) + DEFAULT_DOSE_HOUR
            c1, c2 = st.columns(2)
            d = c1.date_input("Ordered Start Date", od["start"].date())
            t = c2.time_input("Ordered Start Time", od["start"].time())
            od["start"] = datetime.combine(d, t)

    # --- Measured levels --------------------------------------------------
    with st.container(border=True):
        st.header("Measured Vancomycin Levels")
        for i, e in enumerate(ss.level_entries):
            if e["time"] is None:
                e["time"] = sim_start + DEFAULT_LEVEL_HOUR
            e["lvl"] = st.number_input(f"Level {i + 1} (mg/L)", 0.0, 100.0,
                                       float(e["lvl"]), step=0.1, key=f"lvl_v_{e['id']}")
            c1, c2, c3 = st.columns([2, 2, 0.5])
            d = c1.date_input(f"Level date {i + 1}", e["time"].date(), key=f"lvl_d_{e['id']}")
            t = c2.time_input(f"Level time {i + 1}", e["time"].time(), key=f"lvl_t_{e['id']}")
            e["time"] = datetime.combine(d, t)
            c3.markdown("<div style='height:28px'></div>", unsafe_allow_html=True)
            if c3.button("✖", key=f"lvl_del_{e['id']}", help="Remove this level"):
                ss.level_entries.pop(i)
                st.rerun()
            st.divider()

        if st.button("✚ Add Measured Level", key="lvl_add"):
            ss.level_entries.append({"id": str(uuid.uuid4()), "lvl": 15.0, "time": None})
            st.rerun()

    # --- Auto-rewind sim start date to the earliest entered date ----------
    all_dates = [e["time"].date() for e in ss.dose_entries]
    all_dates += [e["time"].date() for e in ss.level_entries]
    if ss.ordered["show"] and ss.ordered["start"]:
        all_dates.append(ss.ordered["start"].date())
    if all_dates:
        earliest = min(all_dates)
        if earliest < sim_start_date:
            ss.sim_start_date = earliest
            st.warning(f"⚠️ **Simulation Start Date auto-adjusted** to "
                       f"{earliest.strftime('%b %d, %Y')} to accommodate earlier input.")
            st.rerun()

    # --- Proceed ----------------------------------------------------------
    st.markdown("<br>", unsafe_allow_html=True)
    if st.button("Proceed to Results & Simulation ➡️", width="stretch", type="primary"):
        go_to("results")


# ===========================================================================
# VIEW 2 — RESULTS & SIMULATION
# ===========================================================================
def render_results_view():
    top_l, top_r = st.columns([1, 1])
    with top_l:
        if st.button("⬅️ Back to Patient Data", width="stretch"):
            go_to("input")
    top_r.caption("Step 2 of 2 — Results & Simulation")

    # --- Rebuild inputs from persisted state ------------------------------
    pt = ss.patient
    muscle_factor = MUSCLE_FACTORS[pt["muscle_choice"]]
    p_info = {"age": pt["age"], "sex": pt["sex"], "weight": pt["weight"],
              "height": pt["height"], "muscle_factor": muscle_factor}

    sim_start = datetime.combine(ss.sim_start_date, datetime.min.time())
    cr_data = [(e["time"], e["val"]) for e in ss.cr_entries]
    cr_func = build_creatinine_function(cr_data=cr_data, patient_params=p_info)

    doses = build_manual_doses([e["dose"] for e in ss.dose_entries],
                               [e["time"] for e in ss.dose_entries], sim_start)
    od = ss.ordered
    if od["show"]:
        max_sim_end = sim_start + timedelta(days=30)   # wide window supports auto-extend
        doses += build_ordered_doses(od["dose"], od["interval"], od["start"],
                                     sim_start, max_sim_end)

    levels = [e["lvl"] for e in ss.level_entries]
    level_times = [e["time"] for e in ss.level_entries]

    # --- Base PK + Bayesian fit to levels ---------------------------------
    params = pk_params_from_patient(pt["age"], pt["sex"], pt["weight"],
                                    pt["height"], cr_func, sim_start,
                                    muscle_factor=muscle_factor)
    pk = VancoPK(params["ke"], params["vd"])

    if levels:
        pk.fit_ke_from_levels(doses, level_times, levels, sim_start,
                              cr_func=cr_func, patient_info=p_info, mode="crcl")
        is_fitted, fit_msg = True, f"Model fitted to {len(levels)} level(s)."
    else:
        is_fitted, fit_msg = False, "Using population PK estimates (no levels entered)."

    # --- Auto simulation duration (>= 5 half-lives, clamped 7-30 d) -------
    eff_ke_crcl = pk.ke * pk.ke_multiplier
    hl_crcl = (np.log(2) / eff_ke_crcl) if eff_ke_crcl > 0 else 24
    hl_kgfr = 0
    if len(cr_data) >= 2:
        _, latest_kgfr = cr_func(cr_data[-1][0])
        if latest_kgfr is not None:
            vd_safe = pk.vd if (pk.vd and pk.vd > 0) else 50.0
            eff_ke_kgfr = ((latest_kgfr * 0.06) / vd_safe) * pk.ke_multiplier
            hl_kgfr = (np.log(2) / eff_ke_kgfr) if eff_ke_kgfr > 0 else 24
    max_hl = max(hl_crcl, hl_kgfr)
    auto_duration = int(np.ceil((5 * max_hl) / 24.0))
    auto_duration = max(7, min(auto_duration, 30))

    duration_days = st.slider(
        "Simulation Duration (Days)", 1, 30, auto_duration,
        help="Defaults to capturing at least 5 half-lives to show steady state (max 30 days).",
    )
    if auto_duration > 7 and duration_days == auto_duration:
        st.info(f"⏳ Auto-extended to **{auto_duration} days** (5 × t½ of ~{max_hl:.1f}h).")

    sim_end = sim_start + timedelta(days=duration_days)

    # --- Final simulation runs (kGFR first, then CrCl — order matters for
    #     the downstream regimen suggestion, which reads pk.ke) ------------
    results_kgfr = None
    if len(cr_data) >= 2:
        results_kgfr = pk.run(doses=doses, duration_days=duration_days, sim_start=sim_start,
                              cr_func=cr_func, patient_info=p_info, mode="kgfr")
    results = pk.run(doses=doses, duration_days=duration_days, sim_start=sim_start,
                     cr_func=cr_func, patient_info=p_info, mode="crcl")

    # --- Try / suggested regimen ------------------------------------------
    with st.container(border=True):
        st.header("Try Regimen / Suggested Regimen")
        show_try = st.checkbox("Show try/suggested regimen on graph", value=False)

        use_kgfr = False
        if results_kgfr is not None:
            use_kgfr = st.checkbox("Use estimated PK parameters from kGFR", value=False)

        if use_kgfr:
            base_ke_kgfr = results_kgfr["ke"] / max(pk.ke_multiplier, 0.01)
            pk_sugg = VancoPK(base_ke_kgfr, results_kgfr["vd"])
            pk_sugg.ke_multiplier = pk.ke_multiplier
            sugg_mode = "kgfr"
        else:
            pk_sugg = pk
            sugg_mode = "crcl"

        sugg_dose, sugg_interval, _ = suggest_regimen(pk_sugg, target_auc=AUC_TARGET,
                                                      patient_info=p_info)
        sugg_sim = pk_sugg.simulate_regimen(sugg_dose, sugg_interval, sim_start, sim_end,
                                            cr_func, p_info, mode=sugg_mode)
        st.markdown(f"**Suggested: {sugg_dose} mg q{sugg_interval}h** "
                    f"(Simulated AUC24 ≈ {sugg_sim['auc24']:.0f})")

        c1, c2 = st.columns(2)
        try_dose = c1.selectbox("Try dose (mg)", DOSE_OPTIONS,
                                index=DOSE_OPTIONS.index(sugg_dose))
        try_interval = c2.selectbox("Try interval (h)", INTERVAL_OPTIONS,
                                    index=INTERVAL_OPTIONS.index(sugg_interval))
        try_results = (pk_sugg.simulate_regimen(try_dose, try_interval, sim_start, sim_end,
                                                cr_func, p_info, mode=sugg_mode)
                       if show_try else None)

    # --- Confidence interval (from fitted multiplier SD) ------------------
    ci_bounds = None
    if levels:
        mult_lo, mult_hi = pk.compute_ci(level=0.5)
        fitted_mult = pk.ke_multiplier
        pk.ke_multiplier = mult_hi
        res_hi = pk.run(doses, duration_days=duration_days, sim_start=sim_start,
                        cr_func=cr_func, patient_info=p_info, mode="crcl")
        pk.ke_multiplier = mult_lo
        res_lo = pk.run(doses, duration_days=duration_days, sim_start=sim_start,
                        cr_func=cr_func, patient_info=p_info, mode="crcl")
        pk.ke_multiplier = fitted_mult
        ci_bounds = (res_lo, res_hi)

    # --- Static Cockcroft-Gault CrCl for plot overlay ---------------------
    static_crcl = None
    if cr_data:
        static_crcl = pk_params_from_patient(
            pt["age"], pt["sex"], pt["weight"], pt["height"],
            cr_func, when=cr_data[-1][0], muscle_factor=muscle_factor)["crcl"]

    # --- Plot -------------------------------------------------------------
    fig = plot_vanco_simulation(sim_start, results, cr_func, levels, level_times,
                                try_results, ci_bounds, static_crcl=static_crcl,
                                results_kgfr=results_kgfr)
    st.plotly_chart(fig, width="stretch")

    (st.info if is_fitted else st.warning)(fit_msg)

    # --- Metrics ----------------------------------------------------------
    def show_metrics(label, res, dose=None, interval=None):
        st.subheader(f"{label} ({dose:.0f} mg q{interval:.0f}h)" if dose and interval else label)
        cols = st.columns(6)
        cols[0].metric("ke (1/h)", f"{res['ke']:.3f}")
        cols[1].metric("Half-life (h)", f"{res['half_life']:.1f}")
        cols[2].metric("Vd (L)", f"{res['vd']:.1f}")
        cols[3].metric("AUC24", f"{res['auc24']:.0f}")
        if dose and interval:
            cpk, ctr = calculate_ss_conc(res["ke"], res["vd"], dose, interval)
            cols[4].metric("Cpkss (mg/L)", f"{cpk:.1f}")
            cols[5].metric("Ctrss (mg/L)", f"{ctr:.1f}")
        else:
            cols[4].metric("Cpkss", "N/A")
            cols[5].metric("Ctrss", "N/A")

        auc = res["auc24"]
        if AUC_LOW <= auc <= AUC_HIGH:
            st.success(f"AUC24 of {auc:.0f} is within target range ({AUC_LOW}-{AUC_HIGH}).")
        elif auc < AUC_LOW:
            st.error(f"AUC24 of {auc:.0f} is below target range (< {AUC_LOW}).")
        else:
            st.error(f"AUC24 of {auc:.0f} is above target range (> {AUC_HIGH}).")

    od_dose = od["dose"] if od["show"] else None
    od_interval = od["interval"] if od["show"] else None

    with st.container(border=True):
        show_metrics("Summary: Ordered Regimen", results, dose=od_dose, interval=od_interval)
    if try_results:
        with st.container(border=True):
            show_metrics("Summary: Try Regimen", try_results, dose=try_dose, interval=try_interval)
    if results_kgfr is not None:
        with st.container(border=True):
            show_metrics("Summary: Kinetic GFR", results_kgfr, dose=od_dose, interval=od_interval)

    st.markdown("<br>", unsafe_allow_html=True)
    if st.button("⬅️ Back to Patient Data & Dosing", width="stretch", key="back_bottom"):
        go_to("input")


# ===========================================================================
# ROUTER
# ===========================================================================
_normalize_entries()
if ss.view == "results":
    render_results_view()
else:
    render_input_view()
