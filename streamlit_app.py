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

import copy
import math
import uuid
from datetime import datetime, timedelta
from zoneinfo import ZoneInfo

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
TRY_START_LEVEL = 15.0            # mg/L; "start when level falls to" option
MAX_SIM_DAYS = 90                 # hard cap for any automatic extension
TRY_START_SIM = "At simulation start (steady-state view)"
TRY_START_AT_LEVEL = f"When level falls to {TRY_START_LEVEL:g} mg/L (from now)"

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
SCHEMA_VERSION = 3   # v3: no default individual dose; key-seeded widgets
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
ss.setdefault("dose_entries", [])           # individual doses are opt-in
ss.setdefault("level_entries", [])          # measured levels are optional
ss.setdefault("ordered", {
    "show": False, "dose": 1000, "interval": 12, "start": None,
})


def go_to(view):
    ss.view = view
    st.rerun()


def local_now():
    """Current wall-clock time in the viewer's time zone, as a naive datetime
    (all app datetimes are naive local times). Uses the time zone the browser
    reports, so a server running on UTC still gives the viewer's local time;
    falls back to the server clock if no time zone is available."""
    try:
        tz = st.context.timezone
        if tz:
            return datetime.now(ZoneInfo(tz)).replace(tzinfo=None)
    except Exception:
        pass
    return datetime.now()


def time_level_falls_to(pk_model, doses, now_h, sim_start, cr_func, p_info, mode, level):
    """Hours after sim_start when the predicted level falls to `level`,
    assuming every dose after now is held.

    Returns (hours, already_below, half_life_h). hours is None if the level
    stays above `level` for the whole 14-day search window. half_life_h is the
    model's half-life at the end of the window (renal function held at the
    last PCr), i.e. the half-life the new regimen will have.
    """
    given = [d for d in doses if d[0] < now_h]          # doses already started
    # Whole days, so the 5-min grid lines up exactly with the main simulations.
    horizon_days = math.ceil((max(now_h, 0.0) + 14 * 24) / 24)
    res = copy.copy(pk_model).run(given, duration_days=horizon_days, sim_start=sim_start,
                                  cr_func=cr_func, patient_info=p_info, mode=mode)
    t_end = res["time"] + (res["time"][1] - res["time"][0])   # conc[i] = value at end of step i
    hl = res["half_life"]
    above = np.nonzero((t_end >= now_h) & (res["conc"] > level))[0]
    if len(above) == 0:
        return max(now_h, 0.0), True, hl
    last = above[-1]            # last moment above target (skips any infusion still running)
    if last == len(t_end) - 1:
        return None, False, hl
    return t_end[last + 1], False, hl


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
def _seed(key, value):
    """Give a keyed widget its starting value from our persisted data.

    Widgets are created with a key and NO value=/index= argument, so the widget
    owns its value in st.session_state[key]. (Passing value= derived from data
    the widget itself updates makes Streamlit rebuild the widget on the next
    run and silently drop every second change.) Streamlit forgets widget keys
    while the other view is shown, so we re-seed from persisted data whenever
    the key is missing.
    """
    if key not in ss:
        ss[key] = value


def render_input_view():
    st.caption("Step 1 of 2 — Patient Data & Dosing")
    if msg := ss.pop("_flash", None):     # one-time notice left by the previous run
        st.warning(msg)

    # --- Simulation start date -------------------------------------------
    if ss.pop("_sync_sim_start", False):          # set by the auto-rewind below
        ss["w_sim_start"] = ss.sim_start_date
    _seed("w_sim_start", ss.sim_start_date)
    sim_start_date = st.date_input("Simulation Start Date", key="w_sim_start")
    ss.sim_start_date = sim_start_date
    sim_start = datetime.combine(sim_start_date, datetime.min.time())

    # --- Patient ----------------------------------------------------------
    pt = ss.patient
    with st.container(border=True):
        st.header("Patient")
        _seed("w_age", pt["age"])
        pt["age"] = st.slider("Age (years)", 17, 100, key="w_age")
        _seed("w_sex", pt["sex"])
        pt["sex"] = st.radio("Sex", ["Male", "Female"], horizontal=True, key="w_sex")
        _seed("w_weight", float(pt["weight"]))
        pt["weight"] = st.slider("Weight (kg)", 30.0, 200.0, step=0.5, key="w_weight")
        _seed("w_height", float(pt["height"]))
        pt["height"] = st.slider("Height (cm)", 140.0, 230.0, step=0.5, key="w_height")

    # --- Plasma creatinine ------------------------------------------------
    with st.container(border=True):
        st.header("Plasma Creatinine")
        _seed("w_muscle", pt["muscle_choice"])
        pt["muscle_choice"] = st.selectbox("Presumed Muscle Mass",
                                           options=list(MUSCLE_FACTORS), key="w_muscle")

        st.subheader("Measured PCr")
        for i, e in enumerate(ss.cr_entries):
            k_val, k_d, k_t = f"cr_val_{e['id']}", f"cr_d_{e['id']}", f"cr_t_{e['id']}"
            _seed(k_val, int(e["val"]))
            _seed(k_d, e["time"].date())
            _seed(k_t, e["time"].time())
            st.markdown(f"**PCr {i + 1} (µmol/L)**")
            e["val"] = st.slider(f"PCr {i + 1} (µmol/L)", 35, 500, step=1,
                                 key=k_val, label_visibility="collapsed")
            c1, c2, c3 = st.columns([2, 2, 0.5])
            d = c1.date_input(f"Date {i + 1}", key=k_d)
            t = c2.time_input(f"Time {i + 1}", key=k_t)
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
        if not ss.dose_entries:
            st.caption("No individual doses entered.")
        for i, e in enumerate(ss.dose_entries):
            if e["time"] is None:
                e["time"] = sim_start + DEFAULT_DOSE_HOUR
            k_v, k_d, k_t = f"md_v_{e['id']}", f"md_d_{e['id']}", f"md_t_{e['id']}"
            _seed(k_v, e["dose"])
            _seed(k_d, e["time"].date())
            _seed(k_t, e["time"].time())
            c1, c2, c3, c4 = st.columns([1.2, 1.2, 1.2, 0.4])
            e["dose"] = c1.selectbox("Dose (mg)", DOSE_OPTIONS, key=k_v)
            d = c2.date_input("Date", key=k_d)
            t = c3.time_input("Time", key=k_t)
            e["time"] = datetime.combine(d, t)
            c4.markdown("<div style='height:28px'></div>", unsafe_allow_html=True)
            if c4.button("✖", key=f"md_del_{e['id']}", help="Remove this dose"):
                ss.dose_entries.pop(i)
                st.rerun()

        if st.button("✚ Add Dose", key="md_add"):
            if ss.dose_entries:
                new_time = ss.dose_entries[-1]["time"] + timedelta(hours=12)
            else:
                new_time = sim_start + DEFAULT_DOSE_HOUR
            ss.dose_entries.append({"id": str(uuid.uuid4()), "dose": 1000, "time": new_time})
            st.rerun()

    # --- Ordered regimen --------------------------------------------------
    with st.container(border=True):
        st.header("Ordered Vancomycin Regimen")
        od = ss.ordered
        _seed("w_od_show", od["show"])
        od["show"] = st.checkbox("Display ordered regimen", key="w_od_show")
        if od["show"]:
            if od["start"] is None:
                od["start"] = sim_start + timedelta(days=1) + DEFAULT_DOSE_HOUR
            _seed("w_od_dose", od["dose"])
            _seed("w_od_int", od["interval"])
            _seed("w_od_d", od["start"].date())
            _seed("w_od_t", od["start"].time())
            od["dose"] = st.selectbox("Ordered dose (mg)", DOSE_OPTIONS, key="w_od_dose")
            od["interval"] = st.selectbox("Interval (h)", INTERVAL_OPTIONS, key="w_od_int")
            c1, c2 = st.columns(2)
            d = c1.date_input("Ordered Start Date", key="w_od_d")
            t = c2.time_input("Ordered Start Time", key="w_od_t")
            od["start"] = datetime.combine(d, t)

    # --- Measured levels --------------------------------------------------
    with st.container(border=True):
        st.header("Measured Vancomycin Levels")
        for i, e in enumerate(ss.level_entries):
            if e["time"] is None:
                e["time"] = sim_start + DEFAULT_LEVEL_HOUR
            k_v, k_d, k_t = f"lvl_v_{e['id']}", f"lvl_d_{e['id']}", f"lvl_t_{e['id']}"
            _seed(k_v, float(e["lvl"]))
            _seed(k_d, e["time"].date())
            _seed(k_t, e["time"].time())
            e["lvl"] = st.number_input(f"Level {i + 1} (mg/L)", 0.0, 100.0, step=0.1, key=k_v)
            c1, c2, c3 = st.columns([2, 2, 0.5])
            d = c1.date_input(f"Level date {i + 1}", key=k_d)
            t = c2.time_input(f"Level time {i + 1}", key=k_t)
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
            ss["_sync_sim_start"] = True      # push the new date into the widget next run
            # Saved and shown on the next run; a warning drawn here would be
            # wiped out immediately by the rerun.
            ss["_flash"] = (f"⚠️ **Simulation Start Date auto-adjusted** to "
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
        # Wide window so the regimen continues through any automatic extension
        # (doses past the end of a simulation have no effect on it).
        max_sim_end = sim_start + timedelta(days=MAX_SIM_DAYS)
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

    slider_days = duration_days
    sim_end = sim_start + timedelta(days=duration_days)

    # --- Final simulation runs (kGFR first, then CrCl — order matters for
    #     the downstream regimen suggestion, which reads pk.ke) ------------
    # AUC24 is averaged over whole ordered-regimen intervals when one is shown
    # (matters for q18/q36/q48/q72h); otherwise it is the final 24 h.
    auc_interval = od["interval"] if od["show"] else None

    def run_main(days):
        res_kgfr = None
        if len(cr_data) >= 2:
            res_kgfr = pk.run(doses=doses, duration_days=days, sim_start=sim_start,
                              cr_func=cr_func, patient_info=p_info, mode="kgfr",
                              auc_interval_h=auc_interval)
        res = pk.run(doses=doses, duration_days=days, sim_start=sim_start,
                     cr_func=cr_func, patient_info=p_info, mode="crcl",
                     auc_interval_h=auc_interval)
        return res_kgfr, res

    results_kgfr, results = run_main(duration_days)

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
        sugg_line = st.empty()        # filled in once the final duration is known

        c1, c2 = st.columns(2)
        try_dose = c1.selectbox("Try dose (mg)", DOSE_OPTIONS,
                                index=DOSE_OPTIONS.index(sugg_dose))
        try_interval = c2.selectbox("Try interval (h)", INTERVAL_OPTIONS,
                                    index=INTERVAL_OPTIONS.index(sugg_interval))

        start_mode = st.radio(
            "Try regimen starts", [TRY_START_SIM, TRY_START_AT_LEVEL], horizontal=True,
            key="try_start_mode",
            help="Simulation start shows the steady state the regimen reaches. "
                 f"The {TRY_START_LEVEL:g} mg/L option uses the current date and time: it holds "
                 "all doses after now and starts the new regimen once the predicted level "
                 f"has fallen to {TRY_START_LEVEL:g} mg/L.")

        try_start_dt, now_h, try_hl = None, None, None
        if start_mode == TRY_START_AT_LEVEL:
            now = local_now()
            now_h = (now - sim_start).total_seconds() / 3600
            t_hit, already, try_hl = time_level_falls_to(pk_sugg, doses, now_h, sim_start,
                                                         cr_func, p_info, sugg_mode,
                                                         TRY_START_LEVEL)
            if t_hit is None:
                st.warning(f"The predicted level does not fall to {TRY_START_LEVEL:g} mg/L "
                           "within 14 days, so no start time can be given.")
            else:
                start_h = math.ceil(t_hit * 4 - 1e-9) / 4      # round up to next 15 min
                try_start_dt = sim_start + timedelta(hours=start_h)
                reason = ("level is already at or below" if already
                          else "predicted level falls to")
                st.info(f"🕒 Start **{try_dose} mg q{try_interval}h** at "
                        f"**{try_start_dt:%b %d, %H:%M}** ({reason} {TRY_START_LEVEL:g} mg/L). "
                        f"Assumes no further doses are given after now "
                        f"({now:%b %d, %H:%M}).")

                # Extend the whole simulation, if needed, so the new regimen runs
                # for 5 half-lives after it starts (i.e. reaches steady state).
                if show_try:
                    needed_days = math.ceil((start_h + 5 * try_hl) / 24 - 1e-9)
                    if needed_days > duration_days:
                        duration_days = min(needed_days, MAX_SIM_DAYS)
                        sim_end = sim_start + timedelta(days=duration_days)
                        results_kgfr, results = run_main(duration_days)
                        if needed_days <= MAX_SIM_DAYS:
                            st.info(f"⏳ Simulation extended from {slider_days} to "
                                    f"**{duration_days} days** so the new regimen runs for "
                                    f"5 half-lives (t½ ≈ {try_hl:.1f} h) after it starts on "
                                    f"{try_start_dt:%b %d, %H:%M}.")
                        else:
                            st.warning(f"⏳ Simulation extended from {slider_days} to "
                                       f"**{MAX_SIM_DAYS} days**, the maximum. Running the new "
                                       f"regimen for 5 half-lives (t½ ≈ {try_hl:.1f} h) after it "
                                       f"starts on {try_start_dt:%b %d, %H:%M} would need "
                                       f"{needed_days} days, so its AUC24 and levels may not be "
                                       f"at steady state.")

        sugg_sim = pk_sugg.simulate_regimen(sugg_dose, sugg_interval, sim_start, sim_end,
                                            cr_func, p_info, mode=sugg_mode)
        sugg_line.markdown(f"**Suggested: {sugg_dose} mg q{sugg_interval}h** "
                           f"(Simulated AUC24 ≈ {sugg_sim['auc24']:.0f})")

        try_results = None
        if show_try:
            if start_mode == TRY_START_SIM:
                try_results = pk_sugg.simulate_regimen(try_dose, try_interval, sim_start,
                                                       sim_end, cr_func, p_info, mode=sugg_mode)
            elif try_start_dt is not None and try_start_dt < sim_end:
                # Actual doses up to now, a hold, then the new regimen from the start time.
                try_doses = [d for d in doses if d[0] < now_h]
                t = (try_start_dt - sim_start).total_seconds() / 3600
                while t < duration_days * 24:
                    try_doses.append((t, try_dose))
                    t += try_interval
                try_results = pk_sugg.run(try_doses, duration_days=duration_days,
                                          sim_start=sim_start, cr_func=cr_func,
                                          patient_info=p_info, mode=sugg_mode,
                                          auc_interval_h=try_interval)

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
        cols[3].metric("AUC24", f"{res['auc24']:.0f}",
                       help="Average 24-h AUC at the end of the simulation, taken over whole "
                            "dosing intervals when the regimen's interval is known.")
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
            try_label = ("Summary: Try Regimen" if try_start_dt is None
                         else f"Summary: Try Regimen from {try_start_dt:%b %d, %H:%M}")
            show_metrics(try_label, try_results, dose=try_dose, interval=try_interval)
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
