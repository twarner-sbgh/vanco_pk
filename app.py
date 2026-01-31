import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from vanco_pk import VancoPK, pk_params_from_patient
from vanco_pk import VancoPK
from creatinine import build_creatinine_function
from dosing import (
    build_manual_doses,
    build_ordered_doses,
    build_try_regimen,
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
sex = st.selectbox("Sex", ["Male", "Female"])

weight = st.slider(
    "Weight (kg)",
    min_value=30.0,
    max_value=200.0,
    value=75.0,
    step=0.1,
)

height = st.slider(
    "Height (cm)",
    min_value=140.0,
    max_value=230.0,
    value=175.0,
    step=0.1,
)

# ---------------------------
# Creatinine
# ---------------------------
st.header("Creatinine")

cr1 = st.slider("Creatinine (µmol/L)", 35, 500, 210)
cr1_time = st.datetime_input("Creatinine time", datetime.now() - timedelta(days=1))

use_second_cr = st.checkbox("Enter second creatinine to simulate trend (to stable value after ~4 days)")

cr2 = None
cr2_time = None

if use_second_cr:
    cr2 = st.slider("Second creatinine (µmol/L)", 35, 500, 210)
    cr2_time = st.datetime_input("Second creatinine time", datetime.now())

use_cr_traj = st.checkbox("Use crude predicted creatinine change (multiple of baseline, reached after ~4 days)")

multiplier = 1.0

if use_cr_traj and not use_second_cr:
    multiplier = st.selectbox(
        "Predicted change",
        [3, 2, 1.5, 1, 0.75, 0.5, 0.25],
        index=3,
        format_func=lambda x: f"{x}x"
    )

cr_func = build_creatinine_function(
    cr1, cr1_time,
    cr2, cr2_time,
    use_cr_traj,
    multiplier
)

# ---------------------------
# Simulation start
# ---------------------------
sim_start_date = st.date_input(
    "Simulation start date",
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
# Build dose list
# ---------------------------
doses = []

# ---------------------------
# Manual doses
# ---------------------------
st.header("Manual doses")

manual_doses = []
manual_times = []

doses += build_manual_doses(
    manual_doses,
    manual_times,
    sim_start
)
for i in range(5):
    if st.checkbox(f"Enable manual dose {i+1}", key=f"md_on_{i}"):
        dose = st.selectbox(
            "Dose (mg)",
            [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500],
            index=3,
            key=f"md_dose_{i}"
        )
        time = st.datetime_input(
            "Time",
            sim_start + timedelta(hours=9, minutes=30),
            key=f"md_time_{i}"
        )

        manual_doses.append(dose)
        manual_times.append(time)

# ---------------------------
# Ordered dose
# ---------------------------
st.header("Ordered dose")

show_ordered_dose = st.checkbox(
    "Display ordered dose regimen",
    value=False,
    key="show_ordered_dose"
)

ordered_dose = None
ordered_interval = None
ordered_start = None
ordered_results = None

if show_ordered_dose:
    ordered_dose = st.selectbox(
        "Ordered dose (mg)",
        [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500],
        index=3,
        key="ordered_dose"
    )

    ordered_interval = st.selectbox(
        "Interval (h)",
        [6, 8, 12, 18, 24, 36, 48, 72],
        index=2,
        key="ordered_interval"
    )

    ordered_start = st.datetime_input(
        "Ordered dose start time",
        sim_start + timedelta(hours=9, minutes=30),
        key="ordered_start"
    )

if (
    show_ordered_dose
    and ordered_start is not None
    and ordered_start < sim_start
):
    st.error("Ordered dose start must be on or after simulation start.")
    st.stop()

if (
    show_ordered_dose
    and ordered_dose is not None
    and ordered_interval is not None
    and ordered_start is not None
):
    doses += build_ordered_doses(
        ordered_dose,
        ordered_interval,
        ordered_start,
        sim_start,
        sim_end
    )

# ---------------------------
# Levels
# ---------------------------
st.header("Measured levels")

levels = []
level_times = []

for i in range(5):
    if st.checkbox(f"Enable level {i+1}", key=f"lvl_on_{i}"):
        lvl = st.number_input(
            "Level (mg/L)",
            0.0, 100.0, 15.0,
            step=0.1, format="%.1f",
            key=f"lvl_val_{i}"
        )
        t = st.datetime_input(
            "Time",
            sim_start + timedelta(hours=9),
            key=f"lvl_time_{i}"
        )
        levels.append(lvl)
        level_times.append(t)

# ---------------------------
# Try regimen
# ---------------------------
st.header("Try regimen")

show_try_regimen = st.checkbox(
    "Show try regimen on graph",
    value=False,
    key="show_try_regimen"
)

# Placeholder (filled after PK sim runs)
try_dose = None
try_interval = None


# ---------------------------
# Run simulation
# ---------------------------
ke, vd = pk_params_from_patient(
    age=age,
    sex=sex,
    weight_kg=weight,
    height_cm=height,
    cr_func=cr_func,
    when=sim_start,
)

pk = VancoPK(ke, vd)

# --- instantiate PK engine ---
pk = VancoPK(ke, vd)

doses = build_manual_doses(manual_doses, manual_times, sim_start)
if show_ordered_dose:
    doses += build_ordered_doses(
        ordered_dose,
        ordered_interval,
        ordered_start,
        sim_start,
        sim_end
    )

results = pk.run(doses)

if show_ordered_dose:
    ordered_results = results  # ordered doses are part of the main simulation

# ---------------------------
# Suggested regimen
# ---------------------------
suggested_dose, suggested_interval, predicted_auc = suggest_regimen(
    pk, target_auc=500
)

st.subheader("Suggested regimen")
st.markdown(
    f"**{suggested_dose} mg IV q{suggested_interval}h**  "
    f"(Predicted AUC24 ≈ {predicted_auc:.0f})"
)

# Prepopulate try regimen
try_dose = st.selectbox(
    "Try dose (mg)",
    [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500],
    index=[250,500,750,1000,1250,1500,1750,2000,2500].index(suggested_dose)
)

try_interval = st.selectbox(
    "Try interval (h)",
    [6, 8, 12, 18, 24, 36, 48, 72],
    index=[6,8,12,18,24,36,48,72].index(suggested_interval)
)

# ===========================
# Plot
# ===========================
fig, ax = plt.subplots(figsize=(14, 6))

# ---- Entered regimen (manual + ordered doses) ----
t_hours = results["time"]
t_dates = [sim_start + timedelta(hours=h) for h in t_hours]
conc = results["conc"]

ax.plot(
    t_dates,
    conc,
    color="blue",
    linewidth=2,
    label="Entered regimen"
)

# ---- Try regimen (optional overlay) ----
try_results = None

if show_try_regimen:
    try_results = pk.simulate_regimen(
        try_dose,
        try_interval,
        sim_start,
        sim_end
    )

if try_results is not None:
    t_try = [sim_start + timedelta(hours=h) for h in try_results["time"]]
    ax.plot(
        t_try,
        try_results["conc"],
        color="green",
        linestyle="--",
        alpha=0.5,
        linewidth=2,
        label="Try regimen"
    )

# ---- Measured levels ----
for lvl, t in zip(levels, level_times):
    ax.scatter(
        t,
        lvl,
        color="#ff5733",
        s=60,
        zorder=5,
        label=None
    )
    ax.text(
        t,
        lvl + 0.8,
        f"{lvl:.1f} mg/L @ {t.strftime('%H:%M')}",
        fontsize=9,
        color="#ff5733",
        ha="center"
    )

# ---- Dose markers & labels ----
ax.autoscale()
y_top = ax.get_ylim()[1] * 0.95

for t_h, dose_mg in doses:
    t_dose = sim_start + timedelta(hours=t_h)
    ax.axvline(
        t_dose,
        color="lightgrey",
        linestyle=":",
        alpha=0.6,
        zorder=0
    )
    ax.text(
        t_dose,
        y_top,
        f"{dose_mg} mg",
        rotation=90,
        fontsize=8,
        color="lightgrey",
        va="top",
        ha="center"
    )

# ---- Axes formatting ----
ax.set_xlabel("Date")
ax.set_ylabel("Vancomycin concentration (mg/L)")
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax.legend()
ax.grid(alpha=0.2)

st.pyplot(fig)

# ---------------------------
# Stats display
# ---------------------------

ke = results["ke"]
t_half = results["half_life"]
vd = results["vd"]
auc24 = results["auc24"]

st.subheader("PK Summary")

col1, col2, col3, col4 = st.columns(4)

col1.metric("ke (1/h)", f"{ke:.3f}")
col2.metric("Half-life (h)", f"{t_half:.1f}")
col3.metric("Vd (L)", f"{vd:.1f}")
col4.metric("AUC24", f"{auc24:.0f}")

if auc24 < 400:
    st.error("AUC24 below target (400–600)")
elif auc24 > 600:
    st.error("AUC24 above target (400–600)")
else:
    st.success("AUC24 within target (400–600)")

# ---------------------------
# Confidence interval
# ---------------------------

if len(levels) >= 1:
    lo, hi = pk.compute_ci(level=0.5)
    st.info(f"50% CI for ke: {lo:.3f} – {hi:.3f}")
else:
    st.warning("Confidence intervals unreliable without measured levels.")
