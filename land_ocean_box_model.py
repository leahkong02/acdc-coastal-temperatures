import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parent / "figures"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 9.0,
    }
)

# ============================================================
# TIME
# ============================================================

day = 86400.0

steps_per_day = 24

dt = day / steps_per_day

nyears = 100

spinup_years = 10

nt = int(nyears * 365 * steps_per_day)

t_days = np.arange(nt) * dt / day

day_of_year = np.mod(t_days, 365.0)

# ============================================================
# PHYSICAL CONSTANTS
# ============================================================

rho_a = 1.25

rho_l = 1000.0

rho_s = 1000.0

cp = 1003.0

cp_s = 1000.0

cp_l = 4182.0

Lv = 2.257e6

sigma = 5.67e-8

eps_s = 0.98

eps_A = 0.75

eps_FT = 0.70

p_s = 1000.0

Rv = 461.0

# ============================================================
# BOX PARAMETERS
# ============================================================

h_BL = 1000.0

C_A = rho_a * cp * h_BL

M_A = rho_a * h_BL

# C_L is calculated each timestep from current soil moisture.

soil_depth = 0.1

porosity = 0.45

Wmax = rho_l * soil_depth * porosity

# ============================================================
# SURFACE RESISTANCES
# ============================================================

r_LO = 250.0

r_SO = 1000.0

r_LS = 250.0

r_SS = 1000.0

# ============================================================
# LAND-OCEAN MIXING
# ============================================================

tau_mix = 2.0 * day

# ============================================================
# FREE-TROPOSPHERIC ENTRAINMENT
# ============================================================

tau_ent_L = 2.0 * day

tau_ent_o = 2.0 * day

We_L = 1.0 / tau_ent_L

We_o = 1.0 / tau_ent_o

# ============================================================
# FIXED OCEAN SST
# ============================================================

To = 290.0

dTo_dt = 0.0

# ============================================================
# NONLINEAR SATURATION SPECIFIC HUMIDITY
# ============================================================

Tc_o = To - 273.15

es_o = 6.11 * 10.0**(7.5 * Tc_o / (237.5 + Tc_o))

qs_o = 0.622 * es_o / (p_s - 0.37 * es_o)

# ============================================================
# STATIONARY STOCHASTIC RADIATIVE FORCING
# ============================================================

Fmean = 160.0

Fnoise_std = 30.0

tau_F = 5.0 * day

TFTmean = 290.0

qFT = 0.003

rng_forcing = np.random.default_rng(123)

rho_F = np.exp(-dt / tau_F)

innovation_std = Fnoise_std * np.sqrt(1.0 - rho_F**2)

Fnoise = np.zeros(nt)

for i in range(nt - 1):
    Fnoise[i + 1] = rho_F * Fnoise[i] + rng_forcing.normal(
        0.0,
        innovation_std,
    )

F = np.maximum(Fmean + Fnoise, 0.0)

theta_FT = np.full(nt, TFTmean)

# ============================================================
# STOCHASTIC PRECIPITATION
# ============================================================

rain_frequency = 0.2

rain_shape = 8.0

rain_scale = 1.0

rain_probability = rain_frequency / steps_per_day

rng = np.random.default_rng(42)

# ============================================================
# STATE VARIABLES
# ============================================================

Ts = np.zeros(nt)

m = np.zeros(nt)

theta_L = np.zeros(nt)

q_L = np.zeros(nt)

theta_o = np.zeros(nt)

q_o = np.zeros(nt)

rain = np.zeros(nt)

Eland = np.zeros(nt)

Eocean = np.zeros(nt)

# ============================================================
# INITIAL CONDITIONS
# ============================================================

Ts[0] = 285.0

m[0] = 0.6

theta_L[0] = 285.0

q_L[0] = 0.006

theta_o[0] = 290.0

q_o[0] = 0.009

# ============================================================
# INTEGRATION
# ============================================================

for i in range(nt - 1):

    # --------------------------------------------------------
    # Saturation specific humidity
    # Tetens saturation vapor pressure and specific humidity
    # --------------------------------------------------------

    Tc_L = Ts[i] - 273.15

    es_L = 6.11 * 10.0**(7.5 * Tc_L / (237.5 + Tc_L))

    qs_L = 0.622 * es_L / (p_s - 0.37 * es_L)

    # --------------------------------------------------------
    # LAND EVAPORATION
    #
    # E_L = rho_a / r_LS * m * [q*(Ts) - q_L]
    # --------------------------------------------------------

    qdiff_L = max(qs_L - q_L[i], 0.0)

    E_L = rho_a * m[i] * qdiff_L / r_LS

    # --------------------------------------------------------
    # OCEAN EVAPORATION
    #
    # E_o = rho_a / r_LO * [q*(To) - q_o]
    # --------------------------------------------------------

    qdiff_o = max(qs_o - q_o[i], 0.0)

    E_o = rho_a * qdiff_o / r_LO

    Eland[i] = E_L

    Eocean[i] = E_o

    # --------------------------------------------------------
    # SENSIBLE HEAT FLUX
    # --------------------------------------------------------

    H_L = rho_a * cp * (Ts[i] - theta_L[i]) / r_SS

    H_o = rho_a * cp * (To - theta_o[i]) / r_SO

    # --------------------------------------------------------
    # LAND LONGWAVE RADIATION
    # --------------------------------------------------------

    OLR_L = eps_s * sigma * Ts[i]**4

    DLR_L = eps_s * eps_A * sigma * theta_L[i]**4

    R_ad = eps_FT * sigma * theta_FT[i]**4

    R_L = eps_A * (R_ad + OLR_L) - 2.0 * eps_A * sigma * theta_L[i]**4

    OLR_o = eps_s * sigma * To**4

    R_o = eps_A * (R_ad + OLR_o) - 2.0 * eps_A * sigma * theta_o[i]**4

    # --------------------------------------------------------
    # STOCHASTIC PRECIPITATION
    # rain_mm is rainfall depth during this timestep
    # --------------------------------------------------------

    does_rain = rng.random()

    if does_rain < rain_probability:
        rain_mm = rng.gamma(rain_shape, rain_scale)
    else:
        rain_mm = 0.0

    rain[i] = rain_mm

    # 1 mm water = 1 kg m^-2
    # Convert event depth to timestep-mean flux
    P_rate = rain_mm / dt

    # ========================================================
    # SIX ODES
    # ========================================================

    # --------------------------------------------------------
    # 1. LAND SURFACE TEMPERATURE
    #
    # C_L dTs/dt =
    # SW - OLR + DLR - L_v E_L - H_L
    # --------------------------------------------------------

    C_L = soil_depth * (
        cp_s * rho_s + cp_l * rho_l * porosity * m[i]
    )

    dTs_dt = (F[i] - OLR_L + DLR_L - Lv * E_L - H_L) / C_L

    # --------------------------------------------------------
    # 2. SOIL MOISTURE
    #
    # dm/dt = (P - E_L) / Wmax
    # --------------------------------------------------------

    dm_dt = (P_rate - E_L) / Wmax

    # --------------------------------------------------------
    # 3. LAND ATMOSPHERIC TEMPERATURE
    #
    # mixing with ocean
    # + sensible heating from land
    # + entrainment from free troposphere
    # --------------------------------------------------------

    dtheta_L_dt = (theta_o[i] - theta_L[i]) / tau_mix

    dtheta_L_dt += H_L / C_A

    dtheta_L_dt += R_L / C_A

    dtheta_L_dt += We_L * (theta_FT[i] - theta_L[i])

    # --------------------------------------------------------
    # 4. LAND ATMOSPHERIC HUMIDITY
    #
    # mixing with ocean
    # + land evaporation
    # + FT entrainment
    # --------------------------------------------------------

    dq_L_dt = (q_o[i] - q_L[i]) / tau_mix

    dq_L_dt += E_L / M_A

    dq_L_dt += We_L * (qFT - q_L[i])

    # --------------------------------------------------------
    # 5. OCEAN ATMOSPHERIC TEMPERATURE
    #
    # mixing with land
    # + sensible heating from ocean
    # + FT entrainment
    # --------------------------------------------------------

    dtheta_o_dt = (theta_L[i] - theta_o[i]) / tau_mix

    dtheta_o_dt += H_o / C_A

    dtheta_o_dt += R_o / C_A

    dtheta_o_dt += We_o * (theta_FT[i] - theta_o[i])

    # --------------------------------------------------------
    # 6. OCEAN ATMOSPHERIC HUMIDITY
    #
    # mixing with land
    # + ocean evaporation
    # + FT entrainment
    # --------------------------------------------------------

    dq_o_dt = (q_L[i] - q_o[i]) / tau_mix

    dq_o_dt += E_o / M_A

    dq_o_dt += We_o * (qFT - q_o[i])

    # ========================================================
    # FORWARD EULER UPDATE
    # ========================================================

    Ts[i + 1] = Ts[i] + dTs_dt * dt

    m[i + 1] = m[i] + dm_dt * dt

    theta_L[i + 1] = theta_L[i] + dtheta_L_dt * dt

    q_L[i + 1] = q_L[i] + dq_L_dt * dt

    theta_o[i + 1] = theta_o[i] + dtheta_o_dt * dt

    q_o[i + 1] = q_o[i] + dq_o_dt * dt

    # ========================================================
    # PHYSICAL BOUNDS
    # ========================================================

    m[i + 1] = np.clip(m[i + 1], 0.0, 1.0)

    q_L[i + 1] = max(q_L[i + 1], 0.0)

    q_o[i + 1] = max(q_o[i + 1], 0.0)

# ============================================================
# SELECT A FIXED CALENDAR-JJA WINDOW AFTER SPIN-UP
# ============================================================

after_spinup = t_days >= spinup_years * 365.0

JJA = (day_of_year >= 151.0) & (day_of_year < 243.0)

summer = after_spinup & JJA

Ts_summer = Ts[summer]

m_summer = m[summer]

theta_L_summer = theta_L[summer]

theta_o_summer = theta_o[summer]

q_L_summer = q_L[summer]

q_o_summer = q_o[summer]

# ============================================================
# SUMMER STATISTICS
# ============================================================

Ts_anom = Ts_summer - np.mean(Ts_summer)

Ts_std = np.std(Ts_summer)

Ts_skew = np.mean(Ts_anom**3) / Ts_std**3

print("Calendar-JJA mean Ts =", np.mean(Ts_summer), "K")

print("Calendar-JJA std Ts =", Ts_std, "K")

print("Calendar-JJA skewness =", Ts_skew)

# ============================================================
# PLOT 1: SUMMER TEMPERATURE PDF
# ============================================================

plt.figure(figsize=(7, 4))

plt.hist(Ts_summer - 273.15, bins=50, density=True)

plt.xlabel("Land surface temperature (°C)")

plt.ylabel("Probability density")

plt.title("Calendar-JJA land temperature under stationary forcing")

plt.tight_layout()

plt.savefig(
    OUTPUT_DIR / "land_ocean_jja_temperature_distribution.png",
    dpi=300,
    bbox_inches="tight",
)

if "agg" not in plt.get_backend().lower():
    plt.show()

plt.close()

# ============================================================
# PLOT 2: SOIL MOISTURE VS SUMMER TEMPERATURE
# ============================================================

plt.figure(figsize=(7, 4))

plt.scatter(
    m_summer,
    Ts_summer - 273.15,
    s=2,
    alpha=0.15,
    rasterized=True,
)

plt.xlabel("Soil moisture")

plt.ylabel("Land surface temperature (°C)")

plt.title("Calendar-JJA soil moisture-temperature relationship")

plt.tight_layout()

plt.savefig(
    OUTPUT_DIR / "land_ocean_jja_soil_moisture_temperature.png",
    dpi=300,
    bbox_inches="tight",
)

if "agg" not in plt.get_backend().lower():
    plt.show()

plt.close()

# ============================================================
# PLOT 3: LAST FIVE YEARS
# ============================================================

last_five_years = t_days >= (nyears - 5) * 365.0

plt.figure(figsize=(9, 4))

plt.plot(t_days[last_five_years] / 365.0, Ts[last_five_years] - 273.15, label="Land surface")

plt.plot(t_days[last_five_years] / 365.0, theta_L[last_five_years] - 273.15, label="Land atmosphere")

plt.plot(t_days[last_five_years] / 365.0, theta_o[last_five_years] - 273.15, label="Ocean atmosphere")

plt.xlabel("Year")

plt.ylabel("Temperature (°C)")

plt.legend()

plt.tight_layout()

plt.savefig(
    OUTPUT_DIR / "land_ocean_last_five_year_temperature_timeseries.png",
    dpi=300,
    bbox_inches="tight",
)

if "agg" not in plt.get_backend().lower():
    plt.show()

plt.close()
