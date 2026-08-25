import numpy as np  # Numerical arrays and calculations
import matplotlib.pyplot as plt  # Plotting
import os
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parent / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

qa_pdf_dir_value = os.environ.get("FIGURE_QA_PDF_DIR")
QA_PDF_DIR = Path(qa_pdf_dir_value).resolve() if qa_pdf_dir_value else None

if QA_PDF_DIR is not None:
    QA_PDF_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 9.0,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)


def save_current_figure(filename):
    plt.savefig(OUTPUT_DIR / filename, dpi=300, bbox_inches="tight")

    if QA_PDF_DIR is not None:
        qa_stem = QA_PDF_DIR / Path(filename).stem

        plt.savefig(qa_stem.with_suffix(".svg"), bbox_inches="tight")

        plt.savefig(
            qa_stem.with_suffix(".pdf"),
            bbox_inches="tight",
        )

        plt.savefig(
            qa_stem.with_suffix(".tiff"),
            dpi=600,
            bbox_inches="tight",
        )

# ============================================================
# TIME

# ============================================================
day = 86400.0  # Seconds in one day
steps_per_day = 20  # Model timesteps per day
dt = day / steps_per_day  # Timestep length [s]
nyears = 100  # Total simulated years
spinup_years = 10  # Discarded spin-up years
nt = int(nyears * 365 * steps_per_day)  # Total number of timesteps
t_days = np.arange(nt) * dt / day  # Model time [days]
day_of_year = np.mod(t_days, 365.0)  # Day within each model year

# ============================================================
# PHYSICAL CONSTANTS
# ============================================================

rho_a = 1.25  # Air density [kg m^-3]
rho_l = 1000.0  # Liquid-water density [kg m^-3]
rho_s = 1000.0  # Dry-soil density [kg m^-3]
cp = 1003.0  # Dry-air heat capacity [J kg^-1 K^-1]
cp_s = 1000.0  # Dry-soil heat capacity [J kg^-1 K^-1]
cp_l = 4182.0  # Liquid-water heat capacity [J kg^-1 K^-1]
Lv = 2.257e6  # Latent heat of vaporization [J kg^-1]
sigma = 5.67e-8  # Stefan-Boltzmann constant [W m^-2 K^-4]
eps_s = 1.  # Surface longwave emissivity
eps_A = 0.67  # ABL effective longwave emissivity/absorptivity
eps_FT = 0.67  # Free-tropospheric effective downward emissivity
p_s = 1000.0  # Surface pressure [hPa]
Rv = 461.0  # Water-vapor gas constant [J kg^-1 K^-1]

# ============================================================
# BOX PARAMETERS

# ============================================================
h_BL = 1000.0  # Boundary-layer depth [m]
C_A = rho_a * cp * h_BL  # Atmospheric heat capacity per unit area [J m^-2 K^-1]
M_A = rho_a * h_BL  # Atmospheric mass per unit area [kg m^-2]
# C_L is calculated each timestep from current soil moisture
soil_depth = 0.1  # Single soil-layer depth [m]
porosity = 0.45  # Maximum volumetric soil-water fraction
Wmax = rho_l * soil_depth * porosity  # Maximum soil-water storage [kg m^-2]

# ============================================================
# SURFACE RESISTANCES

# ============================================================
r_LO = 250.0  # Ocean latent/moisture-transfer resistance [s m^-1]
r_SO = 1000.0  # Ocean sensible-heat-transfer resistance [s m^-1]
r_LS = 250.0
#r_LS = 1000.0  # Land latent/moisture-transfer resistance [s m^-1]
r_SS = 1000.0  # Land sensible-heat-transfer resistance [s m^-1]

# ============================================================
# LAND-OCEAN MIXING

# ============================================================
tau_mix = 2.0 * day  # Land-ocean mixing/relaxation timescale [s]

# ============================================================
# FREE-TROPOSPHERIC ENTRAINMENT

# ============================================================
tau_ent_L = 2.0 * day  # Land BL entrainment timescale [s]
tau_ent_o = 2.0 * day  # Ocean BL entrainment timescale [s]
We_L = 5.0e-6  # Land BL entrainment rate [s^-1]
We_o = 5.0e-6  # Ocean BL entrainment rate [s^-1]

# ============================================================
# FIXED OCEAN SST

# ============================================================
To = 270  # 293.0  # Fixed ocean SST [K]
dTo_dt = 0.0  # Prescribed ocean SST tendency [K s^-1]

# ============================================================
# LINEARIZED SATURATION SPECIFIC HUMIDITY

# # ============================================================
# Tref = 293.0  # Reference temperature for humidity linearization [K]
# Tc_ref = Tref - 273.15  # Reference temperature [degC]
# es_ref = 6.11 * 10.0**(7.5 * Tc_ref / (237.5 + Tc_ref))  # Saturation vapor pressure at Tref [hPa]
# qs_ref = 0.622 * es_ref / (p_s - 0.37 * es_ref)  # Saturation specific humidity at Tref [kg kg^-1]
# gamma_q = qs_ref * Lv / (Rv * Tref**2)  # Linearized Clausius-Clapeyron slope [kg kg^-1 K^-1]
# qs_o = qs_ref + gamma_q * (To - Tref)  # Saturation specific humidity at fixed ocean SST
# ============================================================
# NONLINEAR SATURATION SPECIFIC HUMIDITY
# ============================================================
Tc_o = To - 273.15  # Ocean SST in Celsius
es_o = 6.11 * 10.0**(7.5 * Tc_o / (237.5 + Tc_o))  # Saturation vapor pressure over ocean [hPa]
qs_o = 0.622 * es_o / (p_s - 0.37 * es_o)  # Saturation specific humidity at ocean SST [kg kg^-1]

# ============================================================
# SEASONAL FORCING

# # ============================================================
# Fmean = 160.0  # Annual-mean surface radiative forcing [W m^-2]
# Famp = 80.0  # Seasonal forcing amplitude [W m^-2]
# TFTmean = 290.0  # Annual-mean free-tropospheric temperature [K]
# TFTamp = 4.0  # Seasonal FT temperature amplitude [K]
# qFT = 0.003  # Prescribed free-tropospheric specific humidity [kg kg^-1]
# phase = 2.0 * np.pi * t_days / 365.0  # Annual-cycle phase [rad]
# F = Fmean - Famp * np.cos(phase)  # Seasonal surface radiative forcing [W m^-2]
# theta_FT = TFTmean - TFTamp * np.cos(phase - np.pi / 6.0)  # Seasonal free-tropospheric temperature [K]
# ============================================================
# STOCHASTIC FORCING -- NO ANNUAL CYCLE
# ============================================================

Fmean = 240.0  # Mean surface radiative forcing [W m^-2]
Fnoise_std = 30.0  # Equilibrium standard deviation [W m^-2]
tau_F = 0.001 * day  # Noise decorrelation timescale

TFTmean = 290.0  # Constant free-tropospheric temperature [K]
qFT = 0.003

rng_forcing = np.random.default_rng(1)

rho_F = np.exp(-dt / tau_F)
noise_std = Fnoise_std * np.sqrt(1.0 - rho_F**2)

Fnoise = np.zeros(nt)

for i in range(nt - 1):
    Fnoise[i + 1] = rho_F * Fnoise[i] + rng_forcing.normal(0.0, noise_std)

F = Fmean + Fnoise
F = np.maximum(F, 0.0)

theta_FT = np.full(nt, TFTmean)

# ============================================================
# STOCHASTIC PRECIPITATION

# ============================================================
rain_frequency = 0.2  # Mean rain-event frequency [day^-1]
rain_shape = 8.0  # Gamma shape parameter for event depth
rain_scale = 1.0  # Gamma scale parameter for event depth [mm]
rain_probability = rain_frequency / steps_per_day  # Rain-event probability per timestep
rng = np.random.default_rng(42)  # Reproducible random-number generator

# ============================================================
# STATE VARIABLES

# ============================================================
Ts = np.zeros(nt)  # Land surface temperature [K]
m = np.zeros(nt)  # Normalized soil moisture [0-1]
theta_L = np.zeros(nt)  # Near-surface land-air temperature [K]
q_L = np.zeros(nt)  # Near-surface land-air specific humidity [kg kg^-1]
theta_o = np.zeros(nt)  # Near-surface ocean-air temperature [K]
q_o = np.zeros(nt)  # Near-surface ocean-air specific humidity [kg kg^-1]
rain = np.zeros(nt)  # Rainfall depth per timestep [mm]
Eland = np.zeros(nt)  # Stored land evaporation [kg m^-2 s^-1]
Eocean = np.zeros(nt)  # Stored ocean evaporation [kg m^-2 s^-1]

# ============================================================
# INITIAL CONDITIONS

# ============================================================
Ts[0] = 285.0  # Initial land surface temperature [K]
m[0] = 0.6  # Initial normalized soil moisture
theta_L[0] = 285.0  # Initial land-air temperature [K]
q_L[0] = 0.006  # Initial land-air specific humidity [kg kg^-1]
theta_o[0] = 290.0  # Initial ocean-air temperature [K]
q_o[0] = 0.009  # Initial ocean-air specific humidity [kg kg^-1]

# ============================================================
# INTEGRATION

# ============================================================
for i in range(nt - 1):  # Advance model through time
    # --------------------------------------------------------
    # Saturation specific humidity
    # q*(Ts) = q*_ref + gamma * (Ts - Tref)
    # --------------------------------------------------------
    # qs_L = qs_ref + gamma_q * (Ts[i] - Tref)  # Linearized land saturation specific humidity
    # qs_L = max(qs_L, 0.0)  # Prevent negative saturation humidity
    Tc_L = Ts[i] - 273.15  # Land surface temperature in Celsius
    es_L = 6.11 * 10.0**(7.5 * Tc_L / (237.5 + Tc_L))  # Saturation vapor pressure [hPa]
    qs_L =  0.622 * es_L / (p_s - 0.37 * es_L)  # Nonlinear saturation specific humidity [kg kg^-1]
    # --------------------------------------------------------
    # LAND EVAPORATION
    #
    # E_L = rho_a / r_LS * m * [q*(Ts) - q_L]
    # --------------------------------------------------------
    qdiff_L = max(qs_L - q_L[i], 0.0)  # Positive land surface-air humidity contrast
    E_L = rho_a * m[i] * qdiff_L / r_LS  # Soil-moisture-limited land evaporation
    # --------------------------------------------------------
    # OCEAN EVAPORATION
    #
    # E_o = rho_a / r_LO * [q*(To) - q_o]
    # --------------------------------------------------------
    qdiff_o = max(qs_o - q_o[i], 0.0)  # Positive ocean surface-air humidity contrast
    E_o = rho_a * qdiff_o / r_LO  # Ocean evaporation
    Eland[i] = E_L*1.  # Store land evaporation
    Eocean[i] = E_o*1.  # Store ocean evaporation
    # --------------------------------------------------------
    # SENSIBLE HEAT FLUX
    # --------------------------------------------------------
    H_L = rho_a * cp * (Ts[i] - theta_L[i]) / r_SS  # Land sensible heat flux [W m^-2]
    H_o = rho_a * cp * (To - theta_o[i]) / r_SO  # Ocean sensible heat flux [W m^-2]
    # --------------------------------------------------------
    # LAND LONGWAVE RADIATION
    # --------------------------------------------------------
    OLR_L = eps_s * sigma * Ts[i]**4  # Upward land longwave radiation [W m^-2]
    DLR_L = eps_s * eps_A * sigma * theta_L[i]**4  # Downward ABL longwave absorbed by land surface [W m^-2]
    R_ad = eps_FT * sigma * theta_FT[i]**4  # Downward longwave radiation from free troposphere [W m^-2]
    R_L = eps_A * (R_ad + OLR_L) - 2.0 * eps_A * sigma * theta_L[i]**4  # Net longwave heating of land ABL [W m^-2]
    OLR_o = eps_s * sigma * To**4  # Upward ocean longwave radiation [W m^-2]
    R_o = eps_A * (R_ad + OLR_o) - 2.0 * eps_A * sigma * theta_o[i]**4  # Net longwave heating of ocean ABL [W m^-2]
    # --------------------------------------------------------
    # STOCHASTIC PRECIPITATION
    # rain_mm is rainfall depth during this timestep
    # --------------------------------------------------------
    does_rain = rng.random()  # Random draw for rain occurrence
    if does_rain < rain_probability:  # Check whether a rain event occurs
        rain_mm = rng.gamma(rain_shape, rain_scale)  # Draw event rainfall depth [mm]
    else:  # No rain event
        rain_mm = 0.0  # Set rainfall depth to zero
    rain[i] = rain_mm  # Store timestep rainfall depth [mm]
    # 1 mm water = 1 kg m^-2
    # Convert event depth to timestep-mean flux
    P_rate = rain_mm / dt  # Convert rainfall depth to water mass flux [kg m^-2 s^-1]
    # ========================================================
    # SIX ODES
    # ========================================================
    # --------------------------------------------------------
    # 1. LAND SURFACE TEMPERATURE
    #
    # C_L dTs/dt =
    # SW - OLR + DLR - L_v E_L - H_L
    # --------------------------------------------------------
    C_L = soil_depth * (cp_s * rho_s + cp_l * rho_l * porosity * m[i])  # Soil-moisture-dependent land heat capacity [J m^-2 K^-1]
    dTs_dt = (F[i] - OLR_L + DLR_L - Lv * E_L - H_L) / C_L  # Land surface temperature tendency [K s^-1]
    # --------------------------------------------------------
    # 2. SOIL MOISTURE
    #
    # dm/dt = (P - E_L) / Wmax
    # --------------------------------------------------------
    dm_dt = (P_rate - E_L) / Wmax  # Normalized soil-moisture tendency [s^-1]
    # --------------------------------------------------------
    # 3. LAND ATMOSPHERIC TEMPERATURE
    #
    # mixing with ocean
    # + sensible heating from land
    # + entrainment from free troposphere
    # --------------------------------------------------------
    dtheta_L_dt = (theta_o[i] - theta_L[i]) / tau_mix  # Land-air temperature tendency from land-ocean mixing
    dtheta_L_dt += H_L / C_A  # Add sensible heating from land surface
    dtheta_L_dt += R_L / C_A  # Add longwave radiative heating/cooling of land ABL
    dtheta_L_dt += We_L * (theta_FT[i] - theta_L[i])  # Add FT entrainment to land-air temperature
    # --------------------------------------------------------
    # 4. LAND ATMOSPHERIC HUMIDITY
    #
    # mixing with ocean
    # + land evaporation
    # + FT entrainment
    # --------------------------------------------------------
    dq_L_dt = (q_o[i] - q_L[i]) / tau_mix  # Land-air humidity tendency from land-ocean mixing
    dq_L_dt += E_L / M_A  # Add moistening from land evaporation
    dq_L_dt += We_L * (qFT - q_L[i])  # Add FT entrainment to land-air humidity
    # --------------------------------------------------------
    # 5. OCEAN ATMOSPHERIC TEMPERATURE
    #
    # mixing with land
    # + sensible heating from ocean
    # + FT entrainment
    # --------------------------------------------------------
    dtheta_o_dt = (theta_L[i] - theta_o[i]) / tau_mix  # Ocean-air temperature tendency from land-ocean mixing
    dtheta_o_dt += H_o / C_A  # Add sensible heating from ocean surface
    dtheta_o_dt += R_o / C_A  # Add longwave radiative heating/cooling of ocean ABL
    dtheta_o_dt += We_o * (theta_FT[i] - theta_o[i])  # Add FT entrainment to ocean-air temperature
    # --------------------------------------------------------
    # 6. OCEAN ATMOSPHERIC HUMIDITY
    #
    # mixing with land
    # + ocean evaporation
    # + FT entrainment
    # --------------------------------------------------------
    dq_o_dt = (q_L[i] - q_o[i]) / tau_mix  # Ocean-air humidity tendency from land-ocean mixing
    dq_o_dt += E_o / M_A  # Add moistening from ocean evaporation
    dq_o_dt += We_o * (qFT - q_o[i])  # Add FT entrainment to ocean-air humidity
    # ========================================================
    # FORWARD EULER UPDATE
    # ========================================================
    Ts[i + 1] = Ts[i] + dTs_dt * dt  # Forward-Euler update for land surface temperature
    m[i + 1] = m[i] + dm_dt * dt  # Forward-Euler update for soil moisture
    theta_L[i + 1] = theta_L[i] + dtheta_L_dt * dt  # Forward-Euler update for land-air temperature
    q_L[i + 1] = q_L[i] + dq_L_dt * dt  # Forward-Euler update for land-air humidity
    theta_o[i + 1] = theta_o[i] + dtheta_o_dt * dt  # Forward-Euler update for ocean-air temperature
    q_o[i + 1] = q_o[i] + dq_o_dt * dt  # Forward-Euler update for ocean-air humidity
    # ========================================================
    # PHYSICAL BOUNDS
    # ========================================================
    m[i + 1] = np.clip(m[i + 1], 0.0, 1.0)  # Keep normalized soil moisture within [0,1]
    q_L[i + 1] = max(q_L[i + 1], 0.0)  # Prevent negative land-air humidity
    q_o[i + 1] = max(q_o[i + 1], 0.0)  # Prevent negative ocean-air humidity

# ============================================================
# SELECT THE FULL POST-SPIN-UP PERIOD

# ============================================================
after_spinup = t_days >= spinup_years * 365.0  # Select times after spin-up
Ts_analysis = Ts[after_spinup]  # Post-spin-up land surface temperature
m_analysis = m[after_spinup]  # Post-spin-up normalized soil moisture
theta_L_analysis = theta_L[after_spinup]  # Post-spin-up land-air temperature
theta_o_analysis = theta_o[after_spinup]  # Post-spin-up ocean-air temperature
q_L_analysis = q_L[after_spinup]  # Post-spin-up land-air specific humidity
q_o_analysis = q_o[after_spinup]  # Post-spin-up ocean-air specific humidity

# ============================================================
# POST-SPIN-UP STATISTICS

# ============================================================
Ts_anom = Ts_analysis - np.mean(Ts_analysis)  # Post-spin-up land-temperature anomaly
Ts_std = np.std(Ts_analysis)  # Post-spin-up land-temperature standard deviation
Ts_skew = np.mean(Ts_anom**3) / Ts_std**3  # Post-spin-up land-temperature skewness
print("Post-spin-up sample count =", Ts_analysis.size)
print("Post-spin-up mean Ts =", np.mean(Ts_analysis), "K")
print("Post-spin-up std Ts =", Ts_std, "K")
print("Post-spin-up skewness =", Ts_skew)

# ============================================================
# PLOT 1: POST-SPIN-UP TEMPERATURE PDF

# ============================================================
plt.figure(figsize=(7, 4))  # Create figure
plt.hist(Ts_analysis - 273.15, bins=40, density=True)  # Plot post-spin-up land-temperature PDF
plt.xlabel("Land surface temperature (°C)")  # Label x axis
plt.ylabel("Probability density")  # Label y axis
plt.title("Post-spin-up land temperature under stationary forcing")  # Add plot title
plt.tight_layout()  # Adjust figure layout
save_current_figure("land_ocean_postspinup_temperature_distribution.png")
if "agg" not in plt.get_backend().lower():
    plt.show()
plt.close()

# ============================================================
# PLOT 2: SOIL MOISTURE VS POST-SPIN-UP TEMPERATURE

# ============================================================
plt.figure(figsize=(7, 4))  # Create figure
plt.scatter(m_analysis, Ts_analysis - 273.15, s=2, alpha=0.15, rasterized=True)  # Plot soil moisture versus post-spin-up temperature
plt.xlabel("Soil moisture")  # Label x axis
plt.ylabel("Land surface temperature (°C)")  # Model calculation
plt.title("Post-spin-up soil moisture-temperature relationship")  # Add plot title
plt.tight_layout()  # Adjust figure layout
save_current_figure("land_ocean_postspinup_soil_moisture_temperature.png")
if "agg" not in plt.get_backend().lower():
    plt.show()
plt.close()

# ============================================================
# PLOT 3: LAST FIVE YEARS

# ============================================================
last_five_years = t_days >= (nyears - 5) * 365.0  # Select final five simulated years
plt.figure(figsize=(9, 4))  # Create figure
plt.plot(t_days[last_five_years] / 365.0, Ts[last_five_years] - 273.15, label="Land surface")  # Plot land surface temperature
plt.plot(t_days[last_five_years] / 365.0, theta_L[last_five_years] - 273.15, label="Land atmosphere")  # Plot land-air temperature
plt.plot(t_days[last_five_years] / 365.0, theta_o[last_five_years] - 273.15, label="Ocean atmosphere")  # Plot ocean-air temperature
plt.xlabel("Year")  # Label x axis
plt.ylabel("Temperature (°C)")  # Label y axis
plt.legend()  # Add legend
plt.tight_layout()  # Adjust figure layout
save_current_figure("land_ocean_last_five_year_temperature_timeseries.png")
if "agg" not in plt.get_backend().lower():
    plt.show()
plt.close()
# constant mixing; variable mixing between the ocean and land
