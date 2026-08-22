from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr


coasts = {
    "West (Oregon)": {
        "file_tag": "oregon",
        "pilot": "data/era5_wind_validation/raw/era5_oregon_202407_6hourly.nc",
        "latitude_range": (44.50, 45.50),
        "land_side": "east",
        "color": "#2878B5",
    },
    "East (Maine)": {
        "file_tag": "maine",
        "pilot": "data/era5_wind_validation/raw/era5_maine_sample_20240701T12.nc",
        "latitude_range": (44.00, 45.00),
        "land_side": "west",
        "color": "#E07B39",
    },
}

results = {}

for coast, settings in coasts.items():
    pilot = xr.open_dataset(settings["pilot"])
    lsm = pilot.lsm.isel(valid_time=0).values
    latitudes = pilot.latitude.values
    longitudes = pilot.longitude.values
    lat_min, lat_max = settings["latitude_range"]
    rows = np.where((latitudes >= lat_min) & (latitudes <= lat_max))[0]
    ocean_columns = []
    land_columns = []

    for row in rows:
        ocean_candidates = np.where(lsm[row] <= 0.30)[0]
        land_candidates = np.where(lsm[row] >= 0.70)[0]

        if settings["land_side"] == "east":
            ocean_column = ocean_candidates.max()
            land_column = land_candidates[land_candidates > ocean_column].min()
        else:
            ocean_column = ocean_candidates.min()
            land_column = land_candidates[land_candidates < ocean_column].max()

        ocean_columns.append(ocean_column)
        land_columns.append(land_column)

    ocean_columns = np.asarray(ocean_columns)
    land_columns = np.asarray(land_columns)
    coast_longitudes = 0.5 * (longitudes[ocean_columns] + longitudes[land_columns])
    slope, intercept = np.polyfit(latitudes[rows], coast_longitudes, 1)
    east_per_north = slope * np.cos(np.deg2rad(latitudes[rows].mean()))

    if settings["land_side"] == "east":
        normal = np.array([1.0, -east_per_north])
    else:
        normal = np.array([-1.0, east_per_north])

    normal = normal / np.linalg.norm(normal)
    pair_distance = 111320.0 * np.cos(np.deg2rad(latitudes[rows]))
    pair_distance = pair_distance * np.abs(longitudes[land_columns] - longitudes[ocean_columns])
    cross_shore_distance = np.mean(pair_distance * np.abs(normal[0]))
    daily_frames = []

    for year in range(2015, 2025):
        path = Path(f"data/era5_wind_validation/raw/era5_{settings['file_tag']}_{year}_maysep_6hourly.nc")
        ds = xr.open_dataset(path)
        u10 = ds.u10.values
        v10 = ds.v10.values
        t2m = ds.t2m.values
        u_coast = 0.5 * (u10[:, rows, land_columns] + u10[:, rows, ocean_columns]).mean(axis=1)
        v_coast = 0.5 * (v10[:, rows, land_columns] + v10[:, rows, ocean_columns]).mean(axis=1)
        t_land = t2m[:, rows, land_columns].mean(axis=1)
        t_ocean = t2m[:, rows, ocean_columns].mean(axis=1)
        frame = pd.DataFrame(index=pd.to_datetime(ds.valid_time.values))
        frame["u_normal"] = normal[0] * u_coast + normal[1] * v_coast
        frame["delta_t"] = t_ocean - t_land
        frame["temperature_advection"] = 86400.0 * frame["u_normal"] * frame["delta_t"] / cross_shore_distance
        frame["t_land"] = t_land
        daily = frame.resample("1D").mean()
        daily["tmax_6h"] = frame["t_land"].resample("1D").max()
        daily_frames.append(daily)

    daily = pd.concat(daily_frames)
    calendar_day = daily.index.strftime("%m-%d")
    climatology = daily["tmax_6h"].groupby(calendar_day).mean()
    climatology = climatology.rolling(15, center=True, min_periods=1).mean()
    daily["tmax_anomaly"] = daily["tmax_6h"].values - climatology.loc[calendar_day].values
    temperature = daily["tmax_anomaly"].to_numpy()
    wind = daily["u_normal"].to_numpy()
    advection = daily["temperature_advection"].to_numpy()

    assert len(daily) == 1530
    assert np.isfinite(temperature).all()
    assert np.isfinite(wind).all()
    assert np.isfinite(advection).all()

    temperature_skew = np.mean((temperature - temperature.mean()) ** 3) / np.std(temperature) ** 3
    wind_skew = np.mean((wind - wind.mean()) ** 3) / np.std(wind) ** 3
    advection_skew = np.mean((advection - advection.mean()) ** 3) / np.std(advection) ** 3
    onshore_fraction = np.mean(wind > 0.0)
    results[coast] = {
        "temperature": temperature,
        "wind": wind,
        "advection": advection,
        "temperature_skew": temperature_skew,
        "wind_skew": wind_skew,
        "advection_skew": advection_skew,
        "onshore_fraction": onshore_fraction,
        "cross_shore_distance": cross_shore_distance,
        "color": settings["color"],
    }

    print(
        f"{coast}: n={len(daily)}, dx={cross_shore_distance / 1000.0:.1f} km, "
        f"temperature skew={temperature_skew:+.3f}, wind skew={wind_skew:+.3f}, "
        f"advection skew={advection_skew:+.3f}, onshore={onshore_fraction:.3f}, "
        f"mean advection={advection.mean():+.3f} K/day"
    )

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.75), constrained_layout=True)
all_temperature = np.concatenate([result["temperature"] for result in results.values()])
all_wind = np.concatenate([result["wind"] for result in results.values()])
all_advection = np.concatenate([result["advection"] for result in results.values()])
temperature_bins = np.linspace(all_temperature.min(), all_temperature.max(), 44)
wind_bins = np.linspace(all_wind.min(), all_wind.max(), 44)
advection_bins = np.linspace(all_advection.min(), all_advection.max(), 44)

for coast, result in results.items():
    label = f"{coast}: skew={result['temperature_skew']:+.2f}"
    axes[0].hist(result["temperature"], bins=temperature_bins, density=True, histtype="step", linewidth=1.8, color=result["color"], label=label)
    label = f"{coast}: onshore={100 * result['onshore_fraction']:.1f}%"
    axes[1].hist(result["wind"], bins=wind_bins, density=True, histtype="step", linewidth=1.8, color=result["color"], label=label)
    label = f"{coast}: skew={result['advection_skew']:+.2f}"
    axes[2].hist(result["advection"], bins=advection_bins, density=True, histtype="step", linewidth=1.8, color=result["color"], label=label)

for axis in axes:
    axis.axvline(0.0, color="black", linewidth=0.8, linestyle="--")
    axis.set_ylabel("Probability density")
    axis.legend(fontsize=6)

axes[0].set_xlabel("Land Tmax anomaly (K)\n6-hourly sampled")
axes[0].set_title("a  Temperature")
axes[1].set_xlabel("Coast-normal wind (m s$^{-1}$)\npositive: ocean→land")
axes[1].set_title("b  Wind")
axes[2].set_xlabel("Temperature advection (K day$^{-1}$)\nnegative: cooling")
axes[2].set_title("c  Temperature advection")

output = Path("data/era5_wind_validation/figures/era5_temperature_wind_advection_distributions")
fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
fig.savefig(output.with_suffix(".svg"), bbox_inches="tight")

print(output.with_suffix(".png"))
