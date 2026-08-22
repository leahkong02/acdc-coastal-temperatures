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
    mean_latitude = latitudes[rows].mean()
    east_per_north = slope * np.cos(np.deg2rad(mean_latitude))

    if settings["land_side"] == "east":
        normal = np.array([1.0, -east_per_north])
    else:
        normal = np.array([-1.0, east_per_north])

    normal = normal / np.linalg.norm(normal)
    yearly_frames = []

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
        frame["marine_proxy"] = np.maximum(frame["u_normal"], 0.0) * frame["delta_t"]
        daily = frame.resample("1D").mean()
        daily["year"] = year
        yearly_frames.append(daily)

    daily = pd.concat(yearly_frames)
    yearly = daily.groupby("year")
    yearly_onshore = yearly["u_normal"].apply(lambda values: np.mean(values > 0.0))
    yearly_active = yearly.apply(lambda values: np.mean((values.u_normal > 0.0) & (values.delta_t < 0.0)), include_groups=False)
    onshore = daily.u_normal > 0.0
    active = onshore & (daily.delta_t < 0.0)
    results[coast] = {
        "daily": daily,
        "normal": normal,
        "pairs": len(rows),
        "onshore_mean": yearly_onshore.mean(),
        "onshore_sd": yearly_onshore.std(ddof=1),
        "active_mean": yearly_active.mean(),
        "active_sd": yearly_active.std(ddof=1),
        "onshore_delta_t": daily.loc[onshore, "delta_t"],
        "color": settings["color"],
    }

    print(
        f"{coast}: pairs={len(rows)}, normal=({normal[0]:+.3f},{normal[1]:+.3f}), "
        f"onshore={yearly_onshore.mean():.3f}+/-{yearly_onshore.std(ddof=1):.3f}, "
        f"active_cooling={yearly_active.mean():.3f}+/-{yearly_active.std(ddof=1):.3f}, "
        f"mean_un={daily.u_normal.mean():+.3f} m/s, "
        f"mean_delta_t_onshore={daily.loc[onshore, 'delta_t'].mean():+.3f} K"
    )

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), constrained_layout=True)
all_wind = np.concatenate([result["daily"].u_normal.values for result in results.values()])
all_delta_t = np.concatenate([result["onshore_delta_t"].values for result in results.values()])
wind_bins = np.linspace(all_wind.min(), all_wind.max(), 42)
delta_t_bins = np.linspace(all_delta_t.min(), all_delta_t.max(), 42)

for coast, result in results.items():
    daily = result["daily"]
    label = f"{coast}: onshore {100 * result['onshore_mean']:.1f}±{100 * result['onshore_sd']:.1f}%"
    axes[0].hist(daily.u_normal, bins=wind_bins, density=True, histtype="step", linewidth=2.0, color=result["color"], label=label)
    label = f"{coast}: active {100 * result['active_mean']:.1f}±{100 * result['active_sd']:.1f}%"
    axes[1].hist(result["onshore_delta_t"], bins=delta_t_bins, density=True, histtype="step", linewidth=2.0, color=result["color"], label=label)

axes[0].axvline(0.0, color="black", linewidth=0.9, linestyle="--")
axes[0].set_xlabel("Daily-mean coast-normal wind (m s$^{-1}$)")
axes[0].set_ylabel("Probability density")
axes[0].set_title("a  Positive means ocean→land")
axes[0].legend(fontsize=7)
axes[1].axvline(0.0, color="black", linewidth=0.9, linestyle="--")
axes[1].set_xlabel("$T_{ao}-T_{al}$ on onshore days (K)")
axes[1].set_ylabel("Probability density")
axes[1].set_title("b  Negative means marine cooling")
axes[1].legend(fontsize=7)

output = Path("data/era5_wind_validation/figures/era5_coastal_wind_first_diagnostic")
fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
fig.savefig(output.with_suffix(".svg"), bbox_inches="tight")

print(output.with_suffix(".png"))
