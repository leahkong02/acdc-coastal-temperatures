from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


files = {
    "Oregon": Path("data/era5_wind_validation/raw/era5_oregon_202407_6hourly.nc"),
    "Maine": Path("data/era5_wind_validation/raw/era5_maine_sample_20240701T12.nc"),
}

latitude_ranges = {
    "Oregon": (44.50, 45.50),
    "Maine": (44.00, 45.00),
}

output = Path("data/era5_wind_validation/figures/era5_coastal_grid_selection")
output.parent.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.25), constrained_layout=True)

for ax, (coast, path) in zip(axes, files.items()):
    ds = xr.open_dataset(path)
    lsm = ds.lsm.isel(valid_time=0).values
    latitudes = ds.latitude.values
    longitudes = ds.longitude.values
    lat_min, lat_max = latitude_ranges[coast]
    rows = np.where((latitudes >= lat_min) & (latitudes <= lat_max))[0]
    ocean_points = []
    land_points = []

    for row in rows:
        ocean_columns = np.where(lsm[row] <= 0.30)[0]
        land_columns = np.where(lsm[row] >= 0.70)[0]

        if coast == "Oregon":
            ocean_column = ocean_columns.max()
            land_column = land_columns[land_columns > ocean_column].min()
        else:
            ocean_column = ocean_columns.min()
            land_column = land_columns[land_columns < ocean_column].max()

        ocean_points.append((longitudes[ocean_column], latitudes[row]))
        land_points.append((longitudes[land_column], latitudes[row]))

    ocean_points = np.asarray(ocean_points)
    land_points = np.asarray(land_points)
    coast_longitude = 0.5 * (ocean_points[:, 0] + land_points[:, 0])
    slope, intercept = np.polyfit(ocean_points[:, 1], coast_longitude, 1)
    mean_latitude = ocean_points[:, 1].mean()
    east_per_north = slope * np.cos(np.deg2rad(mean_latitude))

    if coast == "Oregon":
        normal = np.array([1.0, -east_per_north])
    else:
        normal = np.array([-1.0, east_per_north])

    normal = normal / np.linalg.norm(normal)
    fitted_latitudes = np.linspace(latitudes[rows].min(), latitudes[rows].max(), 100)
    fitted_longitudes = slope * fitted_latitudes + intercept
    arrow_longitude = slope * mean_latitude + intercept
    arrow_dx = 0.60 * normal[0] / np.cos(np.deg2rad(mean_latitude))
    arrow_dy = 0.60 * normal[1]

    mesh = ax.pcolormesh(longitudes, latitudes, lsm, vmin=0.0, vmax=1.0, cmap="BrBG", shading="nearest")
    ax.plot(fitted_longitudes, fitted_latitudes, color="black", linewidth=1.2, linestyle="--")
    ax.scatter(ocean_points[:, 0], ocean_points[:, 1], s=24, color="#2878B5", edgecolor="white", linewidth=0.5, label="Ocean cells")
    ax.scatter(land_points[:, 0], land_points[:, 1], s=28, marker="^", color="#C85A3A", edgecolor="white", linewidth=0.5, label="Land cells")
    ax.arrow(arrow_longitude, mean_latitude, arrow_dx, arrow_dy, width=0.018, head_width=0.13, head_length=0.13, color="black", length_includes_head=True)
    ax.text(arrow_longitude + arrow_dx, mean_latitude + arrow_dy, " ocean→land", ha="left" if normal[0] > 0 else "right", va="bottom")
    ax.set_title(f"{coast}: {len(rows)} coastal pairs")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect(1.0 / np.cos(np.deg2rad(mean_latitude)))
    ax.legend(loc="lower left", fontsize=7)
    print(f"{coast}: n_east={normal[0]:+.3f}, n_north={normal[1]:+.3f}, pairs={len(rows)}")

cbar = fig.colorbar(mesh, ax=axes, shrink=0.82, pad=0.02)
cbar.set_label("ERA5 land fraction")

fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
fig.savefig(output.with_suffix(".svg"), bbox_inches="tight")

print(output.with_suffix(".png"))
