from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlretrieve
import json
import tempfile
import time

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


years = range(1985, 2015)
data_directory = Path("data/cpc_northern_hemisphere_skewness/raw")
processed_directory = Path("data/cpc_northern_hemisphere_skewness/processed")
figure_directory = Path("data/cpc_northern_hemisphere_skewness/figures")
data_directory.mkdir(parents=True, exist_ok=True)
processed_directory.mkdir(parents=True, exist_ok=True)
figure_directory.mkdir(parents=True, exist_ok=True)

base_url = "https://psl.noaa.gov/thredds/ncss/grid/Datasets/cpc_global_temp"

longitude_sections = {
    "lon_000_090": (0.25, 89.75),
    "lon_090_180": (90.25, 179.75),
    "lon_180_270": (180.25, 269.75),
    "lon_270_360": (270.25, 359.75),
}

for year in years:
    for section_name, (west, east) in longitude_sections.items():
        output = data_directory / f"cpc_tmax_{year}_maysep_north_of_20n_{section_name}.nc"

        if output.exists():
            with xr.open_dataset(output, engine="h5netcdf") as dataset:
                if dataset.sizes.get("time") != 153:
                    raise RuntimeError(f"Existing file has the wrong time length: {output}")
                if dataset.sizes.get("lat") != 141 or dataset.sizes.get("lon") != 180:
                    raise RuntimeError(f"Existing file has the wrong grid: {output}")
            print(f"skip verified: {output}", flush=True)
            continue

        parameters = {
            "var": "tmax",
            "north": 90,
            "west": west,
            "east": east,
            "south": 20,
            "horizStride": 1,
            "time_start": f"{year}-05-01T00:00:00Z",
            "time_end": f"{year}-09-30T00:00:00Z",
            "timeStride": 1,
            "accept": "netcdf4",
        }
        url = f"{base_url}/tmax.{year}.nc?{urlencode(parameters)}"

        for attempt in range(1, 4):
            with tempfile.NamedTemporaryFile(prefix=f"cpc_tmax_nh_{year}_", suffix=".nc", delete=False) as temporary_file:
                temporary_path = Path(temporary_file.name)

            try:
                print(f"download: {year}, {section_name}, attempt {attempt}", flush=True)
                urlretrieve(url, temporary_path)

                with xr.open_dataset(temporary_path, engine="h5netcdf") as dataset:
                    if dataset.sizes.get("time") != 153:
                        raise RuntimeError(f"Downloaded file has the wrong time length: {temporary_path}")
                    if dataset.sizes.get("lat") != 141 or dataset.sizes.get("lon") != 180:
                        raise RuntimeError(f"Downloaded file has the wrong grid: {temporary_path}")
                    if dataset.tmax.dims != ("time", "lat", "lon"):
                        raise RuntimeError(f"Unexpected dimensions: {dataset.tmax.dims}")

                temporary_path.replace(output)
                break
            except Exception as error:
                print(f"failed: {year}, {section_name}, attempt {attempt}: {error}", flush=True)
                if attempt == 3:
                    raise
                time.sleep(5 * attempt)

        print(f"verified: {output}", flush=True)

year_files = []

for year in years:
    section_files = [
        data_directory / f"cpc_tmax_{year}_maysep_north_of_20n_{section_name}.nc"
        for section_name in longitude_sections
    ]

    if not all(path.exists() for path in section_files):
        raise RuntimeError(f"Missing longitude section for {year}")

    year_files.append(section_files)

climatology_sum = None
climatology_count = None
reference_months = None
reference_days = None
latitudes = None
longitudes = None

for section_files in year_files:
    temperature_sections = []
    longitude_sections_for_year = []

    for path in section_files:
        with xr.open_dataset(path, engine="h5netcdf") as dataset:
            temperature_sections.append(dataset.tmax.load().values.astype(np.float32))
            longitude_sections_for_year.append(dataset.lon.values.copy())

            if latitudes is None:
                latitudes = dataset.lat.values.copy()
            elif not np.array_equal(latitudes, dataset.lat.values):
                raise RuntimeError(f"Latitude mismatch: {path}")

            if reference_months is None:
                months = dataset.time.dt.month.values
                days = dataset.time.dt.day.values
            elif not np.array_equal(reference_months, dataset.time.dt.month.values):
                raise RuntimeError(f"Month mismatch: {path}")

    temperature = np.concatenate(temperature_sections, axis=2)

    if climatology_sum is None:
        climatology_sum = np.zeros(temperature.shape, dtype=np.float64)
        climatology_count = np.zeros(temperature.shape, dtype=np.uint16)
        reference_months = months.copy()
        reference_days = days.copy()
        longitudes = np.concatenate(longitude_sections_for_year)

    if not np.array_equal(months, reference_months) or not np.array_equal(days, reference_days):
        raise RuntimeError(f"Calendar mismatch: {section_files[0]}")

    valid = np.isfinite(temperature)
    climatology_sum += np.where(valid, temperature, 0.0)
    climatology_count += valid

climatology = np.divide(
    climatology_sum,
    climatology_count,
    out=np.full(climatology_sum.shape, np.nan, dtype=np.float32),
    where=climatology_count > 0,
)
del climatology_sum
del climatology_count

window = 31
half_window = window // 2
jja_indices = np.where((reference_months >= 6) & (reference_months <= 8))[0]
smoothed_jja_climatology = np.full(
    (jja_indices.size, latitudes.size, longitudes.size),
    np.nan,
    dtype=np.float32,
)

for position, day_index in enumerate(jja_indices):
    block = climatology[day_index - half_window : day_index + half_window + 1]
    block_valid = np.isfinite(block)
    block_sum = np.where(block_valid, block, 0.0).sum(axis=0, dtype=np.float64)
    block_count = block_valid.sum(axis=0)
    smoothed_jja_climatology[position] = np.divide(
        block_sum,
        block_count,
        out=np.full(block_sum.shape, np.nan, dtype=np.float32),
        where=block_count > 0,
    )

del climatology

sample_count = np.zeros((latitudes.size, longitudes.size), dtype=np.int32)
sum_1 = np.zeros((latitudes.size, longitudes.size), dtype=np.float64)
sum_2 = np.zeros((latitudes.size, longitudes.size), dtype=np.float64)
sum_3 = np.zeros((latitudes.size, longitudes.size), dtype=np.float64)

for section_files in year_files:
    temperature_sections = []

    for path in section_files:
        with xr.open_dataset(path, engine="h5netcdf") as dataset:
            temperature_sections.append(dataset.tmax.isel(time=jja_indices).load().values.astype(np.float32))

    temperature = np.concatenate(temperature_sections, axis=2)

    anomaly = temperature - smoothed_jja_climatology
    valid = np.isfinite(anomaly)
    anomaly = np.where(valid, anomaly, 0.0)
    sample_count += valid.sum(axis=0)
    sum_1 += anomaly.sum(axis=0, dtype=np.float64)
    sum_2 += (anomaly**2).sum(axis=0, dtype=np.float64)
    sum_3 += (anomaly**3).sum(axis=0, dtype=np.float64)

mean = np.divide(sum_1, sample_count, out=np.full_like(sum_1, np.nan), where=sample_count > 0)
raw_moment_2 = np.divide(sum_2, sample_count, out=np.full_like(sum_2, np.nan), where=sample_count > 0)
raw_moment_3 = np.divide(sum_3, sample_count, out=np.full_like(sum_3, np.nan), where=sample_count > 0)
variance = raw_moment_2 - mean**2
third_central_moment = raw_moment_3 - 3.0 * mean * raw_moment_2 + 2.0 * mean**3
skewness = third_central_moment / variance**1.5
minimum_valid_days = int(0.80 * len(years) * 92)
skewness[(sample_count < minimum_valid_days) | (variance <= 0.0)] = np.nan
standard_deviation = np.sqrt(variance)

processed_output = processed_directory / "cpc_northern_hemisphere_jja_tmax_anomaly_skewness_1985_2014.nc"
processed = xr.Dataset(
    data_vars={
        "skewness": (("lat", "lon"), skewness.astype(np.float32)),
        "standard_deviation": (("lat", "lon"), standard_deviation.astype(np.float32)),
        "valid_days": (("lat", "lon"), sample_count),
    },
    coords={"lat": latitudes, "lon": longitudes},
    attrs={
        "title": "Northern Hemisphere JJA daily Tmax anomaly skewness north of 20N",
        "source": "NOAA CPC Global Unified Temperature",
        "period": "1985-2014",
        "season": "June-August",
        "anomaly_method": "Daily climatology smoothed with a 31-day centered running mean",
        "skewness_definition": "mean((x-mean(x))**3) / mean((x-mean(x))**2)**1.5",
        "validity_rule": "At least 80 percent of the 2760 expected JJA days",
    },
)
processed.to_netcdf(processed_output, engine="h5netcdf")

boundary_output = Path("data/natural_earth_110m_countries.geojson")
boundary_output.parent.mkdir(parents=True, exist_ok=True)

if not boundary_output.exists():
    boundary_url = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_110m_admin_0_countries.geojson"
    with tempfile.NamedTemporaryFile(prefix="natural_earth_", suffix=".geojson", delete=False) as temporary_file:
        temporary_boundary = Path(temporary_file.name)
    print("download: Natural Earth country boundaries", flush=True)
    urlretrieve(boundary_url, temporary_boundary)
    temporary_boundary.replace(boundary_output)

with boundary_output.open() as stream:
    boundaries = json.load(stream)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
})

theta = np.deg2rad(longitudes)
radius = 90.0 - latitudes
figure = plt.figure(figsize=(7.2, 7.2))
axis = figure.add_subplot(111, projection="polar")
normalization = mpl.colors.TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
mesh = axis.pcolormesh(
    theta,
    radius,
    skewness,
    cmap="RdBu_r",
    norm=normalization,
    shading="auto",
)
axis.contour(theta, radius, skewness, levels=[0.0], colors="0.25", linewidths=0.35)

for feature in boundaries["features"]:
    geometry = feature["geometry"]
    polygons = geometry["coordinates"] if geometry["type"] == "MultiPolygon" else [geometry["coordinates"]]

    for polygon in polygons:
        for ring in polygon:
            coordinates = np.asarray(ring)

            if coordinates[:, 1].max() < 19.0:
                continue

            breaks = np.where(np.abs(np.diff(coordinates[:, 0])) > 180.0)[0] + 1

            for segment in np.split(coordinates, breaks):
                if segment.shape[0] < 2:
                    continue
                axis.plot(
                    np.deg2rad(segment[:, 0]),
                    90.0 - segment[:, 1],
                    color="0.18",
                    linewidth=0.45,
                )

axis.set_theta_zero_location("N")
axis.set_theta_direction(1)
axis.set_ylim(0.0, 70.0)
axis.set_thetagrids(
    np.arange(0.0, 360.0, 45.0),
    labels=["0°", "45°E", "90°E", "135°E", "180°", "135°W", "90°W", "45°W"],
)
axis.set_rticks([10.0, 30.0, 50.0, 70.0])
axis.set_yticklabels(["80°N", "60°N", "40°N", "20°N"])
axis.set_rlabel_position(18.0)
axis.grid(color="0.45", linewidth=0.35, alpha=0.55)
axis.spines["polar"].set_linewidth(0.8)
axis.set_title("Northern Hemisphere summer daily-maximum temperature anomaly skewness", fontsize=11, pad=22)

colorbar = figure.colorbar(mesh, ax=axis, orientation="horizontal", pad=0.085, fraction=0.055, extend="both")
colorbar.set_label("Skewness of JJA daily Tmax anomalies")
colorbar.set_ticks(np.arange(-1.0, 1.01, 0.2))

figure.text(
    0.5,
    0.018,
    "NOAA CPC 0.5° observations, 1985–2014; 31-day smoothed daily climatology; ≥80% valid JJA days",
    ha="center",
    fontsize=7.0,
    color="0.25",
)
figure.subplots_adjust(left=0.05, right=0.95, top=0.91, bottom=0.13)

figure_output = figure_directory / "cpc_northern_hemisphere_jja_tmax_anomaly_skewness_1985_2014"
figure.savefig(figure_output.with_suffix(".png"), dpi=300, bbox_inches="tight")
plt.close(figure)

finite = np.isfinite(skewness)
print(f"valid grid cells: {finite.sum()}")
print(f"skewness range: {np.nanmin(skewness):+.3f} to {np.nanmax(skewness):+.3f}")
print(f"median skewness: {np.nanmedian(skewness):+.3f}")
print(processed_output)
print(figure_output.with_suffix(".png"))
