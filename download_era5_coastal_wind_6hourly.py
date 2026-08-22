from pathlib import Path

import cdsapi


sites = {
    "oregon": [47.0, -125.5, 44.0, -122.5],
    "new_jersey": [41.0, -75.5, 38.0, -72.5],
    "maine": [46.0, -70.5, 44.0, -66.0],
}

years = range(2015, 2025)
months = ["05", "06", "07", "08", "09"]
days = [f"{day:02d}" for day in range(1, 32)]
times = ["00:00", "06:00", "12:00", "18:00"]
variables = [
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "2m_temperature",
    "sea_surface_temperature",
]

output_directory = Path("data/era5_wind_validation/raw")
output_directory.mkdir(parents=True, exist_ok=True)
client = cdsapi.Client()

pilot_output = output_directory / "era5_maine_sample_20240701T12.nc"

if not pilot_output.exists():
    pilot_request = {
        "product_type": ["reanalysis"],
        "variable": variables + ["land_sea_mask"],
        "year": ["2024"],
        "month": ["07"],
        "day": ["01"],
        "time": ["12:00"],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": sites["maine"],
    }

    print(f"download: {pilot_output}")
    client.retrieve("reanalysis-era5-single-levels", pilot_request, str(pilot_output))

for site, area in sites.items():
    for year in years:
        output = output_directory / f"era5_{site}_{year}_maysep_6hourly.nc"

        if output.exists():
            print(f"skip existing: {output}")
            continue

        request = {
            "product_type": ["reanalysis"],
            "variable": variables,
            "year": [str(year)],
            "month": months,
            "day": days,
            "time": times,
            "data_format": "netcdf",
            "download_format": "unarchived",
            "area": area,
        }

        print(f"download: {output}")
        client.retrieve("reanalysis-era5-single-levels", request, str(output))

print("all downloads complete")
