# Date covariates

One file belongs here; it is not shipped with the repository.

| File | Status | Contents |
|---|---|---|
| `date_covariates.csv` | **built** by `../code/fetch_date_covariates.py` (run from the repository root, after `build/build_core.py`; needs `requests`, `meteostat` and network access) | One row per distinct survey date (32 dates): daily mean, minimum and maximum temperature, relative humidity, precipitation, wind speed and cloud cover for Shenyang from Meteostat (stations ZYTX0 and 54342, coalesced per column), and the daily mean AQI and PM2.5 for Shenyang from the hourly archive at quotsoft.net/air |

`run_all.R` builds the file automatically when it is absent. If the build fails
(no network, or a source no longer serves the dates), analyses 22 and 27 and
Supplementary Tables S16 and S17 are skipped with a message; every other
display item is unaffected.
