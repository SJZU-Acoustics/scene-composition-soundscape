# Analysis 22 (fetch step) - Build the survey-date covariate table:
# Shenyang daily weather (meteostat, stations ZYTX0 and 54342) and daily-mean
# AQI and PM2.5 (quotsoft.net/air hourly china_cities archive, Shenyang
# column). One row per distinct survey date. Run once, from the repository
# root, after build/build_core.py; analyses 22 and 27 join the resulting
# analysis_22_weather_aqi/data/date_covariates.csv (see data/README.md there).
# Requires the packages requests and meteostat (2.x) and network access.

import re
import time
from datetime import datetime, timedelta
from pathlib import Path

import pandas as pd
import requests
import meteostat

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = ANALYSIS_DIR / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

CLEAN_CSV = (ANALYSIS_DIR.parent / "analysis_01_clean_core" / "outcome"
             / "clean_data" / "clean_core_analysis.csv")


def parse_date(x):
    if pd.isna(x):
        return None
    s = re.sub(r"\.0$", "", str(x).strip()).replace(".", "-")
    for fmt in ("%Y-%m-%d", "%Y%m%d"):
        try:
            return datetime.strptime(s, fmt).date()
        except ValueError:
            continue
    return None


core = pd.read_csv(CLEAN_CSV, usecols=["survey_date"], dtype=str)
dates = sorted({d for d in core["survey_date"].map(parse_date) if d})
print(f"{len(dates)} distinct survey dates: {dates[0]} .. {dates[-1]}")

# Weather (meteostat 2.x; ZYTX0 = Shenyang 2.5 km, 54342 = Luoguantun 10 km,
# coalesced per column because prcp/cldc coverage differs by station) ---------
WX_COLS = ["temp", "tmin", "tmax", "rhum", "prcp", "wspd", "cldc"]


def fetch_station(sid):
    ts = meteostat.daily(sid, datetime(dates[0].year, 1, 1),
                         datetime(dates[-1].year, 12, 31))
    df = ts.fetch().reset_index()
    df["date"] = pd.to_datetime(df["time"]).dt.date
    cols = [c for c in WX_COLS if c in df.columns]
    return df[df["date"].isin(dates)][["date"] + cols].set_index("date")


wx = fetch_station("ZYTX0")
wx = wx.combine_first(fetch_station("54342")).reset_index()
print(f"weather rows matched: {len(wx)} / {len(dates)}")
print(wx.notna().sum())

# AQI (quotsoft.net/air) ------------------------------------------------------
aqi_rows = []
for d in dates:
    url = ("https://quotsoft.net/air/data/china_cities_"
           f"{d.strftime('%Y%m%d')}.csv")
    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        df = pd.read_csv(pd.io.common.StringIO(r.text))
        sub = df[df["type"].isin(["AQI", "PM2.5"])]
        rec = {"date": d}
        for t in ("AQI", "PM2.5"):
            vals = pd.to_numeric(
                sub.loc[sub["type"] == t, "沈阳"], errors="coerce").dropna()
            rec["aqi_mean" if t == "AQI" else "pm25_mean"] = (
                vals.mean() if len(vals) else None)
        aqi_rows.append(rec)
    except Exception as e:
        print(f"AQI fetch failed for {d}: {e}")
        aqi_rows.append({"date": d, "aqi_mean": None, "pm25_mean": None})
    time.sleep(1.0)
aqi = pd.DataFrame(aqi_rows)
print(f"AQI rows with data: {aqi['aqi_mean'].notna().sum()} / {len(aqi)}")

out = pd.DataFrame({"date": dates}).merge(wx, on="date", how="left").merge(
    aqi, on="date", how="left")
out.to_csv(DATA_DIR / "date_covariates.csv", index=False)
print(out.to_string())
