import pandas as pd
import numpy as np
from astropy.time import Time


def astro2cart(long: float, cos: float, sin: float, r: float = 6371.0) -> tuple[float, float, float]:
    x = np.cos(long) * cos * r
    y = np.sin(long) * cos * r
    z = sin * r 

    return x, y, z

def utc2tdb(date_str: str) -> Time:
    year, month, day_frac = date_str.split()
    
    year = int(year)
    month = int(month)
    day_whole = int(float(day_frac))
    time_fraction = float(day_frac) - day_whole
    
    seconds = time_fraction * 24 * 60 * 60
    hours = int(seconds / 3600)
    minutes = int((seconds % 3600) / 60)
    seconds = seconds % 60 
    
    iso = f"{year:04d}-{month:02d}-{day_whole:02d}T{hours:02d}:{minutes:02d}:{seconds:09.6f}"
    
    date = Time(iso, format='isot', scale="utc")
    
    return date.tdb
  
def jd2sec(jd_date: float, reference_jd: float = 2458000.5) -> float:
    return (jd_date - reference_jd) * 24 * 60 * 60

def read_table() -> pd.DataFrame:
    colspecs_oumuamua  = [(0, 2), (4, 21), (23, 34), (36, 48), (49, 56), (57, 64), (66, None)]
    col_names_oumuamua = ["note1", "date_utc", "right_ascension", "declination", "magnitude_band", "obs_name", "code"]
    oumuamua_obs  = pd.read_fwf("1I.txt", colspecs=colspecs_oumuamua, names=col_names_oumuamua)
    
    
    colspecs_obs  = [(0, 3), (4, 12), (13, 23), (25, 34), (36, None)]
    col_names_obs = ["code", "longitude", "cos", "sin", "place"]
    observatories = pd.read_fwf("obscodes.txt", colspecs=colspecs_obs, names=col_names_obs)
    
    merged_df = pd.merge(left=oumuamua_obs, right=observatories, how="inner", on="code")
    
    return merged_df

def preproc_table(table: pd.DataFrame) -> pd.DataFrame:
    # observatory coords
    table.drop(["note1", "obs_name", "place", "code"], inplace=True, axis=1)
    
    long = table["longitude"].to_numpy()
    cos  = table["cos"].to_numpy()
    sin  = table["sin"].to_numpy()
    
    carts = np.array([astro2cart(long[i], cos[i], sin[i]) for i in range(len(long))])

    x = carts[:, 0]
    y = carts[:, 1]
    z = carts[:, 2]
    
    table.drop(["cos", "sin", "longitude"], inplace=True, axis=1)

    table["x"] = x
    table["y"] = y
    table["z"] = z
    
    table.dropna(subset=["declination", "right_ascension", "x", "y", "z"], inplace=True)

    # observation time
    date_utc = table["date_utc"].to_numpy()
    
    date_tdb = np.array([jd2sec(utc2tdb(date_utc[i]).jd) for i in range(len(date_utc))])
    
    table.insert(0, "date_tdb", date_tdb)
    table.drop("date_utc", axis=1, inplace=True)

    return table
    
def main() -> None:
    table = read_table()
    table = preproc_table(table)
    
    print(table)


if __name__ == "__main__":
    main()
    