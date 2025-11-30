import pandas as pd
import numpy as np
from astropy.time import Time


def obs2cart(long: float, cos: float, sin: float, r: float = 6371.0) -> tuple[float, float, float]:
    long_rad = np.radians(long)
    
    x = np.cos(long_rad) * cos * r
    y = cos * np.sin(long_rad) * r
    z = sin * r

    return x, y, z

def celestial2cart(ra: float, dec: float) -> tuple[float, float, float]:
    ra_rad  = np.radians(ra)
    dec_rad = np.radians(dec)
    
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    
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

def ra2angle(ra: str) -> float:
    hours, mins, secs = ra.split()
    
    angle = float(hours) + (float(mins) / 60) + (float(secs) / 3600) # to hours
    angle *= 15 # to angle, 1 hour = 15 degree
    
    return angle
    
def dec2angle(dec: str) -> float:
    sign = dec[0]
    sign = -1 if sign == "-" else +1
    
    dec = dec[1:]
    hours, mins, secs = dec.split()
    
    angle = float(hours) + (float(mins) / 60) + (float(secs) / 3600) # to hours
    angle *= 15 # to angle, 1 hour = 15 degree
    
    angle *= sign

    return angle

def band_conversion(magn_band: str) -> float:
    V_BAND_CORRECTIONS = \
    {
        '': -0.8,
        'U': -1.3,
        'B': -0.8,
        'g': -0.35,
        'V': 0.0,
        'r': 0.14,
        'R': 0.4,
        'C': 0.4,
        'W': 0.4,
        'i': 0.32,
        'z': 0.26,
        'I': 0.8,
        'J': 1.2,
        'w': -0.13,
        'y': 0.32,
        'L': 0.2,
        'H': 1.4,
        'K': 1.7,
        'Y': 0.7,
        'G': 0.28,
        'v': 0.0,
        'c': -0.05,
        'o': 0.33,
        'u': 2.5
    }

    parts = magn_band.split()
    
    band = "" if len(parts) == 2 else parts[1]
    magn = float(parts[0])
    
    magn += V_BAND_CORRECTIONS[band]
    
    return magn    

def magn2weight(magn: np.array) -> float:
    min_val = magn.min()
    max_val = magn.max()
    
    min_max_magn = [1.0 - (magn[i] - min_val) / (max_val - min_val) for i in range(len(magn))] # inverted 1.0 = bright, 0.0 = dark
    
    weights = [max(x**0.5, 0.1) for x in min_max_magn] 
    
    return weights

def obsobj2objgeo(obs: np.array, obj: np.array) -> np.array:
    obs_dir = np.array([vec / np.linalg.norm(vec) for vec in obs])
    
    res = np.array([obs_dir[i] + obj[i] for i in range(len(obs))])
    
    return res

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
    
    carts = np.array([obs2cart(long[i], cos[i], sin[i]) for i in range(len(long))])

    x = carts[:, 0]
    y = carts[:, 1]
    z = carts[:, 2]
    
    table.drop(["cos", "sin", "longitude"], inplace=True, axis=1)

    table["x"] = x
    table["y"] = y
    table["z"] = z
    
    table.dropna(subset=["declination", "right_ascension", "x", "y", "z", "magnitude_band"], inplace=True)

    table.rename(columns={"x":"x_obs", "y":"y_obs", "z":"z_obs"}, inplace=True)
    
    # observation time
    date_utc = table["date_utc"].to_numpy()
    
    date_tdb = np.array([jd2sec(utc2tdb(date_utc[i]).jd) for i in range(len(date_utc))])
    
    table.insert(0, "date_tdb", date_tdb)
    table.drop("date_utc", axis=1, inplace=True)
    
    table["date_tdb"] -= min(date_tdb) # shift time 
    
    table.sort_values(by="date_tdb", inplace=True) # order by time
    
    # right ascension and declination
    ra  = table["right_ascension"].to_numpy()
    dec = table["declination"].to_numpy()
    
    ra_angle = np.array([ra2angle(ra[i]) for i in range(len(ra))])
    dec_angle = np.array([dec2angle(dec[i]) for i in range(len(dec))])
    
    table["right_ascension"] = ra_angle
    table["declination"] = dec_angle
    
    cart_coords = np.array([celestial2cart(ra_angle[i], dec_angle[i]) for i in range(len(ra))])

    table["x_obj"] = cart_coords[:, 0]
    table["y_obj"] = cart_coords[:, 1]
    table["z_obj"] = cart_coords[:, 2]
    table.drop(columns=["right_ascension", "declination"], inplace=True)
    
    # magnitude == weight for obs 
    magn_band = table["magnitude_band"].to_numpy()
    magn = np.array([band_conversion(magn_band[i]) for i in range(len(magn_band))])    
    weights = magn2weight(magn)
    table["magnitude_band"] = weights
    table.rename(columns={"magnitude_band": "weight"}, inplace=True)
    
    # obj coords to geocenter
    x_obj, y_obj, z_obj = table["x_obj"].to_numpy(), table["y_obj"].to_numpy(), table["z_obj"].to_numpy()
    obj_coords = np.column_stack((x_obj,y_obj, z_obj))
    
    x_obs, y_obs, z_obs = table["x_obs"].to_numpy(), table["y_obs"].to_numpy(), table["z_obs"].to_numpy()
    obs_coords = np.column_stack((x_obs, y_obs, z_obs))

    obj_geo = obsobj2objgeo(obs_coords, obj_coords)
    
    x, y, z = obj_geo[:, 0], obj_geo[:, 1], obj_geo[:, 2]
    
    table["x_obj"] = x
    table["y_obj"] = y
    table["z_obj"] = z
    
    table.drop(columns=["x_obs", "y_obs",  "z_obs"], inplace=True)
    
    return table
    
def main() -> None:
    table = read_table()
    table = preproc_table(table)
    
    print(table)


if __name__ == "__main__":
    main()
    