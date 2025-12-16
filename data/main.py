import pandas as pd
import numpy as np
from astropy.time import Time


def obs2cart(obs_time: Time, long: float, cos: float, sin: float, r: float = 6371.0) -> tuple[float, float, float]:
    # ITRS coordinates (Earth-fixed)
    long_rad = np.deg2rad(long)
    X = cos * np.cos(long_rad)
    Y = cos * np.sin(long_rad)
    Z = sin
    
    # Greenwich Sidereal Time (radians)
    gst = np.deg2rad(obs_time.sidereal_time('mean', 'greenwich').value * 15.0)
    
    # Rotate from ITRS to GCRS (rotate by -gst around Z-axis)
    cos_gst, sin_gst = np.cos(gst), np.sin(gst)
    x = X * cos_gst + Y * sin_gst
    y = -X * sin_gst + Y * cos_gst
    z = Z
    
    return np.array([x, y, z]) * r

def celestial2cart(ra: float, dec: float) -> tuple[float, float, float]:
    ra_rad  = np.radians(ra)
    dec_rad = np.radians(dec)
    
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    
    return x, y, z

def cart2celestial(x: float, y: float, z: float) -> tuple[float, float]:
    dec_rad = np.arctan2(z, np.sqrt(x**2 + y**2)) # mb wrong -> arcsin(z)
    dec = np.degrees(dec_rad)
    
    ra_rad = np.arctan2(y, x) # [-pi; +pi]
    ra = np.degrees(ra_rad) % 360.0 # [0; 2pi]
    
    return ra, dec

def utc2tdb(date_str: str) -> Time:
    year, month, day_frac = date_str.split()
    
    year = int(year)
    month = int(month)
    day_frac = float(day_frac)
    
    day_int = int(day_frac)
    hour_frac = (day_frac - day_int) * 24.0
    hour = int(hour_frac)
    minute_frac = (hour_frac - hour) * 60.0
    minute = int(minute_frac)
    second = (minute_frac - minute) * 60.0
    
    date_iso = f"{year:04d}-{month:02d}-{day_int:02d}T{hour:02d}:{minute:02d}:{second:06.3f}"
    
    t = Time(date_iso, format='isot', scale='utc')
    return t.tdb
  
def jd2sec(jd_date: float, reference_jd: float = 2458000.5) -> float:
    return (jd_date - reference_jd) * 24.0 * 60.0 * 60.0

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
    
    date_utc = table["date_utc"].to_numpy()
    date_tdb = [utc2tdb(date_utc[i]) for i in range(len(date_utc))]
    
    carts = np.array([obs2cart(date_tdb[i], long[i], cos[i], sin[i]) for i in range(len(long))])

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
    
    date_tdb = np.array([jd2sec(utc2tdb(date_utc[i]).jd, 0.0) for i in range(len(date_utc))])
    
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
    
    # magnitude == weight for obs 
    magn_band = table["magnitude_band"].to_numpy()
    magn = np.array([band_conversion(magn_band[i]) for i in range(len(magn_band))])    
    weights = magn2weight(magn)
    table["magnitude_band"] = weights
    table.rename(columns={"magnitude_band": "weight"}, inplace=True)
    
    return table

def write_table(table: pd.DataFrame) -> None:
    with open("output.txt", mode="w") as file:
        time = table["date_tdb"].to_numpy()
        ra, dec = table["right_ascension"].to_numpy(), table["declination"].to_numpy()  
        x = table["x_obs"].to_numpy()
        y = table["y_obs"].to_numpy()
        z = table["z_obs"].to_numpy()
        
        for i in range(len(time)):
            file.write(f"{time[i]:.10f} {ra[i]:.10f} {dec[i]:.10f} {x[i]:.10f} {y[i]:.10f} {z[i]:.10f}\n")

        file.close()
            
def main() -> None:
    table = read_table()
    table = preproc_table(table)
    write_table(table)
    print(table)


if __name__ == "__main__":
    main()
    