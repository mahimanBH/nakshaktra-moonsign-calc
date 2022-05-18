########################################## ॐ  ##########################################

## Import required packages ##

import numpy as np
from astropy.coordinates import get_body
from datetime import datetime
import pytz
from tzlocal import get_localzone_name
from astropy.time import Time
from astropy.table import Table
from constants import rashi_names, nakshatram_names

def get_lat_lon(ra,dec,inclination):
    """
    Convert (ra,dec) equatorial coordinates of an object to (lat,lon) ecliptic coordinates
    
    Inputs :
    ra = right ascension of object in deg
    dec = declination of object in deg
    inclination = inclination of ecliptic to the celectial equator in deg

    Outputs :
    lat = ecliptic latitude in deg
    lon = ecliptic longitude in deg
    """
    ra = np.deg2rad(ra)
    dec = np.deg2rad(dec)
    inclination = np.deg2rad(inclination)

    beta_coord = np.arcsin(np.sin(dec)*np.cos(inclination) - np.cos(dec)*np.sin(inclination)*np.sin(ra))
    lambda_coord = np.arccos(np.cos(ra)*np.cos(dec)/np.cos(beta_coord))
    if (np.pi <= ra < 2*np.pi):
        lambda_coord = np.pi + np.arccos(-np.cos(ra)*np.cos(dec)/np.cos(beta_coord))
    if (np.pi/2 <= dec < 3*np.pi/2):
        beta_coord = (np.pi/2) + np.arccos(np.sin(dec)*np.cos(inclination) - np.cos(dec)*np.sin(inclination)*np.sin(ra))

    lat = np.rad2deg(beta_coord)
    lon = np.rad2deg(lambda_coord)

    return lat, lon

def calc_nakshatra_tithi(time,filename="nakshatra_at_test_time.pdf",tz_str="Asia/Calcutta",time_format="%Y-%m-%dT%H:%M:%S",inclination=23.4,plot_other_grahas=False):
    """
    Calculates nakshatra and tithi at input time and makes plot of grahas

    Inputs :
    time = time at which to calculate panchanga
    filename = name of file to write the plot of position of grahas in nakshatra and rashi
    tz_str = time zone of the location at which to calculate panchanga
    time_format = format of input time
    inclination = inclination of ecliptic to the celectial equator in deg

    Output :
    Plot saved at filename
    """
    nakshatram_extent = 360/27 ## deg

    start_coord = 23 + 46/60
    nakshatram_coords = np.linspace(start_coord,start_coord+nakshatram_extent*26,27)
    nakshatram_coord_labels = []

    for i in range(len(nakshatram_coords)):
        coord = nakshatram_coords[i]
        if (coord > 360):
            coord = nakshatram_coords[i] - 360
        deg_coord = coord - coord%1
        min_coord = 60*(coord%1)
        nakshatram_coord_labels.append(r"%d$^{{\circ}}$\,%d$^{{\prime}}$"%(deg_coord,min_coord))

    fmt = time_format
    date_str = time

    tz = pytz.timezone(tz_str)
    test_date = datetime.strptime(date_str, fmt)
    local_time = tz.localize(test_date,is_dst=None)
    test_date_utc = local_time.astimezone(pytz.utc)

    test_date_utc_time = Time(test_date_utc.strftime(fmt),format="isot",scale="utc")

    moon_ra = (get_body("moon",test_date_utc_time)).ra.deg
    moon_dec =  (get_body("moon",test_date_utc_time)).dec.deg

    moon_beta, moon_lambda = get_lat_lon(moon_ra,moon_dec,inclination)

    for coord_id, coord in enumerate(nakshatram_coords):
        coord %= 360
        nakshatram_coords[coord_id] = coord
        nakshatram_extent %= 360
        moon_lambda %= 360
        if (90 < coord < 270):
            if ( coord < moon_lambda < coord + nakshatram_extent or coord < moon_lambda + 360 < coord + nakshatram_extent):
                final_nakshatram = nakshatram_names[coord_id]
        else:
            if ( coord < moon_lambda < coord + nakshatram_extent or coord < moon_lambda + 360 < coord + nakshatram_extent):
                final_nakshatram = nakshatram_names[coord_id]

    rashi_start = 23 + 46/60
    rashi_extent = 360/12
    rashi_coords = np.linspace(rashi_start,rashi_start+rashi_extent*11,12)

    for coord_id,coord in enumerate(rashi_coords):
        coord %= 360
        rashi_extent %= 360
        moon_lambda %= 360
        if 77 < coord < 255:
            if coord < moon_lambda < coord + rashi_extent or coord < moon_lambda+360 < coord + rashi_extent:
                final_rashi = rashi_names[coord_id]
        else:
            if coord < moon_lambda < coord + rashi_extent or coord < moon_lambda+360 < coord + rashi_extent:
                final_rashi = rashi_names[coord_id]

    return { "Rashi" :final_rashi, "Nakshaktra" : final_nakshatram}

#local_tz = "US/Central" #"Asia/Calcutta"

date_str = "1999-02-18T02:30:00"

now = datetime.now()
# date_str = now.strftime("%Y-%m-%dT%H:%M:%S")

local_tz = get_localzone_name()

print(calc_nakshatra_tithi(date_str,tz_str=local_tz))
