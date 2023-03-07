import math
import ephem
import pandas as pd
import os
from astropy import units as u
from astropy.time import Time
EPOCH = "2000-01-01 00:00:00"
shifted_epoch = (pd.to_datetime(EPOCH))
# time_since_epoch_ns * u.ns


EARTH_RADIUS = 6378135.0 # WGS72 value; taken from https://geographiclib.sourceforge.io/html/NET/NETGeographicLib_8h_source.html

# CONSTELLATION GENERATION GENERAL CONSTANTS
ECCENTRICITY = 0.0000001  # Circular orbits are zero, but pyephem does not permit 0, so lowest possible value
ARG_OF_PERIGEE_DEGREE = 0.0
PHASE_DIFF = True
EPOCH = "2000-01-01 00:00:00"
MEAN_MOTION_REV_PER_DAY = 15.19  # Altitude ~550 km
ALTITUDE_M = 550000  # Altitude ~550 km
SATELLITE_CONE_RADIUS_M = 940700 # From https://fcc.report/IBFS/SAT-MOD-20181108-00083/1569860.pdf (minimum angle of elevation: 25 deg)
MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))
MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2)) # ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
NUM_ORBS = 72
NUM_SATS_PER_ORB = 22
INCLINATION_DEGREE = 53
def read_tles(gen_data):
    """
    Read a constellation of satellites from the TLES file.

    :param filename_tles:                    Filename of the TLES (typically /path/to/tles.txt)

    :return: Dictionary:{
        "n_orbits": n_orbits,
        "n_sats_per_orbit": n_sats_per_orbit,
        "num_of_all_satellite": n_orbits * n_sats_per_orbit,
        "epoch": epoch,
        "satellites":satellites
        }
    """

    sat_data = os.listdir(gen_data)
    satellites = []
    with open(gen_data + sat_data[0] +'/tles.txt', 'r') as f:
        n_orbits, n_sats_per_orbit = [int(n) for n in f.readline().split()]
        universal_epoch = None

        i = 0
        for tles_line_1 in f:
            tles_line_2 = f.readline()
            tles_line_3 = f.readline()

            # Retrieve name and identifier
            name = tles_line_1
            sid = int(name.split()[1])
            if sid != i:
                raise ValueError("Satellite identifier is not increasing by one each line")
            i += 1

            # Fetch and check the epoch from the TLES data
            # In the TLE, the epoch is given with a Julian data of yyddd.fraction
            # ddd is actually one-based, meaning e.g. 18001 is 1st of January, or 2018-01-01 00:00.
            # As such, to convert it to Astropy Time, we add (ddd - 1) days to it
            # See also: https://www.celestrak.com/columns/v04n03/#FAQ04
            epoch_year = tles_line_2[18:20]
            epoch_day = float(tles_line_2[20:32])
            epoch = Time("20" + epoch_year + "-01-01 00:00:00", scale="tdb") + (epoch_day - 1) * u.day
            if universal_epoch is None:
                universal_epoch = epoch
            if epoch != universal_epoch:
                raise ValueError("The epoch of all TLES must be the same")

            # Finally, store the satellite information
            satellites.append(ephem.readtle(tles_line_1, tles_line_2, tles_line_3))

    return {
        "n_orbits": n_orbits,
        "n_sats_per_orbit": n_sats_per_orbit,
        "num_of_all_satellite": n_orbits * n_sats_per_orbit,
        "epoch": epoch,
        "satellites":satellites
    }


gen_data = './gen_data2/'
sats_info = read_tles(gen_data)
satellites = sats_info['satellites']
# print('%s %s' % (satellite.sublong, satellite.sublat))  # 卫星的星下点
epoch = sats_info['epoch']
time = epoch + 0 * u.minute
d = ephem.Date(str(time))



sat = ephem.EarthSatellite()
sat._epoch = str(epoch)
sat._inc = ephem.degrees(INCLINATION_DEGREE)
sat._e = ECCENTRICITY
sat._raan = ephem.degrees(5)
sat._ap = ARG_OF_PERIGEE_DEGREE
sat._M = ephem.degrees(0)
sat._n = MEAN_MOTION_REV_PER_DAY


# observer = ephem.Observer()
# observer.epoch = str(epoch)
# observer.date = str(d)
# observer.lat = 0
# observer.lon = 5
# observer.elevation = 550000
# observer._inc = ephem.degrees(53)


sat1 = satellites[22]
sat2 = satellites[23]
sat1.compute(d)
sat2.compute(d)
angle = ephem.separation(sat1, sat2)
angle = float(angle)
# a = 1 if angle > 10 else  0
angle1 = ephem.separation(sat, sat1)

# up = 0, down = 1
up_or_down_list = [0] * len(satellites)
#

#
#
#
# for i in range(len(satellites)):
#     base = i - i % 22
#     sat1 = satellites[base]
#     sat2 = satellites[i]
#     sat1.compute(d)
#     sat2.compute(d)
#     print(sat2.M)
#     angle = float(ephem.separation(sat1, sat2))
#     if angle > math.pi/2:
#         up_or_down_list[i] = 1
#
# print(up_or_down_list)

# print(float(satellite.ra))  # 赤经
# print(repr(satellite.ra))
#
# print(satellite.ra / pi * 180.0)
# print(satellite.ra)   # 赤经
# print(satellite.dec)  # 赤纬
# # print(satellite.elevation)
# print(satellite.raan)
# print(satellite.inc)
# print(satellite.n)


# EarthSatellite elements of man-made satellites:
#
# epoch — Reference epoch
# n — Mean motion, in revolutions per day
# inc — Inclination (°)
# raan — Right Ascension of ascending node (°)
# e — Eccentricity
# ap — Argument of perigee at epoch (°)
# M — Mean anomaly from perigee at epoch (°)
# decay — Orbit decay rate in revolutions per day, per day
# drag — Object drag coefficient in per earth radii
# orbit — Integer orbit number of epoch
# print('------------')

