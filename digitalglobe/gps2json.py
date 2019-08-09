# ---------------------------------------------------------------- #
#                                                                  #
# Name:     gps2json.py                                            #
# Author:   Justin Spurbeck                                        #
# Function: Converts raw gps data into JSON format needed          #
#           for input to orbdetpy. UTC time stamps are converted   #
#           to ISO 8601 format and ECEF states are rotated to ECI  #
#           via astropy.                                           #
# Input:    raw_pos_telem2.dat                                     #
# Output:   wv04_gps_input.json                                    #
#                                                                  #
# ---------------------------------------------------------------- #

import os
import warnings
# warnings.filterwarnings("ignore")
import sys
from astropy.coordinates import CartesianRepresentation, CartesianDifferential, GCRS, ITRS
from astropy import units as u
from astropy.time import Time
import time
import numpy as np
import json

if len(sys.argv) < 3:
    print("Usage: python %s input.dat output.json" % sys.argv[0])
    exit()

print("Conversion start : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
output_file = sys.argv[2]
# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

state = []
utc = []

# ---------- Jonathan's test case ---------- #

# r_test = [-3866.513821, -2544.069958, 5225.540181]
# v_test = [-5.67275571, -1.45489028, -4.90255187]
# epoch_test = '2019-02-23T21:08:35'
# epoch_utc = Time(epoch_test, scale='utc', format='isot')
# v = CartesianDifferential(list(np.float_(r_test))*(u.km/u.s))
# r = CartesianRepresentation(list(np.float_(v_test))*u.km, differentials={'s': v})
# r_ecef = ITRS(r, obstime=epoch_utc)
# r_eci_temp = r_ecef.transform_to(GCRS(obstime=epoch_utc))
# v_eci = r_eci_temp.cartesian.differentials['s'].d_xyz
# r_eci = r_eci_temp.cartesian.xyz
#
# print(r_eci)
# print(v_eci)

json_list = []

print('Parsing file and converting to ECI...')

with open(sys.argv[1], "r") as f:

    for line in f:

        temp = line.split()

        epoch_utc = Time(temp[4], scale='utc', format='isot')
        v = CartesianDifferential(list(np.float_(temp[17:20]))*(u.m/u.s))
        r = CartesianRepresentation(list(np.float_(temp[14:17]))*u.m, differentials={'s': v})
        r_ecef = ITRS(r, obstime=epoch_utc)
        r_eci_temp = r_ecef.transform_to(GCRS(obstime=epoch_utc))
        v_eci = r_eci_temp.cartesian.differentials['s'].d_xyz  # this may be off from Jonathan's initial notebook test
        r_eci = r_eci_temp.cartesian.xyz
        utc_ms = temp[4][:-3]  # FIXME -- orbdetpy can handle 6 digits of precision in time so need to remove this

        json_object = {"Time": utc_ms + "Z",  # add this to maintain ISO8601 format for orbdetpy
                       "PositionVelocity": [
                           r_eci[0].value,
                           r_eci[1].value,
                           r_eci[2].value,
                           v_eci[0].value*1000,  # not sure why this changed to km/s
                           v_eci[1].value*1000,
                           v_eci[2].value*1000
                       ]
                       }
        json_list.append(json_object)


print('Outputting JSON file...')

with open(output_file, 'w') as outfile:
    json.dump(json_list, outfile, sort_keys=True, indent=4, ensure_ascii=False)

print("Conversion end : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))




