# plotodet.py - Module to plot OD output.
# Copyright (C) 2018-2019 University of Texas
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Modifications by J. Spurbeck:
# 	1. Forced a linear plot x-axis
#	2. Added RMS calculation output
#	3. Added optional numerical integration of burn span


import sys
import math
import json
import numpy
import dateutil.parser
import matplotlib.pyplot as plt
import datetime
import time
from scipy import integrate


def plot(cfgfile, inpfile, outfile, dv_start, dv_end, interactive=False, filepath=None):
    with open(cfgfile, "r") as fp:
        cfg = json.load(fp)
    with open(inpfile, "r") as fp:
        inp = json.load(fp)
    with open(outfile, "r") as fp:
        out = json.load(fp)["Estimation"]

    key = tuple(cfg["Measurements"].keys())
    dmcrun = (cfg["Estimation"].get("DMCCorrTime", 0.0) > 0.0 and
              cfg["Estimation"].get("DMCSigmaPert", 0.0) > 0.0)

    tstamp, prefit, posfit, inocov, params, estmacc, estmcov = [], [], [], [], [], [], []
    for i, o in zip(inp, out):
        tstamp.append(dateutil.parser.isoparse(i["Time"]))

        if ("PositionVelocity" in key):
            prefit.append([ix - ox for ix, ox in zip(i["PositionVelocity"],
                                                     o["PreFit"]["PositionVelocity"])])
            posfit.append([ix - ox for ix, ox in zip(i["PositionVelocity"],
                                                     o["PostFit"]["PositionVelocity"])])
        else:
            prefit.append([i[key[0]] - o["PreFit"][key[0]][-1],
                           i[key[1]] - o["PreFit"][key[1]][-1]])
            posfit.append([i[key[0]] - o["PostFit"][key[0]][-1],
                           i[key[1]] - o["PostFit"][key[1]][-1]])

        p = []
        for m in range(len(o["InnovationCovariance"])):
            p.append(3.0*numpy.sqrt(o["InnovationCovariance"][m][m]))
        inocov.append(p)

        if (len(o["EstimatedState"]) > 6):
            if (dmcrun):
                params.append(o["EstimatedState"][6:-3])
            else:
                params.append(o["EstimatedState"][6:])

        if (dmcrun):
            estmacc.append(o["EstimatedState"][-3:])

    pre = numpy.array(prefit)
    pos = numpy.array(posfit)
    cov = numpy.array(inocov)
    par = numpy.array(params)
    estmacc = numpy.array(estmacc)
    # tim = [(t - tstamp[0]).total_seconds()/3600 for t in tstamp]
    tim = tstamp

    angles = ("Azimuth", "Elevation", "RightAscension", "Declination")
    if (key[0] in angles and key[1] in angles):
        pre *= 648000.0/math.pi
        pos *= 648000.0/math.pi
        cov *= 648000.0/math.pi
        units = ("arcsec", "arcsec")
    else:
        if ("PositionVelocity" in key):
            units = ("m", "m", "m", "m/s", "m/s", "m/s")
        else:
            units = ("m", "m/s")

    if ("PositionVelocity" in key):
        ylabs = (r"$\Delta x$", r"$\Delta y$", r"$\Delta z$",
                 r"$\Delta v_x$", r"$\Delta v_y$", r"$\Delta v_z$")
        order = (1, 3, 5, 2, 4, 6)
    else:
        ylabs = key

    outfiles = []

    plt.figure(0)
    plt.suptitle("Pre-fit residuals")
    for i in range(pre.shape[-1]):
        if ("PositionVelocity" in key):
            plt.subplot(3, 2, order[i])
        else:
            plt.subplot(2, 1, i + 1)
        plt.plot(tim, pre[:,i], "ob")
        plt.xlabel("Time [utc]")
        plt.ylabel("%s [%s]" % (ylabs[i], units[i]))

    plt.tight_layout(rect = [0, 0.03, 1, 0.95])
    if (filepath is not None):
        outfiles.append(filepath + "_prefit.png")
        plt.savefig(outfiles[-1], format = "png")

    pos_post_fit_rms = 0
    vel_post_fit_rms = 0

    plt.figure(1) 
    plt.suptitle("Post-fit residuals")
    for i in range(pre.shape[-1]):
        if ("PositionVelocity" in key):
            plt.subplot(3, 2, order[i])
        else:
            plt.subplot(2, 1, i + 1)
        plt.plot(tim, pos[:,i], "ob")
        plt.plot(tim, -cov[:,i], "-r")
        plt.plot(tim,  cov[:,i], "-r", label = r"Innov. 3$\sigma$")
        plt.xlabel("Time [utc]")
        plt.ylabel("%s [%s]" % (ylabs[i], units[i]))
        if ("PositionVelocity" not in key):
            plt.ylim(-cov[i,0], cov[i,0])
        if i < 3:
            pos_post_fit_rms = numpy.mean(pos[:, i]**2) + pos_post_fit_rms
        else:
            vel_post_fit_rms = numpy.mean(pos[:, i]**2) + vel_post_fit_rms

    print("Position post-fit 3D residual RMS = " + str(numpy.sqrt(pos_post_fit_rms)) + " meters.")
    print("Velocity post-fit 3D residual RMS = " + str(numpy.sqrt(vel_post_fit_rms)*100) + " cm/s.")

    plt.tight_layout(rect = [0, 0.03, 1, 0.95])
    if (filepath is not None):
        outfiles.append(filepath + "_postfit.png")
        plt.savefig(outfiles[-1], format = "png")

    parnames, parmvals = [], []
    if (cfg["Drag"]["Coefficient"]["Estimation"] == "Estimate"):
        parnames.append(r"$C_D$")

    if (cfg["RadiationPressure"]["Creflection"]["Estimation"] == "Estimate"):
        parnames.append(r"$C_R$")

    for i in range(par.shape[-1]):
        if i == 0:
            plt.figure(2)
            plt.suptitle("Estimated parameters")

        plt.subplot(par.shape[1], 1, i + 1)
        plt.plot(tim, par[:, i], "ob")
        plt.xlabel("Time [utc]")
        plt.ylabel(parnames[i])

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    if filepath is not None:
        outfiles.append(filepath + "_estpar.png")
        plt.savefig(outfiles[-1], format="png")

    if dmcrun:

        dv_start = datetime.datetime.strptime(dv_start, '%Y-%m-%dT%H:%M:%S.%f')
        dv_end = datetime.datetime.strptime(dv_end, '%Y-%m-%dT%H:%M:%S.%f')
        dmc_radial = []
        dmc_along = []
        dmc_cross = []
        dmc_time = []

        plt.figure(3)
        plt.suptitle("DMC estimated accelerations")

        lab = [r"Radial [$\frac{m}{s^2}$]", r"In track [$\frac{m}{s^2}$]",
               r"Cross track [$\frac{m}{s^2}$]"]
        for i in range(3):
            plt.subplot(3, 1, i+1)
            plt.plot(tim, estmacc[:, i], "-b")
            plt.xlabel("Time [utc]")
            plt.ylabel(lab[i])

            # ----- numerically integrate the burn span to get delta-v estimate ----- #

            if dv_start != 0 and dv_end != 0:

                jj = 0

                for t_step in tim:

                    t_step = t_step.strftime('%Y-%m-%dT%H:%M:%S.%f')
                    t_step = datetime.datetime.strptime(t_step, '%Y-%m-%dT%H:%M:%S.%f')

                    if dv_start <= t_step <= dv_end:
                        if i == 0:
                            dmc_radial.append(abs(estmacc[jj, i]))
                            dmc_time.append(time.mktime(t_step.timetuple()))
                        elif i == 1:
                            dmc_along.append(abs(estmacc[jj, i]))
                        elif i == 2:
                            dmc_cross.append(abs(estmacc[jj, i]))

                    jj += 1

            dv_radial = integrate.cumtrapz(dmc_radial, dmc_time)
            dv_along = integrate.cumtrapz(dmc_along, dmc_time)
            dv_cross = integrate.cumtrapz(dmc_cross, dmc_time)
            print(dv_radial[-1])
            print(dv_along[-1])
            print(dv_cross[-1])
            print(math.sqrt(dv_radial[-1]**2 + dv_along[-1]**2 + dv_cross[-1]**2))

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        if filepath is not None:
            outfiles.append(filepath + "_estacc.png")
            plt.savefig(outfiles[-1], format="png")

    if interactive:
        plt.show()

    plt.close("all")
    return outfiles


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python %s config_file measurement_file output_file" % sys.argv[0])
        exit()

    if len(sys.argv) == 4:
        plot(sys.argv[1], sys.argv[2], sys.argv[3], 0, 0, interactive=True)
    elif len(sys.argv) == 6:
        plot(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], interactive=True)
