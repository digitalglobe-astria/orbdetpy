/*
 * Measurements.java - Functions to parse OD measurement files.
 * Copyright (C) 2018-2019 University of Texas
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.astria;

import com.google.gson.Gson;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.estimation.measurements.AbstractMeasurement;
import org.orekit.estimation.measurements.AngularAzEl;
import org.orekit.estimation.measurements.AngularRaDec;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.estimation.measurements.ObservableSatellite;
import org.orekit.estimation.measurements.ObservedMeasurement;
import org.orekit.estimation.measurements.PV;
import org.orekit.estimation.measurements.Range;
import org.orekit.estimation.measurements.RangeRate;
import org.orekit.estimation.measurements.modifiers.Bias;
import org.orekit.frames.Frame;
import org.orekit.frames.Transform;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.utils.PVCoordinates;

public class Measurements
{
    class Kepler
    {
	double SMA;
	double Ecc;
	double Inc;
	double RAAN;
	double ArgP;
	double MeanAnom;

	public Kepler(double a, double e, double i, double W, double o, double M)
	{
	    this.SMA = a;
	    this.Ecc = e;
	    this.Inc = i;
	    this.RAAN = W;
	    this.ArgP = o;
	    this.MeanAnom = M;
	}
    }

    class Equinoctial
    {
	double SMA;
	double Ex;
	double Ey;
	double Hx;
	double Hy;
	double Lm;

	public Equinoctial(double a, double ex, double ey, double hx, double hy, double lm)
	{
	    this.SMA = a;
	    this.Ex = ex;
	    this.Ey = ey;
	    this.Hx = hx;
	    this.Hy = hy;
	    this.Lm = lm;
	}
    }

    class State
    {
	double[] Cartesian;
	Kepler Kepler;
	Equinoctial Equinoctial;
    }

    class Measurement
    {
	String Time;
	String Station;
	Double Azimuth;
	Double Elevation;
	Double Range;
	Double RangeRate;
	Double RightAscension;
	Double Declination;
	Double[] PositionVelocity;

	public Measurement()
	{
	}

	public Measurement(Measurement src)
	{
	    Time = src.Time;
	    Station = src.Station;
	    Azimuth = src.Azimuth;
	    Elevation = src.Elevation;
	    Range = src.Range;
	    RangeRate = src.RangeRate;
	    RightAscension = src.RightAscension;
	    Declination = src.Declination;
	    PositionVelocity = src.PositionVelocity;
	}
    }

    class SimulatedMeasurement extends Measurement
    {
	State TrueState;
	Double AtmDensity;
	double[] AccGravity;
	double[] AccDrag;
	double[] AccOceanTides;
	double[] AccSolidTides;
	double[] AccThirdBodies;
	double[] AccRadiationPressure;
	double[] AccThrust;
	double[] StationState;

	public SimulatedMeasurement()
	{
	}

	public SimulatedMeasurement(SimulatedMeasurement src)
	{
	    super(src);
	    TrueState = src.TrueState;
	    AtmDensity = src.AtmDensity;
	    AccGravity = src.AccGravity;
	    AccDrag = src.AccDrag;
	    AccOceanTides = src.AccOceanTides;
	    AccSolidTides = src.AccSolidTides;
	    AccThirdBodies = src.AccThirdBodies;
	    AccRadiationPressure = src.AccRadiationPressure;
	    AccThrust = src.AccThrust;
	    StationState = src.StationState;
	}
    }

    Measurement[] rawmeas;
    ArrayList<ObservedMeasurement<?>> measobjs;

    public static Measurements loadJSON(Settings odcfg, String json)
    {
	Measurements m = new Measurements();
	m.rawmeas = new Gson().fromJson(json, Measurement[].class);
	m.getMeasurementObjects(odcfg);
	return(m);
    }

    private void getMeasurementObjects(Settings odcfg)
    {
	Settings.Measurement cazim = odcfg.cfgMeasurements.get("Azimuth");
	Settings.Measurement celev = odcfg.cfgMeasurements.get("Elevation");
	Settings.Measurement crigh = odcfg.cfgMeasurements.get("RightAscension");
	Settings.Measurement cdecl = odcfg.cfgMeasurements.get("Declination");
	Settings.Measurement crang = odcfg.cfgMeasurements.get("Range");
	Settings.Measurement crrat = odcfg.cfgMeasurements.get("RangeRate");
	Settings.Measurement cposi = odcfg.cfgMeasurements.get("PositionVelocity");

	ArrayList<Measurement> tempraw = new ArrayList<Measurement>(rawmeas.length);
	for (Measurement m: rawmeas)
	{
	    if (m.Station != null || m.PositionVelocity != null)
		tempraw.add(m);
	}
	rawmeas = tempraw.toArray(new Measurement[0]);

	measobjs = new ArrayList<ObservedMeasurement<?>>();
	for (Measurement m: rawmeas)
	{
	    GroundStation gs = null;
	    Settings.Station jsn = null;
	    if (m.Station != null)
	    {
		gs = odcfg.stations.get(m.Station);
		jsn = odcfg.cfgStations.get(m.Station);
	    }

	    AbsoluteDate time = new AbsoluteDate(DateTimeComponents.parseDateTime(m.Time), DataManager.utcscale);
	    if (m.Azimuth != null && cazim != null && celev != null)
	    {
		AngularAzEl obs = new AngularAzEl(gs, time, new double[]{m.Azimuth, m.Elevation},
						  new double[]{cazim.Error[0], celev.Error[0]},
						  new double[]{1.0, 1.0}, new ObservableSatellite(0));
		if (jsn.AzimuthBias != 0.0 || jsn.ElevationBias != 0.0)
		    obs.addModifier(new Bias<AngularAzEl>(
					new String[] {"Az", "El"}, new double[] {jsn.AzimuthBias, jsn.ElevationBias},
					new double[] {1.0, 1.0}, new double[] {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY},
					new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}));
		measobjs.add(obs);
	    }

	    if (m.RightAscension != null && crigh != null && cdecl != null)
	    {
		AngularRaDec obs = new AngularRaDec(gs, DataManager.getFrame("EME2000"), time,
						    new double[]{m.RightAscension, m.Declination},
						    new double[]{crigh.Error[0], cdecl.Error[0]},
						    new double[]{1.0, 1.0}, new ObservableSatellite(0));
		if (jsn.RightAscensionBias != 0.0 || jsn.DeclinationBias != 0.0)
		    obs.addModifier(new Bias<AngularRaDec>(
					new String[] {"RA", "Dec"}, new double[] {jsn.RightAscensionBias, jsn.DeclinationBias},
					new double[] {1.0, 1.0}, new double[] {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY},
					new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}));
		measobjs.add(obs);
	    }

	    if (m.Range != null && crang != null)
	    {
		Range obs = new Range(gs, crang.TwoWay, time, m.Range, crang.Error[0], 1.0, new ObservableSatellite(0));
		if (jsn.RangeBias != 0.0)
		    obs.addModifier(new Bias<Range>(
					new String[] {"Range"}, new double[] {jsn.RangeBias}, new double[] {1.0},
					new double[] {Double.NEGATIVE_INFINITY}, new double[] {Double.POSITIVE_INFINITY}));
		measobjs.add(obs);
	    }

	    if (m.RangeRate != null && crrat != null)
	    {
		RangeRate obs = new RangeRate(gs, time, m.RangeRate, crrat.Error[0], 1.0, crrat.TwoWay, new ObservableSatellite(0));
		if (jsn.RangeRateBias != 0.0)
		    obs.addModifier(new Bias<RangeRate>(
					new String[] {"RangeRate"}, new double[] {jsn.RangeRateBias}, new double[] {1.0},
					new double[] {Double.NEGATIVE_INFINITY}, new double[] {Double.POSITIVE_INFINITY}));
		measobjs.add(obs);
	    }

	    if (m.PositionVelocity != null && cposi != null)
	    {
		Double[] X = m.PositionVelocity;
		if (cposi.ReferenceFrame != null)
		{
		    Frame fromframe = DataManager.getFrame(cposi.ReferenceFrame);
		    Transform xfm = fromframe.getTransformTo(odcfg.propframe, time);
		    PVCoordinates frompv = new PVCoordinates(new Vector3D(X[0], X[1], X[2]), new Vector3D(X[3], X[4], X[5]));
		    PVCoordinates topv = xfm.transformPVCoordinates(frompv);
		    Vector3D p = topv.getPosition();
		    Vector3D v = topv.getVelocity();
		    X = new Double[]{p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ()};
		}

		PV obs = new PV(time, new Vector3D(X[0], X[1], X[2]), new Vector3D(X[3], X[4], X[5]),
				cposi.Error, 1.0, new ObservableSatellite(0));
		if (jsn != null && jsn.PositionVelocityBias != null)
		    obs.addModifier(new Bias<PV>(
					new String[] {"x", "y", "z", "Vx", "Vy", "Vz"}, jsn.PositionVelocityBias,
					new double[] {1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
					new double[] {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY,
						      Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY},
					new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY,
						      Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}));
		measobjs.add(obs);
	    }
	}
    }
}
