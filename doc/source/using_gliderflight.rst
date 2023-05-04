Using the GliderFlight module
=============================

How to use the GliderFlight module is probably best done with a worked
example.


First we need some data to work with. Let's say we have some data
files from a Slocum glider. Typically, we would work with the
high-density data files, and would need to have access to the
engineering data, which are stored in files with the dbd extension,
and science data, which are stored in files with the ebd extension. In
this example, we use ``dbdreader`` to read the glider files. The
Python module dbdreader can be installed from PyPi using ``pip3
install --user dbdreader`` or install from source from `github
<https://github.com/smerckel/dbdreader>`_

::
   
   import numpy as np

   import dbdreader
   dbd = dbdreader.MultiDBD(pattern = "/path/to/data/*.[de]bd")


Now we have a handle to the data files, we need to extract the
required parameters. From the engineering data, we need the pitch and
the buoyancy drive. From the science data, we need the CTD
parameters. Instead of relying on the time of publishing, we take the
ctd sampling time *sci_ctd41cp_timestamp*. Occasionally, the CTD
fields contain data that have not been sampled, and default values are
returned. These default values can be detected from time stamps to be
zero, for example. After reading the values, we simply remove odd ones.

::

   tmp = dbd.get_sync("sci_ctd41cp_timestamp", "sci_water_temp", "sci_water_cond", "sci_water_pressure", "m_pitch", "m_ballast_pumped")
   _, tctd, T, C, P, pitch, buoyancy_change = np.compress(tmp[1]>0, tmp, axis=1)


Since version 0.4.0 of ``dbdreader`` we can also MultiDBD's method
get_CTD_sync():

::
   
   tctd, T, C, P, pitch, buoyancy_change = dbd.get_CTD_sync("m_pitch", "m_ballast_pumped")


In the example, we used the sensor name *m_ballast_pumped*, which is
appropriate for a shallow glider. When data from a deep glider were
used, the name should be replaced by *m_de_oil_vol*.

One of the input parameters to the glider flight model is the in-situ
density. So let's compute that one first. For this, we use the Gibbs
Seawater module, instllable using pip (for example, ``pip3
install --user gsw``), or from source from
https://github.com/TEOS-10/GSW-python. Also, we need latitude and
longitude information. This information could be retrieved from the
glider parameters ``m_gps_lat`` and ``m_gps_lon``, condensed into a
single scalar, as an array of same length as ``T``. 

::

   import gsw
   # C is given in S/m, and P in bar
   SP = gsw.SP_from_C(C*10, T, P*10)
   SA = gsw.SA_from_SP(SP, P*10, lon, lat)
   rho = gsw.rho_t_exact(SA, T, P*10)


Now we have density, we can pack all the required data into a dictionary::

  data = dict(time = tctd, pressure = P, pitch = m_pitch, buoyancy_change=buoyancy_change, density=density)


Most likely, you would want to use that data to calibrate some model
coefficients that change from depolyment to deployment, such as the
glider volume, and drag coefficient. To that end, we create an
instance from the SteadyStateCalibrate class. and populate it with the data we have got.

::

   import gliderflight

   gm = gliderflight.SteadyStateCalibrate(rho0=1024)
   gm.set_input_data(**data)
   # or, alternatively
   # gm.set_input_data(tctd, P, pitch, buoyancy_change, rho)

When we calibrate a steady-state model, we don't want to include data
points were we know the steady-state model is invalid, such as around
the transitions from down to up casts, and near the surface. Assuming
that we have dive profiles down to 100 m, may discard the first 20 m,
and the last 20 m of the dives when optimising the model against the
observed pressure rate.

::

   condition = np.logical_or(P*10<20, P*10>80)
   gm.OR(condition)

Before we can start calibrating the model, we need to set glider and
deployment specific model coefficients. Let's say we weighted the
glider and found its mass to be :math:`m_g=70` kg. (If the mass of the
glider is not know, it can be guessed. The calibration method will
adjust the volume such that the glider *density* is correct. An error
in the volume will have a small effect on the buoyancy force
calculated. As long as the mass (or the volume) is correct withing a
few percent, the errors involved are negligible.

So, we will set the mass and the volume. Also, we will set the
parasite drag to a realistic value. ::

  gm.define(mg=70.000, Vg=70e-3, Cd0=0.16)
  if not gm.undefined_parameters():
      print("We still have undefined parameters...")
      print(gm.undefined_parameters())

Now we're good to run the calibration and store the results in
calibration_result. ::

  calibration_result = gm.calibrate("Vg", "Cd0")

upon which we should have a dictionary with the keys "Vg" and "Cd0" and their optimised values. The coefficients are also updated in the model itself. So, ``gm.Cd0`` would return the same value as reported in the dictionary.

To the the glider flight results, such as angle of attack and incident
water velocity, it is not necessary to solve the model again. So, all these parameters are accessible via ``gm.modelresult`` or using the properties ``t``, ``ug``, ``wg``, ``alpha`` and ``U``.

So, we could now plot the incident water velocity as function of time::

  import matplotlib.pyplot as plt

  f, ax = plt.subplots(1,1)

  ax.plot(gm.t, gm.U, label='Incident water velocity')
  ax.set_xlabel('time (s)')
  ax.set_ylabel('U m s$^{-1}$')
  ax.legend()
  
