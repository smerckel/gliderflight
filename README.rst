|PyPI version| |Docs badge| |License|

GliderFlight for Slocum ocean gliders
=====================================

Synopsis
--------

Gliderflight is a python module to calibrate a model that predicts the
glider flight through water. The model results can be used to estimate
the speed through water, a parameter which is required to compute
turbulent dissipation rates from temperature microstructure or shear
probe data, collected with a turbulence profiler mounted on top of an
ocean glider.

Background
------------

The dissipation rate of turbulent kinetic energy is a parameter that
plays a key role in many physical and biogeo chemical processes in
oceans and coastal seas. However, direct oceanic measurements of
turbulence are relatively scarce, as most observations stem from
free-falling profilers, operated from seagoing vessels.


An emerging alternative to ship-based profiling is the use of ocean
gliders with mounted turbulence profilers.  A required parameter in
the processing of microstructure shear and temperature measurements is
the speed of flow past the sensors. This speed can be measured
directly with additional sensor, such as an electromagnetic current
meter or mounted acoustic Doppler current profiler, but often gliders
are not equipped with additional velocity sensors. Alternatively, a
glider flight model can be used to estimate the speed through
water. Such a model is described in the paper *A dynamic flight model
for Slocum gliders and implications for turbulence microstructure
measurements* [merckelbach2019]_. This Python
model implements the steady-state and dynamic glider flight models,
described therein.

Documentation
-------------

Documentation of this software package can be found at 
https://gliderflight.readthedocs.io/en/latest/

Steady-state model
------------------

The steady-state model implemented, considers a horizontal and
vertical force balance. Vertical forces are a balance between
buoyancy, gravity and the vertical components of the lift and drag
forces. The horizontal force balance consists of the horizontal
components of the lift and drag forces only. These two equations can
be solved for the angle of attack and the speed through water,
determining the flight at any instance of time.

Input to the model comes from parameters measured by the glider, such
as the measured pitch angle (m_pitch), buoyancy change
(m_ballast_pumped or m_de_oil_vol) and the in-situ
density. Furthermore, the model requires the specification of a number
of coefficients:

* mg: mass of the glider (kg)
* Vg: volume of the glider (m³)
* Cd0: parasite drag coefficient
* epsilon: compressibility of the hull (1/Pa)
* ah: lift angle coefficient due to the hull (1/rad)
* Cd1: induced drag coefficient (1/rad²)

Using the depth-rate from the pressure sensor as only model
constraint, the mass (or glider volume) and the parasite drag
coefficient can be determined. To determine the lift angle coefficient
requires an additional constraint that contains a horizontal velocity
component. Details of this procedure are given in [merckelbach2019]_.

Dynamic model
-------------
In addition to a steady-state model, this code also implements a
dynamic model, that is, including the inertial terms. Since this model
needs to be integrated, for which the Runge-Kutta method is used, it
is more computational expensive. The dynamic model produces more
accurate results when forcing conditions change rapidly, such as when
crossing a sharp pycnocline or during the transition from dive to
climb. Apart from the mathematical model underlying, the interfaces to
both models are the same.

Model calibration and data masking
----------------------------------

To calibrate a model, either steady-state or dynamic, we may wish not
to include all the data in the evaluation of the cost-function. To
that end, data can be masked. The Calibrate class provides boolean
operators to do this:

* OR()
* AND()
* NAND()

By default a mask set to False for all data. To mask data for which a
condition evaluates to True, the OR() method should be used. For
example, ::

   gm = SteadyStateCalibrate(rho0=1024)
   gm.set_input_data(datadict)
   
   condition = depth<10
   gm.OR(condition)
   

which would exclude all data points for which the depth is less than
10 m from the evaluation of the cost-function.

A truth table:

+------+----------+----+-----+----+
| mask | conditon | OR | AND |NAND|
+------+----------+----+-----+----+
|  0   |    0     |  0 |  0  | 1  |
+------+----------+----+-----+----+
|  1   |    0     |  1 |  0  | 1  |
+------+----------+----+-----+----+
|  1   |    1     |  1 |  1  | 0  |
+------+----------+----+-----+----+
|  0   |    1     |  1 |  0  | 1  |
+------+----------+----+-----+----+


Example
-------

An example to calibrate a model::

   # create a dictionary with the data

   data = dict(time=t, pressure=P, pitch=pitch, buoyancy_change=deltaV)

   gm = SteadyStateCalibrate()
   # we have to define mass and volume at the minimum
   gm.define(mg=70, Vg=70)

   gm.set_input_data(data)

   # mask all data below 10 m
   gm.OR(pressure*10<10)
   # mask all data exceeding 60 m
   gm.OR(pressure*10>60)

   result = gm.calibrate("mg", "Cd0")
   
   print("Calibrated parameters:")
   for k,v in result.items():
       print("{}: {}".format(k,v)

   # Instead of printing the parameters from the results, we could also
   # get them from the corresponding attributes: print("Cd0:", gm.Cd0).

   print("Cd0:", gm.Cd0)

   # We also don't need to run the model again either. The model output
   # is also accessible from attributes:
   #
   # gm.t # time
   # gm.U # incident velocity
   # gm.alpha # angle of attack
   # gm.ug    # horizontal speed
   # gm.wg    # vertical speed
   # gm.w     # vertical water velocity
   
   # if we want to run a model with a given set of parameters

   fm = DynamicGLiderModel(dt=1, rho0=1024, k1=0.02, k2=0.92)
   # copy the settings from the steady state model
   fm.copy_settings(gm)

   solution = fm.solve(data)
   
   # solution is now a named tuple, according to the definition:
   # Modelresult = namedtuple("Modelresult", "t u w U alpha pitch ww")


How to cite
-----------
When you publish results that were obtained with this software, please use the
following citation:

|   Merckelbach, L., A. Berger, G. Krahmann, M. Dengler, and J. Carpenter, 2019: A
|            dynamic flight model for Slocum gliders and implications for turbulence
|            microstructure measurements. J. Atmos. Oceanic Technol., 36(2),
|            281-296, doi:10.1175/JTECH-D-18-0168.1.


Copyright information
---------------------
Copyright (c) 2018, 2019 Helmholtz Zentrum Geesthacht, Germany
                   Lucas Merckelbach, lucas.merckelbach@hzg.de

Software is licensed under the MIT licence.

References
----------
.. [merckelbach2019] Merckelbach, L., A. Berger, G. Krahmann, M. Dengler, and J. Carpenter, 2019: A
   dynamic flight model for Slocum gliders and implications for
   turbulence microstructure measurements. J. Atmos. Oceanic
   Technol. 36(2), 281-296, doi:10.1175/JTECH-D-18-0168.1

.. |PyPI version| image:: https://badgen.net/pypi/v/gliderflight
   :target: https://pypi.org/project/gliderflight
.. |Docs badge| image:: https://readthedocs.org/projects/dbdreader/badge/?version=latest
   :target: https://gliderflight.readthedocs.io/en/latest/
.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://github.com/smerckel/gliderflight/blob/master/LICENSE 

   
