Glidertrim
==========

Glidertrim is a simple commandline utility, that, given one or more
pairs of dbd/ebd files, and a target density, computes the optimum
weight change. The typical application is that during the deployment
a test dive is made, and the resulting dbd/edb files are analysed
quickly, producing an objective estimate of how much a glider is
overweight or underweight, given a target density.

Synopsis
--------

**glidertrim <glidername> <dbd file> [dbd file] [dbd file]**

glidername is the name or identifier of the glider. It is used to
write settings in the configuration file ($HOME/.glidertrimrc), so
that the settings can be prepared prior to deployment and recorded for
later use.

dbd_file is a path to a dbd file and can included wildcards. It is
important that the dbd filename is provided only, but the matching ebd
file is present in the same directory the dbd file resides in.

Description
-----------
Glidertrim takes a glider id and one or more dbd (with matching ebd)
files as input, and given some configuration settings which the user
can change, calculates the ideal weight change.

An example is given below::

  $ glidertrim comet comet-2018-136-00-000.dbd
  comet-2018-136-00-000.dbd found. Ok

  Enter value for target_density (kg/m^3)        (current: 1026.000000): 
  Enter value for mg (kg)                        (current: 69.500000): 
  Enter value for Vg (m^3)                       (current: 0.065000): 
  Enter value for minlimitdepth (m)              (current: 3.000000): 
  Enter value for maxlimitdepth (m)              (current: 55.000000): 
  Enter value for cond_a (m/S)                   (current: 1.000000): 
  Enter value for cond_b (-)                     (current: 0.000000): 
  Enter value for buoyancy_engine (shallow|deep) (current: deep): 
  Enter value for latitude (decimal deg)         (current: 54.000000): 
  Enter value for longitude (decimal deg)        (current: 8.000000): 
  Enter value for calibrate_epsilon (yes|no)     (current: no): 
  Error: 1.1079233e-01  -  Cd0=0.1500  Vg=0.0650
  Error: 1.0586534e-01  -  Cd0=0.1575  Vg=0.0650
  :
  :
  Error: 1.7395719e-03  -  Cd0=0.3100  Vg=0.0678
  Error: 1.7395722e-03  -  Cd0=0.3100  Vg=0.0678
  Error: 1.7395718e-03  -  Cd0=0.3099  Vg=0.0678

  
  Drag coefficient Cd    : 0.309885 (-)
  Glider volume Vg       : 0.067768 (m^3)
  Glider compressibility : 5.000000 (*e-10)
  Glider density         : 1025.555502 (kg/m^3)
  Weight change          : 30.122814 (g)

  Estimated pitch relationship:
  tan(pitch) = T1 * buoyancy(m^3) + T2 * battpos(m) + T3 (tan(pitch0)) + T4 P (kbar)
  T1 : 1784.9588582927518
  T2 : -23.879872610869665
  T3 : 0.09352561010068428
  T4 : 0
  Press enter to exit


A table with configurable parameters is shown below:

+------------------------------+--------------------------------------------------+
|Parameter (unit)              | Description                                      |
+==============================+==================================================+
|target_density (kg/m^3)       | target density                                   |
+------------------------------+--------------------------------------------------+
|mg (kg)                       | (measured) mass of glider                        |
+------------------------------+--------------------------------------------------+
|Vg (m^3)                      | estimate of glider volume                        |
+------------------------------+--------------------------------------------------+
|minlimitdepth (m)             | minimum depth allowed in optimisation routine    |
+------------------------------+--------------------------------------------------+
|maxlimitdepth (m)             | maximum depth allowed in optimisation routine    |
+------------------------------+--------------------------------------------------+
|cond_a (m/S)                  | scaling factor for correcting conductivity       |
+------------------------------+--------------------------------------------------+
|cond_b (-)                    | offset for correcting conductivity               |
+------------------------------+--------------------------------------------------+
|buoyancy_engine (shallow|deep)| sets buoyancy engine used, deep or shallow       |
+------------------------------+--------------------------------------------------+
|latitude (decimal deg)        | latitude of experiment  (used for density)       |
+------------------------------+--------------------------------------------------+
|longitude (decimal deg)       | longitude of experiment (idem)                   |
+------------------------------+--------------------------------------------------+
|calibrate_epsilon (yes|no)    | whether or not to calibrate for compressibility  |
+------------------------------+--------------------------------------------------+

The user is presented with default values and has the option to change
the value or simply press enter, which retains the defailt
value. After entering the last configuration parameter, the utility
optimises for the parasite drag coefficient Cd0, the glider volume Vg,
and, if calibrate_epsilon is set to "yes", the compressibility.

After the optimised values are found, they are displayed. Among the
results returned are the glider density and the required weight change
to match the glider's density to the target density.

The results are also shown graphically. The left panel shows the
vertical profiles of water velocity (raw and filtered), the glider
vertical velocity and the glider speed through water. The right panel
shows the in-situ density profile, the target density, and the actual
glider density (accounting for the compressibility). The dashed lines
refer to the glider's density with increments of 50 g of weight
change.


Estimated pitch relationship
----------------------------

The pitch the glider assumes during diving and climbing depends on the
total torque exerted on the glider. Components influening the torque
balance are the mass, the so-called h-moment (distance between the
centres of buoyancy and gravity), the buoyancy drive, pitch battery
position and pressure.  Based on a linear regression model the
contributions to the pitch by the buoyancy change, battery position,
and pressure are estimated. The intended purpose is for glider flight
simulations to use the correct pitch, when pitch is not set directly,
but through a fixed battery position, for example.

The relationship for the pitch is given by

.. math::
   \tan(\phi) = T_1 \cdot V_b + T_2 \cdot b_p + T_3 + T_4 \cdot P,

where

:math:`V_b` is the buoyancy change in m\ :sup:`3`

:math:`b_p` is the battery position in m, and

:math:`P` is the pressure in kbar.
