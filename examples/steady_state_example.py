# part of glider flight model for Slocum ocean gliders
#
# Licensed under MIT license
#
# Author: Lucas Merckelbach
#         lucas.merckelbach@hzg.de
#
# Sept. 2018
#
# A python script showing the use of the gliderflight model.
#
# It calibrates and runs a steady-state model and then the same for
# the dynamic model
#
# This script uses a subset data that were also used in the JTECH paper.

import numpy as np
import sys
sys.path.insert(0,'..')
import gliderflight


# load some measurements. Data extracted from a glider mission and processed DVL data.
tctd, P, pitch, buoyancy_drive, density, U_adcp = np.load("../data/gliderflight_data.npy")

# Setup a steady-state model to calibrate
GM = gliderflight.SteadyStateCalibrate(rho0=1004)
GM.define(ah=3.8, Cd1=10.5)
GM.set_input_data(time=tctd, pressure=P, pitch=pitch, buoyancy_change=buoyancy_drive, density=density)
GM.OR(P*10<10)
GM.OR(P*10>35)
GM.define(Cd0=0.15, mg=70, Vg=70/1004)
calibration_result = GM.calibrate("Cd0", "mg", verbose=True)
