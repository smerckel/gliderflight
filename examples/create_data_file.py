# script used to create the data set gliderflight_data.npy
import numpy as np

import dbdreader
try:
    # use my fast private gsw implementation if available
    from fast_gsw import rho as gsw_rho
except ImportError:
    # or fall back on the community one, installable via pip, for example.
    import gsw
    
    def gsw_rho(C, T, P, lon, lat):
        SP = gsw.SP_from_C(C, T, P)
        # This particular data set was collected in the Baltic.
        SA = gsw.SA_from_SP_Baltic(SP, lon, lat)
        # in-situ density
        rho = gsw.density.rho_t_exact(SA, T, P)
        return rho
    
dbd = dbdreader.MultiDBD(pattern="/home/lucas/gliderdata/subex2016/hd/comet-2016-174-05-00?.[de]bd")

_tmp = dbd.get_sync("sci_ctd41cp_timestamp", "sci_water_pressure m_pitch m_ballast_pumped sci_water_cond sci_water_temp".split())
_, tctd, P, pitch, buoyancy_drive, C, T = _tmp.compress(_tmp[1]>0, axis=1)

tilt_correction_factor=0.86
pitch_offset=-0.004,
pitch = pitch*tilt_correction_factor + pitch_offset
density = gsw_rho(C*10, T, P*10, 15, 54)
U_adcp = np.interp(tctd, *np.load("../data/U_kalman_datasetIII.npy"))
U_adcp[np.where(np.logical_or(P*10<12, pitch>0))[0]] = np.nan
np.save("../data/gliderflight_data.npy", (tctd, P, pitch, buoyancy_drive, density, U_adcp))

