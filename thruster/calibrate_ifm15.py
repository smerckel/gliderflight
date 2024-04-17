import sys
sys.path.insert(0, '..')

import matplotlib.pyplot as plt

import dbdreader
import gliderflight
import fast_gsw

path = "/home/lucas/gliderdata/caboverde_201911/hd/ifm15-2019-333-04-00?.dbd"
dbd = dbdreader.MultiDBD(path, complement_files=True)

tctd, C, T, P, pitch, buoyancy_drive, thruster_setting = dbd.get_CTD_sync("m_pitch", "m_de_oil_vol", "c_thruster_on")
insitu_density = fast_gsw.rho(C*10, T, P*10, -20.6, 17.74)

gm = gliderflight.SteadyStateCalibrate(rho0=1025)
gm.define(Cd0=0.1, mg=65, Vg=65e-3)
gm.set_input_data(tctd, P, pitch, buoyancy_drive, insitu_density)
gm.define(**gm.calibrate('Cd0', 'Vg'))
mr = gm.solve()



# We have some files where the thruster value was increased step wise. Select them here.
path = "/home/lucas/gliderdata/caboverde_201911/hd/ifm15-2019-338-00-*.dbd"
ps = dbdreader.DBDPatternSelect(date_format="%d %m %Y %H:%M")
fns = ps.select(pattern=path, from_date="5 12 2019 13:07", until_date="5 12 2019 17:32")

dbd_thruster = dbdreader.MultiDBD(fns, complement_files=True)
tctd, C, T, P, pitch, buoyancy_drive, thruster_setting = dbd_thruster.get_CTD_sync("m_pitch", "m_de_oil_vol", "m_thruster_power")
insitu_density = fast_gsw.rho(C*10, T, P*10, -20.6, 17.74)

from gliderflight import SteadyStateCalibrate
from scipy.optimize import fsolve as so_fsolve

class SteadyStateCalibrateThruster(SteadyStateCalibrate):
    def __init__(self, rho0=None):
        super().__init__(rho0=rho0)

    def compute_Ft(self):
        data = self.input_data
        pitch = data['pitch']
        rho = data['density']
        try:
            dhdt = data['dhdt']
        except KeyError:
            dhdt = self.compute_dhdt(data['time'], data['pressure'])
            data['dhdt'] = dhdt
            
        pressure, Vbp = self.convert_pressure_Vbp_to_SI(data['pressure'],
                                                        data['buoyancy_change'])
        FB, Fg = self.compute_FB_and_Fg(pressure, rho, Vbp)
        return FB-Fg

gmt = SteadyStateCalibrateThruster(rho0=1025)
gmt.copy_settings(gm)
gmt.set_input_data(tctd, P, pitch, buoyancy_drive, insitu_density)
mrt = gmt.solve()

FB = gmt.compute_Ft()

def equations(p, pitch, w, FB, gmt):
    alpha, U, Ft = p
    q =  0.5 * gmt.RHO0 * U**2 * gmt.S
    FL = q * alpha * (gmt.aw + gmt.ah)
    FD = q * (gmt.Cd0 + gmt.Cd1*alpha**2)
    eq0 = np.sin(pitch + alpha) * U - w
    eq1 = FB + np.sin(pitch) * Ft - np.sin(pitch + alpha) * FD - np.cos(pitch + alpha) * FL
    eq2 =      np.cos(pitch) * Ft - np.cos(pitch + alpha) * FD + np.sin(pitch + alpha) * FL
    return (eq0, eq1, eq2)


alpha = []
U = []
Ft = []
for i, (_pitch, _w, _Fb) in enumerate(zip(pitch, gmt.input_data['dhdt'], FB)):
    r = so_fsolve(equations, (_pitch/10, 0.2, 1), args=(_pitch, _w, _Fb, gmt))
    alpha.append(r[0])
    U.append(r[1])
    Ft.append(r[2])
Ft= np.array(Ft)
alpha= np.array(alpha)
U= np.array(U)

# make a histogram and find the four spikes
n, xb = np.histogram(thruster_setting, 40)
n_limit = 300
idx = np.where(n>n_limit)[0]
bl = xb[idx]
br = xb[idx+1]

Ft_mean = []
thruster_power = []
for _bl, _br in zip(bl, br):
    condition = np.logical_and(thruster_setting>=_bl,
                               thruster_setting<_br)
    Ft_mean.append(Ft.compress(condition).mean())
    thruster_power.append(thruster_setting.compress(condition).mean())
Ft_mean = np.array(Ft_mean)
thruster_power= np.array(thruster_power)


fig, ax = plt.subplots(3,1)
ax[0].plot(tctd-tctd[0], P*10)
ax[0].set_xlabel("time (s)")
ax[0].set_ylabel("depth (m)")

ax[1].plot(tctd-tctd[0], thruster_setting)
ax[1].set_xlabel("time (s)")
ax[1].set_ylabel("m_thruster_power (W)")

ax[2].plot(thruster_setting, Ft, 'o',alpha=0.1)
ax[2].plot(thruster_power, Ft_mean, 'ro')
ax[2].set_xlabel("m_thruster_power (W)")
ax[2].set_ylabel("Fthruster (N)")
ax[0].set_xlim(0, 16200)
ax[1].set_xlim(0, 16200)
ax[0].set_ylim(220, 0)
ax[1].set_ylim(0, 5.5)
ax[2].set_xlim(0, 5.5)
ax[2].set_ylim(0, 10)

x = np.log(thruster_power)
y = np.log(Ft_mean)
a, b = np.polyfit(x, y, 1)
tp = np.linspace(0, 6)
ax[2].plot(tp, np.exp(b) * tp**a, '-')
