#Copyright (c) 2018 Helmholtz Zentrum Geesthacht, Germany
#                   Lucas Merckelbach, lucas.merckelbach@hzg.de

import configparser
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot
from matplotlib.dates import epoch2num, DateFormatter

import gsw

import dbdreader
import gliderflight


DEFAULTS = dict(target_density=1025, mg=65, Vg=65e-3, minlimitdepth=10, maxlimitdepth=30,
                cond_a=1.0, cond_b=0., buoyancy_engine='shallow',latitude=54, longitude=8,
                calibrate_epsilon='no')
class UI(object):
    CONFIGFILE=os.path.join(os.environ['HOME'],'.glidertrimrc')    
    FLOAT_OPTIONS = "target_density mg Vg minlimitdepth maxlimitdepth cond_a cond_b latitude longitude".split()
    STR_OPTIONS = ["buoyancy_engine", "calibrate_epsilon"]
    UNITS = dict(mg="kg", Vg="m^3", minlimitdepth="m", maxlimitdepth="m", cond_a="m/S",
                 cond_b="-", buoyancy_engine="shallow|deep", latitude='decimal deg',
                 longitude='decimal deg', calibrate_epsilon='yes|no',
                 target_density='kg/m^3')
    
    def __init__(self):
        self.configparser=configparser.ConfigParser()
        self.settings = {}

    def get_glidername_and_filenames(self):
        if len(sys.argv)<3:
            print('''
    Provide:
            glider name or glider ID
            the path(s) of one or more dbd files. Make sure that the corresponding
                  ebd files are in the same directory.
    ''')
            sys.exit(1)
        glidername = sys.argv[1]
        fns=sys.argv[2:]
        for fn in fns:
            if os.path.exists(fn):
                print("%s found. Ok"%(fn))
            else:
                print("%s not found. Exiting."%(fn))
                sys.exit(2)
        fns.sort()
        return glidername, fns
    
    def read_config(self, glidername):
        for k,v in DEFAULTS.items():
            self.settings[k] = v
            
        self.configparser.read(UI.CONFIGFILE)
        if self.configparser.has_section(glidername):
            options=self.configparser.options(glidername)

            for option in options:
                if option in UI.FLOAT_OPTIONS:
                    self.settings[option] = self.configparser.getfloat(glidername,option)
                elif option in UI.STR_OPTIONS:
                    self.settings[option] = self.configparser.get(glidername,option)
    
    def write_config(self, glidername):
        if not self.configparser.has_section(glidername):
            self.configparser.add_section(glidername)
        for k, v in self.settings.items():
            if k in UI.FLOAT_OPTIONS:
                self.configparser.set(glidername, k, v.__str__())
            else:
                self.configparser.set(glidername, k, v)
        with open(UI.CONFIGFILE,'w') as fp:
            self.configparser.write(fp)

    def configure(self):
        print()
        for k,v in self.settings.items():
            if k in UI.FLOAT_OPTIONS:
                s="%s (%s)"%(k, UI.UNITS[k])
                mesg = "Enter value for %-30s (current: %f): "%(s,v)
                ans = input(mesg)
                if ans == '':
                    continue
                else:
                    self.settings[k]=float(ans)
            else:
                if k == 'calibrate_epsilon' and self.settings['buoyancy_engine']=='shallow':
                    self.settings[k]='no'
                    continue
                s="%s (%s)"%(k, UI.UNITS[k])
                mesg = "Enter value for %-30s (current: %s): "%(s,v)
                ans = input(mesg)
                if ans == '':
                    continue
                else:
                    self.settings[k] = ans

class Model(object):
    def __init__(self, fns, settings):
        self.fns = fns
        self.settings = settings
        self.model = gliderflight.SteadyStateCalibrate()

    def read_dbds(self, fns, settings):
        if settings['buoyancy_engine']=='shallow':
            buoyancy_variable='m_ballast_pumped'
        elif settings['buoyancy_engine']=='deep':
            buoyancy_variable='m_de_oil_vol'
        else:
            raise ValueError('Unknown buoyancy engine. Accepted values: (shallow|deep)')
        dbd = dbdreader.MultiDBD(filenames=fns, include_paired=True)
        tmp = dbd.get_sync("sci_ctd41cp_timestamp",["sci_water_pressure",
                                                    "sci_water_cond",
                                                    "sci_water_temp",
                                                    "m_pitch",
                                                    "m_battpos",
                                                    buoyancy_variable])
        _, tctd, P, C, T, pitch, battpos, buoyancy_change = tmp.compress(tmp[3]>0.01, axis=1)
        SP = gsw.SP_from_C(C*10, T, P*10)
        SA = gsw.SA_from_SP(SP, P*10, settings['longitude'], settings['latitude'])
        density = gsw.rho_t_exact(SA, T, P*10)
        
        #density = fast_gsw.rho(C*10, T, P*10, settings['longitude'], settings['latitude'])

        data = dict(time=tctd, pressure=P, pitch=pitch,
                    buoyancy_change=buoyancy_change, battpos=battpos, density=density)
        self.model.input_data = data
        mask=np.zeros(tctd.shape,'int').astype(bool) # use all data by default
        self.model.set_mask(mask)
        self.model.OR(P*10<settings['minlimitdepth'])
        self.model.OR(P*10>settings['maxlimitdepth'])
        
        
    def calibrate(self, settings):
        self.model.define(ah=3.8, Cd1=10.5, Cd0=0.15, mg=settings['mg'], Vg=settings['Vg'])
        cal_params = "Cd0 Vg".split()
        if settings['calibrate_epsilon']=='yes':
            cal_params.append('epsilon')
        results = self.model.calibrate(*cal_params, verbose=True)
        return results

    def estimate_pitch_relation(self,settings):
        b=self.model.input_data['buoyancy_change']*1e-6
        p=self.model.input_data["pitch"]
        bp=self.model.input_data["battpos"]*2.56e-2
        P=self.model.input_data['pressure']*1e-3
        if settings['calibrate_epsilon']=='yes':
            x=np.vstack([b,bp,np.ones(bp.shape),P])
        else:
            x=np.vstack([b,bp,np.ones(bp.shape)])
        condition=np.logical_and(abs(p)>11*np.pi/180.,np.isfinite(b))
        x = x.compress(condition,axis=1)
        p = p.compress(condition)
        coefs = np.linalg.lstsq(x.T,np.tan(p), rcond=None)[0]
        return coefs

    def report(self, results, settings):
        print()
        print()
        print("Drag coefficient Cd    : %f (-)"%(self.model.Cd0))
        print("Glider volume Vg       : %f (m^3)"%(self.model.Vg))
        print("Glider compressibility : %f (*e-10)"%(self.model.epsilon))
        print("Glider density         : %f (kg/m^3)"%((self.model.mg/self.model.Vg)))
        print("Weight change          : %f (g)"%((settings['target_density']-(self.model.mg/self.model.Vg))*self.model.Vg*1000.))
        #
        print()
        print("Estimated pitch relationship:")
        coefs = self.estimate_pitch_relation(settings)
        print("tan(pitch) = T1 * buoyancy(m^3) + T2 * battpos(m) + T3 (tan(pitch0)) + T4 P (kbar)")
        print("T1 :",coefs[0])
        print("T2 :",coefs[1])
        print("T3 :",coefs[2])
        if len(coefs)==4:
            print("T4 :",coefs[3])
        else:
            print("T4 : 0")

    def graphic_report(self, results, settings):
        f, ax = pyplot.subplots(1,2,sharey=True)
        mr = self.model.modelresult
        z = self.model.input_data['pressure']*10
        P = self.model.input_data['pressure']
        zi = np.linspace(z.min(), z.max())
        Pi = zi/10.
        
        ax[0].plot(mr.ww, z, label='w_water', color='C0', alpha=0.5)
        ax[0].plot(np.convolve(mr.ww, np.ones(15)/15, 'same'),
                   z, label='w_water', color='C5')
        ax[0].plot(mr.w, z, label='w_glider', color='C1')
        ax[0].plot(mr.U, z, label='U_glider', color='C3')
        ax[1].plot(zi*0+settings['target_density'], zi, label='Target density')
        ax[1].plot(self.model.input_data['density'], z, label='Density')
        rho_g = self.model.mg/(self.model.Vg*(1-self.model.epsilon*1e-10*Pi*1e5))
        ax[1].plot(rho_g, zi, label='Glider density', color='C2')
        for i in range(1, 5):
            rho_g = (50e-3*i+self.model.mg)/(self.model.Vg*(1-self.model.epsilon*1e-10*Pi*1e5))
            ax[1].plot(rho_g, zi, color='C2', ls='--', alpha=0.5)
            rho_g = (-50e-3*i+self.model.mg)/(self.model.Vg*(1-self.model.epsilon*1e-10*Pi*1e5))
            ax[1].plot(rho_g, zi, color='C2', ls='--', alpha=0.5)
        ax[0].set_ylabel('Depth (m)')
        ax[0].set_xlabel('Velocity (m/s)')
        ax[1].set_xlabel('Density (kg/m^3)')
        ax[0].legend()
        ax[1].legend()
        ax[0].invert_yaxis()
        f.show()

def main():
    ui = UI()
    glidername, fns = ui.get_glidername_and_filenames()
    ui.read_config(glidername)
    ui.configure()
    ui.write_config(glidername)

    model = Model(fns, ui.settings)
    model.read_dbds(fns, ui.settings)
    results = model.calibrate(ui.settings)
    model.report(results, ui.settings)
    model.graphic_report(results, ui.settings)
    input("Press enter to exit")
