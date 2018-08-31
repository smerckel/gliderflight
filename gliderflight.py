from collections import namedtuple

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve, fmin as fminimize


Modelresult = namedtuple("Modelresult", "t u w U alpha pitch ww")


class ModelParameterError(BaseException):
    pass

class ModelParameters(object):
    '''Configuration class for glider model parameters


    This class defines the configuration parameters of the glider
    model. The class is meant to be subclassed from model
    implementations.

    
    Methods defined in this class:

    * define(): define or set a parameter
    * show_settings(): prints the current settings
    * copy_settings(): updates model parameter settings from another model
    * undefined_parameters(): returns which parameters have not been set yet.
    * cd1_estimate(): estimates the induced drag coefficient
    * aw_estimate: estimates the lift coefficient due to the wings using a parameterisation

    '''
    

    def __init__(self, parameterised_parameters_dict):
        ''' Constructor

        :param parameterised_parameters_dict: a dictionary of parameters that should be computed rather then being set explicitly
        :type parameterised_parameters_dict: dictionary

        '''
        
        self.parameters = 'Cd0 Cd1 S eOsborne AR Omega mg epsilon Vg Cd1_hull aw ah'.split()
        self.parametersAoa = 'Cd0 Cd1 eOsborne AR Omega S Cd1_hull ah aw'.split()
        self.parameterised_parameters = parameterised_parameters_dict
        self._parameterised_parameters = list(self.parameterised_parameters)
       
        for i in self.parameters:
            self.__dict__[i]=None
        self.__aoa_parameter_changed=True
        #set some default values for "constant" parameters
        self.define(AR=7,Omega=43*np.pi/180.,S=0.1,ah=3.8,eOsborne=0.8,Cd1_hull=9.7,epsilon=5)
        # note epsilon is given without e-10 factor.
        self.modelresult = None
        
    def show_settings(self):
        '''Print model parameters '''
        
        p = list(self.parameters)
        p.sort()
        for _p in p:
            if self.__dict__[_p]:
                print("{:20s} : {:.3g}".format(_p, self.__dict__[_p]))
                
    def copy_settings(self, other):
        ''' Copy model parameters 

        :param other: an other instance of this class (or subclassing this class)
        :type other: ModelParameters

        :example:

        dynamic_model.copy_settings(steady_state_model)

        '''
        p = list(self.parameters)
        for _p in p:
            self.__dict__[_p] = other.__dict__[_p]

    def get_settings(self):
        ''' Get model settings 

        :return: a dictionary with the current parameter setting
        :rtype: dict

        '''
        settings = dict((p, self.__dict__[p]) for p in list(self.parameters))
        return settings
    
    def define(self,**kw):
        ''' Define or set one or more glider configuration parameters.

        :param **kw: keywords with parameter name and values
        :type **kw: keyword assignment

        :example:

           glidermodel.define(Cd0=0.24)
        
           glidermodel.define(Vg=50e-3, mg=60)

        '''
        for k,v in kw.items():
            self.__dict__[k] = v
            if v is None: # None values to set parameter to be computed, if computable.
                if k not in self._parameterised_parameters:
                    self._parameterised_parameters.append(k)
                else:
                    raise ValueError("Parameter requires a finite value.")
            else:
                # remove from parameterised parameter list if not None and computable
                if k in self._parameterised_parameters:
                    self._parameterised_parameters.remove(k)
            if k in self.parametersAoa:
                self.__aoa_parameter_changed=True
        # set parameterised values or update them:        
        for k,v in self.parameterised_parameters.items():
            if k in self._parameterised_parameters:
                try:
                    self.__dict__[k]=v()
                except ModelParameterError:
                    pass
                
    def has_aoa_parameter_changed(self):
        ''' test whether any of the parameters that appear in the angle of attack estimate have changed

        :return: test result
        :rtype: bool

        '''
        rv=self.__aoa_parameter_changed
        return rv

    def undefined_parameters(self):
        ''' Returns undefined parameters

        :return: list of undefined parametrs
        :rtype: list
        '''
        return [i for i in self.parameters if self.__dict__[i]==None and i not in self._parameterised_parameters]

    def cd1Estimate(self):
        ''' Parameterisation for Cd1 

        :return: Value for Cd1
        :rtype: float
        '''
        # check whether all parameters are defined:
        for i in 'aw eOsborne AR Cd1_hull'.split():
            if self.__dict__[i]==None:
                raise ModelParameterError()
        Cd1_w = self.aw**2/np.pi/self.eOsborne/self.AR
        Cd1_hull = self.Cd1_hull
        return Cd1_w + Cd1_hull

    def awEstimate(self):
        ''' Parameterisation for aw 

        :return: Value for aw
        :rtype: float
        '''

        # check whether all parameters are defined:
        for i in 'AR Omega'.split():
            if self.__dict__[i]==None:
                raise ModelParameterError()
        return 2*np.pi*self.AR/(2+np.sqrt(self.AR**2*(1+np.tan(self.Omega)**2)+4))



    
class GliderModel(object):
    ''' Common glider model class

    This class, meant to be subclassed, implements the physical glider model description

    '''
    RHO0=1024
    G = 9.81

    def __init__(self, rho0=None):
        self.RHO0 = rho0 or GliderModel.RHO0


    def convert_pressure_Vbp_to_SI(self, m_water_pressure, m_de_oil_vol):
        ''' converts units of glider sensor data into SI data 

        :param m_water_pressure: water pressure in bar
        :type m_water_pressure: array or float
        :param m_de_oil_vol: buoyancy change reported by glider in cc
        :type m_de_oil_vol: array or float
        :return: converted values
        :rtype: tuple of numpy arrays
        '''
        pressure=m_water_pressure*1e5 # Pa
        Vbp=m_de_oil_vol*1e-6 # m^3
        return pressure, Vbp

    def compute_FB_and_Fg(self, pressure, rho, Vbp):
        ''' Computes the vertical forces FB and Fg

        :param pressure: pressure (Pa)
        :param rho: in-situ density (kg m$^{-3}$)
        :param Vbp: volume of buoyancy change (m$^{-3}$)
        :type pressure: numpy array or float
        :type rho: numpy array or float
        :type Vbp: numpy array or float
        :return: Computed buoyancy and gravity forces
        :rtype: tuple length 2 of numpy arrays or floats
        '''
        g=self.G
        #
        FB=g*rho*(self.Vg*(1.-self.epsilon*1e-10*pressure) + Vbp)
        Fg=self.mg*g
        return FB, Fg
        

    def stall_factor(self, alpha, **kwds):
        return 1.

    def compute_dhdt(self, time, pressure):
        ''' Compute the depth rate from the pressure

        :param time: time (s)
        :param pressure: pressure (Pa)
        :type time: numpy array or float
        :type pressure: numpy array of float
        :return: depth rate
        :rtype: numpy array or float

        .. note::
            The density used to convert pressure into depth is given by self.RHO0
        '''
        return -np.gradient(pressure*1e5)/np.gradient(time)/self.RHO0/self.G
    
    def compute_lift_and_drag(self,alpha, U, rho):
        ''' Compute lift and drag forces

        :param alpha: angle of attack (rad)
        :param U: incident water velocity
        :param rho: in-situ density
        :type alpha: numpy array or float
        :type U: numpy array or float
        :type rho: numpy array or float
        :return: dynamic pressure, lift force and drag force
        :rtype: tuple of length 3 with numpy arrays or floats
        '''
        q = 0.5 * rho * self.S * U**2
        L = q * (self.aw + self.ah)*alpha*self.stall_factor(alpha)
        D = q * (self.Cd0 + self.Cd1*alpha**2)
        return q, L, D

    # providing easy access to useful computational results as attributes
    def __get_modelresult(self, p):
        if self.modelresult:
            i = self.modelresult._fields.index(p)
            return self.modelresult[i]

    @property        
    def t(self):
        ''' time (s)'''
        return self.__get_modelresult('t')
    @property        
    def U(self):
        ''' incident water velocity (m s$^{-1}$) '''
        return self.__get_modelresult('U')
    @property        
    def alpha(self):
        '''angle of attack (rad)'''
        return self.__get_modelresult('alpha')
    @property        
    def pitch(self):
        '''pitch angle (rad) '''
        return self.__get_modelresult('pitch')
    @property
    def wg(self):
        '''vertical velocity of glider relative to surface (m s$^{-1}$)'''
        return self.__get_modelresult('w')
    @property        
    def w(self):
        '''vertical water velocity (m s$^{-1}$)'''
        return self.__get_modelresult('ww')

    
class SteadyStateGliderModel(ModelParameters, GliderModel):
    ''' Steady-state implementation 

    This class inherits from ModelParameters and GliderModel. The physcis are
    provided by GliderModel. Interacting with ModelParameters is done through methods
    provided by ModelParameters.

    The only method provided by this class that is of interest to the user is solve().
    The input to solve is a dictionary with time, pressure, pitch, buoyancy change density.

    Methods inherited from ModelParameters can be used to define/set model coefficients, and 
    to copy settings from another model instance.
    
    After solving the model results are available as properties (t, U, wg, w, alpha)

    :example:
    
    >>>gm = SteadyStateGliderModel(rho0=1024)
    >>>gm.define(mg=70)
    >>>gm.define(Vg=68, Cd0=0.15)
    >>>gm.solve(dict(time=tctd, pressure=P, pitch=pitch, buoyancy_change=buoyancy_drive, density=density))
    >>>print(gm.U)

    '''

    
    def __init__(self, rho0=None):
        '''Constructor
        
        :param rho0: background in-situ density
        :type rho0: float
        '''
        ModelParameters.__init__(self, dict(aw=self.awEstimate, Cd1=self.cd1Estimate))
        GliderModel.__init__(self, rho0=rho0)
        self.pitch_i = np.linspace(0.05*np.pi/180, 60*np.pi/180,100)
        
    def reset(self):
        ''' Resets angle of attack interpolation function '''
        self.ifun = None

    def model_fun(self, x,  m_pitch):
        ''' implicit function of the angle of attack '''
        alpha = x
        equation = (self.Cd0 +self.Cd1*alpha**2)/((self.ah+self.aw)*np.tan(alpha+m_pitch))-alpha
        return equation

        
    def solve_for_angle_of_attack(self, pitch):
        ''' Solves for the angle of attack

        :param pitch: pitch
        :type pitch: float or numpy array
        :return: angle of attack
        :rtype: float or numpy array
        
        .. note::
           This method uses an interpolating function. If any parameter on which this calculation 
           depends, changes, the interpolating function is recomputed. Whether any of these parameters
           is changed, is tracked by the define() method. 
        '''
        
        if self.has_aoa_parameter_changed() or self.ifun is None:
            aoas = np.zeros_like(self.pitch_i)
            for i, _pitch in enumerate(self.pitch_i):
                tmp = fsolve(self.model_fun,
                             _pitch/50,
                             args=(_pitch,),
                             full_output=1)
                num_result, result, err, mesg = tmp
                aoas[i] = num_result
            self.ifun = interp1d(self.pitch_i, aoas)
        #
        idx = np.where(np.logical_and(np.abs(pitch)>=self.pitch_i.min(),
                                      np.abs(pitch)<=self.pitch_i.max()))[0]
        aoa = np.zeros_like(pitch)
        aoa[idx] = self.ifun(abs(pitch[idx]))*np.sign(pitch[idx])
        return aoa
        
    def solve_model(self, rho, FB, pitch, Fg):
        ''' Solves first for angle of attack and then incident velocity'''
        alpha = self.solve_for_angle_of_attack(pitch)
        
        q=(FB-Fg)*np.sin(pitch+alpha)/(self.Cd0+self.Cd1*alpha**2)
        Usquared=2.*q/rho/self.S
        jdx = np.where(Usquared<0)[0]
        Usquared[jdx]=0
        U = Usquared**0.5
        return alpha, U

    def solve(self, data):
        ''' Solve the model

        :param data: environment data
        :type data: dict
        :return: model result
        :rytpe: Modelresult (named tuple with t, u, w, U, alpha,  pitch, ww fields)

        The data supplied should contain at least time, pressure, pitch, buoyancy_change and
        density, as reported by the glider. Depth rate (dhdt) will be added if not already present.
        Other data are ignored.

        :example:
        
        >>> gm = SteadyStateGliderModel()
        >>> gm.define(mg=70, Vg=68)
        >>> data = dict(time=time, pressure=P, pitch=pitch, buoyancy_change=vb, density=rho)
        >>> gm.solve(data)
        >>> plot(gm.U)
        '''
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
        alpha, U = self.solve_model(rho, FB, pitch, Fg)
        wg = np.sin(pitch+alpha)*U
        ug = np.cos(pitch+alpha)*U
        ww = dhdt - wg
        q, L,D = self.compute_lift_and_drag(alpha, U, rho)
        self.modelresult = Modelresult(data["time"], ug, wg, U, alpha, pitch, ww)  
        return self.modelresult
    

class DynamicGliderModel(ModelParameters, GliderModel):
    ''' Dynamic glider model implementation

    This class inherits from ModelParameters and GliderModel. The physcis are
    provided by GliderModel. Interacting with ModelParameters is done through methods
    provided by ModelParameters.

    The only method provided by this class that is of interest to the user is solve().
    The input to solve is a dictionary with time, pressure, pitch, buoyancy change density.

    Methods inherited from ModelParameters can be used to define/set model coefficients, and 
    to copy settings from another model instance.
    
    After solving the model results are available as properties (t, U, wg, w, alpha)

    The dynamic model solves the force balances including the intertial forces by numerical integration
    using a Runge-Kutta scheme. The inertial terms include the added mass terms. The relevant parameters
    can be set when creating an instance of this class.

    **Added mass**

    Added mass terms are specified by the coefficients k1 and k2, which refer to the added mass terms
    along the principle glider axis (k1) and vertically perpendicular (k2), where k1 and k2 are given
    as fraction of the glider mass mg. 
    :example:
    
    >>>dm = DynamicGliderModel(rho0=1024, k1=0.2, k2=0.92, mg=70)
    >>>dm.define(mg=70)
    >>>dm.define(Vg=68, Cd0=0.15)
    >>>dm.solve(dict(time=tctd, pressure=P, pitch=pitch, buoyancy_change=buoyancy_drive, density=density))
    >>>print(dm.U)

    '''

    def __init__(self, dt=None, rho0=None, k1=0.20, k2=0.92, alpha_linear=90, alpha_stall=90,
                 max_depth_considered_surface=0.5):
        ''' Constructor

        :param dt: time step (s)
        :param rho0: background density (kg m$^{-3}$)
        :param k1: added mass fraction in longitudinal direction
        :param k2: added mass fraction perpendicular to longitudinal direction
        :param alpha_linear: angle in degree, up to where lift force is linear with angle of attack
        :param alpha_stall: angle in degree where stall occurs
        :param max_depth_considered_surface: depths shallower than this are considered at the surface and (u,v)=0
        :type dt: float
        :type rho0: float
        :type k1: float
        :type k2: float
        :type alpha_linear: float
        :type alpha_stall: float
        :type max_depth_considered_surface: float

        '''
        ModelParameters.__init__(self, dict(aw=self.awEstimate, Cd1=self.cd1Estimate))
        GliderModel.__init__(self, rho0=rho0)
        self.k1 = k1
        self.k2 = k2 
        self.alpha_linear = alpha_linear*np.pi/180 # 90 degrees : disable
        self.alpha_stall = 90*np.pi/180 # 90 degrees: disable
        self.dt = dt
        self.max_depth_considered_surface = max_depth_considered_surface
        
    def stall_factor(self, alpha):
        alpha_abs = abs(alpha)
        if alpha_abs<self.alpha_linear:
            r = 1.
        else:
            d_alpha = alpha_abs - self.alpha_linear
            tau_alpha = self.alpha_stall - self.alpha_linear
            if tau_alpha<1e-9:
                r = 1
            else:
                r = np.exp(-d_alpha/tau_alpha)
        return r
    
    def __compute_k(self, _u, _w, _rho, _pitch, _FBg, _m11, _m12, _m21, _m22, h):
        U = (_u**2 + _w**2)**(0.5)
        alpha = np.arctan2(_w, _u) - _pitch
        q, L, D = self.compute_lift_and_drag(alpha, U, _rho)
        Fx = -np.cos(_pitch + alpha) * D + np.sin(_pitch + alpha)*L
        Fy = -np.cos(_pitch + alpha) * L - np.sin(_pitch + alpha)*D + _FBg
        k_u = h * (_m11*Fx + _m12*Fy)
        k_w = h * (_m21*Fx + _m22*Fy)
        return k_u, k_w
    
    def RK4(self, h, M, FBg, pitch, rho, at_surface, u, w):
        ''' Runge-Kutta integration method '''
        N = u.shape[0]
        data = np.array([FBg, pitch, rho, at_surface, *M]).T
        # compute input data on mid points h/2
        data_mid = (data[:-1] + data[1:])*0.5
        
        for i in range(N-1):
            _u, _w = u[i], w[i]
            _FBg, _pitch, _rho, _at_surface, _m11, _m12, _m21, _m22 = data[i]
            
            if _at_surface:
                u[i+1] = w[i+1] = 0
            else:
                # stage 1
                k1_u, k1_w = self.__compute_k(_u, _w, _rho, _pitch, _FBg, _m11, _m12, _m21, _m22, h)

                # stage 2
                _FBg, _pitch, _rho, _at_surface, _m11, _m12, _m21, _m22 = data_mid[i]
                __u = _u + k1_u*0.5
                __w = _w + k1_w*0.5
                k2_u, k2_w = self.__compute_k(__u, __w, _rho, _pitch, _FBg, _m11, _m12, _m21, _m22, h)
                
                # stage 3
                __u = _u + k2_u*0.5
                __w = _w + k2_w*0.5
                k3_u, k3_w = self.__compute_k(__u, __w, _rho, _pitch, _FBg, _m11, _m12, _m21, _m22, h)

                #stage 4
                _FBg, _pitch, _rho, _at_surface, _m11, _m12, _m21, _m22 = data[i+1]
                __u = _u + k3_u
                __w = _w + k3_w
                k4_u, k4_w = self.__compute_k(__u, __w, _rho, _pitch, _FBg, _m11, _m12, _m21, _m22, h)
                
                # Euler
                #u[i+1] = _u + k1_u
                #w[i+1] = _w + k1_w
                u[i+1] = _u + (k1_u + 2*k2_u + 2*k3_u + k4_u)/6
                w[i+1] = _w + (k1_w + 2*k2_w + 2*k3_w + k4_w)/6
                
                    
    def compute_inverted_mass_matrix(self, pitch):
        ''' Computes the inverse of the mass matrix

        not to be called directly

        '''
        m11 = self.k1*self.mg
        m22 = self.k2*self.mg
        #
        C2 = np.cos(pitch)**2
        CS = np.cos(pitch)*np.sin(pitch)
        denom = self.mg**2 + (m22+m11)*self.mg + m11*m22
        M11 = ((m22-m11)*C2 + m11 + self.mg) / denom
        M12 = ((m22-m11)*CS) / denom
        M21 = M12
        M22 = (-(m22-m11)*C2 + m22 + self.mg) / denom
        return M11, M12, M21, M22

    
    def integrate(self, data, dt=None):
        ''' integrate system

        not to be called directly
        '''
        tm = data['time']
        pitch = data['pitch']
        m_water_pressure = data['pressure']
        rho = data['density']
        m_de_oil_vol = data['buoyancy_change']
        dhdt = data['dhdt']
        pressure, Vbp = self.convert_pressure_Vbp_to_SI(m_water_pressure, m_de_oil_vol)

        h = dt or 2.
        ti = np.arange(tm.min(), tm.max()+h, h)
        pitch = np.interp(ti, tm, pitch)
        pressure = np.interp(ti, tm, pressure)
        rho = np.interp(ti, tm, rho)
        Vbp = np.interp(ti, tm, Vbp)
        dhdt = np.interp(ti, tm, dhdt)
        
        FB, Fg = self.compute_FB_and_Fg(pressure, rho, Vbp)
        M = self.compute_inverted_mass_matrix(pitch)

        u = np.zeros_like(FB)
        w = np.zeros_like(FB)
        threshold = self.max_depth_considered_surface * 1e4 # in Pa.
        self.RK4(h, M, FB-Fg, pitch, rho, pressure<threshold, u, w)

        gamma = np.arctan2(w, u)
        alpha = gamma - pitch
        Umag = (u**2 + w**2)**(0.5)
        ur = np.cos(pitch)*u + np.sin(pitch)*w
        wr = -np.sin(pitch)*u + np.cos(pitch)*w
        return Modelresult(ti, u, w, Umag, alpha, pitch, dhdt-w)

    def solve(self, data):
        ''' Solve the model

        :param data: environment data
        :type data: dict
        :return: model result
        :rytpe: Modelresult (named tuple with t, u, w, U, alpha,  pitch, ww fields)

        The data supplied should contain at least time, pressure, pitch, buoyancy_change and
        density, as reported by the glider. Depth rate (dhdt) will be added if not already present.
        Other data are ignored.

        :example:
        
        >>> dm = DynamciGliderModel(dt=1, k1=0.2, k2=0.98, rho0=1024)
        >>> dm.define(mg=70, Vg=68)
        >>> data = dict(time=time, pressure=P, pitch=pitch, buoyancy_change=vb, density=rho)
        >>> dm.solve(data)
        >>> plot(dm.U)
        '''

        pitch = data['pitch']
        try:
            dhdt = data['dhdt']
        except KeyError:
            dhdt = self.compute_dhdt(data['time'], data['pressure'])
            data['dhdt'] = dhdt

        integration_result = self.integrate(data, self.dt)
        self.modelresult = Modelresult(*[np.interp(data['time'],
                                                   integration_result.t, v) for v in integration_result])
        return self.modelresult

    
class Calibrate(object):
    ''' Generic class providing the calibration machinery for steady-state and dynamic flight models 

        Basic steps are to
  
        * set input data using the set_inputdata() method,
        * mask data that should not take part in the minimisation routine (data near the dive apices, 
          at the start of the dive, pycnoclines, etc.
        * run the calibration.
       

        Logical operators are used to mask data:

        OR, AND, and NAND are implemented, and work on the *mask*

    '''
    def __init__(self):
        self.mask=None
        
    def set_mask(self,mask):
        '''Set a mask 

        Masks those data that should not be used to calibrate.
        
        :param mask: mask
        :type mask: bool or numpy array
        

        .. note::
        If already set ones (after set_input_data(), then mask can be
        True or False to set all elements in mask.
        '''

        if not self.mask is None and type(mask)==bool: # the mask has been set before and True of False is set
            if mask:
                self.mask = np.ones_like(self.mask).astype(bool)
            else:
                self.mask = np.zeros_like(self.mask).astype(bool)
        else:
            if mask.dtype==bool:
                self.mask = mask
            else:
                self.mask = mask.astype(bool)

        
    def AND(self,mask):
        ''' Logical AND
        
        The new mask is the intersection (AND) of the existing mask and supplied mask

        :param mask: mask
        :type mask: bool or numpy array
        '''
        self.mask &= self.__ensure_booltype(mask)

    def OR(self,mask):
        ''' Logical OR

        The new mask is the union (OR) of the existing mask and the supplied mask.
        :param mask: mask
        :type mask: bool or numpy array
        '''
        self.mask |= self.__ensure_booltype(mask)

    def NAND(self,mask):
        ''' Logical NAND

        The new mask is the inverted intersection of the existing mask and the supplied mask.
        :param mask: mask
        :type mask: bool or numpy array
        '''
        self.mask &= ~self.__ensure_booltype(mask)
        
        
    def set_input_data(self,time ,pressure, pitch, buoyancy_change, density,
                       dhdt=None, u_relative=None, w_relative=None, U_relative=None):
        ''' Sets the input data
        time 
        pressure
        pitch
        buoyancy_change
        in-situ density
        and optionally u_relative and w_relative
        
        :param time: time (s)
        :param pressure: pressure (Pa)
        :param pitch: pitch (rad)
        :param buoyancy_change: buoyancy change reported by the glider (cc)
        :param density: in-situ density (kg m${-3}$)
        :param dhdt: depth rate m s$^{-1}$ (optional, if not given it is computed from pressure)
        :param u_relative: measured horizontal speed m s$^{-1}$ (optional)
        :param w_relative: measured vertical speed m s$^{-1}$ (optional)
        :param U_relative: measured speed through water m s$^{-1}$ (optional)
        :type time: float or numpy array
        :type pressure: float or numpy array
        :type pitch: float or numpy array
        :type buoyancy_change: float or numpy array
        :type density: float or numpy array
        :type dhdt: float or numpy array
        :type u_relative: float or numpy array
        :type w_relative: float or numpy array
        :type U_relative: float or numpy array

        .. note::
        A mask is automatically created (including all data) when this method is called.
        '''
        dhdt = dhdt or self.compute_dhdt(time, pressure)
        self.input_data = dict(time=time ,pressure=pressure, pitch=pitch, buoyancy_change=buoyancy_change,
                               density=density, dhdt=dhdt, u_relative=u_relative, w_relative=w_relative,
                               U_relative=U_relative)
        mask=np.zeros(time.shape,'int').astype(bool) # use all data by default
        self.set_mask(mask)
        
                        
    def cost_function(self,x,parameters, constraints, weights, verbose):
        ''' Cost-function used to optimise parameters

        This method first sets the parameters which are to be optimised for, and then
        computes the glider flight. A "cost" is computed from relatively weighted constraints.

        :param x: values of the parameters to be varied
        :param parameters: parameter names
        :param constraints: names of measured velocities against which glider flight is evaluated. These must be present in the dictionary supplied by the set_input_data() method.
        :param weights: if more than one constraint is provided, weights sets their relative importance.
        :param verbose: allows to print intermediate results during optimising
        :type x: numpy array
        :type parameters: tuple/list of strings
        :type constraints: tuple/list of string, or single string when using one constraint
        :type weights: tuple/list of floats
        :type verbose: bool
        :return: RMS value of exposed measurements (not masked)
        :rtype: float
        '''
        # set the data
        kwds=dict([(_p,_x) for _p,_x in zip(parameters,x)])
        self.define(**kwds) # ensures that aoa is reinitialised if necessary
        self.solve(self.input_data)
        U = self.modelresult.U
        alpha = self.modelresult.alpha
        weights = weights or np.ones_like(constraints, float)/len(constraints)

        w_mod = U*np.sin(self.input_data['pitch']+alpha)
        u_mod = U*np.cos(self.input_data['pitch']+alpha)

        error = 0
        
        for constraint, weight in zip(constraints, weights):
            if constraint == 'dhdt':
                error += weight * (w_mod - self.input_data['dhdt'])**2
            elif constraint == 'w_relative':
                error += weight * (w_mod - self.input_data['w_relative'])**2
            elif constraint == 'u_relative':
                error += weight * (u_mod - self.input_data['u_relative'])**2
            elif constraint == 'U_relative':
                error += weight * (U - self.input_data['U_relative'])**2
            else:
                raise NotImplementedError()
        error = error.compress(~self.mask)
        mse = np.nanmean(error)
        if verbose:
            s="  ".join(["{:s}={:.4f}".format(k,v) for k,v in kwds.items()])
            print("Error: {:.7e}  -  {}".format(mse, s))
        return mse

    def calibrate(self,*p, constraints = ('dhdt',), weights=None, verbose=False):
        ''' Calibrate model

        Given one or more model coefficients and specifications of measurements to use, this 
        method calibrates the model, using the self.cost_funtion() method. The interface is flexible
        so any parameter that is used in the model description can be optimised for. Also the velocity
        component or combination of components can be set.

        :param *p: variable length of parameter names to be optimised.
        :param constraints: names of measured velocities against which glider flight is evaluated. These must be present in the dictionary supplied by the set_input_data() method.
        :param weights: if more than one constraint is provided, weights sets their relative importance.
        :param verbose: allows to print intermediate results during optimising
        :type *p: string
        :type constraints: tuple/list of string, or single string when using one constraint
        :type weights: tuple/list of floats
        :type verbose: bool
        :return: the result of the optimisation routine
        :rtype: dict

        :example:
        
        >>> # calibrating for mass and drag coefficient (implicitly using depth-rate) and printing  
        >>> # intermediate results (mainly for debugging/progress monitoring)
        >>> results = gm.calibrate("mg", "Cd0", verbose=True)
        >>> print(results)
            {'mg': 70.00131, 'Cd0':0.145343}
        >>>
        >>> # Also calibrating the lift coefficient using measured incident water velocity
        >>> results = gm.calibrate("mg", "Cd0", "ah", constraints=('dhdt', 'U_relative'), weights=(0.5, 0.5), verbose = True)
        >>> print(results)
            {'mg': 70.00131, 'Cd0':0.145343, 'ah':3.78787}
        
        .. note::
        The default measurement to evaluate the model against is the depth rate dhdt. If not specified when
        setting the input data using the set_input_data() method, it is computed automatically. Other velocity
        components that are to be used to calibrate the model have to be set specifically.
        '''
        constraints = self.__ensure_iterable(constraints)
        x0=[self.__dict__[i] for i in p]
        args = (p, constraints, weights, verbose)
        R=fminimize(self.cost_function,x0,args=args,disp=False)
        rv=dict([(k,v) for k,v in zip(p,R)])
        return rv

    def __ensure_booltype(self, mask):
        if mask.dtype==bool:
            return mask
        else:
            return mask.astype(bool)

    def __ensure_iterable(self, c):
        if isinstance(c, str):
            return (c,)
        else:
            return c





class SteadyStateCalibrate(SteadyStateGliderModel, Calibrate):
    ''' Steady-state glider flight model, with calibration interface '''
    
    def __init__(self, rho0=None):
        SteadyStateGliderModel.__init__(self, rho0)
        Calibrate.__init__(self)


class DynamicCalibrate(DynamicGliderModel, Calibrate):
    ''' Dynamic glider flight model, with calibration interface '''
    def __init__(self, rho0=None, k1=0.02, k2=0.92, dt = None):
        DynamicGliderModel.__init__(self, dt=dt, rho0=rho0, k1=k1, k2=k2)
        Calibrate.__init__(self)
        


