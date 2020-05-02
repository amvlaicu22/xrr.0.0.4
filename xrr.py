import wx
import sys

import time
import numpy as np
import scipy
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import differential_evolution

import matplotlib.pyplot as plt
import CromerMan
from CromerMan import Cromer
from where import where

class _PhysConst:
    _wl_A = 1.54     #//A Cu Ka1
    _hh = 4.135667662E-15 #// eV*s
    _cc = 299792458 #// m/s
    _hc = 1.23984193 #// eV*um
    _keV = _hc/_wl_A*10.0
    _re = 2.8179403E-15      #//m
    _Na = 6.0221E23              #// 1/mol
        
    def __init__(self):
        self._keV = self.wlA2keV(self._wl_A)
        
    def wlA2keV(self, wlA):
        self._wl_A = wlA
        self._keV = self._hc/wlA*10.0
        return self._keV
        
    def keV2wlA(self, keV):
        self._keV = keV
        self._wl_A = self._hc/keV*10.0
        return self._wl_A

class _XLayer:
    _name = ""
    _comp = ""
    _dd = 0  # thickness A
    _sg = 0  # roughness A
    _rh = 0  # density g/cm^3
    
    _hd = ['Name', 'Composition', 'd [A]', 'sg [A]', 'rho [g/cm^3]', 'Mm [g/mol]', 'Vm [A^3]']
    
    AtomS = []
    AtomN = []
    AtomM = []
    AtomV = []
    _Mm = 0  # Molas mass
    _Vm = 0  # Molar volume
    
    Q = []   # q vector
    SF = []   # scattering factor(q) for the layer
    SFa = []  # scattering factor(q) for each atom type
    SLD = []  # SLD
    K = []   # wave-vector k
    F = []   # fresnel
    Phi = []  # phase
    Atn = []  # attenuation
    FA = []  # fresnel with attenuation
    Bt = []   # betha
    R = []   # reflectivity (complex)
    R2 = []
    
    _re = _PhysConst()._re
    _Na = _PhysConst()._Na
    _hc = _PhysConst()._hc
    
    _wl_A = 1.54
    _keV = 8.9
    
    
    def __init__(self, name="", comp="", d=1.0, s=0.0, r=1.0, M=1.0, V=1.0):
        self._name = name
        self._comp = comp
        self._dd = d  # thickness     A
        self._sg = s  # roughness     A
        self._rh = r  # density       g/cm^3
        self._Mm = M  # molar Mass    g
        self._Vm = V  # molar Volume  m^3
        
        self.AtomS = []
        self.AtomN = []
        self.AtomM = []
        self.AtomV = []
                
        self._re = _PhysConst()._re
        self._Na = _PhysConst()._Na
        self._hc = _PhysConst()._hc
        
        self._wlA = 1.54
        self.f_wlA(self._wlA)
        
        self.f_comp2chem()

        
    def f_qwaves(self, qq):
        # qq = Q vector
        nq = len(qq)
        self.Q = qq.copy()

        self.SF   = np.array(qq, dtype=np.complex_)   # scattering factor for layer
        self.SFa = []
        #https://periodictable.readthedocs.io/en/latest/api/xsf.html
        self.SLD  = np.array(qq, dtype=np.complex_)   # SLD scattering length density
        self.K   = np.array(qq, dtype=np.complex_)    # wave-vector k
        self.F   = np.array(qq, dtype=np.complex_)   # Fresnel
        self.Phi  = np.array(qq, dtype=np.complex_)  # phase
        self.Atn  = np.array(qq, dtype=np.complex_)     # attenuation
        self.FA  = np.array(qq, dtype=np.complex_)      #Fresnel w. attn
        self.Bt   = np.array(qq, dtype=np.complex_)    # betha
        self.R   = np.array(qq, dtype=np.complex_)    # reflectivity (complex)
        self.R2 = np.array(qq)                          #refl. int.
        
        
    def f_keV(self, keV):
        self._keV = keV
        self._wlA = self._hc/self._keV*10.0
        
    def f_wlA(self, wlA):
        self._wlA = wlA
        self._keV = self._hc/self._wlA*10.0
        
    def f_SF(self, qq, keV):
        # ~ input is "Atom name", Q = 2*pi/d=4*pi*sin(theta)/Lambda [A^-1], and energy in keV
        # ~ returns Complex f0+fp + i* fpp   in [electron units]        
        self.f_keV(keV)
        for atm, AN in zip(self.AtomS, self.AtomN):
            f0 = Cromer()._f0Q(atm, qq)
            fp = Cromer()._fpE(atm, keV)
            asfwv = f0 + fp[0] + fp[1]*1j 
            self.SFa.append(asfwv)
            self.SF += asfwv*AN
            
    def f_SLD(self):
        self.SLD = self._re*1.0E-14*(self._Na*self._rh/self._Mm)*self.SF
    
    def f_comp2chem(self):
        self.AtomS = []
        self.AtomN = []
        if self._comp == "":
            return
        chemlist = blank_strip(self._comp).split(';')
        self.AtomS = chemlist[0::2]
        atomNs = np.asarray(chemlist[1::2])
        self.AtomN = atomNs.astype(np.float)
        self.ntot = sum(i for i in self.AtomN)
        AtomM = []
        for atom in self.AtomS:
            mm = _MolarMass().mm(atom)
            AtomM.append( mm )
        self.AtomM = np.asarray(AtomM)
        self._Mm = sum(i for i in self.AtomM*self.AtomN)
        Vm_cm3 = self._Mm/self._rh/self._Na #//(g/mol)/(g/cm^3)/(1/mol) = cm^3
        Vm_A3  = Vm_cm3*1.0e24 #//A^3
        self._Vm = Vm_A3
        
    
    def __header__(self):
        return    "Name    | Chemical Composition | d_[A]| s_[A]|[g/cm^3]|[g/mol]| [A^3]"
                 
    def __line__(self):
        return    "_____________________________________________________________________"
        
        
    def __str__(self):
        t =  f"{self._name:5s}  {self._comp:12s}  {self._dd:4.1f}  {self._sg:4.1f}  {self._rh:7.4f}  "
        t += f"{self._Mm:5.1f}  {self._Vm:7.1f}"
        return t


    
########################################################################
########################################################################
########################################################################

class XRR():
    
    _a0 = 1.0 #//scale  
    _b0_m = 1.0 # bgnd mantisse
    _b0_e = -6  # bgnd exponent
    _b0 = _b0_m*10.0**_b0_e
    # ~ self.b0 = self.b0_m*np.exp(self.b0_e*np.ln(10.0)) #//backgroung zero level
    _2tht0 = 0.2
    
    _wl_A = 1.54     #//A Cu Ka1
    _keV = _PhysConst().wlA2keV(_wl_A)
    
    LL = [] # layer array
    
    mx = []     # measured data
    my = []
    xx = []     # truncated for xx < _2tht0
    yy = []
    rr0 = []    # reflectivity calculated
    rr1 = []    # with scale _a0 and background _b0

    qq = []     # q = 4 pi sin(tht)/wl
    qq4 = []    # q^4
    yq4 = []    # yy * q^4
    r0q4 = []   # rr0 * q^4
    r1q4 = []   # rr1 * q^4
    ryq4 = []   # fit rezidue (yq4 - r1q4)
    
    fitn = 0
    fiterr = 0
    fitkeys = [ 'ab', 'dd', 'rh', 'sg']
    fitmode = 0
    abort = 0
    fitmax = 100
    
    popt = []
    pcov = []
    perr = []
    ferr = []
    
    whrn = 3
        
    def __init__(self):
        super().__init__()
        where("", self.whrn)
        self.fitkeys = [ 'ab', 'dd', 'rh', 'sg']
        self.fitmode = 0
        self.fitmax = 100
        self.abort = 0
        self.fitn = 0
        
        self._a0 = 1.0 #//scale  
        self._b0_m = 1.0 # bgnd mantisse
        self._b0_e = -6  # bgnd exponent
        self._b0 = self._b0_m*10.0**self._b0_e
        
        self._keV = _PhysConst().wlA2keV(self._wl_A)
        
        self.xrr_str()
        # ~ self.data_load('xrr_test.dat')
        self.data = self.data_str
        self.data_parse()   # ~ >>  mx, my "measured data": delete data points < _2tht0
        self.data_init()   # >> xx, yy, qq, qq4, yq4, rr0, rr1, r0q4, r1q4, ryq4
        
        # ~ self.layer_load('xrr_test.xrr')
        # ~ self.layer_save('xrr_test.xrr')
        self.layer = self.layer_str
        self.layer_parse() 
        self.layer_init()
        self.layer_print()
        self.layer_profile()
        
        self.xrr()
        # ~ self.fit()
        self.fitmax = 0


    
    def layer_add(self, layer = _XLayer()):
        self.LL.append(layer)
    
    def layer_save(self, fname): 
        fileh = open(fname,"w")
        header = '# name , chem_comp , thickness_[A] , roughness_[A] , density_[g/cm^3] \n'
        fileh.write(header)
        # ~ print(header)
        for L in self.LL:
            line = L._name + ' , ' + L._comp + ' , ' + str(L._dd) + ' , ' + str(L._sg) + ' , ' + str(L._rh) + ' \n' 
            fileh.write(line)
            # ~ print(line)
        fileh.close()
        
        
    def layer_load(self, fname):
        self.file_layer = fname
        with open(fname) as fileh:  
            layer = fileh.readlines()
            fileh.close()
        self.layer = ""
        for l in layer:
            self.layer += l
            
    def layer_parse(self):
        where("", self.whrn)
        self.LL = []
        for line in self.layer.splitlines():
            if line[0] != '#' :
                # ~ print(line)
                [name, comp, s_dd, s_sg, s_rh] = line.split(',')
                self.layer_add(_XLayer(name,  comp, d=float(s_dd),  s=float(s_sg), r=float(s_rh) ))
        return
    
    def layer_print(self):
        print( _XLayer().__line__() )
        print( _XLayer().__header__() )
        print( _XLayer().__line__() )
        for  L in self.LL: 
            print( L )
        print( _XLayer().__line__() )
        
    
    def layer_init(self):
        where("", self.whrn)
        for  L in self.LL: 
            L.f_comp2chem()
            L.f_qwaves(self.qq)
            L.f_SF(self.qq, self._keV)
            L.f_SLD()
        
        
    def data_load(self, fname):
        self.file_data = fname
        with open(fname) as fileh:  
            data = fileh.readlines()
            fileh.close()
        self.data = ""
        for d in data:
            self.data += d

            
    def data_parse(self):
        where("", self.whrn)
        xx, yy = [], []
        for line in self.data.splitlines():
            vx, vy = line.split()
            xx.append(float(vx))
            yy.append(float(vy))
        self.mx = np.asarray(xx)
        self.my = np.asarray(yy)
        
    def data_init(self):
        where("", self.whrn)
        self.xx = self.mx.copy()
        self.yy = self.my.copy()

        boolArr = ( self.xx < self._2tht0 )
        indxArr = np.argwhere( boolArr )
        self.xx = np.delete(self.xx, indxArr)
        self.yy = np.delete(self.yy, indxArr)

        vv =self.yy.max()
        self.yy /= vv 
        self.my /= vv
        
        self.qq = self.xx.copy()
        self.qq = 4.0*np.pi*np.sin(self.xx/2.0*np.pi/180.)/self._wl_A
        
        self.qq4 = self.qq.copy()
        self.qq4 = self.qq**4   #   qq*qq*qq*qq
        
        self.yq4 = self.yy.copy()
        self.yq4 = self.yy * self.qq4
        
        self.rr0 = self.yy.copy()  # without scale and background
        self.rr0 = self.yy*0.8
                
        self.rr1 = self.yy.copy()  # = a0*refl0+b0
        self._a0=1.0
        self._b0=1.0e-6
        self.rr1 = self.rr0*self._a0 + self._b0
        
        self.r0q4 = self.yy.copy()   # = refl0*qq4
        self.r0q4 = self.rr0*self.qq4
        
        self.r1q4 = self.yy.copy()   # = refl1*qq4
        self.r1q4 = self.rr1*self.qq4
        
        self.ryq4 = np.zeros(len(self.qq)) # fit residue
        self.ryq4 = np.abs(self.yq4)/10.0
        
        
    def plot(self):
        fig, axs = plt.subplots(2, 2)
        ax1 = axs[0][0]
        ax2 = axs[1][0]
        ax3 = axs[0][1]
        ax4 = axs[1][1]
        # Adjust the scaling factor to fit your legend text completely outside the plot
        # (smaller value results in more space being made for the legend)
        fig.subplots_adjust(right=0.98)
        ax = 1.02
        ay = 1.1
        drag = True
        self.plot1(ax1, ax, ay, drag)
        self.plot2(ax2, ax, ay, drag)
        self.plot3(ax3, ax, ay, drag)
        self.plot4(ax4, ax, ay, drag)
        plt.tight_layout()
        plt.show()
        return

        # ~ ax = pg.figure.add_subplot(111)
        # ~ ax = pg.figure.gca(projection = '3d') 
        
        # ~ axes = plt.gca()
        # ~ axes.set_xlim([xmin,xmax])
        # ~ axes.set_ylim([ymin,ymax])

        # ~ ticks = np.arange(0., 8.1, 2.)
        # ~ tickLb = ['%1.2f' % tick for tick in ticks] # list comprehension to get all tick labels...
        # ~ ax1.xaxis.set_ticks(ticks)
        # ~ ax2.xaxis.set_ticks(ticks)
        # ~ ax2.xaxis.set_ticklabels(tickLb)
        # ~ ax2.set_xlim(ax1.get_xlim())
        # ~ plot(x, y, 'go--', linewidth=2, markersize=12)
        # ~ plot(x, y, color='green', marker='o', linestyle='dashed', linewidth=2, markersize=12)
        
    def plot1(self, ax1, ax, ay, drag):
        ax1.plot(self.mx, self.my, label='$meas.$', c='k')
        ax1.plot(self.xx, self.yy, label='$meas.''$', c='g', linewidth=1.5)
        ax1.plot(self.xx, self.rr1, label='$calc.$', c='r', linewidth=1.5)
        ax1.legend(loc=2, bbox_to_anchor=(ax, ay), ncol=1, fontsize=8, title='').draggable(drag)
        ax1.grid(True)
        ax1.set_xlabel(r'$2 \theta \ [deg]$')  # raw string r'jfkejrf' or use '$\\theta \\rho$
        ax1.set_ylabel(r'$I/I_{max}$')
        ax1.set_yscale('log')

        
    def plot2(self, ax2, ax, ay, drag):
        ax2.plot(self.qq, self.yq4, label='meas.', c='g', linewidth=1.5)
        ax2.plot(self.qq, self.r1q4, label='calc.', c='r', linewidth=1.5)
        for L in self.LL:
            ax2.plot(self.qq, abs(L.R2*L.Q**4), label=L._name, linewidth=0.75)
        ax2.legend(loc=2, bbox_to_anchor=(ax, ay), ncol=1, fontsize=8, title=r'$R^2Q^4$').draggable(drag)  
        ax2.grid(True)
        ax2.set_yscale('log')
        ax2.set_xlabel(r'$Q=4 \pi \ sin(\theta)/ \lambda \ [1/ \AA]$')
        ax2.set_ylabel(r'$int*\ Q^4$')
        ax2.set_ylim([min(self.yq4)/5., max(self.yq4)*5.])

        
        
    def plot3(self, ax3, ax, ay, drag):
        ax3.plot(self.p_z, self.p_rh, c='r', label=r'$profile $')
        ax3.plot(self.p_z, self.p_rL, c='k', label=r'$layer $')
        ax3.legend(loc=2, bbox_to_anchor=(ax, ay), ncol=1, fontsize=8, title=r'$\rho \ [g/cm^3]$').draggable(drag)
        ax3.grid(True)
        ax3.set_ylabel(r'$\rho \ [g/cm^3]$')
        ax3.set_xlabel(r'$z \ [\AA]$')

    def plot4(self, ax4, ax, ay, drag):
        for a, p_c in zip(self.atoms, self.atomp):
            ax4.plot(self.p_z, p_c, label=a)
        ax4.legend(loc=2, bbox_to_anchor=(ax, ay), ncol=1, fontsize=8, title=r'$atoms$').draggable(drag)
        ax4.grid(True)
        ax4.set_ylabel(r'$comp. index$')
        ax4.set_xlabel(r'$z \ [\AA]$')

       
    def layer_to_fitparam(self):
        where("", self.whrn)
        self.pfit = []
        self.pkey = []
        def keyadd( key, idx): 
            self.pkey.append( [ key, idx, key+'['+str(idx)+']' ] )
        
        for key in self.fitkeys:
            if key == 'ab':
                self.pfit.append(self._a0)
                self.pfit.append(self._b0)
                keyadd('a0', -1)
                keyadd('b0', -1)
            if key == 'dd':
                l = 1
                for L in self.LL[1:-1]:
                    self.pfit.append( abs(L._dd))
                    keyadd('dd', l)
                    l = l+1
            if key == 'rh':
                l = 0
                for L in self.LL[:-1]:
                    self.pfit.append( abs(L._rh))
                    keyadd('rh', l)
                    l = l+1
            if key == 'sg':
                l = 0
                for L in self.LL[:-1]:
                    self.pfit.append( abs(L._sg))
                    keyadd('sg', l)
                    l = l+1
        # ~ for k, p in zip( self.pkey, self.pfit):
            # ~ print(k[2],'=', p)
        return self.pfit
        
                    
    def layer_from_fitparam(self, pw):
        # ~ where(f"{self.fitn}", 3)
        N=0
        for key in self.fitkeys:
            if key == 'ab':
                self._a0 = pw[N+0]    
                self._b0 = pw[N+1]
                N = N+2
            if key == 'dd':
                l=0
                for L in self.LL[1:-1]:
                    L._dd = abs(pw[N+l])
                    l=l+1
                N = N + len(self.LL[1:-1])
            if key == 'rh':
                l=0
                for L in self.LL[:-1]:
                    L._rh = abs(pw[N+l])
                    l=l+1
                N = N + len(self.LL[:-1])
            if key == 'sg':
                l=0
                for L in self.LL[:-1]:
                    L._sg = abs(pw[N+l])
                    l=l+1
                N = N + len(self.LL[:-1])
                
                
                
                
    def set_fitkeys(self, fitkeys=[ 'ab', 'dd', 'rh', 'sg']):
    #===================================================
        self.fitkeys = fitkeys
    
    def print_kwargs(self, **kwargs):
    #================================
        for key, value in kwargs.items():
            print("{} = {}, ".format(key, value), end='')
        print()
        
    def print_dictio(self, dictionary):
    #================================
        for key, value in dictionary.items():
            print(f"{key} = {value}")
            # ~ print("________________________________________________")
                
    def post_event(self, **kw):
    #=================================
        # ~ d = {'fitn':self.fitn, 'err':"{0:.1e}".format(self.fiterr)}
        d = {'fitn':self.fitn, 'err':f"{self.fiterr:.2e}" }
        d.update(kw)    # add more key:values from **kw
        # ~ m = {0:'L-M', 1:'D-E'}
        # ~ n = {'@':'XRR', 'fit':m[self.fitmode]}
        # ~ d.update(n)
        where(f"\n>>{d}",1)
    

    def fit(self):
    #========================
        where("", self.whrn)
        # ~ print('\t XRR.fit()')

        # ~ def zero(): self.r1q4[:] = 0.1
        def zero(): self.r1q4 = self.yq4
        #============================
        def xrrp(qq, *pw):
        #===================
            # ~ qq=4*pi*sin(tht)/wl
            # ~ pw=fit parameter wave
                
            self.layer_from_fitparam(pw)    
            self.r1q4 = self.xrr()  
            self.fitn = self.fitn + 1
            yc = self.r1q4
            ym = self.yq4
            self.fiterr = np.sqrt(sum((ym - yc)**2) / len(ym))
            self.ferr.append(self.fiterr)
            print(".", end='', flush=True)
        
            if not (self.fitn % update ):
                self.post_event(key=1, msg='fitting' )

            if (self.fitmax !=0):
                if (self.fitn > self.fitmax): 
                    zero()
                    self.post_event(key=1, msg='n > nmax')
                    #--------------------------------------
            if (self.abort ): 
                zero()
                self.post_event(key=-1, msg='aborting')
                #------------------------------------
            return self.r1q4 
        #def_end

        def xrrd(p, *data):
        #===================
            if not (self.fitn % update):
                self.post_event(key=1,  msg='fitting')
            self.fitn =self.fitn +1
            x, y = data
            self.layer_from_fitparam(p)
            self.xrr()   # >> xrr.r1q4 
            yc = self.r1q4
            ym = self.yq4
            self.fiterr = np.sqrt(sum((ym - yc)**2) / len(ym))
            self.ferr.append(self.fiterr)
            return self.fiterr
        #def_end
        
        class Callback(object):
            def __init__(self):
                self.nit = 0
            def __call__(self, par, convergence=0):
                # ~ print(self.nit, par, flush=True)
                # ~ print(self.nit, flush=True)
                self.nit += 1

        
        self.fitn = 0
        self.nit = 0
        self.abort = 0
        self.ferr = []
        update = 20
        maxiter = 5
        
        # ~ self.layer_print()
        self.layer_to_fitparam()
            
        if self.fitmode is 0:
        #====================
            
            print('fit Nonlinear Least-Squares ', self.fitkeys)  
            self.r1q4 = xrrp(self.qq, *self.pfit)

            # try using full_output
            res =  scipy.optimize.curve_fit( xrrp, self.qq, self.yq4, 
                p0 = self.pfit, full_output=1)
            popt, pcov, infodict, errmsg, ier = res
            print('_________________________________________________________________')
            print('errmsg =', errmsg)
            # ~ print('ier =', ier)
            # ~ print('infodict =')
            # ~ self.print_dictio(infodict)
                
            # ~ normal output fit
            #####################
            # ~ res =  scipy.optimize.curve_fit( xrrp, 
                # ~ self.qq, self.yq4, p0 = self.pfit)
            # ~ popt, pcov = res
                
            self.popt, self.pcov = popt, pcov
            self.perr = np.sqrt(np.diag(pcov))
            print('perr = ', self.perr)

            
        if self.fitmode is 1:
        #=====================
            print('fit differential evolution :', self.fitkeys)

            bounds = []
            v = 0.25
            for p in self.pfit:
                bounds.append( (p*(1-v), p*(1+v)) )
            # ~ print(bounds)

            x = [v for v in self.qq]
            y = [v for v in self.yq4]
            p = [v for v in self.pfit]

            args = (x, y)
            xrrd(p, *args) # check calling 

            cback = Callback()
            result = differential_evolution(xrrd, bounds, args=args, 
                disp=True, callback=cback, tol=0.1, maxiter=maxiter)
            #https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html
            # ~ from scipy.optimize import differential_evolution
            # ~ scipy.optimize.differential_evolution(func, bounds, args=(), 
            # ~ strategy='best1bin', maxiter=1000, popsize=15, tol=0.01, mutation=(0.5, 1), 
            # ~ recombination=0.7, seed=None, callback=None, disp=False, polish=True, 
            # ~ init='latinhypercube', atol=0, updating='immediate', workers=1, constraints=())
            
            self.popt = result.x
            print(result)
        
        ###############################################
        self.layer_from_fitparam(self.popt)
        self.xrr()
        self.ryq4 = self.yq4 - self.r1q4
        self.post_event(key=0, msg='finished')


    def xrr(self):
        # ~ where("",3)
        nL = len(self.LL)
        nq = len(self.qq)
        nf = nL -1

        for L in self.LL: L.f_SLD()     # scattering length density
        SLDN = self.LL[nf].SLD          # sld in last outer layer
        for L in self.LL:               # K wave vector
            L.K = np.sqrt(self.qq*self.qq/4.0 - 4.0*np.pi*(L.SLD-SLDN))
        
        self.LL[ 0].F[:] = 1.0 + 0j 
        self.LL[ 0].FA[:]= 1.0 + 0j 
        self.LL[ 0].R[:] = 0.0 + 0j 
        self.LL[ 0].Bt[:] = 1.0 + 0j
        self.LL[nf].Bt[:] = 1.0 + 0j
        self.LL[ 0].R2[:] = 0.0 
        
        R0   = self.LL[0].R
        K0   = self.LL[0].K
        sg0  = self.LL[0]._sg
        
        for L in self.LL[1:]:
            L.F     = (L.K - K0)/(L.K + K0)     #Fresnel
            L.Phi   = L.K * K0 * sg0 * sg0      #phase
            L.Atn   = np.exp(-L.Phi)            #attenuation
            L.FA    = L.F * L.Atn               #Fresnel w. attn.
            
            L.Bt    = np.exp(-2.0*1j*L.K*abs(L._dd)) #betha
            
            L.R     = L.Bt*(R0 + L.FA)/(1.0 + R0*L.FA)  #refl.
            L.R2    = abs(L.R)**2                       #refl. int.
            
            sg0     = L._sg
            K0      = L.K
            R0      = L.R
            
        self.rr0 = self.LL[nf].R2       # refl. int in last outer layer
        self.rr1 = abs(self._a0)*self.rr0+abs(self._b0)     # with scale and bgnd
        
        self.r0q4 = self.qq4*self.rr0
        self.r1q4 = self.qq4*self.rr1
        
        # ~ self.r0q4 = np.log(self.qq4*self.rr0)
        # ~ self.r1q4 = np.log(self.qq4*self.rr1)

        return self.r1q4
        
        
    def layer_profile(self):
        where("",self.whrn)
        # ~ print('\t XRR.layer_profile(), fitn=', self.fitn)

        from mpmath import erf
        def step0(x0, y1, y2, sg, x):
            return (y2 - y1)*(1. + erf( (x - x0)/sg) )/2. 
    
        def step(x0, y1, y2, sg, x):
            return np.array( [ (y2 - y1)*(1. + erf( (xi - x0)/sg) )/2. for xi in x ] )

        t0 = time.time()
        L_dd = np.asarray([ L._dd for L in self.LL ])
        L_sg = np.asarray([ L._sg for L in self.LL ])
        L_rh = np.asarray([ L._rh for L in self.LL ])
        
        nl = len(L_dd)

        L_dt = L_dd.copy()
        # ~ L_dt = [ 0 for i in L_dt ]
        L_dt[:] = 0
        dtt = 0
        for i in range(nl):
            dtt += L_dd[i]
            L_dt[i] = dtt
        dtt += 20 + 20 #//20 A at the substrate and top, each

        dx = 2.0
        nn = int(dtt/dx)
        p_z = np.array(nn)
        p_z = np.array([ -20.0 + x*dx for x in  range(nn) ])

        # ~ //density profile: (y) with sigma; (b) square
        p_rh = np.array(nn)
        p_rL = np.array(nn)
        st   = np.array(nn)
        
        p_rh = np.array( [ L_rh[0] for i in range(nn) ] )
        p_rL = np.array( [ L_rh[0] for i in range(nn) ] )

    
        # ~ print('profile step density')
        # ~ t0 = time.time()
        for i in range(nl-1):
            st = [ step0(L_dt[i],L_rh[i],L_rh[i+1], L_sg[i], x) for x in p_z]
            p_rh = p_rh + st
            st = [ step0(L_dt[i],L_rh[i],L_rh[i+1], 0.1, x) for x in p_z ]
            p_rL = p_rL + st
        # ~ t1 = time.time()
        # ~ print( t1-t0 )
        
        # for each layer type get composition in each layer
        atoms = []
        for L in self.LL:
            for a in L.AtomS:
                if a not in atoms:
                    atoms.append(a)
        
        atomc = []
        for a in atoms:
            LC = []
            for L in self.LL:
                if a in L.AtomS:
                    p = L.AtomS.index(a)
                    c = L.AtomN[p]
                else:
                    c = 0.0
                LC.append(c)
            atomc.append(LC)
        
        # ~ print('profile step atoms')
        # ~ t0 = time.time()
        atomp = []
        k=0
        for a in atoms:
            PA = np.array(nn)
            PA = [ atomc[k][0] for x in range(nn) ]
            for i in range(nl-1):
                PA = PA + step( L_dt[i], atomc[k][i], atomc[k][i+1], L_sg[i], p_z) 
            atomp.append(PA)
            k=k+1
        t1 = time.time()
        # ~ print( t1-t0 )
        
        self.p_z = p_z
        self.p_rh = p_rh
        self.p_rL = p_rL
        self.atoms = atoms
        self.atomc = atomc
        self.atomp = atomp
        
    def xrr_str(self):
        where("",self.whrn)
        # ~ print('\t XRR.xrr_str()')
        self.layer_str = \
"""# name , chem_comp , thickness_[A] , roughness_[A] , density_[g/cm^3] 
glass    , Si;1;O;2           ,  0.00 ,  6.13 ,  3.5200 
NoName   , Si;0.1;O;0.5;C;0.4 , 12.77 ,  3.43 ,  0.0030 
W        , W;1                , 42.36 ,  3.87 , 20.1600 
W-Si-O   , W;0.6;Si;0.3;O;0.1 ,  7.82 ,  3.90 , 10.1700 
NoName   , Si;0.1;O;0.5;C;0.4 ,  7.47 ,  5.59 ,  0.0014 
Si       , Si;1               , 62.96 , 35.62 ,  1.4800 
SiO2     , Si;1;O;2           , 14.99 ,  4.82 ,  3.3600 
air      , N;1.5;O;0.5        ,  0.00 ,  0.00 ,  0.0012"""

        self.data_str = \
"""	-0.1	0.003532061
	-0.09	0.007885087
	-0.08	0.0171938
	-0.07	0.03694645
	-0.06	0.1086672
	-0.05	0.3739794
	-0.04	1.098235
	-0.03	1.89013
	-0.02	2.849063
	-0.01	3.775216
	-3.4694e-18	4.602186
	0.01	5.048631
	0.02	5.247839
	0.03	5.384246
	0.04	5.415466
	0.05	5.381364
	0.06	5.316523
	0.07	4.900577
	0.08	4.182397
	0.09	3.282901
	0.1	2.395894
	0.11	1.595941
	0.12	1.031448
	0.13	0.7582853
	0.14	0.6469861
	0.15	0.6774736
	0.16	0.7258646
	0.17	0.6982829
	0.18	0.7445005
	0.19	0.7650576
	0.2	0.8182998
	0.21	0.8497479
	0.22	0.8887008
	0.23	0.9513689
	0.24	0.965646
	0.25	1
	0.26	1.017183
	0.27	1.117111
	0.28	1.151945
	0.29	1.194777
	0.3	1.185819
	0.31	1.228626
	0.32	1.267411
	0.33	1.28134
	0.34	1.32853
	0.35	1.402618
	0.36	1.441523
	0.37	1.443804
	0.38	1.517291
	0.39	1.552354
	0.4	1.536383
	0.41	1.568084
	0.42	1.624039
	0.43	1.631004
	0.44	1.688881
	0.45	1.723703
	0.46	1.707973
	0.47	1.705812
	0.48	1.719861
	0.49	1.633165
	0.5	1.68792
	0.51	1.67195
	0.52	1.644332
	0.53	1.677714
	0.54	1.625
	0.55	1.694645
	0.56	1.699688
	0.57	1.672671
	0.58	1.712056
	0.59	1.707733
	0.6	1.707973
	0.61	1.725865
	0.62	1.692724
	0.63	1.740634
	0.64	1.721542
	0.65	1.732829
	0.66	1.713017
	0.67	1.720821
	0.68	1.713257
	0.69	1.701489
	0.7	1.684918
	0.71	1.667987
	0.72	1.662464
	0.73	1.600745
	0.74	1.606388
	0.75	1.579251
	0.76	1.565442
	0.77	1.537584
	0.78	1.512968
	0.79	1.508646
	0.8	1.489914
	0.81	1.463857
	0.82	1.457973
	0.83	1.376681
	0.84	1.403818
	0.85	1.398535
	0.86	1.351105
	0.87	1.345941
	0.88	1.338737
	0.89	1.361912
	0.9	1.318684
	0.91	1.307517
	0.92	1.268612
	0.93	1.292747
	0.94	1.292027
	0.95	1.307997
	0.96	1.309918
	0.97	1.257925
	0.98	1.270293
	0.99	1.258646
	1	1.277738
	1.01	1.257445
	1.02	1.250961
	1.03	1.232949
	1.04	1.224784
	1.05	1.209294
	1.06	1.20305
	1.07	1.15993
	1.08	1.15365
	1.09	1.15945
	1.1	1.168408
	1.11	1.131136
	1.12	1.112992
	1.13	1.07452
	1.14	1.064842
	1.15	1.067027
	1.16	0.995161
	1.17	1.000733
	1.18	0.9603266
	1.19	0.9654059
	1.2	0.9274136
	1.21	0.8928195
	1.22	0.887488
	1.23	0.8681316
	1.24	0.8158742
	1.25	0.8042628
	1.26	0.7561119
	1.27	0.7469141
	1.28	0.7391691
	1.29	0.6852186
	1.3	0.65933
	1.31	0.6220701
	1.32	0.6363473
	1.33	0.5923031
	1.34	0.5267292
	1.35	0.539073
	1.36	0.521902
	1.37	0.4904419
	1.38	0.4846302
	1.39	0.4386648
	1.4	0.4096302
	1.41	0.3934198
	1.42	0.3915106
	1.43	0.3828771
	1.44	0.3728386
	1.45	0.3605788
	1.46	0.3466378
	1.47	0.3302714
	1.48	0.3150696
	1.49	0.3018012
	1.5	0.2875
	1.51	0.2779899
	1.52	0.265514
	1.53	0.2547911
	1.54	0.2441283
	1.55	0.2340178
	1.56	0.2240394
	1.57	0.2157901
	1.58	0.207805
	1.59	0.2007685
	1.6	0.1940442
	1.61	0.1860951
	1.62	0.1798391
	1.63	0.1748439
	1.64	0.1685519
	1.65	0.1620557
	1.66	0.1579131
	1.67	0.151537
	1.68	0.1480908
	1.69	0.1430956
	1.7	0.1385327
	1.71	0.1350504
	1.72	0.1289265
	1.73	0.1254203
	1.74	0.1215178
	1.75	0.1178771
	1.76	0.1135123
	1.77	0.1089481
	1.78	0.1049376
	1.79	0.1021122
	1.8	0.0973499
	1.81	0.09299472
	1.82	0.08916546
	1.83	0.08630044
	1.84	0.08163905
	1.85	0.07876441
	1.86	0.07330572
	1.87	0.06960615
	1.88	0.06744476
	1.89	0.06244357
	1.9	0.05877281
	1.91	0.05504323
	1.92	0.05147815
	1.93	0.04814601
	1.94	0.04549112
	1.95	0.04056196
	1.96	0.03802954
	1.97	0.03700408
	1.98	0.03427954
	1.99	0.03154899
	2	0.0285999
	2.01	0.02606268
	2.02	0.02338617
	2.03	0.02082493
	2.04	0.01872719
	2.05	0.01672671
	2.06	0.01460735
	2.07	0.01260927
	2.08	0.01093132
	2.09	0.009238113
	2.1	0.007698487
	2.11	0.006527738
	2.12	0.00532721
	2.13	0.004326849
	2.14	0.003528218
	2.15	0.00279707
	2.16	0.002151297
	2.17	0.001607109
	2.18	0.001155163
	2.19	0.0008041427
	2.2	0.0005538064
	2.21	0.0003866235
	2.22	0.0003041427
	2.23	0.0003017772
	2.24	0.0003552954
	2.25	0.0004663545
	2.26	0.0006335015
	2.27	0.0008208814
	2.28	0.001061924
	2.29	0.001322767
	2.3	0.001597022
	2.31	0.001866234
	2.32	0.002154299
	2.33	0.00245257
	2.34	0.002752762
	2.35	0.003024016
	2.36	0.003306436
	2.37	0.003553554
	2.38	0.003790706
	2.39	0.004032181
	2.4	0.004234991
	2.41	0.004450528
	2.42	0.004637968
	2.43	0.004791187
	2.44	0.004930716
	2.45	0.00510939
	2.46	0.00518684
	2.47	0.00532829
	2.48	0.005455331
	2.49	0.005535663
	2.5	0.005612752
	2.51	0.005682517
	2.52	0.005801993
	2.53	0.005884487
	2.54	0.005965778
	2.55	0.006075168
	2.56	0.006153699
	2.57	0.006259246
	2.58	0.00635879
	2.59	0.006533982
	2.6	0.006649136
	2.61	0.006823487
	2.62	0.006991955
	2.63	0.007070005
	2.64	0.007276177
	2.65	0.007556196
	2.66	0.007701729
	2.67	0.007929035
	2.68	0.008231868
	2.69	0.008498319
	2.7	0.008675072
	2.71	0.008989794
	2.72	0.009277859
	2.73	0.00957709
	2.74	0.009845221
	2.75	0.01019549
	2.76	0.01045245
	2.77	0.01066066
	2.78	0.01093588
	2.79	0.01125048
	2.8	0.01157553
	2.81	0.0118152
	2.82	0.01215898
	2.83	0.01246038
	2.84	0.012494
	2.85	0.01276057
	2.86	0.01303915
	2.87	0.01317483
	2.88	0.01324328
	2.89	0.0134294
	2.9	0.01353026
	2.91	0.01373679
	2.92	0.01362992
	2.93	0.01383405
	2.94	0.01387608
	2.95	0.01374039
	2.96	0.01379083
	2.97	0.01364794
	2.98	0.01365634
	2.99	0.01354587
	3	0.01340538
	3.01	0.01319164
	3.02	0.01309198
	3.03	0.01277378
	3.04	0.01244717
	3.05	0.01229227
	3.06	0.01208213
	3.07	0.01184762
	3.08	0.01151873
	3.09	0.0112153
	3.1	0.01087872
	3.11	0.0106513
	3.12	0.01012848
	3.13	0.009892412
	3.14	0.009528939
	3.15	0.009289986
	3.16	0.008916186
	3.17	0.008478386
	3.18	0.008242675
	3.19	0.007957374
	3.2	0.007565322
	3.21	0.007255284
	3.22	0.006978026
	3.23	0.006664746
	3.24	0.006311839
	3.25	0.00615694
	3.26	0.005768732
	3.27	0.005384486
	3.28	0.005216259
	3.29	0.00501537
	3.3	0.004853266
	3.31	0.004650936
	3.32	0.004476465
	3.33	0.00429611
	3.34	0.004117916
	3.35	0.003955812
	3.36	0.003780379
	3.37	0.003621638
	3.38	0.003507925
	3.39	0.003354707
	3.4	0.003247839
	3.41	0.003137848
	3.42	0.00304707
	3.43	0.002962416
	3.44	0.002867916
	3.45	0.002803314
	3.46	0.002727065
	3.47	0.002660783
	3.48	0.002581532
	3.49	0.002525216
	3.5	0.002488473
	3.51	0.002443924
	3.52	0.002398175
	3.53	0.002326369
	3.54	0.002278819
	3.55	0.002240034
	3.56	0.002188521
	3.57	0.002159222
	3.58	0.002105548
	3.59	0.002053794
	3.6	0.002004683
	3.61	0.001961095
	3.62	0.00189013
	3.63	0.001839938
	3.64	0.001793228
	3.65	0.001740274
	3.66	0.001673511
	3.67	0.001625961
	3.68	0.001552594
	3.69	0.001490394
	3.7	0.001416547
	3.71	0.00135987
	3.72	0.001279179
	3.73	0.001215058
	3.74	0.00115598
	3.75	0.001081208
	3.76	0.001006532
	3.77	0.0009410543
	3.78	0.0008854347
	3.79	0.0008129443
	3.8	0.0007350265
	3.81	0.0006831412
	3.82	0.000626525
	3.83	0.0005677594
	3.84	0.0005177354
	3.85	0.0004695005
	3.86	0.0004158862
	3.87	0.0003813281
	3.88	0.0003346422
	3.89	0.0002985831
	3.9	0.0002637488
	3.91	0.0002373919
	3.92	0.0002143732
	3.93	0.0001867435
	3.94	0.0001700168
	3.95	0.0001567843
	3.96	0.0001552353
	3.97	0.0001491715
	3.98	0.0001423871
	3.99	0.0001493636
	4	0.0001517771
	4.01	0.0001570605
	4.02	0.0001745197
	4.03	0.0001814962
	4.04	0.0002064842
	4.05	0.0002160543
	4.06	0.0002292387
	4.07	0.0002507085
	4.08	0.0002822166
	4.09	0.0003047791
	4.1	0.0003283982
	4.11	0.0003519693
	4.12	0.0003785423
	4.13	0.0004023055
	4.14	0.0004209534
	4.15	0.0004469861
	4.16	0.0004683237
	4.17	0.0004925312
	4.18	0.0005127762
	4.19	0.0005295509
	4.2	0.0005607829
	4.21	0.0005617795
	4.22	0.0005748199
	4.23	0.0005957493
	4.24	0.0006090658
	4.25	0.0006191403
	4.26	0.0006248799
	4.27	0.0006289385
	4.28	0.0006319885
	4.29	0.0006403338
	4.3	0.0006393732
	4.31	0.0006386047
	4.32	0.0006415226
	4.33	0.0006439362
	4.34	0.0006308958
	4.35	0.0006060519
	4.36	0.0006174952
	4.37	0.000595245
	4.38	0.00059055
	4.39	0.0005838977
	4.4	0.0005615034
	4.41	0.0005636528
	4.42	0.000548463
	4.43	0.0005331052
	4.44	0.0005244837
	4.45	0.0005177354
	4.46	0.0005005043
	4.47	0.0004912464
	4.48	0.000479719
	4.49	0.000471926
	4.5	0.0004583814
	4.51	0.0004478026
	4.52	0.0004478987
	4.53	0.0004283382
	4.54	0.0004142411
	4.55	0.000423367
	4.56	0.0004122839
	4.57	0.0004069044
	4.58	0.000412464
	4.59	0.0004051753
	4.6	0.0004136048
	4.61	0.0004154299
	4.62	0.0004029419
	4.63	0.0004057157
	4.64	0.0004068564
	4.65	0.0004113233
	4.66	0.0004206292
	4.67	0.0004248679
	4.68	0.0004252762
	4.69	0.0004226345
	4.7	0.0004381364
	4.71	0.0004373559
	4.72	0.0004534102
	4.73	0.000453134
	4.74	0.0004570125
	4.75	0.0004628962
	4.76	0.0004668588
	4.77	0.0004701369
	4.78	0.0004712416
	4.79	0.0004819525
	4.8	0.0004902498
	4.81	0.0004938041
	4.82	0.0004864193
	4.83	0.000490658
	4.84	0.000485963
	4.85	0.0004826369
	4.86	0.0004811239
	4.87	0.0004824568
	4.88	0.000478074
	4.89	0.0004708214
	4.9	0.0004679515
	4.91	0.0004579731
	4.92	0.0004552834
	4.93	0.0004432397
	4.94	0.0004367195
	4.95	0.0004244597
	4.96	0.0004151177
	4.97	0.0004080932
	4.98	0.0003910903
	4.99	0.0003768612
	5	0.0003743996
	5.01	0.0003600384
	5.02	0.0003475024
	5.03	0.0003316763
	5.04	0.0003216018
	5.05	0.0003080692
	5.06	0.0002910543
	5.07	0.0002755163
	5.08	0.0002670269
	5.09	0.0002554515
	5.1	0.000246878
	5.11	0.0002253122
	5.12	0.000217159
	5.13	0.0002090778
	5.14	0.0001932157
	5.15	0.0001798631
	5.16	0.0001758886
	5.17	0.0001644092
	5.18	0.0001596662
	5.19	0.0001482229
	5.2	0.0001414745
	5.21	0.0001362272
	5.22	0.0001270173
	5.23	0.0001211335
	5.24	0.0001211335
	5.25	0.0001199976
	5.26	0.0001099676
	5.27	0.0001122466
	5.28	0.0001082805
	5.29	0.0001066847
	5.3	0.0001049063
	5.31	0.0001044513
	5.32	0.000111609
	5.33	0.0001058646
	5.34	0.0001097394
	5.35	0.0001109246
	5.36	0.0001117459
	5.37	0.0001160315
	5.38	0.000117536
	5.39	0.0001221422
	5.4	0.0001275216
	5.41	0.000129659
	5.42	0.0001277978
	5.43	0.0001237872
	5.44	0.0001306196
	5.45	0.0001356868
	5.46	0.0001376441
	5.47	0.0001369116
	5.48	0.0001400096
	5.49	0.0001416066
	5.5	0.0001444765
	5.51	0.0001451129
	5.52	0.0001483069
	5.53	0.0001416547
	5.54	0.0001460255
	5.55	0.0001425673
	5.56	0.0001472142
	5.57	0.0001440202
	5.58	0.0001428794
	5.59	0.0001408742
	5.6	0.0001380043
	5.61	0.0001393732
	5.62	0.000134354
	5.63	0.0001365034
	5.64	0.0001338136
	5.65	0.0001298511
	5.66	0.000131172
	5.67	0.0001212776
	5.68	0.0001153014
	5.69	0.0001172623
	5.7	0.0001162584
	5.71	0.0001116547
	5.72	0.0001145269
	5.73	0.0001120641
	5.74	0.0001110615
	5.75	0.0001015778
	5.76	9.916187e-05
	5.77	0.0001000288
	5.78	9.943565e-05
	5.79	9.419308e-05
	5.8	9.524135e-05
	5.81	9.738353e-05
	5.82	9.054515e-05
	5.83	8.85843e-05
	5.84	8.963377e-05
	5.85	9.122959e-05
	5.86	9.049952e-05
	5.87	9.040826e-05
	5.88	8.812921e-05
	5.89	8.995197e-05
	5.9	9.255164e-05
	5.91	9.227786e-05
	5.92	9.031701e-05
	5.93	9.100145e-05
	5.94	9.18672e-05
	5.95	9.163906e-05
	5.96	9.847863e-05
	5.97	9.633526e-05
	5.98	0.0001005295
	5.99	9.893372e-05
	6	0.0001037668
	6.01	0.0001062284
	6.02	0.0001043144
	6.03	0.0001093288
	6.04	0.0001041775
	6.05	0.0001064565
	6.06	0.000111244
	6.07	0.0001126573
	6.08	0.0001100132
	6.09	0.0001176717
	6.1	0.0001123835
	6.11	0.0001147538
	6.12	0.0001166691
	6.13	0.0001189037
	6.14	0.000116806
	6.15	0.0001173991
	6.16	0.0001177173
	6.17	0.0001137056
	6.18	0.0001116547
	6.19	0.0001164409
	6.2	0.000112201
	6.21	0.0001117459
	6.22	0.000110377
	6.23	0.0001058646
	6.24	0.0001054995
	6.25	0.0001031736
	6.26	0.0001024904
	6.27	0.000102171
	6.28	9.943565e-05
	6.29	9.40562e-05
	6.3	9.401057e-05
	6.31	8.822046e-05
	6.32	8.416187e-05
	6.33	8.70341e-05
	6.34	8.060639e-05
	6.35	7.932998e-05
	6.36	7.581893e-05
	6.37	7.021134e-05
	6.38	6.647335e-05
	6.39	6.679155e-05
	6.4	6.268852e-05
	6.41	5.926874e-05
	6.42	5.954251e-05
	6.43	5.575889e-05
	6.44	5.042387e-05
	6.45	4.846422e-05
	6.46	4.800793e-05
	6.47	4.727906e-05
	6.48	4.31292e-05
	6.49	4.0622e-05
	6.5	4.071326e-05
	6.51	3.761287e-05
	6.52	3.268972e-05
	6.53	3.706653e-05
	6.54	3.305355e-05
	6.55	3.21878e-05
	6.56	3.132205e-05
	6.57	2.981748e-05
	6.58	2.954371e-05
	6.59	2.726345e-05
	6.6	2.945245e-05
	6.61	2.48475e-05
	6.62	2.630644e-05
	6.63	2.480187e-05
	6.64	2.571326e-05
	6.65	2.352546e-05
	6.66	2.434558e-05
	6.67	2.475624e-05
	6.68	2.34342e-05
	6.69	2.471062e-05
	6.7	2.347983e-05
	6.71	2.252282e-05
	6.72	2.161023e-05
	6.73	2.165586e-05
	6.74	2.142772e-05
	6.75	2.179275e-05
	6.76	2.357108e-05
	6.77	2.23403e-05
	6.78	2.179275e-05
	6.79	2.060759e-05
	6.8	2.229467e-05
	6.81	2.192964e-05
	6.82	1.996878e-05
	6.83	2.037944e-05
	6.84	1.946806e-05
	6.85	1.782661e-05
	6.86	1.846422e-05
	6.87	1.923991e-05
	6.88	1.782661e-05
	6.89	1.604827e-05
	6.9	1.814601e-05
	6.91	1.668708e-05
	6.92	1.613953e-05
	6.93	1.632205e-05
	6.94	1.481748e-05
	6.95	1.472622e-05
	6.96	1.344981e-05
	6.97	1.42243e-05
	6.98	1.221902e-05
	6.99	1.372358e-05
	7	1.144356e-05
	7.01	1.308478e-05
	7.02	1.253722e-05
	7.03	1.180824e-05
	7.04	1.121554e-05
	7.05	1.162584e-05
	7.06	1.121554e-05
	7.07	1.016691e-05
	7.08	1.240154e-05
	7.09	9.847863e-06
	7.1	1.15347e-05
	7.11	1.17171e-05
	7.12	1.24916e-05
	7.13	1.144356e-05
	7.14	1.162584e-05
	7.15	1.126117e-05
	7.16	1.199064e-05
	7.17	1.367796e-05
	7.18	1.326729e-05
	7.19	1.404179e-05
	7.2	1.472622e-05
	7.21	1.449808e-05
	7.22	1.618516e-05
	7.23	1.591138e-05
	7.24	1.668708e-05
	7.25	1.910303e-05
	7.26	1.778098e-05
	7.27	1.933117e-05
	7.28	2.006004e-05
	7.29	2.320605e-05
	7.3	2.006004e-05
	7.31	2.083574e-05
	7.32	2.165586e-05
	7.33	2.192964e-05
	7.34	2.352546e-05
	7.35	2.585015e-05
	7.36	2.685399e-05
	7.37	2.316042e-05
	7.38	2.457373e-05
	7.39	2.648895e-05
	7.4	2.872238e-05
	7.41	2.885927e-05
	7.42	2.749159e-05
	7.43	2.813041e-05
	7.44	2.863112e-05
	7.45	2.963497e-05
	7.46	2.963497e-05
	7.47	2.990874e-05
	7.48	2.726345e-05
	7.49	2.972623e-05
	7.5	2.808477e-05
	7.51	3.036383e-05
	7.52	2.849544e-05
	7.53	2.917867e-05
	7.54	3.031821e-05
	7.55	2.790226e-05
	7.56	2.803915e-05
	7.57	3.009006e-05
	7.58	2.785663e-05
	7.59	2.708093e-05
	7.6	2.790226e-05
	7.61	2.708093e-05
	7.62	2.5622e-05
	7.63	2.658021e-05
	7.64	2.753722e-05
	7.65	2.616955e-05
	7.66	2.557637e-05
	7.67	2.771974e-05
	7.68	2.475624e-05
	7.69	2.708093e-05
	7.7	2.493876e-05
	7.71	2.548631e-05
	7.72	2.516691e-05
	7.73	2.475624e-05
	7.74	2.466499e-05
	7.75	2.183838e-05
	7.76	2.329731e-05
	7.77	2.302353e-05
	7.78	2.393612e-05
	7.79	2.384486e-05
	7.8	2.42087e-05
	7.81	2.31148e-05
	7.82	2.074448e-05
	7.83	2.147334e-05
	7.84	2.197527e-05
	7.85	2.060759e-05
	7.86	2.274976e-05
	7.87	1.978626e-05
	7.88	2.119957e-05
	7.89	2.028819e-05
	7.9	1.873799e-05
	7.91	1.864674e-05
	7.92	2.151897e-05
	7.93	1.983189e-05
	7.94	1.869236e-05
	7.95	1.983189e-05
	7.96	1.960495e-05
	7.97	2.006004e-05
	7.98	1.791787e-05
	7.99	1.750721e-05
	8	1.723343e-05
	8.01	1.741595e-05
	8.02	1.60939e-05
	8.03	1.910303e-05
	8.04	1.800913e-05
	8.05	1.586575e-05
	8.06	1.691403e-05
	8.07	1.732469e-05
	8.08	1.559198e-05
	8.09	1.586575e-05
	8.1	1.417867e-05
	8.11	1.741595e-05
	8.12	1.518252e-05
	8.13	1.509126e-05
	8.14	1.540946e-05
	8.15	1.436119e-05
	8.16	1.449808e-05
	8.17	1.20365e-05
	8.18	1.331292e-05
	8.19	1.303914e-05
	8.2	1.144356e-05
	8.21	1.185387e-05
	8.22	1.11244e-05
	8.23	1.15347e-05
	8.24	9.756604e-06
	8.25	9.437441e-06
	8.26	9.300673e-06
	8.27	9.209534e-06
	8.28	8.571205e-06
	8.29	8.480068e-06
	8.3	8.84474e-06
	8.31	7.613833e-06
	8.32	7.112272e-06
	8.33	7.932998e-06
	8.34	6.70197e-06
	8.35	6.474063e-06
	8.36	6.10927e-06
	8.37	6.291667e-06
	8.38	5.835735e-06
	8.39	5.653338e-06
	8.4	4.923872e-06
	8.41	4.878362e-06
	8.42	4.194404e-06
	8.43	4.46794e-06
	8.44	4.376801e-06
	8.45	4.650336e-06
	8.46	3.237032e-06
	8.47	3.100264e-06
	8.48	3.237032e-06
	8.49	3.510567e-06
	8.5	2.826729e-06
	8.51	3.556196e-06
	8.52	2.553194e-06
	8.53	2.644332e-06
	8.54	2.416306e-06
	8.55	2.644332e-06
	8.56	2.598703e-06
	8.57	2.097262e-06
	8.58	1.732469e-06
	8.59	1.823727e-06
	8.6	2.006004e-06
	8.61	2.325168e-06
	8.62	2.142771e-06
	8.63	2.461936e-06
	8.64	2.188401e-06
	8.65	2.188401e-06
	8.66	1.823727e-06
	8.67	2.23403e-06
	8.68	2.006004e-06
	8.69	2.006004e-06
	8.7	1.960495e-06
	8.71	1.732469e-06
	8.72	1.641331e-06
	8.73	1.778098e-06
	8.74	1.823727e-06
	8.75	1.778098e-06
	8.76	1.869236e-06
	8.77	2.097262e-06
	8.78	2.006004e-06
	8.79	2.507565e-06
	8.8	1.960495e-06
	8.81	2.23403e-06
	8.82	3.009006e-06
	8.83	2.735471e-06
	8.84	2.051633e-06
	8.85	2.416306e-06
	8.86	3.054635e-06
	8.87	2.553194e-06
	8.88	2.7811e-06
	8.89	2.826729e-06
	8.9	2.006004e-06
	8.91	2.461936e-06
	8.92	2.416306e-06
	8.93	2.735471e-06
	8.94	3.556196e-06
	8.95	3.145773e-06
	8.96	2.917868e-06
	8.97	3.601705e-06
	8.98	3.054635e-06
	8.99	2.7811e-06
	9	3.100264e-06
	9.01	2.917868e-06
	9.02	3.829731e-06
	9.03	3.191402e-06
	9.04	3.784102e-06
	9.05	3.784102e-06
	9.06	4.148895e-06
	9.07	3.966499e-06
	9.08	3.87524e-06
	9.09	4.012128e-06
	9.1	4.923872e-06
	9.11	4.832733e-06
	9.12	4.012128e-06
	9.13	4.103266e-06
	9.14	4.103266e-06
	9.15	3.784102e-06
	9.16	3.692964e-06
	9.17	4.148895e-06
	9.18	4.42243e-06
	9.19	4.148895e-06
	9.2	5.151897e-06
	9.21	4.695965e-06
	9.22	4.878362e-06
	9.23	5.379803e-06
	9.24	5.471062e-06
	9.25	4.46794e-06
	9.26	4.559198e-06
	9.27	4.9695e-06
	9.28	4.46794e-06
	9.29	5.425432e-06
	9.3	4.741595e-06
	9.31	4.923872e-06
	9.32	5.653338e-06
	9.33	5.060639e-06
	9.34	5.425432e-06
	9.35	5.243036e-06
	9.36	5.607829e-06
	9.37	5.379803e-06
	9.38	4.741595e-06
	9.39	4.923872e-06
	9.4	5.379803e-06
	9.41	4.787104e-06
	9.42	5.01513e-06
	9.43	5.197407e-06
	9.44	4.741595e-06
	9.45	5.471062e-06
	9.46	5.379803e-06
	9.47	4.787104e-06
	9.48	5.01513e-06
	9.49	4.832733e-06
	9.5	4.559198e-06
	9.51	5.516571e-06
	9.52	5.607829e-06
	9.53	4.148895e-06
	9.54	6.337296e-06
	9.55	4.741595e-06
	9.56	4.194404e-06
	9.57	5.106268e-06
	9.58	5.151897e-06
	9.59	5.243036e-06
	9.6	4.559198e-06
	9.61	4.513569e-06
	9.62	5.288665e-06
	9.63	4.832733e-06
	9.64	4.103266e-06
	9.65	4.787104e-06
	9.66	3.829731e-06
	9.67	4.331172e-06
	9.68	3.920869e-06
	9.69	3.829731e-06
	9.7	3.829731e-06
	9.71	3.692964e-06
	9.72	4.240034e-06
	9.73	4.787104e-06
	9.74	4.42243e-06
	9.75	3.100264e-06
	9.76	4.194404e-06
	9.77	3.464937e-06
	9.78	3.556196e-06
	9.79	3.556196e-06
	9.8	4.012128e-06
	9.81	4.057637e-06
	9.82	4.103266e-06
	9.83	3.464937e-06
	9.84	3.738473e-06
	9.85	3.920869e-06
	9.86	4.194404e-06
	9.87	3.647334e-06
	9.88	2.917868e-06
	9.89	3.510567e-06
	9.9	2.963497e-06
	9.91	2.917868e-06
	9.92	2.325168e-06
	9.93	3.464937e-06
	9.94	2.553194e-06
	9.95	3.510567e-06
	9.96	2.826729e-06
	9.97	2.689962e-06
	9.98	2.188401e-06
	9.99	2.416306e-06
	10	2.188401e-06
	10.01	2.7811e-06
	10.02	1.732469e-06
	10.03	2.461936e-06
	10.04	2.279539e-06
	10.05	2.051633e-06
	10.06	1.869236e-06
	10.07	2.006004e-06
	10.08	1.823727e-06
	10.09	2.051633e-06
	10.1	1.641331e-06
	10.11	1.641331e-06
	10.12	1.322166e-06
	10.13	1.322166e-06
	10.14	1.550072e-06
	10.15	1.778098e-06
	10.16	1.68684e-06
	10.17	1.231028e-06
	10.18	1.458934e-06
	10.19	1.231028e-06
	10.2	1.778098e-06
	10.21	1.641331e-06
	10.22	1.139794e-06
	10.23	1.0942e-06
	10.24	1.504563e-06
	10.25	1.550072e-06
	10.26	1.231028e-06
	10.27	1.504563e-06
	10.28	1.322166e-06
	10.29	1.003014e-06
	10.3	1.276537e-06
	10.31	1.458934e-06
	10.32	1.413305e-06
	10.33	9.118395e-07
	10.34	1.0942e-06
	10.35	1.322166e-06
	10.36	1.413305e-06
	10.37	1.322166e-06
	10.38	1.367795e-06
	10.39	1.413305e-06
	10.4	1.68684e-06
	10.41	1.003014e-06
	10.42	1.0942e-06
	10.43	1.0942e-06
	10.44	7.750601e-07
	10.45	1.185387e-06
	10.46	8.662464e-07
	10.47	1.0942e-06
	10.48	1.276537e-06
	10.49	1.413305e-06
	10.5	1.276537e-06
	10.51	1.231028e-06
	10.52	1.413305e-06
	10.53	1.276537e-06
	10.54	9.118395e-07
	10.55	5.471061e-07
	10.56	1.231028e-06
	10.57	7.750601e-07
	10.58	1.048607e-06
	10.59	8.662464e-07
	10.6	7.750601e-07
	10.61	1.231028e-06
	10.62	7.294669e-07
	10.63	1.048607e-06
	10.64	1.367795e-06
	10.65	1.0942e-06
	10.66	9.574208e-07
	10.67	1.185387e-06
	10.68	1.322166e-06
	10.69	1.048607e-06
	10.7	1.413305e-06
	10.71	1.504563e-06
	10.72	1.0942e-06
	10.73	1.139794e-06
	10.74	1.595701e-06
	10.75	1.048607e-06
	10.76	1.185387e-06
	10.77	1.458934e-06
	10.78	1.185387e-06
	10.79	1.0942e-06
	10.8	1.68684e-06
	10.81	1.276537e-06
	10.82	1.595701e-06
	10.83	1.276537e-06
	10.84	1.231028e-06
	10.85	1.185387e-06
	10.86	1.276537e-06
	10.87	1.048607e-06
	10.88	1.458934e-06
	10.89	1.048607e-06
	10.9	1.732469e-06
	10.91	1.231028e-06
	10.92	1.960495e-06
	10.93	1.960495e-06
	10.94	1.732469e-06
	10.95	1.732469e-06
	10.96	1.504563e-06
	10.97	1.413305e-06
	10.98	2.051633e-06
	10.99	1.367795e-06
	11	1.595701e-06
	11.01	1.458934e-06
	11.02	2.325168e-06
	11.03	1.276537e-06
	11.04	1.185387e-06
	11.05	2.142771e-06
	11.06	2.325168e-06
	11.07	1.68684e-06
	11.08	1.231028e-06
	11.09	1.504563e-06
	11.1	1.595701e-06
	11.11	1.641331e-06
	11.12	1.185387e-06
	11.13	1.458934e-06
	11.14	1.778098e-06
	11.15	1.823727e-06
	11.16	1.914866e-06
	11.17	1.367795e-06
	11.18	1.367795e-06
	11.19	1.367795e-06
	11.2	1.595701e-06
	11.21	1.641331e-06
	11.22	2.416306e-06
	11.23	2.097262e-06
	11.24	1.458934e-06
	11.25	1.367795e-06
	11.26	1.732469e-06
	11.27	1.641331e-06
	11.28	1.823727e-06
	11.29	1.504563e-06
	11.3	1.276537e-06
	11.31	1.641331e-06
	11.32	9.118395e-07
	11.33	1.778098e-06
	11.34	1.413305e-06
	11.35	1.048607e-06
	11.36	1.732469e-06
	11.37	1.550072e-06
	11.38	1.550072e-06
	11.39	1.413305e-06
	11.4	1.595701e-06
	11.41	1.550072e-06
	11.42	1.595701e-06
	11.43	1.276537e-06
	11.44	1.595701e-06
	11.45	1.68684e-06
	11.46	1.914866e-06
	11.47	1.458934e-06
	11.48	1.732469e-06
	11.49	2.188401e-06
	11.5	1.231028e-06
	11.51	1.458934e-06
	11.52	1.276537e-06
	11.53	1.413305e-06
	11.54	1.276537e-06
	11.55	1.550072e-06
	11.56	1.778098e-06
	11.57	9.574208e-07
	11.58	1.68684e-06
	11.59	9.118395e-07
	11.6	9.574208e-07
	11.61	1.413305e-06
	11.62	1.0942e-06
	11.63	1.550072e-06
	11.64	1.504563e-06
	11.65	1.550072e-06
	11.66	1.413305e-06
	11.67	1.367795e-06
	11.68	1.185387e-06
	11.69	1.048607e-06
	11.7	1.322166e-06
	11.71	1.185387e-06
	11.72	1.0942e-06
	11.73	1.322166e-06
	11.74	1.139794e-06
	11.75	1.139794e-06
	11.76	1.185387e-06
	11.77	1.322166e-06
	11.78	9.118395e-07
	11.79	1.048607e-06
	11.8	1.869236e-06
	11.81	1.185387e-06
	11.82	1.413305e-06
	11.83	1.322166e-06
	11.84	8.206533e-07
	11.85	1.732469e-06
	11.86	1.139794e-06
	11.87	8.662464e-07
	11.88	1.139794e-06
	11.89	1.276537e-06
	11.9	1.048607e-06
	11.91	1.550072e-06
	11.92	9.574208e-07
	11.93	1.185387e-06
	11.94	1.550072e-06
	11.95	1.504563e-06
	11.96	1.048607e-06
	11.97	1.003014e-06
	11.98	1.413305e-06
	11.99	1.185387e-06
	12	1.322166e-06
	12.01	8.662464e-07
	12.02	1.504563e-06
	12.03	1.139794e-06
	12.04	9.574208e-07
	12.05	5.926873e-07
	12.06	1.139794e-06
	12.07	1.276537e-06
	12.08	1.413305e-06
	12.09	1.048607e-06
	12.1	7.750601e-07
	12.11	9.118395e-07
	12.12	1.048607e-06
	12.13	1.413305e-06
	12.14	8.662464e-07
	12.15	7.294669e-07
	12.16	9.118395e-07
	12.17	1.139794e-06
	12.18	9.574208e-07
	12.19	1.0942e-06
	12.2	1.367795e-06
	12.21	1.276537e-06
	12.22	8.662464e-07
	12.23	8.206533e-07
	12.24	9.574208e-07
	12.25	7.750601e-07
	12.26	1.550072e-06
	12.27	8.206533e-07
	12.28	1.367795e-06
	12.29	1.139794e-06
	12.3	1.68684e-06
	12.31	1.185387e-06
	12.32	1.367795e-06
	12.33	1.048607e-06
	12.34	9.118395e-07
	12.35	1.185387e-06
	12.36	8.206533e-07
	12.37	1.003014e-06
	12.38	9.574208e-07
	12.39	1.185387e-06
	12.4	1.231028e-06
	12.41	1.504563e-06
	12.42	8.662464e-07
	12.43	1.231028e-06
	12.44	1.185387e-06
	12.45	9.118395e-07
	12.46	1.139794e-06
	12.47	8.662464e-07
	12.48	1.003014e-06
	12.49	6.838737e-07
	12.5	1.003014e-06
	12.51	6.382805e-07
	12.52	1.276537e-06
	12.53	1.048607e-06
	12.54	1.0942e-06
	12.55	8.206533e-07
	12.56	1.139794e-06
	12.57	1.139794e-06
	12.58	5.01513e-07
	12.59	1.595701e-06
	12.6	1.139794e-06
	12.61	9.574208e-07
	12.62	1.003014e-06
	12.63	1.139794e-06
	12.64	1.139794e-06
	12.65	1.048607e-06
	12.66	1.276537e-06
	12.67	1.367795e-06
	12.68	1.048607e-06
	12.69	1.139794e-06
	12.7	1.276537e-06
	12.71	1.185387e-06
	12.72	8.662464e-07
	12.73	1.276537e-06
	12.74	1.139794e-06
	12.75	1.003014e-06
	12.76	9.118395e-07
	12.77	1.276537e-06
	12.78	1.231028e-06
	12.79	1.322166e-06
	12.8	1.595701e-06
	12.81	1.367795e-06
	12.82	1.367795e-06
	12.83	1.003014e-06
	12.84	1.139794e-06
	12.85	9.574208e-07
	12.86	1.276537e-06
	12.87	1.231028e-06
	12.88	1.276537e-06
	12.89	1.139794e-06
	12.9	1.0942e-06
	12.91	1.504563e-06
	12.92	8.662464e-07
	12.93	1.276537e-06
	12.94	1.0942e-06
	12.95	1.0942e-06
	12.96	1.458934e-06
	12.97	1.0942e-06
	12.98	1.185387e-06
	12.99	6.382805e-07
	13	1.367795e-06
	13.01	1.048607e-06
	13.02	1.322166e-06
	13.03	1.231028e-06
	13.04	1.139794e-06
	13.05	1.367795e-06
	13.06	1.231028e-06
	13.07	1.504563e-06
	13.08	7.294669e-07
	13.09	9.574208e-07
	13.1	9.574208e-07
	13.11	1.641331e-06
	13.12	1.0942e-06
	13.13	1.231028e-06
	13.14	1.139794e-06
	13.15	1.0942e-06
	13.16	1.550072e-06
	13.17	1.139794e-06
	13.18	7.750601e-07
	13.19	1.413305e-06
	13.2	1.003014e-06
	13.21	1.185387e-06
	13.22	1.231028e-06
	13.23	1.231028e-06
	13.24	1.003014e-06
	13.25	1.0942e-06
	13.26	1.185387e-06
	13.27	1.003014e-06
	13.28	1.185387e-06
	13.29	1.322166e-06
	13.3	1.276537e-06
	13.31	1.413305e-06
	13.32	1.231028e-06
	13.33	1.458934e-06
	13.34	1.550072e-06
	13.35	1.185387e-06
	13.36	9.574208e-07
	13.37	1.048607e-06
	13.38	1.68684e-06
	13.39	1.869236e-06
	13.4	1.550072e-06
	13.41	1.276537e-06
	13.42	1.048607e-06
	13.43	1.413305e-06
	13.44	1.68684e-06
	13.45	1.68684e-06
	13.46	9.118395e-07
	13.47	1.413305e-06
	13.48	1.139794e-06
	13.49	1.550072e-06
	13.5	1.003014e-06
	13.51	1.322166e-06
	13.52	1.185387e-06
	13.53	9.118395e-07
	13.54	1.504563e-06
	13.55	8.206533e-07
	13.56	1.003014e-06
	13.57	1.276537e-06
	13.58	1.139794e-06
	13.59	1.048607e-06
	13.6	1.550072e-06
	13.61	1.003014e-06
	13.62	1.595701e-06
	13.63	1.139794e-06
	13.64	1.139794e-06
	13.65	1.048607e-06
	13.66	1.0942e-06
	13.67	9.574208e-07
	13.68	1.139794e-06
	13.69	1.322166e-06
	13.7	1.185387e-06
	13.71	1.231028e-06
	13.72	9.574208e-07
	13.73	1.68684e-06
	13.74	1.0942e-06
	13.75	8.662464e-07
	13.76	9.118395e-07
	13.77	1.68684e-06
	13.78	1.550072e-06
	13.79	1.0942e-06
	13.8	1.550072e-06
	13.81	1.048607e-06
	13.82	1.139794e-06
	13.83	9.118395e-07
	13.84	1.367795e-06
	13.85	1.185387e-06
	13.86	9.574208e-07
	13.87	1.458934e-06
	13.88	9.118395e-07
	13.89	1.322166e-06
	13.9	1.185387e-06
	13.91	1.048607e-06
	13.92	1.458934e-06
	13.93	1.185387e-06
	13.94	1.367795e-06
	13.95	1.322166e-06
	13.96	1.185387e-06
	13.97	1.231028e-06
	13.98	1.0942e-06
	13.99	1.003014e-06
	14	1.276537e-06"""
        



def blank_strip(string):
    bstring = ""
    for c in string:
        if c not in (' ', '\t', '\n'):
            bstring += c
    return bstring

def peak_gauss(x, x0, w, A):    # w = FWHM
    # f(x=x0) = A
    # S(f)[-inf, +inf) = A*w*np.sqrt(2.0*np.pi)
    # FWHM = w = v*2.0*np.sqrt(2.0*np.ln(2.0))
    # v = (w/2.0/np.sqrt(2.0*np.log(2.0))
    # return A*(np.exp((-0.5)*(((x-x0)/v)**2)))
    return A*np.exp(-2.0*np.log(2.0)*((x-x0)/w)**2)
    
def peak_func(x, *p):
    y = 0.0; 
    n = len(p)/3
    if n != 0 : 
        for i in range(n) : 
            y += peak_gauss(x, *p[i*3:i*3+3])   # *p :== expand list p[0], p[1], p[2]
    return y

def val2idx(array, value):
    # ~ import numpy as np
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    # ~ idx = array.index(value)
    return idx  #, array[idx]
    

class _MolarMass:
    MM, MZ, MS, MN = [], [], [], []
    TT = """0;nullum;--;0.000
    1;hydrogen;H;2.000
    2;helium;He;4.002
    3;lithium;Li;6.938
    4;beryllium;Be;9.012
    5;boron;B;10.806
    6;carbon;C;12.0096
    7;nitrogen;N;14.006
    8;oxygen;O;15.999
    9;fluorine;F;18.998
    10;neon;Ne;20.1797
    11;sodium;Na;22.989
    12;magnesium;Mg;24.3050
    13;aluminium;Al;26.981
    14;silicon;Si;28.084
    15;phosphorus;P;30.973
    16;sulfur;S;32.059
    17;chlorine;Cl;35.446
    18;argon;Ar;39.9481
    19;potassium;K;39.0983
    20;calcium;Ca;40.078
    21;scandium;Sc;44.955
    22;titanium;Ti;47.867
    23;vanadium;V;50.9415
    24;chromium;Cr;51.9961
    25;manganese;Mn;54.938
    26;iron;Fe;55.845
    27;cobalt;Co;58.933
    28;nickel;Ni;58.6934
    29;copper;Cu;63.546
    30;zinc;Zn;65.38
    31;gallium;Ga;69.723
    32;germanium;Ge;72.63
    33;arsenic;As;74.921
    34;selenium;Se;78.96
    35;bromine;Br;79.904
    36;krypton;Kr;83.798
    37;rubidium;Rb;85.4678
    38;strontium;Sr;87.62
    39;yttrium;Y;88.905
    40;zirconium;Zr;91.224
    41;niobium;Nb;92.906
    42;molybdenum;Mo;95.96
    43;technetium;Tc;98.5
    44;ruthenium;Ru;101.07
    45;rhodium;Rh;102.905
    46;palladium;Pd;106.42
    47;silver;Ag;107.8682
    48;cadmium;Cd;112.411
    49;indium;In;114.818
    50;tin;Sn;118.710
    51;antimony;Sb;121.760
    52;tellurium;Te;127.60
    53;iodine;I;126.904
    54;xenon;Xe;131.293
    55;caesium;Cs;132.905
    56;barium;Ba;137.327
    57;lanthanum;La;138.905
    58;cerium;Ce;140.116
    59;praseodymium;Pr;140.907
    60;neodymium;Nd;144.242
    61;promethium;Pm;147.3
    62;samarium;Sm;150.36
    63;europium;Eu;151.964
    64;gadolinium;Gd;157.25
    65;terbium;Tb;158.925
    66;dysprosium;Dy;162.500
    67;holmium;Ho;164.930
    68;erbium;Er;167.259
    69;thulium;Tm;168.934
    70;ytterbium;Yb;173.054
    71;lutetium;Lu;174.9668
    72;hafnium;Hf;178.49
    73;tantalum;Ta;180.947
    74;tungsten;W;183.84
    75;rhenium;Re;186.207
    76;osmium;Os;190.23
    77;iridium;Ir;192.217
    78;platinum;Pt;195.084
    79;gold;Au;196.966
    80;mercury;Hg;200.59
    81;thallium;Tl;204.382
    82;lead;Pb;207.2
    83;bismuth;Bi;208.980
    84;polonium;Po;208.98
    85;astatine;At;210
    86;radon;Rn;210
    87;francium;Fr;212
    88;radium;Ra;226
    89;actinium;Ac;225
    90;thorium;Th;232.038
    91;protactinium;Pa;231.035
    92;uranium;U;238.028"""

    def __init__(self):
        for L in self.TT.splitlines():
            [Z, Name, Symb, Mmass] = L.split(';')
            self.MM.append(float(Mmass))
            self.MZ.append(int(Z))
            self.MS.append(Symb)
            self.MN.append(Name)
        
    def mm(self, atom):
        # ~ for i, s in enumerate(self.Symb):
            # ~ if s == atom:
                # ~ return self.Mm[i]
        p = self.MS.index(atom)
        MM = self.MM[p]
        # ~ print '_MolarMass.mm(',atom,')=', MM
        return MM

    
    
def scipy_optimize_minimize_test(self):
    import pylab as pl
    def f(x,*p):
        a1,a2,a3,a4,a5 = p
        return a1*np.abs(x-a2)**a3 * np.exp(-a4 * np.abs(x)**a5)

    points = 500
    xlim = 4.
    x = (np.random.rand(points) * 2. - 1.) * xlim
    x.sort()
    p0 = [1.4, -0.8, 1.1, 1.2, 2.2]
    y = f(x,*p0)
    y_noise = y + np.random.randn(points) * .05
    # mean squared error wrt. noisy data as a function of the parameters
    x = self.qq
    err = lambda p: np.mean((f(x,*p)-y_noise)**2)
    
    # bounded optimization using scipy.minimize
    p_init = [1., -1., .5, .5, 2.]
    print(p_init)
    print(p_bounds)
    p_opt = scipy.optimize.minimize(
        err, # minimize wrt to the noisy data
        p_init, 
        bounds=p_bounds,     #[(None,None),(-1,1),(None,None),(0,None),(None,None)], # set the bounds
        method="L-BFGS-B" # this method supports bounds
    ).x
    
    print(p_opt)
    # plot everything
    pl.scatter(x, y_noise, alpha=.2, label="f + noise")
    # ~ pl.plot(x, y, c='#000000', lw=2., label="f")
    pl.plot(x, f(x,*p_opt) ,'-', c='r', lw=2., label="fitted f")

    pl.xlabel("x")
    pl.ylabel("f(x)")
    pl.legend(loc="best")
    # ~ pl.xlim([-xlim*1.01,xlim*1.01])
    pl.yscale('log')

    pl.show()

def step_test():
    n = 100
    stepx = np.array(n)
    stepy = np.array(n)
    stepx = np.asarray( [-20. + x*50.0/n for x in range(n)] )
    stepy = np.array([ 5 for sx in stepx])
    stepy = stepy + step(-10., 5, 10., 3., stepx) 
    stepy = stepy + step( 0., 10, 3, 3., stepx) 
    stepy = stepy + step( 10., 3, 7, 3., stepx) 
    stepy = stepy + step( 20., 7, 0, 3., stepx) 

    fig, axs = plt.subplots(1)
    axs.plot(stepx, stepy, label='step')
    axs.grid(True)
    plt.show()
    
    
    
def fit_lmfit_test(self):
    from lmfit import minimize, Parameters, Parameter, report_fit
    import numpy as np

    xdata = np.array([177.,180.,183.,187.,189.,190.,196.,197.,201.,202.,203.,204.,206.,218.,225.,231.,234.,
          252.,262.,266.,267.,268.,277.,286.,303.])

    ydata = np.array([0.81,0.74,0.78,0.75,0.77,0.81,0.73,0.76,0.71,0.74,0.81,0.71,0.74,0.71,
          0.72,0.69,0.75,0.59,0.61,0.63,0.64,0.63,0.35,0.27,0.26])

    def fit_fc(params, x, data):

        n1 = params['n1'].value
        n2 = params['n2'].value
        n3 = params['n3'].value
        n4 = params['n4'].value
        n5 = params['n5'].value

        model = n1 + (n2 * x + n3) * 1./ (1. + np.exp(n4 * (n5 - x)))

        return model - data #that's what you want to minimize

    # create a set of Parameters
    # 'value' is the initial condition
    # 'min' and 'max' define your boundaries
    params = Parameters()
    params.add('n1', value= 0.2*1.5, min=0.2, max=0.8, vary=True)
    params.add('n2', value= -0.005*1.5, min=-0.3, max=10**(-10), vary=True)
    params.add('n3', value= 1.0766*1.5, min=-1000., max=1000., vary=True)
    params.add('n4', value= -0.36379*1.5, min=-1000., max=1000., vary=True)
    params.add('n5', value= 274.0*1.5, min=0., max=1000., vary=True)

    # do fit, here with leastsq model
    result = minimize(fit_fc, params, args=(xdata, ydata), method='leastsq')

    # write error report
    report_fit(params)

    xplot = np.linspace(min(xdata), max(xdata), 1000)
    yplot = result.params['n1'] + (result.params['n2'] * xplot + result.params['n3']) * \
                                  1./ (1. + np.exp(result.params['n4'] * (result.params['n5'] - xplot)))
    #plot results
    try:
        import pylab
        pylab.plot(xdata, ydata, 'k+')
        pylab.plot(xplot, yplot, 'r')
        pylab.show()
    except:
        pass
        
def flatten(self):
    
    # ~ Given a list of lists l,

    flat_list = [item for sublist in l for item in sublist]

    # ~ which means:

    flat_list = []
    for sublist in l:
        for item in sublist:
            flat_list.append(item)

    # ~ is faster than the shortcuts posted so far. (l is the list to flatten.)

    # ~ Here is the corresponding function:

    flatten = lambda l: [item for sublist in l for item in sublist]

    # ~ As evidence, you can use the timeit module in the standard library:

    # ~ $ python -mtimeit -s'l=[[1,2,3],[4,5,6], [7], [8,9]]*99' '[item for sublist in l for item in sublist]'
    # ~ 10000 loops, best of 3: 143 usec per loop
    # ~ $ python -mtimeit -s'l=[[1,2,3],[4,5,6], [7], [8,9]]*99' 'sum(l, [])'
    # ~ 1000 loops, best of 3: 969 usec per loop
    # ~ $ python -mtimeit -s'l=[[1,2,3],[4,5,6], [7], [8,9]]*99' 'reduce(lambda x,y: x+y,l)'
    # ~ 1000 loops, best of 3: 1.1 msec per loop
    
def scipy_optimize_differential_evolution_test():
    '''
    This is kinda straightforward solution which shows the idea, 
    also code isn`t very pythonic but for simplicity i think its 
    good enough. Ok as example we want to fit equation of a kind 
    y = ax^2 + bx + c 
    to a data obtained from equation 
    y = x^2. 
    It obvious that parameter a = 1 and b,c should equal to 0. 
    Since differential evolution algorithm finds minimum of a 
    function we want to find a minimum of a root mean square deviation 
    (again, for simplicity) of analytic solution of general equation 
    (y = ax^2 + bx + c) with given parameters (providing some initial 
    guess) vs "experimental" data. So, to the code:
    '''
    
        
    from scipy.optimize import differential_evolution
    def func(parameters, *data):
        #we have 3 parameters which will be passed as parameters and
        #"experimental" x,y which will be passed as data
        a,b,c = parameters
        x,y = data
        result = 0
        for i in range(len(x)):
            result += (a*x[i]**2 + b*x[i]+ c - y[i])**2
        return result**0.5

    #initial guess for variation of parameters
    #             a            b            c
    bounds = [(1.5, 0.5), (-0.3, 0.3), (0.1, -0.1)]
    #producing "experimental" data 
    x = [i for i in range(100)]
    y = [v**2 for v in x]
    #packing "experimental" data into args
    args = (x,y)
    result = differential_evolution(func, bounds, args=args)
    print(result.x)
    
    
if __name__ == "__main__":
   
    Xrr = XRR()
    Xrr.fitmax = 100
    Xrr.fit()
    Xrr.layer_print()
    Xrr.layer_profile()
    Xrr.fitmax = 0
    Xrr.plot()
    
    
   
    
