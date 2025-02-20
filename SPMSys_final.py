import numpy as np
from sympy import symbols, O, expand
import matplotlib.pyplot as plt
import matplotlib as mpl
import reader

#Earth atrributes
MEarth = 5.972 * 10**24                #mass (kg)
REarth = 6.378 * 10**6                 #radius (m)
stiffnessEarth = 10**11                #stiffness (Pa)
viscosityEarth = 10**17                #viscosity (Pa*s)
densityEarth = 5514                    #density (kg/(m**3))

#Moon attributes
MLuna = 0.012 * 5.972 * 10**24         #mass (kg)
RLuna = 1.7374 * 10**6                 #radius (m)
stiffnessLuna = 6 * 10**10             #stiffness (Pa)
viscosityLuna = 10**17                 #viscosity (Pa*s)
densityLuna = 3340                     #density (kg/(m**3))

#Sun attributes
MSun = 2 * 10**30                      #mass (kg)

#Moon orbit attributes
aLuna = 384399000                      #major axis (m)
eccLuna = 0.055                        #eccentricity

#Earth orbit attributes
aEarth = 1 * 1.495978707 * 10**11      #major axis (m)
eccEarth = 0.017                       #eccentricity

class SelfOtherSystem():
    G = 6.6743 * 10**(-11)             #gravitational konstant (N*(m*2)/(kg**2))
    eccSymbol = symbols('ecc')

    def __init__(self,
                 MSelf=MEarth, RSelf = REarth, stiffnessSelf = stiffnessEarth, viscositySelf = viscosityEarth,
                 MOther = MLuna,
                 a = aLuna, ecc = eccLuna,
                 deg = 1, pow = 2):
        self.MSelf = MSelf
        self.RSelf = RSelf
        self.stiffnessSelf = stiffnessSelf
        self.viscositySelf = viscositySelf
        self.densitySelf = MSelf/(4/3 * np.pi * RSelf**3)

        self.MOther = MOther

        self.a = a
        self.ecc = ecc

        self.deg = deg
        self.pow = pow

        # self.rotSpeedSelf  = 1e-6

        self.qVector = np.arange(-1 * self.deg, self.deg + 1,1) #DisturbFuncM - values of (l - 2p + q)    

    def MeanMotion(self):
        return np.sqrt(self.G * (self.MSelf + self.MOther)/self.a**3)
    
    def EccFuncSq(self,
                  p = 0, ecc=0):

        G_20q = np.array((
        # e^0, e^1, e^2, e^3, e^4
        [0, 0, 0, 0, 0],            #q = -2
        [0, -1/2, 0, 1/16, 0],      #q = -1
        [1, 0, -5/2, 0, 13/16],     #q =  0
        [0, 7/2, 0, -123/16, 0],    #q =  1
        [0, 0, 17/2, 0, -115/6]))    #q =  2

        G_21q = np.array((
        # e^0, e^1, e^2, e^3, e^4
        [0, 0, 9/4, 0, 7/4],            #q = -2
        [0, 3/2, 0, 27/16, 0],      #q = -1
        [1, 0, 3/2, 0, 15/8],     #q =  0
        [0, 3/2, 0, 27/16, 0],    #q =  1
        [0, 0, 9/4, 0, 7/4]))    #q =  2

        q_0 = 2 #row with coefs for q = 0

        G2 = []
        for i in range(q_0 - self.deg, q_0 + self.deg + 1):
            if p == 0:
                G2.append(self.BinomialExpansion(G_20q[i]).subs(self.eccSymbol, ecc))
            if p == 1:
                G2.append(self.BinomialExpansion(G_21q[i]).subs(self.eccSymbol, ecc))
    
        return np.float64(np.array(G2))
    
    def BinomialExpansion(self, coefs):
        coefsSum = 0
        for i in range(len(coefs)):
            coefsSum += coefs[i] * self.eccSymbol**i

        return (expand(coefsSum**2)+O(self.eccSymbol**(self.pow+1))).removeO()
    
    def QualFunc(self,
                 rotSpeedSelf, norb,
                 l = 2, m = 2, p = 0,
                 returnLoveNum = False):
        imLoveNum = []
        for q in range(-1* self.deg, self.deg + 1):
            fouTidMode = (l - 2 * p + q) * norb - m * rotSpeedSelf    #dělá to samé co funkce FouTidMode z DU2
            if self.viscositySelf == 0 or fouTidMode == 0:
                #imLoveNum.append(3/2)
                imLoveNum.append(0)
            else:
                compliance = 1/self.stiffnessSelf - 1j/(self.viscositySelf * fouTidMode)
                loveNum = 3/2 * 1/(1 + 19/(2 * compliance * self.densitySelf**2 * self.RSelf**2 * self.G * 4*np.pi/3))
                if returnLoveNum == True:
                    return loveNum
                imLoveNum.append(-np.imag(loveNum))
                #print(compliance, loveNum, imLoveNum[-1])
        
        return imLoveNum
    
    #indexy l,m,p = 2,2,0
    #3/2 = F_220{i=0}^2 (inclination function)
    def MeanTorque(self,
                   rotSpeedSelf, norb, a, ecc):
        T = self.G * self.MOther**2 / a * (self.RSelf / a)**5 * 3/2 * np.sum(self.EccFuncSq(ecc=ecc) * self.QualFunc(rotSpeedSelf, norb))
        return T
    
    def ImpulseTheorem(self,
                       rotSpeedSelf, norb, a, ecc):
        C = 2 / 5 * self.MSelf * self.RSelf**2
        T = self.MeanTorque(rotSpeedSelf, norb, a, ecc)
        return T / C

    @classmethod
    def SwitchBodies(cls,
                     Sys,
                     ROther = RLuna, stiffnessOther = stiffnessLuna, viscosityOther = viscosityLuna):
        return cls(MSelf=Sys.MOther, RSelf=ROther, stiffnessSelf=stiffnessOther, viscositySelf=viscosityOther, MOther=Sys.MSelf, a=Sys.a, ecc=Sys.ecc, deg=Sys.deg, pow=Sys.pow)

    @classmethod
    def ExoplanetSys(cls,
                     name,
                     stiffnessExop = stiffnessEarth, viscosityExop = viscosityEarth,
                     deg = 1, pow = 2):
        RExop = float(reader.radii[np.where(reader.names == name)] * REarth)
        MExop = float(reader.masses[np.where(reader.names == name)] * MEarth)
        eccExop = float(reader.eorbs[np.where(reader.names == name)])
        aExop = float(reader.aorbs[np.where(reader.names == name)] * 1.496 * 10**11)

        MExos = float(reader.mstars[np.where(reader.names == name)] * MSun)
        
        return cls(MSelf = MExop, RSelf = RExop, stiffnessSelf = stiffnessExop, viscositySelf = viscosityExop, ecc = eccExop, a = aExop, MOther = MExos, deg = deg, pow = pow)
    
    @classmethod
    def Exoplanet(cls,
                  name,
                  MExom = MLuna,
                  aExom = aLuna, eccExom = eccLuna,
                  stiffnessExop = stiffnessEarth, viscosityExop = viscosityEarth,
                  deg = 1, pow = 2):
        RExop = float(reader.radii[np.where(reader.names == name)] * REarth)
        MExop = float(reader.masses[np.where(reader.names == name)] * MEarth)
        # eccExop = float(reader.eorbs[np.where(reader.names == name)])
        # aExop = float(reader.aorbs[np.where(reader.names == name)] * 1.496 * 10**11)

        # MExos = float(reader.mstars[np.where(reader.names == name)] * MSun)
        
        return cls(MSelf = MExop, RSelf = RExop, stiffnessSelf = stiffnessExop, viscositySelf = viscosityExop, ecc = eccExom, a = aExom, MOther = MExom, deg = deg, pow = pow)

def MakeSecTorqueMMoonRatioGraph(Sys,
                                 xmin = 0, xmax = 2, steps = 500,
                                 measurements = 5, MMoonRatioMin = 0.01, MMoonRatioMax = 0.1):
    plt.rcParams['text.usetex'] = True

    fig, ax = plt.subplots()
    cmap = mpl.colormaps['viridis_r']
    x = np.linspace(xmin,xmax, num = steps)
    plt.xlim(xmin,xmax)
    y1 = [None] * len(x)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm = plt.Normalize(MMoonRatioMin, MMoonRatioMax))
    cbar = fig.colorbar(sm, ax=ax, ticks = np.linspace(MMoonRatioMin, MMoonRatioMax,measurements))

    plt.xlabel(r'$r$')
    plt.ylabel(r'$\langle \mathcal{T} \rangle $')

    cbar.set_label(r'$m_m/m_p$')

    for i in range(measurements):
        Sys.MOther = Sys.MSelf * (MMoonRatioMin + (i * (MMoonRatioMax - MMoonRatioMin)/(measurements - 1)))

        for j in range(len(x)):
            y1[j] = Sys.MeanTorque(x[j] * Sys.MeanMotion(), Sys.MeanMotion(), Sys.a, Sys.ecc)
        print(i+1,'/',measurements)

        ax.plot(x,y1,'-',c=cmap(i/measurements))

    ax.grid()
    # plt.legend()
    plt.show()

    plt.rcParams['text.usetex'] = False

studied = ['HD 100777 b', 'GJ 3470 b', 'LHS 1815 b', 'GJ 1252 b', 'HIP 39017 b', 'Kepler-22 b', 'HD 88986 b']

#exoplanet attributes
exopName = studied[5]                  #studied exoplanet name
viscosityExop = 10**17                 #exoplanet viscosity (Pa*s)

MMoonRatio = 0.1                       #exomoon/exoplanet ratio
RExom = RLuna                          #exomoon radius (m)
aExom = aLuna                          #exomoon orbit major axis (m) 
eccExom = eccLuna                      #exomoon orbit eccentricity

Exop = SelfOtherSystem.ExoplanetSys(exopName, viscosityExop = viscosityExop)

# MakeSecTorqueMMoonRatioGraph(Sys = Exop, measurements=5)
