import SPMSys_final as SPMSys
import reader
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

class SatelliteBoundaries():

    def __init__(self,
                 MPlanet = SPMSys.MEarth, RPlanet = SPMSys.REarth, # stiffness = SPMSys.stiffnessEarth, viscosity = SPMSys.viscosityEarth, density = SPMSys.densityEarth,
                 MMoon = SPMSys.MLuna, RMoon = SPMSys.RLuna,
                 MStar = SPMSys.MSun,
                 aMoon = SPMSys.aLuna, eccMoon = SPMSys.eccLuna,
                 aPlanet = SPMSys.aEarth, eccPlanet = SPMSys.eccEarth):

        self.MPlanet = MPlanet
        self.RPlanet = RPlanet
        # self.stiffness = stiffness
        # self.viscosity = viscosity
        # self.density = density

        self.MMoon = MMoon
        self.RMoon = RMoon

        self.MStar = MStar

        self.aMoon = aMoon
        self.eccMoon = eccMoon

        self.aPlanet = aPlanet
        self.eccPlanet = eccPlanet
    
    def HillSphere(self):
        return self.aPlanet * (1 - self.eccPlanet) * (self.MPlanet / (3 * self.MStar))**(1/3)

    #reduced Hill sphere for a circular orbit (eccMoon = 0):
    def RedHillSphere(self,
                      retrograde = False):
        if retrograde == 0:
            B = 0.49
            b = 1.03
        else:
            B = 0.93
            b = 1.08

        return B * self.aPlanet * (1 - b * self.eccPlanet) * (self.MPlanet / (3 * self.MStar))**(1/3)
    
    def RocheRad(self,
                 A = 2.2):
        return A * self.RMoon * (self.MPlanet / self.MMoon)**(1/3)
    
    def maxEcc(self):
        return min(1 - self.RocheRad() / self.aMoon, self.RedHillSphere() / self.aMoon - 1)
    
    @classmethod
    def ExoplanetSys(cls,
                     name,
                     MExom = SPMSys.MLuna, RExom = SPMSys.RLuna,
                     aExom = SPMSys.aLuna, eccExom = SPMSys.eccLuna):
        RExop = float(reader.radii[np.where(reader.names == name)] * SPMSys.REarth)
        MExop = float(reader.masses[np.where(reader.names == name)] * SPMSys.MEarth)
        eccExop = float(reader.eorbs[np.where(reader.names == name)])
        aExop = float(reader.aorbs[np.where(reader.names == name)] * 1.496 * 10**11)

        MExos = float(reader.mstars[np.where(reader.names == name)] * SPMSys.MSun)
        
        return cls(MPlanet=MExop, RPlanet=RExop, MMoon=MExom, RMoon=RExom, MStar=MExos, aMoon=aExom, eccMoon=eccExom, aPlanet=aExop,  eccPlanet=eccExop)

def MakeRedHillGraph():
    noOfPlanets = len(reader.radii)
    print('noOfPlanets: ',noOfPlanets)
    vals = np.empty((noOfPlanets, 2))

    j = 0
    
    for i in range(noOfPlanets):
        Sys = SatelliteBoundaries(MPlanet=reader.masses[i] * SPMSys.MEarth, MStar=reader.mstars[i] * SPMSys.MSun, aPlanet=reader.aorbs[i] * 1.496 * 10**11, eccPlanet=reader.eorbs[i])
        vals[i,:] = [reader.radii[i] * SPMSys.REarth, Sys.RedHillSphere()]

        if vals[i,1] / (7 * 10**7) <=5:
            j+=1
        print(vals[i,:])
    
    print('no. of points: ', j)

    plt.plot(vals[:,0] / (7 * 10**7),vals[:,1] / (7 * 10**7),'.')
    plt.plot([0,5],[0,5],'-')

    plt.rcParams['text.usetex'] = True
    plt.xlabel(r'$R_p \ (r_{jup})$')
    plt.ylabel(r'$r´_H \ (r_{jup})$')

    plt.axis('square')
    plt.xlim(0,5)
    plt.ylim(0,5)

    plt.grid()
    plt.show()
    plt.rcParams['text.usetex'] = False

def MakeMaxEccGraph(Sys):

    amin = Sys.RocheRad()
    amax = Sys.RedHillSphere()

    x = np.linspace(amin, amax, num = 100)
    y = np.empty(len(x))

    for i in range(len(x)):
        Sys.aMoon = x[i]
        y[i] = Sys.maxEcc()

    plt.plot(x,y,'.')

    plt.xlabel('a')
    plt.ylabel('e_max')

    plt.grid()
    plt.show()

#max ecc graf, ale se započítáním vývoje dráhy v čase
#TODO závislost této f-ce na pozorovaném systému
def MakeDynMaxEccGraph(exopName, MMoonRatio):
    f = open('vals.txt')
    lines = f.readlines()
    cmap = mpl.colormaps['coolwarm']

    Sys = SatelliteBoundaries.ExoplanetSys(exopName, MExom= MMoonRatio * reader.masses[np.where(reader.names == exopName)] * SPMSys.MEarth)

    # print(Sys.aPlanet)

    fig, ax = plt.subplots()
    plt.rcParams['text.usetex'] = True

    sm = plt.cm.ScalarMappable(cmap = cmap, norm = plt.Normalize(lines[0].split(',')[2], lines[-1].split(',')[2]))
    cbar = fig.colorbar(sm, ax=ax, ticks = (np.linspace(float(lines[0].split(',')[2]), float(lines[-1].split(',')[2]), 10)))
    cbar.set_label(r'$\dot \theta$')

    amin = Sys.RocheRad()
    amax = Sys.RedHillSphere()
    print(amin, amax)
    x = np.linspace(amin, amax, num = 100)
    y = np.empty(len(x))

    for j in range(len(x)):
        Sys.aMoon = x[j]
        y[j] = Sys.maxEcc()

    plt.plot(x,y,'-', color = 'g')

    for i in range(len(lines)):
        # Sys = SatelliteBoundaries(aMoon=(lines[i].split(','))[0], eccMoon=(lines[i].split(',')[1])

        # Sys = SatelliteBoundaries(aPlanet=float(lines[i].split(',')[0]), eccPlanet=float(lines[i].split(',')[1]))
        print('a:',Sys.aPlanet,'ecc:',Sys.eccPlanet)

        # x = np.linspace(amin, amax, num = 100)
        # y = np.empty(len(x))

        # for j in range(len(x)):
        #     Sys.aMoon = x[j]
        #     y[j] = Sys.maxEcc()

        plt.plot(float(lines[i].split(',')[0]), float(lines[i].split(',')[1]), 'o', c=cmap(float(lines[i].split(',')[2])/(float(lines[-1].split(',')[2]) - float(lines[0].split(',')[2]))))

        print(i+1,'/',len(lines))

    plt.xlabel('a')
    plt.ylabel('e_max')

    plt.grid()
    plt.show()
    plt.rcParams['text.usetex'] = False

def MakeHillRocheGraph():
    noOfPlanets = len(reader.radii)
    print('noOfPlanets: ',noOfPlanets)
    vals = np.empty((noOfPlanets, 2))

    PToMMass = 0.01

    for i in range(noOfPlanets):
        Sys = SatelliteBoundaries(MPlanet=reader.masses[i] * SPMSys.MEarth, MStar=reader.mstars[i] * SPMSys.MSun, aPlanet=reader.aorbs[i] * 1.496 * 10**11, eccPlanet=reader.eorbs[i], MMoon=reader.masses[i]*PToMMass, RMoon= (3/4/np.pi * reader.masses[i] * PToMMass / SPMSys.densityLuna)**(1/3))
        vals[i,:] = [Sys.RocheRad(), Sys.RedHillSphere()]

    plt.plot(vals[:,0] / (7 * 10**7),vals[:,1] / (7 * 10**7),'.')

    plt.plot([0,5],[0,5],'-')

    plt.rcParams['text.usetex'] = True
    plt.xlabel(r'$r_R \ (r_{jup})$')
    plt.ylabel(r'$r´_H \ (r_{jup})$')

    plt.axis('square')

    plt.ylim(0,5)
    plt.xlim(0,5)
    plt.grid()
    plt.show()

    plt.rcParams['text.usetex'] = False

# MakeDynMaxEccGraph(SPMSys.exopName, SPMSys.MMoonRatio)

# MakeDynMaxEccGraph(SPMSys.OrbitEvo(), SPMSys.RLuna/SPMSys.REarth)


# MakeHillRocheGraph()
# MakeRedHillGraph()
# MakeMaxEccGraph(SELSys)
