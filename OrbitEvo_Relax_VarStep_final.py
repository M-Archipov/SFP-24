import SPMSys_final as SPMSys
import SatelliteBoundaries_final as SB
import numpy as np
import matplotlib.pyplot as plt

class OrbitEvo(SPMSys.SelfOtherSystem):

    def DisturbFuncM(self,
                rotSpeedSelf, norb, a, ecc):
        #independent on m,p,q:
        independent = norb**2 * a**2 * self.MOther / self.MSelf * (self.RSelf / a)**5

        #lmpq = 220q
        RM_220q = - 3/4 * (2 + self.qVector) * self.EccFuncSq(ecc=ecc) * self.QualFunc(rotSpeedSelf, norb)

        #lmpq = 201q
        RM_201q = - 1/4 * self.qVector * self.EccFuncSq(p=1, ecc=ecc) * self.QualFunc(rotSpeedSelf, norb, m=0, p=1)

        result = independent * np.sum(RM_220q + RM_201q)
        return result 

    def DisturbFuncOmega(self,
                     rotSpeedSelf, norb, a, ecc):
        #independent on m,p,q:
        independent = norb**2 * a**2 * self.MOther / self.MSelf * (self.RSelf / a)**5

        #lmpq = 220q
        RM_220q = - 3/2 * self.EccFuncSq(ecc=ecc) * self.QualFunc(rotSpeedSelf, norb)

        result = independent * np.sum(RM_220q)
        return result 
    
    # should not use self.a and self.ecc - RK4 requires intermediate values
    def SMAxisEvoFunc(self,
                      rotSpeedSelf, norb, a, ecc):
        return 2/(norb * a) * self.DisturbFuncM(rotSpeedSelf, norb, a, ecc)

    def EccEvoFunc(self,
                   rotSpeedSelf, norb, a, ecc):
        output = ((1-ecc**2)/(norb * a**2 * ecc) * self.DisturbFuncM(rotSpeedSelf, norb, a, ecc)) - (
            (np.sqrt(1-ecc**2))/(norb * a**2 * ecc) * self.DisturbFuncOmega(rotSpeedSelf, norb, a, ecc))
        return output
    
    def ComboFunc(self,
                  rotSpeedSelf, a, ecc, regime='despin'):
        
        if regime=='despin':
            rotSpeedSelf_temp = rotSpeedSelf
            a_temp   = a
            ecc_temp = ecc
        elif regime=='relax':
            rotSpeedSelf_temp = rotSpeedSelf
            a_temp   = self.a
            ecc_temp = self.ecc
        else:
            rotSpeedSelf_temp = self.rotSpeedSelf
            a_temp   = a
            ecc_temp = ecc

        norb = np.sqrt(self.G * (self.MSelf + self.MOther)/a_temp**3)
        y = np.array((self.ImpulseTheorem(rotSpeedSelf_temp, norb, a_temp, ecc_temp), 
                      self.SMAxisEvoFunc(rotSpeedSelf_temp, norb, a_temp, ecc_temp), 
                      self.EccEvoFunc(rotSpeedSelf_temp, norb, a_temp, ecc_temp)))
        return y

#step, time in years
#longstep = longstepModifier * step
def MakeTwoBodyRotSpeedEvoGraph(MainBody, MainBodySORatio = 2, RevolvingBodySORatio = 2,
                                step = 1000, time = 1e9, longstepModifier = 25, err=1e-6):

    step = step * 86400 * 365   #convert step to seconds
    RevolvingBody = OrbitEvo.SwitchBodies(MainBody)

    def RKSolverVec(func, beginState, step = step, regime='despin'):
        k1 = func(*beginState, regime=regime)
        k2 = func(*(beginState + 0.5 * step * k1), regime=regime)
        k3 = func(*(beginState + 0.5 * step * k2), regime=regime)
        k4 = func(*(beginState + step * k3), regime=regime)

        end = step * (1/6 * (k1 + k4) + 1/3 * (k2 + k3))
        return end

    beginState = np.array((MainBodySORatio * MainBody.MeanMotion(), RevolvingBodySORatio * RevolvingBody.MeanMotion(), MainBody.a, MainBody.ecc))
    print(beginState)

    y0,y1,y2,y3 = [[0, beginState[0] / MainBody.MeanMotion()]],[[0, beginState[1] / MainBody.MeanMotion()]],[[0, beginState[2]]],[[0, beginState[3]]]
    # y0, y1, y2, y3 = np.empty((time,2)), np.empty((time,2)),np.empty((time,2)),np.empty((time,2))

    y0_aux = [[0, beginState[0] / MainBody.MeanMotion()]]

    # x = np.arange(time) * (step / (86400 * 365))

    f = open('vals.txt','w')
    f.write('')
    f.close()
    f = open('vals.txt','a')

    records = 2000
    logTime = 0

    initial_despinning = True
    last_step_small = True
    try:
        while y0[-1][0] < time:
            longstep = 0
            if len(y0) > 1:
                print(np.abs(1 - y0_aux[-1][1]/y0_aux[-2][1]))
            if len(y0) > 1 and (np.abs(1 - y0_aux[-1][1]/y0_aux[-2][1]) < err) and last_step_small:

                print('TIMESKIP')
                beginState[0] = y0_aux[-1][1] * MainBody.MeanMotion()

                # Relaxed spin state => take a long step
                longstep = step * longstepModifier
                beginStatePM = RKSolverVec(MainBody.ComboFunc, np.array((beginState[0], beginState[2], beginState[3])), step = longstep, regime='longstep')
                beginStateMP = RKSolverVec(RevolvingBody.ComboFunc, np.array((beginState[1], beginState[2], beginState[3])), step = longstep, regime='longstep')
                beginStatePM[0] = 0

                beginState[0] += beginStatePM[0]
                beginState[1] += beginStateMP[0]
                beginState[2] += beginStatePM[1] + beginStateMP[1]
                beginState[3] += beginStatePM[2] + beginStateMP[2]

                MainBody.rotSpeedSelf, RevolvingBody.rotSpeedSelf = beginState[0], beginState[1]
                MainBody.a, RevolvingBody.a = beginState[2], beginState[2]
                MainBody.ecc, RevolvingBody.ecc = beginState[3], beginState[3]
                beginState[0] = y0_aux[-1][1] * MainBody.MeanMotion()

                y0.append([y0[-1][0] + (longstep / 86400 / 365), beginState[0] / MainBody.MeanMotion()])
                y1.append([y1[-1][0] + (longstep / 86400 / 365), beginState[1] / RevolvingBody.MeanMotion()])
                y2.append([y2[-1][0] + (longstep / 86400 / 365),  beginState[2]])
                y3.append([y3[-1][0] + (longstep / 86400 / 365),  beginState[3]])

                y0_aux.append([y0[-1][0], beginState[0] / MainBody.MeanMotion()])

                initial_despinning = False
                last_step_small = False
            else:
                if initial_despinning:
                    print('INITIAL DESPINNING')

                    # Unrelaxed spin state => take a small step
                    beginStatePM = RKSolverVec(MainBody.ComboFunc, np.array((beginState[0], beginState[2], beginState[3])), regime='despin')
                    beginStateMP = RKSolverVec(RevolvingBody.ComboFunc, np.array((beginState[1], beginState[2], beginState[3])), regime='despin')

                    last_step_small = True

                    beginState[0] += beginStatePM[0]
                    beginState[1] += beginStateMP[0]
                    beginState[2] += beginStatePM[1] + beginStateMP[1]
                    beginState[3] += beginStatePM[2] + beginStateMP[2]

                    MainBody.rotSpeedSelf, RevolvingBody.rotSpeedSelf = beginState[0], beginState[1]
                    MainBody.a, RevolvingBody.a = beginState[2], beginState[2]
                    MainBody.ecc, RevolvingBody.ecc = beginState[3], beginState[3]
                    y0.append([y0[-1][0] + (step / 86400 / 365), beginState[0] / MainBody.MeanMotion()])
                    y1.append([y1[-1][0] + (step / 86400 / 365), beginState[1] / RevolvingBody.MeanMotion()])
                    y2.append([y2[-1][0] + (step / 86400 / 365),  beginState[2]])
                    y3.append([y3[-1][0] + (step / 86400 / 365),  beginState[3]])

                    y0_aux.append([y0[-1][0], y0[-1][1]])

                else:
                    print('RELAXATION')

                    # Unrelaxed spin state => take a small step
                    beginStatePM = RKSolverVec(MainBody.ComboFunc, np.array((beginState[0], beginState[2], beginState[3])), regime='relax')
                    beginStateMP = RKSolverVec(RevolvingBody.ComboFunc, np.array((beginState[1], beginState[2], beginState[3])), regime='relax')

                    last_step_small = True

                    beginState[0] += beginStatePM[0]
                    #beginState[1] += beginStateMP[0]

                    y0_aux.append([y0_aux[-1][0], beginState[0] / MainBody.MeanMotion()])

            print(len(y0)-1, ' : ', y0[-1][1],y1[-1][1],y2[-1][1],y3[-1][1]) #PMSys.MeanMotion())
            print('years left:', time - y0[-1][0])

            # if y1[-1][1]<0.98:
            #     print("ERROR")
            #     break


            if y0[-1][0] > logTime:
                string = str(y2[-1][1]) + ',' + str(y3[-1][1]) + ',' + str(y0[-1][0]) + '\n'
                f.write(string)
                logTime += time/records
    except KeyboardInterrupt:
        pass

    f.close()

    #graph
    plt.rcParams['text.usetex'] = True

    fig, ax1 = plt.subplots()
    ax1.plot([y[0] for y in y0],[y[1] for y in y0],'#ff0000', label='rotSpeedPlanet')
    ax1.plot([y[0] for y in y1],[y[1] for y in y1],'#ffaaaa', label='rotSpeedMoon', marker = '.')
    ax1.set_xlabel(r"$t \ (yrs)$")
    ax1.set_ylabel(r"$\dot{\theta}$", color='r')
    ax1.tick_params('y', colors='r')

    ax2 = ax1.twinx()
    ax2.plot([y[0] for y in y2],[y[1] for y in y2],'g', label='a')
    ax2.set_ylabel('a', color='g')
    # ax2.tick_params('y', colors='g')

    ax3 = ax1.twinx()
    ax3.plot([y[0] for y in y3],[y[1] for y in y3],'b', label='ecc')
    ax3.set_ylabel('e', color='b')
    # ax3.tick_params('y', colors='b')
    ax3.spines['right'].set_position(('outward', 60))
    current_values = plt.gca().get_yticks()
    plt.gca().set_yticklabels(['{:.2e}'.format(x) for x in current_values])

    fig.tight_layout()
    plt.grid()
    plt.show()
    plt.rcParams['text.usetex'] = False

SPSys = OrbitEvo.ExoplanetSys(SPMSys.exopName, viscosityExop = SPMSys.viscosityExop)
SEESys = SB.SatelliteBoundaries.ExoplanetSys(SPMSys.exopName, MExom=SPSys.MSelf * SPMSys.MMoonRatio, eccExom= 0.2)

# PMSys = OrbitEvo(MSelf=PSSys.MSelf, RSelf=PSSys.RSelf, stiffnessSelf=PSSys.stiffnessSelf, viscositySelf=PSSys.viscositySelf, MOther=PSSys.MSelf * SPMSys.MMoonRatio, a = SPMSys.aLuna / SPMSys.REarth * PSSys.RSelf, ecc = SPMSys.eccLuna)
PMSys = OrbitEvo(MSelf=SPSys.MSelf, RSelf=SPSys.RSelf, stiffnessSelf=SPSys.stiffnessSelf, viscositySelf=SPSys.viscositySelf, MOther=SEESys.MMoon, a = SEESys.RedHillSphere()*0.90, ecc = SEESys.eccMoon)

MakeTwoBodyRotSpeedEvoGraph(PMSys, step = 100, longstepModifier= 100, time = 5e8, err = 1e-5)
