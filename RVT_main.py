from __future__ import print_function
import numpy as np
import sympy
import sys
import math

if sys.platform.lower().startswith("win", 0, 3):
    sys.path.append("C:\\Users\\ADMIN\\Documents\\FreeFlyer\\FreeFlyer 7.6.0.54542 (64-Bit)\\Runtime API\\python\\examples\\src")
    sys.path.append("C:\\Users\\ADMIN\\Documents\\FreeFlyer\\FreeFlyer 7.6.0.54542 (64-Bit)\\Runtime API\\python\\src")
else:
    sys.path.append("../../../src")
    sys.path.append("../../src")

# aisolutions.freeflyer.runtimeapi can be found in "Runtime API" folder
# run a search for "Runtime API" if needed. Add path to library.
try:
    from ExampleUtilities.ExampleUtilities import ExampleUtilities
    from array import array
    from aisolutions.freeflyer.runtimeapi.RuntimeApiEngine import RuntimeApiEngine
    from aisolutions.freeflyer.runtimeapi.RuntimeApiEngine import FFTimeSpan
    from aisolutions.freeflyer.runtimeapi.ConsoleOutputProcessingMethod import ConsoleOutputProcessingMethod
    from aisolutions.freeflyer.runtimeapi.RuntimeApiException import RuntimeApiException
except ImportError:
    print("Import error!")

# constants
rE = 6378.137
mu = 398600.436
J2 = 0.0010826
ecc = 0.000755155
#alt_init = 6378.137+550
#inc_init = 80
raan = 302.271831900499               # RAAN of LTAN 13:30 for SSO
argprg = 90            # Argument of Perigee
ta = 0                 # True Anomaly
propType = "J2Mean" #Ephemeris, Norad, TwoBody, RK45, J2Mean, RK78, RK89, Cowell
stepSize = 10          # F.F Scenario Step Size
sensorHeight = 10

#Repeat Ground Track Method
def find_revpday(alt_rev, inc):
    # alt_rev = alt0 +rE
    # J2 effects
    i = 0
    while i < 10000:
        n = np.sqrt(mu/alt_rev)/alt_rev  # Mean Motion

        f = 3 * n * J2 * rE * rE / (2 * alt_rev ** 2 * (1 - ecc ** 2) ** 2)
        c = np.cos(np.deg2rad(inc))
        s = np.sin(np.deg2rad(inc))
        wDot = -f * (2.5 * s ** 2 - 2)
        WDot = -f * c
        MDot = -f * np.sqrt(1 - ecc ** 2) * (1.5 * s ** 2 - 1)

        # satellite orbit rate with J2 perturbation
        wS = n + wDot + MDot

        # earth rotation rate with respect to orbit plane that drifts at WDot
        wE = (2 * np.pi / 86400) + 2 * np.pi/365.25/86400 - WDot
        rat = wS / wE
        k = sympy.nsimplify(rat, tolerance=  0.00106)
        numer = k.numerator()
        denom = k.denominator()
        if 5 < denom < 31:
            i = 10000
        elif i%2 == 1:
            alt_rev = alt_rev + i/10
            i = i + 1
        else:
            alt_rev = alt_rev - i/10
            i = i + 1
    revperday = numer/denom
    SMA = round(alt_rev, 3)
    return SMA


def RGT(revperday, inc):
    n = revperday*(2 * np.pi / 86400)
    c = np.cos(np.deg2rad(inc))
    s = np.sin(np.deg2rad(inc))
    SMA = (mu**(1/3))*((1/n)**(2/3))
    perigeeRadius = (1 - ecc**2)

    for i in range(1, 13):
        semiparameter = SMA*(1-ecc**2)
        RAANPert = -1.5*n*J2*((rE/semiparameter)**2)*c
        argPeriGPert = 0.75*n*J2*((rE/semiparameter)**2)*(4-(5*s**2))
        MAPert = 0.75*n*J2*((rE/semiparameter)**2)*np.sqrt(1-ecc**2)*(2-(3*s**2))
        n = revperday*((2 * np.pi / 86400)-RAANPert) - (MAPert - argPeriGPert)
        SMA = (mu ** (1 / 3)) * ((1 / n) ** (2 / 3))
    return SMA

#Calculate the inclination from given altitude
def SSO(alt0):
    #rad = -(alt0/12352)**(7/2)
    rad = -(2*alt0**(7/2)*(1.991063853*10**(-7))*(1-ecc**2)**2)/(3*rE**2*J2*np.sqrt(mu))
    inc_new = math.acos(rad)*180/np.pi
    return inc_new


class RevisitTime:
    """Run its create_and_run_engines method"""

    def __init__(self):
        pass

    def create_and_run_engine(self, alt, inc, numberSats, numberPlanes, lookAngle, raanSpreading, period):
        """Creates and runs engine"""
        #Check orbit
        alt_output = alt
        alt = alt + rE
        if inc == -1:
            inc1 = SSO(alt)
            inc = inc1
        else:
            pass
            #alt_new = find_revpday(alt, inc)
            #alt = alt_new

        #Check Number of Planes
        if numberSats == 1:
            numberPlanes = 1
        elif numberSats == 2:
            if numberPlanes == 4:
                numberPlanes = 2
        else:
            pass

        #print(ExampleUtilities.get_examples_path())
        mission_plan_path = ExampleUtilities.combine_paths(ExampleUtilities.get_examples_path()
                                                           ,"AnalysisRevisitTime_v2_ex.MissionPlan")
        ExampleUtilities.set_working_directory_to_program_directory()
        # Get path to runtime library
        ff_install_dir = ExampleUtilities.get_freeflyer_install_directory()
        try:
            with RuntimeApiEngine(ff_install_dir, consoleOutputProcessingMethod=
            ConsoleOutputProcessingMethod.RedirectToRuntimeApi,
                                  windowedOutputMode=None) as engine:
                #Load a Mission Plan into the engine for use
                #print("Load the Mission Plan.")
                engine.loadMissionPlanFromFile(mission_plan_path)
                #Prepare the engine to execute the statements contained in the Mission Plan
                #print("Prepare to execute statements.")
                engine.prepareMissionPlan()
                #Run Script 1
                #print("Run to the 'Description state' label.")
                engine.executeUntilApiLabel("Description state")
                #print("Run to the 'Set state' label.")
                engine.executeUntilApiLabel("Set state")
                engine.setExpressionVariable('sensorAngle', lookAngle)
                engine.setExpressionVariable('h', alt)
                #print("Run to the 'SetSensor state' label.")
                engine.executeUntilApiLabel("SetSensor state")
                #Formation Setting
                engine.setExpressionString("propType", propType)
                engine.setExpressionVariable('scFormation.Count', numberSats)
                engine.setExpressionVariable('scFormation.ViewAsGroup', 0)
                #ER = engine.getExpressionVariable('Earth.Radius')
                engine.setExpressionVariable('period', period)
                # calculate sensor angle
                lookAngle = engine.getExpressionVariable('look_ang_range')
                #print("Run to the 'GenerateSat state' label.")
                engine.executeUntilApiLabel("GenerateSat state")
                numberSatPerPlane = numberSats // numberPlanes

                for j in range(0, numberPlanes):
                    for i in range(0, numberSatPerPlane):
                        engine.setExpressionVariable("scFormation[" + str(i + numberSatPerPlane * j) + "].A", alt)
                        engine.setExpressionVariable("scFormation[" + str(i + numberSatPerPlane * j) + '].I', inc)

                        raant = raan + raanSpreading / numberPlanes * j
                        engine.setExpressionVariable("scFormation[" + str(i + numberSatPerPlane * j) + "].RAAN", raant)

                        engine.setExpressionVariable("scFormation[" + str(i + numberSatPerPlane * j) + "].E", ecc)
                        ta1 = ta + 360 / numberSatPerPlane * i + 360 / numberSats * j
                        engine.setExpressionVariable("scFormation[" + str(i + numberSatPerPlane * j) + "].TA", ta1)
                        engine.setExpressionVariable("scFormation[" + str(i + numberSatPerPlane * j) + "].W", argprg)
                        engine.setExpressionVariable("scFormation[" + str(i + numberSatPerPlane * j) +
                                                     "].Sensors[0].MaskType", 3)
                        engine.setExpressionArray("scFormation[" + str(i + numberSatPerPlane * j) +
                                                  "].Sensors[0].RectangularHalfAngles", [sensorHeight / 2, lookAngle])
                        engine.setExpressionArray("scFormation[" + str(i + numberSatPerPlane * j) +
                                                  "].Sensors[0].BoresightRotationSeq", [3, 1, 2])
                        engine.setExpressionArray("scFormation[" + str(i + numberSatPerPlane * j) +
                                                  "].Sensors[0].BoresightAngles", [0, 0, 0])
                        engine.setExpressionTimeSpan("scFormation[" + str(i + numberSatPerPlane * j) +
                                    "].Propagator.StepSize", FFTimeSpan.fromWholeSecondsAndNanoseconds(stepSize, 0))


                #Run Script2
                #print("Run to the 'PointGroup state' label.")
                engine.executeUntilApiLabel("PointGroup state")
                #print("Run to the 'SimRevisit state' label.")
                engine.executeUntilApiLabel("SimRevisit state")

                #Date Export
                pointrevisitavg = np.array(engine.getExpressionArray("PointRevisit[NumOfPoints:NumOfPoints*2-1]"))
                pointrevisitnum = np.array(engine.getExpressionArray("PointRevisit[0:NumOfPoints-1]"))
                PointRevisitMax = np.array(engine.getExpressionArray('PointRevisit[NumOfPoints*4:NumOfPoints*5-1]'))
                PointCoverage = np.array(engine.getExpressionArray('PointCoverage'))
                PointCoverage_0 = PointCoverage[PointCoverage == 0]
                PointCoverage_1 = PointCoverage[PointCoverage == 1]
                Percentage_Coverage = round(len(PointCoverage_1)/(len(PointCoverage_0)+len(PointCoverage_1))*100,1)
                #Exclude all zeros
                pointrevisitnum_no0 = pointrevisitnum[pointrevisitnum > 0]
                pointrevisitavg_no0 = pointrevisitavg[pointrevisitavg > 0]
                PointRevisitMax_no0 = PointRevisitMax[PointRevisitMax > 0]
                novisit = 0
                if len(pointrevisitnum) - len(pointrevisitnum_no0) > 1:
                    novisit = 1
                    print("There are points where Satellites cannot visit!")
                #Prevent errors for the cases which satellites never visit any point
                if len(pointrevisitavg_no0) == 0:
                    pointrevisitavg_no0 = [0]
                else:
                    pass

                avgVar = np.mean(pointrevisitavg_no0 * 24)
                #minVar = min(pointrevisitavg_no0 * 24 * 3600)
                maxVar = max(pointrevisitavg_no0 * 24)
                maxmax = max(PointRevisitMax_no0 * 24)
                #print("Clean up the Mission Plan.")
                engine.cleanupMissionPlan()

        except RuntimeApiException as exp:
            print(exp.message)

        return alt_output, inc, numberSats, numberPlanes, lookAngle-3.9, avgVar, maxmax, maxVar, novisit, Percentage_Coverage


#if __name__ == '__main__':
#    print(RevisitTime.create_and_run_engine(1, 500, 44.2, 16, 4, 0, 360, 1))
