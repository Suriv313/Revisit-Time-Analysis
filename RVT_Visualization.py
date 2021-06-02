from __future__ import print_function
import numpy as np
import sys
import math
import matplotlib.pyplot as plt


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
raan = 302.271831900499 # RAAN of LTAN 13:30 for SSO
argprg = 90            # Argument of Perigee
ta = 0                 # True Anomaly
propType = "J2Mean" #Ephemeris, Norad, TwoBody, RK45, J2Mean, RK78, RK89, Cowell
stepSize = 10          # F.F Scenario Step Size
sensorHeight = 10

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
                pointrevisitLong = np.array(engine.getExpressionArray("PointGroup1.PointLongitude"))
                pointrevisitLat = np.array(engine.getExpressionArray("PointGroup1.PointLatitude"))
                pointrevisitavg_no0 = pointrevisitavg[pointrevisitavg > 0]
                maxVar = max(pointrevisitavg_no0 * 24 * 3600)
                maxVars = max(pointrevisitavg_no0)
                minVar = min(pointrevisitavg_no0 * 24 * 3600)
                minVars = min(pointrevisitavg_no0)
                maxindx = np.where(pointrevisitavg == maxVars)
                minindx = np.where(pointrevisitavg == minVars)
                maxlat = pointrevisitLat[maxindx]
                maxlong = pointrevisitLong[maxindx]
                minlat = pointrevisitLat[minindx]
                minlong = pointrevisitLong[minindx]
                #Figure Part
                fig, ax = plt.subplots()
                img = plt.imread("map3.jpg")
                ax.imshow(img, extent = [110, 150, 24, 50])
                cp = plt.scatter(pointrevisitLong, pointrevisitLat, c = pointrevisitavg*24, cmap = 'jet', vmin=0, vmax=20)
                plt.text(maxlong[0]+0.3, maxlat[0]+0.3, 'MAX = {}'.format(round(maxVar/3600)), color = 'ghostwhite')
                plt.axhline(maxlat[0], 0, 1, color='whitesmoke', linestyle='--', linewidth='1')
                plt.axvline(maxlong[0], 0, 1, color='whitesmoke', linestyle=':', linewidth='1')
                plt.text(minlong[0]+0.3, minlat[0]+0.3, 'MIN = {}'.format(round(minVar/3600)), color = 'ghostwhite')
                plt.axhline(minlat[0], 0, 1, color='white', linestyle='--', linewidth='1')
                plt.axvline(minlong[0], 0, 1, color='white', linestyle=':', linewidth='1')
                plt.title('Revisit Time Analysis')
                plt.xlabel('Longitude(˚)')
                plt.ylabel('Latitude(˚)')
                cb1 = plt.colorbar(cp, orientation='horizontal')
                cb1.set_label('Average Revisit Time')
                plt.show()
                engine.cleanupMissionPlan()

        except RuntimeApiException as exp:
            print(exp.message)

VisualRGT = RevisitTime()
VisualRGT.create_and_run_engine(800+rE, 40, 16, 4, 0, 360, 3)
