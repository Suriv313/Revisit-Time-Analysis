from itertools import product
import numpy as np
import RVT_main
import multiprocessing
import csv

#Initial Values
period = [60]
rE = 6378.137 #Earth Radius

alt_init = [500,888] #Initial Altitude
inc_init = [-1, 42, 90] #Initial Inclination, -1 represents SSO(Sun Sychronized Oribit
num_sat = [1,2,4,8,16,32] #Number of Satellites
num_plane = [1,2,4] #Number of Planes
look_ang = [0,10,20] #Sensor Angle
constellation = [360] #RAAN Spreading

#Create Input Arrays
items = [alt_init, inc_init, num_sat, num_plane, look_ang, constellation, period]
dff = list(product(*items)) #Full Factorial Design
result = np.zeros((len(dff), 10))
Out = np.zeros((len(dff), 10))
RVTime = RVT_main.RevisitTime()

#Run Engine
def RunEngine(i):
   result[i,:] = RVTime.create_and_run_engine(dff[i][0], dff[i][1], dff[i][2], dff[i][3], dff[i][4], dff[i][5],
                                              dff[i][6])
   print(i+1,'/',len(dff))
   print(dff[i])
   Out = result[i,:]
   return Out


if __name__ == '__main__':
    pool = multiprocessing.Pool(processes=6)
    Output = pool.map(RunEngine, range(0, len(dff), 1))
    with open('Output.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        header_list = ['Altitude', 'Inclination', 'NumberSats','NumberPlanes', 'SenAng', 'AvgHour', 'maxmaxHour', 'MaxHour', 'novisit', 'Coverage Percentage']
        writer.writerow(header_list)
        writer.writerows(Output)


