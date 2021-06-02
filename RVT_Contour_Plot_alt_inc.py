import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
def Draw_Contour(NumSats, SenAng, typeVar):
    #Loading Result File from runRTV.py
    Output = pd.read_csv('Output.csv')

    #Sorting and Filtering the results
    Output = Output.sort_values(by='Altitude',axis=0)
    condition_Satnum = (Output['NumberSats'] == NumSats)
    condition_SenAng = (Output['SenAng'] == SenAng)
    condition_Inc47 = (Output['Inclination'] == 47)
    condition_Inc90 = (Output['Inclination'] == 90)
    condition_IncSSO = (Output['Inclination'] >= 93)

    Output_all_condition47 = Output[condition_Satnum & condition_SenAng & condition_Inc47]
    Output_all_condition90 = Output[condition_Satnum & condition_SenAng & condition_Inc90]
    Output_all_conditionSSO = Output[condition_Satnum & condition_SenAng & condition_IncSSO]


    #Define X, Y coordinate for contour plot
    alt = Output_all_condition47['Altitude']
    alt_list = np.array(alt.values.tolist())
    #inc = [42, 80]
    inc = [47, 90, 98]
    X, Y = np.meshgrid(alt_list, inc)

    #Define Z coordinate for contour plot
    Var47 = Output_all_condition47[typeVar]
    Var90 = Output_all_condition90[typeVar]
    VarSSO = Output_all_conditionSSO[typeVar]

    Var47_list = Var47.values.tolist()
    Var90_list = Var90.values.tolist()
    VarSSO_list = VarSSO.values.tolist()

    #Var22 = np.array([Var42_list, Var80_list])
    Var22 = np.array([Var47_list, Var90_list, VarSSO_list])

    #Plot Contour & Labelling
    cp = plt.contourf(X, Y, Var22)#, levels = 400)#, vmin = 0, vmax = 48)#, cmap = "PuBu")
    kp = plt.contour(X, Y, Var22)
    plt.colorbar(cp)
    plt.clabel(kp, colors = 'white')
    plt.xlabel('Altitude(km)')
    plt.ylabel('Inclination(˚)')
    plt.axis([480, 900, 45, 100])
    plt.title("Average of ART when Number of Sat = {Te_t} Sensor Angle = {n_t}"
              .format(Pe_t = typeVar, Te_t = NumSats, n_t = SenAng))

#Plotting Contour graphs in one page
#a = [1, 4, 16, 32]
#b = [10, 20, 30]
#k = 1
fig = plt.figure()
Draw_Contour(8,0,'avgHour')
"""
for i in a:
    for j in b:
        fig.add_subplot(4,3, k)
        Draw_Contour(i, j, 'Average') #Put 'avgVar' or 'maxVar' or 'minVar' for Average/Maximum/Minimum revisit time
        k = k + 1
plt.subplots_adjust(left=0.125,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.2,
                    hspace=0.5)
"""
plt.show()
