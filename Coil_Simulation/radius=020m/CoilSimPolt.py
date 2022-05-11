# -*- coding: utf-8 -*-
"""
Created on Fri May  6 16:41:31 2022

@author: admin
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
max_pulse = []
velocity  = []
for j in range(1):
   name = '_thetapi2'
   data      = pd.read_csv(r'D:\python 代码\CoilSim\CoilSimulation_For_Charged_Or_Magnetic_Moment_Particles\Coil_Simulation\radius=020m\data\MMP_data'+name)
   #print(data)
   parameter = pd.read_csv(r'D:\python 代码\CoilSim\CoilSimulation_For_Charged_Or_Magnetic_Moment_Particles\Coil_Simulation\radius=020m\data\MMP_parameter'+name)
   theta     = np.arccos(parameter.iloc[1,3])
   t    = data.iloc[0,:]
#print(t)
   U    = data.iloc[1,:]
   print(U)
   index = []
   for i in range(len(U)):
     if abs(U[i]) > 10**-50:
        index.append(i)
#index = []
#for j in range(len(t)):
 #   if t[j]> 3.2 and t[j] < 3.4:
  #      index.append(j)
   plt.figure(figsize = (10,10),dpi= 100)
   U1 = U[index]
   t1 = t[index]
   max_pulse.append(max(np.abs(U1[1:])))
   velocity.append(theta)
   plt.plot(t1[1:],U1[1:],label = 'velocity:'+str(np.array(parameter.iloc[0,1:]))+'\n'+'source location:'\
         +str(np.array(parameter.iloc[1,1:]))+'\n'+'beta:'+str(np.array(parameter.iloc[2,1]))+'\n'+\
             'direction of M:'+str(np.array(parameter.iloc[3,1:])))
   plt.scatter(t1[1:],U1[1:],s=2,color = 'red')
   plt.title('Wave shape for a M_M particle',fontsize=18)
   plt.xlabel('time(μs)',fontsize=18)
   plt.ylabel('U(V)',fontsize=18)
   plt.legend(fontsize=10,loc = 'upper right')
   plt.savefig(r'C:\Users\admin\Desktop\MMP_pulse'+name+format('.png'))
  # plt.close()
print(max_pulse)
print(velocity)
print(len(velocity))
plt.figure(figsize = (10,10),dpi = 100)
plt.plot(velocity,np.abs(max_pulse))
plt.xlabel('theta(arc)')
plt.ylabel('max U(V)')
plt.title('max U as a function of theta')
plt.savefig(r'C:\Users\admin\Desktop\MMP_U_theta')