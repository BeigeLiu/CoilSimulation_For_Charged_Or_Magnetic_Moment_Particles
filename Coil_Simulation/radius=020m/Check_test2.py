# -*- coding: utf-8 -*-
"""
Created on Mon May  9 14:42:18 2022

@author: admin
"""
import numpy as np
import matplotlib.pyplot as plt
pi       = np.pi
SpeedOfLight = 2.998*10**8
v        = SpeedOfLight*np.array([1,0,0])
r0       = -np.array([1,0,0])
t        = np.linalg.norm(r0)/np.linalg.norm(v)
Constant = 3*1.6*10**-26
R        = lambda r,theta:np.array([r*np.cos(theta),r*np.sin(theta),0])
vector_r = lambda r,theta,t: R(r,theta)-v*t-r0
scaler_r = lambda r,theta,t: np.linalg.norm(vector_r(r,theta,t))
v_dot_r  = lambda r,theta,t: np.dot(v,vector_r(r,theta,t))
v_cross_r= lambda r,theta,t:np.dot(
                                             np.cross(v,vector_r(r,theta,t)),
                                             np.array([0,0,1]))
first_com= lambda r,theta,t:(1-beta**2)*pow(scaler_r(r,theta,t),2)
second_com =lambda r,theta,t:(np.dot(v,vector_r(r,theta,t))/SpeedOfLight)**2
function = lambda r,theta,t,v,beta: r*Constant*(1-beta**2)*v_cross_r(r,theta,t)*\
                                    v_dot_r(r,theta,t)/(
                                        first_com(r,theta,t)+second_com(r,theta,t))**2.5
function_no_relativity = lambda r,theta,t,v: r*Constant*np.dot(
                                             np.cross(v,vector_r(r,theta,t)),
                                             np.array([0,0,1]))*                \
                             np.dot(v,vector_r(r,theta,t))/                     \
                             (scaler_r(r,theta,t))**5
betas = []
func  = []
theta = pi/2 + 0.000001
func2 = []
for i in range(1000):
    beta = 0.9999 + (i)/10000000
    v    = beta*SpeedOfLight*np.array([1,0,0])
    t        = np.linalg.norm(r0)/np.linalg.norm(v)
    func.append(function(0.01,theta,t,v,beta))
    func2.append(function_no_relativity(0.01,theta,t,v))
    betas.append(beta)
plt.figure(figsize=(10,10),dpi = 100)
plt.plot(betas,np.log10(np.abs(func)),label = 'relativity')
plt.plot(betas,np.log10(np.abs(func2)),label = 'none relativity')
plt.xlabel('beta',fontsize = 18)
plt.ylabel('intergrate function',fontsize = 18)
plt.legend(fontsize = 18)
plt.savefig(r'C:\Users\admin\Desktop\inter_func')