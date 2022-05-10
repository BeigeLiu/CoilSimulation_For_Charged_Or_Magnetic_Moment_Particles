# -*- coding: utf-8 -*-
"""
Created on Mon May  9 11:23:50 2022

@author: admin
"""
import numpy as np
pi       = np.pi
SpeedOfLight = 3*10**8
beta     = 0.99999999
v        = beta*SpeedOfLight*np.array([1,0,0])
r0       = -np.array([1,0,0])
t        = np.linalg.norm(r0)/np.linalg.norm(v)
Constant = 3*1.6*10**-26
R        = lambda r,theta:np.array([r*np.cos(theta),r*np.sin(theta),0])
vector_r = lambda r,theta,t: R(r,theta)-v*t-r0
scaler_r = lambda r,theta,t: np.linalg.norm(vector_r(r,theta,t))
function = lambda r,theta,t: r*Constant*(1-beta**2)*np.dot(
                                             np.cross(v,vector_r(r,theta,t)),
                                             np.array([0,0,1]))*                \
                             np.dot(v,vector_r(r,theta,t))/                     \
                             ((1-beta**2)*pow(scaler_r(r,theta,t),2) +                         \
                              (np.dot(v,vector_r(r,theta,t))/SpeedOfLight)**2)**2.5
function_no_relativity = lambda r,theta,t: r*Constant*np.dot(
                                             np.cross(v,vector_r(r,theta,t)),
                                             np.array([0,0,1]))*                \
                             np.dot(v,vector_r(r,theta,t))/                     \
                             (scaler_r(r,theta,t))**5
r = 0.01
theta = pi/2
print(function(0.01,pi/2,t))
print(function_no_relativity(0.01,pi/2,t))
dot = lambda r,theta,t:np.dot(
                                             np.cross(v,vector_r(r,theta,t)),
                                             np.array([0,0,1]))
print('cross_z:',dot(r,theta,t))
print('dot', np.dot(v,vector_r(r,theta,t)))
print('fenmu first:',(1-beta**2)*pow(scaler_r(r,theta,t),2))
print('fenmu second:',(np.dot(v,vector_r(r,theta,t))/SpeedOfLight)**2)