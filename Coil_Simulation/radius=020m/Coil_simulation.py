import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import integrate

pi       = np.pi
SpeedOfLight = 2.998e8
Constant = 3*1.6*10**-26

def integration(func,x_bound,y_bound,interval = 10000):
    #interval_r     = interval
    interval_theta = 500
    inte           = 0
    #theta = []
    #integr = []
    for i in range(interval_theta):
        func1 = lambda r: func(r,y_bound[0]+i*(y_bound[1]-y_bound[0])/interval_theta)
        #theta.append(y_bound[0]+i*(y_bound[1]-y_bound[0])/interval_theta)
        #integr.append(integrate.romberg(func1,x_bound[0],x_bound[1]))
        inte += integrate.romberg(func1,x_bound[0],x_bound[1])
    return inte
def diff(B,r,theta,t,stride):
    if False:
       print('B:',B)
       print('stride:',stride)
       print('f(stride+t)',B(r,theta,t+stride))
       print('f(t)',B(r,theta,t))
       print('Δf',B(r,theta,t+stride)-B(r,theta,t))
    return (B(r,theta,t+stride)-B(r,theta,t))/stride
def EMF(r0,v,radius):
    beta        = np.linalg.norm(v)/(SpeedOfLight)
    epoch       = 1000
    stride      = np.linalg.norm(r0)/(0.5*epoch*np.linalg.norm(v))
    def integrated_function(r0,v,r,theta,t,Plot_details = False):
                ####function definition:
                R           = lambda r,theta:   np.array([r*np.cos(theta),r*np.sin(theta),0])
                vector_r    = lambda r,theta,t: R(r,theta)-v*t-r0
                scaler_r    = lambda r,theta,t: np.linalg.norm(vector_r(r,theta,t))
                numerator   = lambda r,theta,t: (1-beta**2)*np.cross(vector_r(r,theta,t),v)
                denominator = lambda r,theta,t: ((1-beta**2)*(scaler_r(r,theta,t))**2 + (np.dot(vector_r(r,theta,t),v)/SpeedOfLight)**2)**2.5
                B           = lambda r,theta,t: Constant*numerator(r,theta,t)/denominator(r,theta,t)
                B_z         = lambda r,theta,t: np.dot(B(r,theta,t),np.array([0,0,1]))
                func        = lambda r,theta,t: r*diff(B_z,r,theta,t,stride = stride*1e-10)##stride = deviation stride
                function    = lambda r,theta,t: r*Constant*(1-beta**2)*np.dot(
                                             np.cross(v,vector_r(r,theta,t)),
                                             np.array([0,0,1]))*                \
                             np.dot(v,vector_r(r,theta,t))/                     \
                             ((1-beta**2)*pow(scaler_r(r,theta,t),2) +                         \
                              (np.dot(v,vector_r(r,theta,t))/SpeedOfLight)**2)**2.5
                ##for debug:
                if Plot_details:
                       print('Consider loction',-r0 -v*t + R(r,theta))
                       print('theta:',theta)
                       print('r:',r)
                       print('Particle location',-r0-v*t)
                       print('beta:',beta)
                       print('v',v)
                       print('R:',R(r,theta))
                       print('vector r',vector_r(r,theta,t))
                       print('scaler r',scaler_r(r,theta,t))
                       print('numerator',numerator(r,theta,t))
                       print('denominator',denominator(r,theta,t))
                       print('B',B(r,theta,t))
                       print('B_z',B_z(r,theta,t))
                       print('func',func(r,theta,t))
                       print('function:',function(r,theta,t))
                return func(r,theta,t)
    t = 0
    time = []
    out  = []
    Plot_details = False
    for i in range(epoch):
        t = t + stride
        fun1 = lambda r,theta:integrated_function(r0,v,r,theta,t,Plot_details)
        inte = integrate.dblquad(fun1,0,radius,0,2*np.pi)[0]
        #inte = integration(fun1,[0,radius],[0,2*np.pi])
        print('Particle location',-r0-v*t)
        print('integrate value',inte)
        if abs(inte) > 10**-5:
            Plot_details = False
        else:
            Plot_details = False
        time.append(10**6*t)  ## in μs
        out.append(inte)
    print('stride'+str(stride*10**6)+'μs')
    print('end location:',r0+v*stride*epoch,'m')
    return time,out




theta = 0
r = 0.1
ma = []
y  = 0.025
z  = 0.001
velocity = []
beta = 0.99
for k in range(1):
    theta = np.pi/2
    v     = beta*SpeedOfLight*np.array([np.sin(theta),0,np.cos(theta)])
    source = -r*np.array([np.sin(theta),0,np.cos(theta)])+np.array([0,y,z])
    t,out = EMF(source,v,0.05)
    data = pd.DataFrame([t,out])
    print('direction',theta/np.pi)
    print('source',source)
    print('velocity:',v)
    data.to_csv(r'/home/dachuang/lbg/CoilSimulation/data_thetapi2')
    #parameters = [v,source,beta]
    parameter = pd.DataFrame([v,source,[beta]])
    parameter.to_csv(r'/home/dachuang/lbg/CoilSimulation/parameter_thetapi2')
