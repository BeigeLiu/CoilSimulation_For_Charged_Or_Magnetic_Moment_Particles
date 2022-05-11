from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
SpeedOfLight = 2.997*1e8 #in m/s
def diff(B,r,theta,t,stride):
    if False:
       print('B:',B)
       print('stride:',stride)
       print('f(stride+t)',B(r,theta,t+stride))
       print('f(t)',B(r,theta,t))
       print('Δf',(B(r,theta,t+stride)-B(r,theta,t))/stride)
    return (B(r,theta,t+stride)-B(r,theta,t))/stride

def integration(func,x_bound,y_bound,interval = 10000):
    interval_r     = interval
    interval_theta = 100
    inte           = 0
    for i in range(interval_theta):
        func1 = lambda r: func(r,y_bound[0]+i*(y_bound[1]-y_bound[0])/interval_theta)
        inte += integrate.romberg(func1,x_bound[0],x_bound[1])
    return inte
def MMF(r0,v,radius,M0):
    mu       = 1.2566*1e-6  ## mu
    beta     = np.linalg.norm(v)/(SpeedOfLight)
    epoch    = 2000
    stride   = np.linalg.norm(r0)/(0.5*epoch*np.linalg.norm(v))
    ####function definition:
    def integrated_function(r0,v,r,theta,t,Plot_details = False):
             R        = lambda r,theta:np.array([r*np.cos(theta),r*np.sin(theta),0])
             vec_r    = lambda r,theta,t:-r0-v*t + R(r,theta)
             sca_r    = lambda r,theta,t:np.linalg.norm(vec_r(r,theta,t))
             first_n  = lambda r,theta,t:3/4*((1-beta**2)**1.5*np.dot(M0,vec_r(r,theta,t)))*vec_r(theta,r,t)
             first_d  = lambda r,theta,t:((1-beta**2)*(sca_r(r,theta,t))**2+(np.dot(vec_r(r,theta,t),v)/SpeedOfLight)**2)**2.5
             second_n = lambda r,theta,t:1/4*(1-beta**2)*M0
             second_d = lambda r,theta,t:((1-beta**2)*(sca_r(r,theta,t))**2 + (np.dot(vec_r(r,theta,t),v)/SpeedOfLight)**2)**1.5
             B        = lambda r,theta,t:mu*(first_n(r,theta,t)/first_d(r,theta,t) - second_n(r,theta,t)/second_d(r,theta,t))
             B_z      = lambda r,theta,t:np.dot(np.array([0,0,1]),B(r,theta,t))
             func     = lambda r,theta,t:r*diff(B_z,r,theta,t,stride = stride*1e-10)##stride = deviation stride
           ##for debug:
             if  Plot_details:
                   print('Consider loction',-r0 -v*t + R(r,theta))
                   print('theta',theta)
                   print('r',r)
                   print('t',t)
                   print('Particle location',-r0-v*t)
                   print('beta:',beta)
                   print('v',v)
                   print('R:',R(r,theta))
                   print('vector r',vec_r(r,theta,t))
                   print('scaler r',sca_r(r,theta,t))
                   print('first n',first_n(r,theta,t))
                   print('first d',first_d(r,theta,t))
                   print('second_n',second_n(r,theta,t))
                   print('second_d',second_d(r,theta,t))
                   print('B',B(r,theta,t))
                   print('B_z',B_z(r,theta,t))
                   print('func',func(r,theta,t))
             return func(r,theta,t)
    t = 0
    time = []
    out  = []
    for i in range(epoch):
        t = t + stride
        fun1 = lambda r,theta:integrated_function(r0,v,r,theta,t,Plot_details = False)
        inte = integrate.dblquad(fun1,0,radius,0,2*np.pi)[0]
        #inte = integration(fun1,[0,radius],[0,2*np.pi])
        print('integrate value:',inte)
        print('Particle location',-r0-v*t)
        time.append(10**6*t)  ## in μs
        out.append(inte)
    print('stride'+str(stride*10**6)+'μs')
    print('end location:',r0+v*stride*epoch,'m')
    return time,out

import pandas as pd



theta = 0
r = 0.1
ma = []
y  = 0.0
z  = 0.0
velocity = []
beta = 0.99
scaler_M  = 9.284e-24
for k in range(1):
    theta    = np.pi/2
    print(k)
    v        = beta*SpeedOfLight*np.array([np.sin(theta),0,np.cos(theta)])
    vector_M = np.array([0,1,0])
    source   = -r*np.array([np.sin(theta),0,np.cos(theta)])+np.array([0,y,z])
    t,out    = MMF(source,v,0.05,scaler_M*vector_M)
    data     = pd.DataFrame([t,out])
    print('direction',theta/np.pi)
    print('source',source)
    print('velocity:',v)
    data.to_csv(r'/home/dachuang/lbg/CoilSimulation/MMP_data_thetapi2')
    #parameters = [v,source,beta]
    parameter = pd.DataFrame([v,source,[beta],vector_M])
    parameter.to_csv(r'/home/dachuang/lbg/CoilSimulation/MMP_parameter_thetapi2')
