import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import integrate

pi       = np.pi
SpeedOfLight = 2.998e8
Constant = 3*1.6*10**-26
def diff(B,r,theta,t,stride):
    #print('B:',B)
    #print('f(x+t)',B(r,theta,t+stride))
    #print('f(t)',B(r,theta,t))
    #print('Δf',B(t,theta,t+stride)-B(r,theta,t))
    return (B(r,theta,t+stride)-B(r,theta,t))/stride

def integration(func,x_bound,y_bound,interval = 10000):
    interval_r     = interval
    interval_theta = 100
    inte           = 0
    for i in range(interval_theta):
        func1 = lambda r: func(r,y_bound[0]+(y_bound[1]-y_bound[0])/interval_theta)
        inte += integrate.romberg(func1,x_bound[0],x_bound[1])
    return inte
def EMF(r0,v,radius):
    beta = np.linalg.norm(v)/(SpeedOfLight)
    epoch = 1000
    stride = np.linalg.norm(r0)/(0.5*epoch*np.linalg.norm(v))
    ####function definition:
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
    ##check formula:
    #print('B_func:',B)
    #print('B',B(0.01,np.pi/2,r0[0]/v[0]))
    #print('分母',((np.linalg.norm((-r0+v*(r0[0]/v[0])+R(0.01,np.pi/2)))\
     #                                                  )**2*(1-beta**2)
     #     +(np.dot(v,-r0+v*(r0[0]/v[0])+R(0,np.pi/2)))/(3*10**8))**1.5)
    #print('分子',C*(1-beta**2)*np.dot(np.cross(v,-r0+R(0.01,np.pi/2)),np.array([0,0,1])))
    #print('diff ',func(0.01,np.pi/2,r0[0]/v[0]))
    #return 0
    t = 0
    time = []
    out  = []
    for i in range(epoch):
        t = t + stride
        fun1 = lambda theta,r:function(r,theta,t)
        #inte = integrate.dblquad(fun1,0,2*np.pi,0,radius)[0]
        inte = integration(fun1,[0,radius],[0,2*np.pi])
        print('inetegrate value',inte)
        time.append(10**6*t)  ## in μs
        out.append(inte)
    print('stride'+str(stride*10**6)+'μs')
    print('end location:',r0+v*stride*epoch,'m')
    return time,out




theta = 0
r = 0.3
ma = []
y  = 0
z  = 0.01
velocity = []
beta = 0.9
for k in range(1):
    theta = np.pi/2
    v     = beta*SpeedOfLight*np.array([np.sin(theta),0,np.cos(theta)])
    source = -r*np.array([np.sin(theta),0,np.cos(theta)])+np.array([0,y,z])
    t,out = EMF(source,v,0.02)
    data = pd.DataFrame([t,out])
    print('direction',theta/np.pi)
    print('source',source)
    print('velocity:',v)
    data.to_csv(r'/home/dachuang/lbg/CoilSimulation/data_beta09')
    #parameters = [v,source,beta]
    parameter = pd.DataFrame([v,source,[beta]])
    parameter.to_csv(r'/home/dachuang/lbg/CoilSimulation/parameter_beta09')
