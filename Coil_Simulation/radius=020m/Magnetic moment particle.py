from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
C2 = 1
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
def MMF(r0,v,radius,M0):
    C  = 9.284e-24  ## miu
    beta = np.linalg.norm(v)/(3*10**8)
    epoch = 3000
    stride = np.linalg.norm(r0)/(0.5*epoch*np.linalg.norm(v))
    ####function definition:
    R = lambda r,theta:np.array([r*np.cos(theta),r*np.sin(theta),0])
    func = lambda r,theta,t:r*diff(B,r,theta,t,stride = stride/100)
    if False:
       a = np.arange(0,10000)*(epoch*stride)/10000
       b = []
       diff_b = []
       for i in range(a.shape[0]):
        b.append(C*func(0.001,np.pi/2,a[i]))
        diff_b.append(C*func(0.001,np.pi/2,a[i]))
       plt.figure(figsize = (10,10),dpi = 100)
       plt.plot(a,b,label = 'location = [0,0.001,0]m')
       plt.xlabel('t(μs)',fontsize=18)
       plt.ylabel('dB(T)/dt',fontsize=18)
       plt.legend(fontsize=18,loc = 'upper right')
       plt.title('integrate value as a function of t',fontsize=18)
    #plt.plot(a,diff_b)
       plt.savefig(r'/home/dachuang/lbg/CoilSimulation/integrate function-t')
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
        fun1 = lambda theta,r:func(r,theta,t)
        #inte = integrate.dblquad(fun1,0,2*np.pi,0,radius)[0]
        inte = integration(fun1,[0,radius],[0,2*np.pi])
        print('inetegrate value',C*inte)
        time.append(10**6*t)  ## in μs
        out.append(C*inte)
    print('stride'+str(stride*10**6)+'μs')
    print('end location:',r0+v*stride*epoch,'m')
    return time,out

import pandas as pd



theta = 0
r = 0.5
ma = []
y  = 0.001
z  = 0.001
velocity = []
beta = 0.9
scaler_M  = 9.284e-24
for k in range(1):
    theta    = np.pi/6
    v        = beta*3*10**8*np.array([np.sin(theta),0,np.cos(theta)])
    vector_M = scaler_M*np.array([0,0,1])
    source   = -r*np.array([np.sin(theta),0,np.cos(theta)])+np.array([0,y,z])
    t,out    = MMF(source,v,0.02)
    data     = pd.DataFrame([t,out])
    print('direction',theta/np.pi)
    print('source',source)
    print('velocity:',v)
    data.to_csv(r'/home/dachuang/lbg/CoilSimulation/data_theta'+str(k))
    #parameters = [v,source,beta]
    parameter = pd.DataFrame([v,source,[beta]])
    parameter.to_csv(r'/home/dachuang/lbg/CoilSimulation/parameter_theta'+str(k))
