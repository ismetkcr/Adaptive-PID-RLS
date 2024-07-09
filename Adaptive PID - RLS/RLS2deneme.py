# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 13:45:49 2024

@author: ismt
"""

import numpy as np
from controlutils import  Plotter,  lowpass_filter_single_value
from RLScontrollerclass import RLScontroller
from Modelclass import TF
from controlutils import Subplotter



dt = 0.5 # bunu farklı yerden vermeyi düşüneceğim..
#--------------construct model with TF-----------
dt = 0.5 
Kp = 0.06018
Tau = 70.09
delay = 3
KU = 350
PU = 8


yprev, yprev2, uprev, uprev2 = 0, 0, 0, 0
modelTF = TF(yprev, uprev)

#KC, TI, TD = 44, 30, 1

#wprev = np.array([[77, (77/9.43) , 77*1.42 ]]).T

wprev = np.array([[0, 0, 0]]).T
Pprev = 10000 * np.eye(3)
mu = 0.05
v = 0.1
t = 0
tfinal = 2700

ft = int(tfinal/dt)
r = 7

yrec = np.array([yprev])
rrec = np.array([0])
trec = np.zeros((1,1))
urec = uprev


ct = RLScontroller(wprev, Pprev, mu, v, dt)
eprev, eprev2, eprev3 = 0, 0, 0

x_array = None
for k in range(ft):
    if k == 0:
        r = 7
    t += dt
    if t % 900 == 0:
        if r<7.2:
            r = 8
        elif r>7.2:
            r = 7
    
    rf, x_array = lowpass_filter_single_value(r, x_array)
    ct.controller_output(eprev)
    ct.trackvalues(track= True)
    ureal = ct.ureal
    #TF output
    modelTF.evaluate(ureal)
    ynew = modelTF.ynew
    error = rf - ynew
    
    
    
    yprev = ynew   
    uprev = ureal                     
    eprev = error
    
    
         
    trec = np.append(trec, t)
    yrec = np.append(yrec, ynew)
    rrec = np.append(rrec, r)
    urec = np.append(urec, ureal)
    

plotter = Plotter(xlabel = 'time, [s]',
                  ylabel = 'out',
                  title = 'APIDRLS')
plotter.addplot(trec, yrec, label = 'processout')
plotter.addplot(trec, rrec, label = 'setpoint')

plotter.show()




ct.weightmatrix = np.array(ct.weightmatrix)
KCrec = ct.weightmatrix[:,0]
KIrec = ct.weightmatrix[:,1]
KDrec = ct.weightmatrix[:,2]
plotter2 = Subplotter(4)
plotter2.addsubplot([trec], [KCrec], title='KC', xlabel='Time', ylabel='KC')
plotter2.addsubplot([trec], [KIrec], title='KI', xlabel='Time', ylabel='KI')
plotter2.addsubplot([trec], [KDrec], title='KD', xlabel='Time', ylabel='KD')
plotter2.show()


















""" denemeler içindi açma buraları boşver...

# normalizedError = (mu * v *error) / (np.sqrt(mu**2 + v**2 * error**2))
# data_matrisprev = np.array([[eprev2, sprev2, (eprev2 - eprev3)]])
# data_matris = np.array([[eprev, sprev, (eprev -  eprev2)]])
# Lmatrix = 1 + (data_matrisprev @ Pprev @ data_matrisprev.T)    
# Lmatrixinv = np.linalg.inv(Lmatrix)
# Rregmatrix = Pprev @ data_matrisprev.T @ data_matrisprev @ Pprev
# Pnew = Pprev - (Lmatrixinv * Rregmatrix)                
# Nmatrix = normalizedError - uprev + (data_matrisprev @ wprev)
# Kmatrix = Pnew @ data_matrisprev.T        
# wnew = wprev + (Nmatrix * Kmatrix)
# unew = data_matris @ wnew
# ureal = max(umin, min(unew, umax))
        
# #update integral data..
#  if unew < umin:
#     s = float((s + normalizedError) - (np.sign(s) * (umin - unew)))
#     #s =sprev
#    # print("min", s)
          
#  elif umin <= unew and unew <= umax:
#     s = (s + normalizedError)
       
#  elif unew > umax:                     
#     s = float((s + normalizedError) - (np.sign(sprev) * (unew - umax)))
#     #print("max", s)
#     #s = sprev
    
#  eprev3 = eprev2
#  eprev2 = eprev
#  eprev = normalizedError        
#  sprev2 = sprev
#  sprev = s       
#  uprevreal = ureal 
#  uprev = unew
#  yprev = ynew
#  wprev = wnew
"""








