# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 19:53:59 2024

@author: ismt
"""
import numpy as np
class RLScontroller(object):
    # w0 - initial estimate used to initialize the estimator
    # P0 - initial estimation error covariance matrix
    # R  - covariance matrix of the measurement noise
    def __init__(self, w0, P0, mu, v, dt):
        self.wprev=w0
        self.Pprev=P0
        self.dt = dt        
        self.I = np.eye(len(self.wprev))                           
        self.weightmatrix=[]
        self.weightmatrix.append(self.wprev)
        self.eprev, self.eprev2, self.eprev3 = 0, 0, 0
        self.s = 0
        self.sprev = 0
        self.sprev2 = 0
        self.uprev = 0
        self.ise = 0
        self.iae = 0
        self.umax = 150
        self.umin = 0
        self.mu = mu
        self.v = v
        self.estimationErrorCovarianceMatrices=[]
        self.estimationErrorCovarianceMatrices.append(self.Pprev)
        self.errors=[]
        self.controlerrormatrix = []
        
    def controller_output(self,error,normalize = 'sigmoid'):
        self.control_error = error
        if normalize == 'sigmoid':
            self.normalizedError = (self.mu * self.v *error) / (np.sqrt(self.mu**2 + self.v**2 * error**2))
        elif normalize == 'tanh':
            self.normalizedError = self.mu * np.tanh(self.v*error / self.mu)        
        
        #construct data matris..
        data_matrisprev = np.array([[self.eprev2, self.sprev2, (self.eprev2 - self.eprev3)]])
        #print("prev",data_matrisprev)
        data_matris = np.array([[self.eprev, self.sprev, (self.eprev -  self.eprev2)]])
        #print("now",data_matris)
        #calculate unew
        Lmatrix = 1 + (data_matrisprev @ (self.Pprev @ data_matrisprev.T))
        #print(Lmatrix)
        Lmatrixinv = np.linalg.inv(Lmatrix)
        #print(Lmatrixinv)
        Rregmatrix = self.Pprev @ data_matrisprev.T @ data_matrisprev @ self.Pprev
        #print(Rregmatrix)
        self.Pnew = self.Pprev - (Lmatrixinv * Rregmatrix)                
        Nmatrix = self.normalizedError - self.uprev + (data_matrisprev @ self.wprev)
        Kmatrix = self.Pnew @ data_matrisprev.T        
        self.wnew = self.wprev + (Nmatrix * Kmatrix)
        
               
        self.unew = data_matris @ self.wnew
        self.ureal = max(self.umin, min(self.unew, self.umax))        
        #update integral data..
        if self.unew < self.umin:
            self.s = float((self.s + self.normalizedError) - (np.sign(self.s) * (self.umin - self.unew)))           
        elif self.umin <= self.unew and self.unew <= self.umax:
            self.s = (self.s + self.normalizedError)               
        elif self.unew > self.umax:                        
            self.s = float((self.s + self.normalizedError) - (np.sign(self.sprev) * (self.unew - self.umax)))
                
        self.wprev = self.wnew
        self.eprev3 = self.eprev2
        self.eprev2 = self.eprev
        self.eprev = self.normalizedError        
        self.sprev2 = self.sprev
        self.sprev = self.s       
        self.uprev = self.unew
        
        
    def trackvalues(self, track=False):
        if track == True:
            self.weightmatrix.append(self.wnew)
            self.estimationErrorCovarianceMatrices.append(self.Pnew)
            self.errors.append(self.normalizedError)
            self.ise +=  self.control_error**2*self.dt
            self.iae +=  abs(self.control_error)*self.dt
            
        
            
        
        
        
        
        
        
        
    