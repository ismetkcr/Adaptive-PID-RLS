
import numpy as np

class TF:
    def __init__(self, yprev, uprev):
        
        self.Kp = 0.1285
        self.Tau = 26.54        
        self.dt = 0.5
        self.delay = int(7/self.dt)
        self.uk_delay = np.zeros(self.delay)
        self.yprev = yprev    
            
    def evaluate(self, unew):
        a = -1/self.Tau
        b = self.Kp/self.Tau
        
        self.ynew = self.yprev + self.dt * (a * self.yprev + b * self.uk_delay[-1])
        self.uk_delay = np.roll(self.uk_delay, 1)
        self.uk_delay[0] = unew
        self.yprev = self.ynew
        
        return self.ynew
    
    def dispparameters(self):
        print('Kp = %.f, Tau = %.f, delay = %.f' % (self.Kp, self.Tau, self.delay))
        

class ParametricModel:
    def __init__(self, yprev, uprev):
       self.yprev = yprev
       self.uprev = uprev
       self.y_k2 = self.yprev
       self.y_k1 = self.yprev
       
       self.u_k2 = self.uprev
       self.u_k1 = self.uprev
       #self.Parameters = np.array([9.91642589e-01,  6.55843000e-04, -7.82225102e-05,  3.36225486e-09])
       
       # self.Parameters = np.array([1.03325254e+00, -3.47668962e-02,  6.94744301e-06, -1.80377910e-06,
       #      4.95126009e-06])
       #self.Parameters = np.array([9.99087371e-01, 3.97501415e-04])
       # self.Parameters = np.array([1.00360140e+00, -5.51842567e-03,  3.80484799e-06, -7.68663375e-08,
       #         3.98200482e-06])
       self.Parameters = np.array([1.03325254e+00, -3.47668962e-02,  6.94744301e-06, -1.80377910e-06,
                4.95126009e-06])
    def output(self, unew):
        #self.model_matrices = np.array([self.y_k1, self.y_k1**2, self.u_k1,  self.u_k2**4])
        # self.model_matrices = np.array([self.y_k1, self.y_k2, self.u_k1, self.u_k2, 
        #                         self.y_k1**2, self.y_k1*self.y_k2, self.y_k1*self.u_k1, 
        #                         self.y_k2**2, self.y_k2*self.u_k1, self.y_k2*self.u_k2])
        #self.model_matrices = np.array([self.y_k1, self.u_k1, self.y_k2**2, self.u_k2**3])
        #self.model_matrices = np.array([self.y_k1, self.y_k2**2, self.uk1, self.u_k2**3, ])
        self.model_matrices = np.array([self.y_k1, self.y_k2, self.y_k2*self.y_k1**2,  self.u_k1**2, self.u_k2**2])

        #self.model_matrices = np.array([np.array([self., self.y_k2, self.y_k2*self.y_k1**2,  self.u_k1**2, self.u_k2**2])])
        self.ynew = self.model_matrices @ self.Parameters.T
        
        self.y_k2 = self.y_k1
        self.y_k1 = self.ynew
        
        self.u_k2 = self.u_k1
        self.u_k1 = unew
        
        return self.ynew