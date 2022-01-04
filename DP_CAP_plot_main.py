#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 13:57:14 2021

@author: vinicius fiocco sciuti
"""


import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field

class Tensor:
    
    """ Creates a tensor object with apropriate 3x3 size. It Starts with 3x3 zeros,
    if only one valu is given, it is considered to be sigma11"""
    
    def __init__(self, values = None):
        
        self.tensor = np.zeros(shape=[3,3])
        
        if not isinstance(values, np.ndarray):
        
            values = np.array(values)   
        
        if values.size > 1:
            
            r,c = values.shape
                
            self.tensor[0:r,0:c] = values
        
        else:
            
            self.tensor[0,0] = values
            
        self.I1 = self.Inv_1()
        self.I2 = self.Inv_2()
        self.I3 = self.Inv_3()
        self.p = self.calc_p()
        self.q = self.calc_q()
                
    def Inv_1(self):
        
        self.I1 = np.trace(self.tensor)
        return self.I1        
                
    def Inv_2(self):
        
        self.I2 = (np.trace(self.tensor)**2  - np.trace(self.tensor**2))/2
        return self.I2
    
    def Inv_3(self):
        
        self.I3 = np.linalg.det(self.tensor)
        return self.I3

    def calc_p(self):

        self.p = self.Inv_1()/3          
        return self.p      
    
    def calc_q(self):
        
        self.q = np.sqrt(((3.0/2.0)*( np.trace(self.tensor**2) - 
                                     (1.0/3.0)*np.trace(self.tensor)**2)))        
        return self.q
                
@dataclass
class DP_cap: 
    
    d: float = 10
    beta: float = 15
    R: float = 0.01
    pb: float = 2
    alpha: float = 0.01
    k: float = 1
    hardening: list[float] = field(default_factory=list)
        
    pa = (pb - R*d)/(1 + R*np.tan(np.deg2rad(beta)))
    
    br =  np.deg2rad(beta)   

class Plot_DP_cap:
    
    def __init__(self, material):
        
        # self.p = np.linspace([0, material.pb])
        
        self.dp_surf(material)
        
        self.cap_surf(material)
        
        
    
    def dp_surf(self,material):
        
        p = np.array([-material.d/material.br, material.pa])
        
        DP_s = p*np.tan(material.br)+material.d  
        
        
        plt.plot(p, DP_s)
        # return DP_s
        
    def cap_surf(self,material):
        
        p = np.linspace(material.pa,material.pb,20)
        
        cap_s = (
            (1+material.alpha -material.alpha/np.cos(material.br))/material.R) * np.sqrt(
                ((material.R*(material.d + material.pa*np.tan(material.br)))**2)-(p-material.pa)**2)
                
        cte = ((material.R*(material.d + material.pa*np.tan(material.br)))**2)
        deltap = (p-material.pa)**2
        cte = np.ones(deltap.shape)*cte
        
        # plt.figure()
        # plt.plot(cte,'-r')
        # plt.plot(deltap,'-b')
                                                      
        print(cap_s)
        print(p)
        plt.plot(p, cap_s,'-or')  
        return cap_s
    
    def trans_surf(self,material):
        
        trans_s = (1-(material.alpha/np.cos(material.br))) * (material.d + material.pa*
                                    np.tan(material.br)) + np.sqrt(material.alpha**2 * 
                                    (material.d + material.pa * np.tan(np.deg2rad(material.br)))**2 -
                                    (self.p - material.pa**2))
        return trans_s    
      

if __name__ == '__main__':

    print('ok')    
    
    t = Tensor(np.array([[2,0],[0,0]]))
   
    material = DP_cap(0.1, 50, 1, 2.5715, 1.0e-6, 1, [0, 5])
    
    
    # ptfe = DP_cap()
    
    Plot_DP_cap(material)
    
    