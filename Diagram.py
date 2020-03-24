# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 11:26:12 2018

@author: Austin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
class Diagram:
    def __init__(self, Dgm, globalmaxdeath = None, infinitedeath=-1, inf_policy="keep"):
        #@param Dgm: A persistence Diagram, i.e. an n by 2 array
        #@param globalmaxdeath: The maximum possible death value for the persistence process, e.g. in the case of 8-bit images, the max death is 255. Leave as None if there is no such value
        #@param infinitedeath: The value that represents an infinite death value, e.g. for Perseus that value is -1, which is default for this package
        self.globalmaxdeath = globalmaxdeath
        self.infinitedeath = infinitedeath
        self.Birth = np.array(Dgm)[:,0]
        self.Death = np.array(Dgm)[:,1]
        if inf_policy=="remove":
            self.Birth = self.Birth[self.Death != self.infinitedeath]
            self.Death = self.Death[self.Death != self.infinitedeath]
        elif inf_policy=="keep":
            if self.globalmaxdeath is None:
                self.Death[self.Death==self.infinitedeath] =np.max(self.Death)
            else:
                self.Death[self.Death==self.infinitedeath] = self.globalmaxdeath
        self.shape = Dgm.shape
        self.diagram = np.stack([self.Birth, self.Death], axis = 1)

    def plot(self, ptsize = 3):
        #Produces a plot of the Diagram with infinite death marked with red diamonds
        Birth = self.Birth
        Death = self.Death
        if self.globalmaxdeath is None:
            finBirth = Birth[Death!=self.infinitedeath]
            finDeath = Death[Death!=self.infinitedeath]
            infBirth =Birth[Death==self.infinitedeath]
            infDeath =Death[Death==self.infinitedeath] + np.max(self.Death)
        else:
            finBirth = Birth[Death!=self.infinitedeath]
            finDeath = Death[Death!=self.infinitedeath]
            infBirth =Birth[Death==self.infinitedeath]
            infDeath = Death[Death==self.infinitedeath] + self.globalmaxdeath
        plt.scatter(finBirth, finDeath, s=ptsize)
        plt.scatter(infBirth, infDeath, marker = "D", c='red')
        plt.plot([np.min(self.Birth) - 1, np.max(self.Death) + 1], [np.min(self.Birth) - 1, np.max(self.Death) + 1], c="green")
        plt.show()        
    def Betticurve(self, meshstart, meshstop, num_in_mesh):
        #Produces the Betti curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath] +2 + self.globalmaxdeath
        bins = np.linspace(meshstart, meshstop, num_in_mesh)
        centers = (bins[1:]+bins[:-1])/2
        tmp = np.zeros([self.shape[0], num_in_mesh])
        FUN = np.ones(self.shape[0])
        for i in range(self.shape[0]):
            x = np.array([Birth[i],Death[i]])
            res =np.where(np.digitize(bins, x, right=False)==1)[0]
            if len(res) !=0:
                tmp[i, res[0]:res[len(res)-1]+1] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve
    
    def landscape_at_t(self, t, k):
        #computes Bubenik's PEristence landscapes at a particular value
        #@param t: value to compute
        #@param k: level of landscape
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath] +2 + self.globalmaxdeath
        tmp= np.stack([Birth, Death], axis=1)[(t>=self.Birth) &(t<self.Death)]
        tmp1 = np.stack([t-tmp[:,0], tmp[:,1]-t], axis = 1)
        return np.sort(np.concatenate((np.min(tmp1, axis = 1),np.zeros(k+1))))[::-1][k]

    def landscape(self, k, meshstart, meshstop, numberinmesh):
        #Produces the Persistence Landscape of the diagram
        #@param k: level of the landscape
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        L = np.array([])
        x = np.linspace(meshstart, meshstop, numberinmesh)
        for t in x:
            L = np.append(L, X.PClandscape_at_t(t, k))
        return L
    def lifecurve(self, meshstart, meshstop, num_in_mesh):
        #Produces the life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath]+2 + self.globalmaxdeath
        bins = np.linspace(meshstart, meshstop, num_in_mesh)
        #centers = (bins[1:]+bins[:-1])/2
        tmp = np.zeros([self.shape[0], num_in_mesh])
        FUN = Death - Birth
        for i in range(self.shape[0]):
            x = np.array([Birth[i],Death[i]])
            res =np.where(np.digitize(bins, x, right=False)==1)[0]
            if len(res) !=0:
                tmp[i, res[0]:res[len(res)-1]+1] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve
    def midlifecurve(self, meshstart, meshstop, num_in_mesh):
        #Produces the midlife curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath]+2 + self.globalmaxdeath
        bins = np.linspace(meshstart, meshstop, num_in_mesh)
        centers = (bins[1:]+bins[:-1])/2
        tmp = np.zeros([self.shape[0], num_in_mesh])
        FUN = 0.5*(Birth+Death)
        for i in range(self.shape[0]):
            x = np.array([Birth[i],Death[i]])
            res =np.where(np.digitize(bins, x, right=False)==1)[0]
            if len(res) !=0:
                tmp[i, res[0]:res[len(res)-1]+1] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve/2
    def multilifecurve(self, meshstart, meshstop, num_in_mesh):
        #Produces the multiplicative life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath]+ np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath] +2 + self.globalmaxdeath
        bins = np.linspace(meshstart, meshstop, num_in_mesh)
        centers = (bins[1:]+bins[:-1])/2
        tmp = np.zeros([self.shape[0], num_in_mesh])
        FUN = Death/Birth
        for i in range(self.shape[0]):
            x = np.array([Birth[i],Death[i]])
            res =np.where(np.digitize(bins, x, right=False)==1)[0]
            if len(res) !=0:
                tmp[i, res[0]:res[len(res)-1]+1] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve
    def stabilizedlifecurve(self, meshstart, meshstop, num_in_mesh):
        #Produces the life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        curve = self.lifecurve(meshstart,meshstop,num_in_mesh)/self.totallife()
        return curve
    def stabilizedmidlifecurve(self, meshstart, meshstop, num_in_mesh):
        #Produces the midlife curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        curve = self.midlifecurve(meshstart,meshstop,num_in_mesh)/self.totalmidlife()
        return curve
    def stabilizedmultilifecurve(self, meshstart, meshstop, num_in_mesh):
        #Produces the multiplicative life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        curve = self.multilifecurve(meshstart,meshstop,num_in_mesh)/self.totalmultilife()
        return curve
    def lifeentropycurve(self, meshstart, meshstop, num_in_mesh):
        #Produces the life entropy curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        tmp = self.stabilizedlifecurve(meshstart,meshstop,num_in_mesh)
        curve = -1*tmp*np.log(tmp)
        return curve
    def multilifeentropycurve(self, meshstart, meshstop, num_in_mesh):
        #Produces the multiplicative life entropy curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        tmp = self.stabilizedmidlifecurve(meshstart,meshstop,num_in_mesh)
        curve = -1*tmp*np.log(tmp)
        return curve
    def midlifeentropycurve(self, meshstart, meshstop, num_in_mesh):
        #Produces the midlife entropy curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        tmp = self.stabilizedmultilifecurve(meshstart,meshstop,num_in_mesh)
        curve = -1*tmp*np.log(tmp)
        return curve
    def custom_curve_at_t(self,fun,stat,t):
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath == -1:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath] +2 + self.globalmaxdeath
        tmpBirth= Birth[(t>=self.Birth) &(t<self.Death)]
        tmpDeath =Death[(t>=self.Birth) &(t<self.Death)]
        values = []
        for i in range(tmpBirth.shape[0]):
            values.append(fun(self,tmpBirth[i], tmpDeath[i], t))
        return stat(values)
    def custom_curve(self,fun, stat, meshstart, meshstop, numberinmesh):
        L = np.array([])
        x = np.linspace(meshstart, meshstop, numberinmesh)
        for t in x:
            L = np.append(L, self.custom_curve_at_t(fun, stat, t))
        return L
    def totallife(self):
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath] +2 + self.globalmaxdeath
        return np.sum(Death-Birth)
    def totalmidlife(self):
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath] +2 + self.globalmaxdeath
        return np.sum(Death+Birth)/2
    def totalmultilife(self):
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath] +2 + self.globalmaxdeath
        return np.sum(Death/Birth)
    def smooth_life(self, meshstart, meshstop, num_in_mesh, spread = 1):
        #Produces the gaussian life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath]+2 + self.globalmaxdeath
        T = np.linspace(meshstart, meshstop, num_in_mesh)

        Life = (Death - Birth)/np.max(Death-Birth)
        B = np.matmul(np.array(Birth).reshape(-1,1),np.ones([1,num_in_mesh]))
        De = np.matmul(np.array(Death).reshape(-1,1),np.ones([1,num_in_mesh]))
        L = np.matmul(np.array(Life).reshape(-1,1),np.ones([1,num_in_mesh]))
        return np.sum(L*stats.norm.cdf((T-B)/spread)*(1-stats.norm.cdf((T-De)/spread)), axis=0)
    def smooth_life_derivative(self, meshstart, meshstop, num_in_mesh, spread = 1):
        #Produces the gaussian life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath]+2 + self.globalmaxdeath
        T = np.linspace(meshstart, meshstop, num_in_mesh)

        Life = (Death - Birth)/np.max(Death-Birth)
        B = np.matmul(np.array(Birth).reshape(-1,1),np.ones([1,num_in_mesh]))
        De = np.matmul(np.array(Death).reshape(-1,1),np.ones([1,num_in_mesh]))
        L = np.matmul(np.array(Life).reshape(-1,1),np.ones([1,num_in_mesh]))
        return np.sum(L*(stats.norm.pdf((T-B)/spread) - stats.norm.pdf((T-B)/spread)*(stats.norm.cdf((T-De)/spread))-stats.norm.pdf((T-De)/spread)*(stats.norm.cdf((T-B)/spread))), axis=0)
    def smooth_Betti(self, meshstart, meshstop, num_in_mesh,spread=1):
        #Produces the gaussian life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath]+2 + self.globalmaxdeath
        T = np.linspace(meshstart, meshstop, num_in_mesh)

        Life = np.ones(Birth.shape[0])/Birth.shape[0]
        B = np.matmul(np.array(Birth).reshape(-1,1),np.ones([1,num_in_mesh]))
        De = np.matmul(np.array(Death).reshape(-1,1),np.ones([1,num_in_mesh]))
        L = np.matmul(np.array(Life).reshape(-1,1),np.ones([1,num_in_mesh]))
        return np.sum(L*stats.norm.cdf((T-B)/spread)*(1-stats.norm.cdf((T-De)/spread)), axis=0)
    def smooth_midlife(self, meshstart, meshstop, num_in_mesh,spread=1):
        #Produces the gaussian life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath]+2 + self.globalmaxdeath
        T = np.linspace(meshstart, meshstop, num_in_mesh)

        Life = (Birth+Death)/np.max(Birth+Death)
        B = np.matmul(np.array(Birth).reshape(-1,1),np.ones([1,num_in_mesh]))
        De = np.matmul(np.array(Death).reshape(-1,1),np.ones([1,num_in_mesh]))
        L = np.matmul(np.array(Life).reshape(-1,1),np.ones([1,num_in_mesh]))
        return np.sum(L*stats.norm.cdf((T-B)/spread)*(1-stats.norm.cdf((T-De)/spread)), axis=0)
    def smooth_midlife_derivative(self, meshstart, meshstop, num_in_mesh, spread = 1):
        #Produces the gaussian life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath]+2 + self.globalmaxdeath
        T = np.linspace(meshstart, meshstop, num_in_mesh)

        Life = (Death + Birth)/np.max(Death+Birth)
        B = np.matmul(np.array(Birth).reshape(-1,1),np.ones([1,num_in_mesh]))
        De = np.matmul(np.array(Death).reshape(-1,1),np.ones([1,num_in_mesh]))
        L = np.matmul(np.array(Life).reshape(-1,1),np.ones([1,num_in_mesh]))
        return np.sum(L*(stats.norm.pdf((T-B)/spread) - stats.norm.pdf((T-B)/spread)*(stats.norm.cdf((T-De)/spread))-stats.norm.pdf((T-De)/spread)*(stats.norm.cdf((T-B)/spread))), axis=0)
    def smooth_multilife(self, meshstart, meshstop, num_in_mesh,spread=1):
        #Produces the gaussian life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath]+2 + self.globalmaxdeath
        T = np.linspace(meshstart, meshstop, num_in_mesh)

        Life = (Death/Birth)/np.max(Death/Birth)
        B = np.matmul(np.array(Birth).reshape(-1,1),np.ones([1,num_in_mesh]))
        De = np.matmul(np.array(Death).reshape(-1,1),np.ones([1,num_in_mesh]))
        L = np.matmul(np.array(Life).reshape(-1,1),np.ones([1,num_in_mesh]))
        return np.sum(L*stats.norm.cdf((T-B)/spread)*(1-stats.norm.cdf((T-De)/spread)), axis=0)
    def smooth_lifeentropy(self, meshstart, meshstop, num_in_mesh,spread=1):
        #Produces the gaussian life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath]+2 + self.globalmaxdeath
        T = np.linspace(meshstart, meshstop, num_in_mesh)

        Life = (Death-Birth)/np.max(Death-Birth)
        Life = Life*np.log(Life)
        B = np.matmul(np.array(Birth).reshape(-1,1),np.ones([1,num_in_mesh]))
        De = np.matmul(np.array(Death).reshape(-1,1),np.ones([1,num_in_mesh]))
        L = np.matmul(np.array(Life).reshape(-1,1),np.ones([1,num_in_mesh]))
        return np.sum(L*stats.norm.cdf((T-B)/spread)*(1-stats.norm.cdf((T-De)/spread)), axis=0)
    def smooth_life_sum(self, meshstart, meshstop, num_in_mesh):
        #Produces the gaussian life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath]+2 + self.globalmaxdeath
        T = np.linspace(meshstart, meshstop, num_in_mesh)

        Life = (Death - Birth)/np.sum(Death-Birth)
        B = np.matmul(np.array(Birth).reshape(-1,1),np.ones([1,num_in_mesh]))
        De = np.matmul(np.array(Death).reshape(-1,1),np.ones([1,num_in_mesh]))
        L = np.matmul(np.array(Life).reshape(-1,1),np.ones([1,num_in_mesh]))
        return np.sum(L*stats.norm.cdf((T-B))*(1-stats.norm.cdf((T-De))), axis=0)
    def smooth_midlife_sum(self, meshstart, meshstop, num_in_mesh):
        #Produces the gaussian life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
            if self.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(self.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                Death[Death==self.infinitedeath] = Death[Death==self.infinitedeath]+2 + self.globalmaxdeath
        T = np.linspace(meshstart, meshstop, num_in_mesh)

        Life = (Death + Birth)/np.sum(Death+Birth)
        B = np.matmul(np.array(Birth).reshape(-1,1),np.ones([1,num_in_mesh]))
        De = np.matmul(np.array(Death).reshape(-1,1),np.ones([1,num_in_mesh]))
        L = np.matmul(np.array(Life).reshape(-1,1),np.ones([1,num_in_mesh]))
        return np.sum(L*stats.norm.cdf((T-B))*(1-stats.norm.cdf((T-De))), axis=0)