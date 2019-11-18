# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 11:26:12 2018

@author: Austin
"""

import numpy as np
import matplotlib.pyplot as plt

class Diagram:
    def __init__(self, Dgm, globalmaxdeath = None, infinitedeath=-1):
        #@param Dgm: A persistence Diagram, i.e. an n by 2 array
        #@param globalmaxdeath: The maximum possible death value for the persistence process, e.g. in the case of 8-bit images, the max death is 255. Leave as None if there is no such value
        #@param infinitedeath: The value that represents an infinite death value, e.g. for Perseus that value is -1, which is default for this package
        self.Birth = np.array(Dgm)[:,0]
        self.Death = np.array(Dgm)[:,1]
        self.globalmaxdeath = globalmaxdeath
        self.infinitedeath = infinitedeath
        self.shape = Dgm.shape
        self.diagram = np.stack([self.Birth, self.Death], axis = 1)

    def plot(self, ptsize = 3):
        #Produces a plot of the Diagram with infinite death marked with red diamonds
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath ==-1:
            if self.globalmaxdeath is None:
                finBirth = Birth[Death>=0]
                finDeath = Death[Death>=0]
                infBirth = Birth[Death<0]
                infDeath =Death[Death<0] + np.max(self.Death) + 2
            else:
                finBirth = Birth[Death>=0]
                finDeath = Death[Death>=0]
                infBirth = Birth[Death<0]
                infDeath = Death[Death<0] +2 + self.globalmaxdeath
        else:
            if self.globalmaxdeath is None:
                finBirth = Birth[Death!=self.infinitedeath]
                finDeath = Death[Death!=self.infinitedeath]
                infBirth =Birth[Death==self.infinitedeath]
                infDeath =Death[Death==self.infinitedeath] + np.max(self.Death) + 2
            else:
                finBirth = Birth[Death!=self.infinitedeath]
                finDeath = Death[Death!=self.infinitedeath]
                infBirth =Birth[Death==self.infinitedeath]
                infDeath = Death[Death==self.infinitedeath] +2 + self.globalmaxdeath
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
        return curve
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
    def lifeentropycurve(self, meshstart, meshstop, num_in_mesh):
        #Produces the life entropy curve of the diagram
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
                Death[Death==self.infinitedeath] =Death[Death==self.infinitedeath]+2 + self.globalmaxdeath
        bins = np.linspace(meshstart, meshstop, num_in_mesh)
        centers = (bins[1:]+bins[:-1])/2
        tmp = np.zeros([self.shape[0], num_in_mesh])
        FUN = -(Death - Birth)/np.sum(Death -Birth)*np.log((Death - Birth)/np.sum(Death -Birth))
        for i in range(self.shape[0]):
            x = np.array([Birth[i],Death[i]])
            res =np.where(np.digitize(bins, x, right=False)==1)[0]
            if len(res) !=0:
                tmp[i, res[0]:res[len(res)-1]+1] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve
    def multilifeentropycurve(self, meshstart, meshstop, num_in_mesh):
        #Produces the multiplicative life entropy curve of the diagram
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
        FUN = -(Death/Birth)/np.sum(Death/Birth)*np.log((Death /Birth)/np.sum(Death /Birth))
        for i in range(self.shape[0]):
            x = np.array([Birth[i],Death[i]])
            res =np.where(np.digitize(bins, x, right=False)==1)[0]
            if len(res) !=0:
                tmp[i, res[0]:res[len(res)-1]+1] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve
    def midlifeentropycurve(self, meshstart, meshstop, num_in_mesh):
        #Produces the midlife entropy curve of the diagram
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
        FUN = -(Death + Birth)/np.sum(Death + Birth)*np.log((Death + Birth)/np.sum(Death -Birth))
        for i in range(self.shape[0]):
            x = np.array([Birth[i],Death[i]])
            res =np.where(np.digitize(bins, x, right=False)==1)[0]
            if len(res) !=0:
                tmp[i, res[0]:res[len(res)-1]+1] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve
    def custom_curve_at_t(Dgm,fun,stat,t):
        Birth = Dgm.Birth
        Death = Dgm.Death
        if Dgm.infinitedeath == -1:
            if Dgm.globalmaxdeath is None:
                Death[Death<0] =Death[Death<0] + np.max(Dgm.Death) + 2
            else:
                Death[Death<0] = Death[Death<0] +2 + Dgm.globalmaxdeath
        else:
            if Dgm.globalmaxdeath is None:
                Death[Death==Dgm.infinitedeath] =Death[Death==Dgm.infinitedeath] + np.max(Dgm.Death) + 2
            else:
                Death[Death==Dgm.infinitedeath] = Death[Death==Dgm.infinitedeath] +2 + Dgm.globalmaxdeath
        tmpBirth= Birth[(t>=Dgm.Birth) &(t<Dgm.Death)]
        tmpDeath =Death[(t>=Dgm.Birth) &(t<Dgm.Death)]
        values = []
        for i in range(tmpBirth.shape[0]):
            values.append(fun(Dgm,tmpBirth[i], tmpDeath[i], t))
        return stat(values)
    def custom_curve(Dgm,fun, stat, meshstart, meshstop, numberinmesh):
        L = np.array([])
        x = np.linspace(meshstart, meshstop, numberinmesh)
        for t in x:
            L = np.append(L, custom_curve_at_t(Dgm, fun, stat, t))
        return L