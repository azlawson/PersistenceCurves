# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 11:26:12 2018

@author: Austin
"""

import numpy as np
import matplotlib.pyplot as plt

class Diagram:
    def __init__(self, Dgm, globalmaxdeath = None, infinitedeath=None):
        self.Birth = np.array(Dgm)[:,0]
        self.Death = np.array(Dgm)[:,1]
        self.globalmaxdeath = globalmaxdeath
        self.infinitedeath = infinitedeath
        self.shape = Dgm.shape
        self.diagram = np.stack([self.Birth, self.Death], axis = 1)

    def plot(self, ptsize = 3):
        Birth = self.Birth
        Death = self.Death
        if self.infinitedeath is None:
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
            res =np.digitize(x, centers)
            tmp[i, res[0]:res[1]] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve
    
    def landscape_at_t(self, t, k):
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
        L = np.array([])
        x = np.linspace(meshstart, meshstop, numberinmesh)
        for t in x:
            L = np.append(L, X.PClandscape_at_t(t, k))
        return L
    def lifecurve(self, meshstart, meshstop, num_in_mesh):
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
        FUN = Death - Birth
        for i in range(self.shape[0]):
            x = np.array([Birth[i],Death[i]])
            res =np.digitize(x, centers)
            tmp[i, res[0]:res[1]] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve
    def birthcurve(self, meshstart, meshstop, num_in_mesh):
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
        FUN = Birth
        for i in range(self.shape[0]):
            x = np.array([Birth[i],Death[i]])
            res =np.digitize(x, centers)
            tmp[i, res[0]:res[1]] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve
    def deathcurve(self, meshstart, meshstop, num_in_mesh):
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
        FUN = Death
        for i in range(self.shape[0]):
            x = np.array([Birth[i],Death[i]])
            res =np.digitize(x, centers)
            tmp[i, res[0]:res[1]] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve
    def midlifecurve(self, meshstart, meshstop, num_in_mesh):
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
            res =np.digitize(x, centers)
            tmp[i, res[0]:res[1]] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve
    def multilifecurve(self, meshstart, meshstop, num_in_mesh):
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
            res =np.digitize(x, centers)
            tmp[i, res[0]:res[1]] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve
    def lifeentropycurve(self, meshstart, meshstop, num_in_mesh):
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
            res =np.digitize(x, centers)
            tmp[i, res[0]:res[1]] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve
    def multilifeentropycurve(self, meshstart, meshstop, num_in_mesh):
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
            res =np.digitize(x, centers)
            tmp[i, res[0]:res[1]] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve
    def midlifeentropycurve(self, meshstart, meshstop, num_in_mesh):
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
            res =np.digitize(x, centers)
            tmp[i, res[0]:res[1]] = FUN[i]
        curve = tmp.sum(axis = 0)
        return curve
    def customcurve_at_t(self, t, fun, stat):
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
        tmpBirth= Birth[(t>=self.Birth) &(t<self.Death)]
        tmpDeath =Death[(t>=self.Birth) &(t<self.Death)]
        if tmpBirth.shape[0]==0:
            return 0
        tmpT = t+np.zeros(tmpBirth.shape[0])
        f = np.vectorize(fun, otypes=[np.float])
        tmp2 = f(tmpBirth, tmpDeath, tmpT)
        return stat(tmp2)

    def customcurve(self, fun, stat, meshstart, meshstop, numberinmesh):
        L = np.array([])
        x = np.linspace(meshstart, meshstop, numberinmesh)
        for t in x:
            L = np.append(L, self.customcurve_at_t(t, fun, stat))
        return L