# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 11:26:12 2018
@author: Austin Lawson (azlawson@uncg.edu)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
class Diagram:
    def __init__(self, Dgm, globalmaxdeath = None, infinitedeath=-1, inf_policy="keep"):
        """
        Transforms a diagram (n by 2 NumPy array or Pandas DataFrame) into the Diagram class.
        
        Parameters
        Dgm: A persistence diagram, i.e. an n by 2 array or pandas dataframe
        
        globalmaxdeath: The maximum possible death value for the persistence process, used to replace infinite death value.
        
        infinitedeath: The value that represents an infinite death value, e.g. for Perseus that value is -1, which is default for this package
        
        inf_policy: 'keep' or 'remove'. If 'keep', infinite death values will be replaced globalmaxdeath value if set otherwise the max death value of the diagram
        """
        self.globalmaxdeath = globalmaxdeath
        self.infinitedeath = infinitedeath
        self.Birth = np.array(Dgm)[:,0]
        self.Death = np.array(Dgm)[:,1]
        if inf_policy=="remove":
            self.Birth = self.Birth[self.Death != self.infinitedeath]
            self.Death = self.Death[self.Death != self.infinitedeath]
        elif inf_policy=="keep":
            if self.globalmaxdeath is None:
                self.Death[self.Death==self.infinitedeath] =np.max(self.Death[self.Death!=self.infinitedeath])
            else:
                self.Death[self.Death==self.infinitedeath] = self.globalmaxdeath
        self.shape = Dgm.shape
        self.diagram = np.stack([self.Birth, self.Death], axis = 1)

    def plot(self, xlim = None,ylim=None, ptsize = 3):
        """
        Produces a plot of the Diagram
        """
        if xlim == None:
            (self.Birth.min(),self.Death.max())
        if ylim == None:
            (self.Birth.min(),self.Death.max())
        Birth = self.Birth
        Death = self.Death
        if self.globalmaxdeath is None:
            plt.scatter(Birth, Death, s=ptsize)
            plt.plot(xlim, ylim, c="green")
        else:
            plt.scatter(Birth, finDeath, s=ptsize)
            plt.scatter(Birth[Death==globalmaxdeath], Death[Death==globalmaxdeath], marker = "D", c='red')
            plt.plot(xlim,ylim, c="green")
    def Betticurve(self, meshstart, meshstop, num_in_mesh):
        """
        Produces the Betti curve of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of Betti curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        Birth = self.Birth
        Death = self.Death
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

    def landscape2(self, k, meshstart, meshstop, numberinmesh):
        """
        Produces k-th persistence landscape (http://www.jmlr.org/papers/volume16/bubenik15a/bubenik15a.pdf) of the diagram
        
        Parameters:
        k: level of the landscape
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        Outputs the k-th persistence landscape values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        Birth = self.Birth
        Death = self.Death
        def landscape_at_t(Birth,Death, t, k):
            tmp= np.stack([Birth, Death], axis=1)[(t>=Birth) &(t<Death)]
            tmp1 = np.stack([t-tmp[:,0], tmp[:,1]-t], axis = 1)
            return np.sort(np.concatenate((np.min(tmp1, axis = 1),np.zeros(k+1))))[::-1][k]
        L = np.array([])
        x = np.linspace(meshstart, meshstop, numberinmesh)
        for t in x:
            L = np.append(L, landscape_at_t(Birth,Death,t, k))
        return L
    def landscape(self, k, meshstart, meshstop, num_in_mesh):
        Birth = self.Birth.reshape([self.shape[0],1])
        Death = self.Death.reshape([self.shape[0],1])
        T = np.matmul(np.linspace(meshstart, meshstop, num_in_mesh).reshape([num_in_mesh,1]),np.ones([1,self.shape[0]])).T
        tmpB = T-Birth
        tmpD = Death-T
        tri=np.min(np.stack([tmpB,tmpD]),axis=0)
        land=np.sort(tri,axis=0)[self.shape[0]-k-1]
        land[land<0]=0
        return land
        
    def lifecurve(self, meshstart, meshstop, num_in_mesh):
        """
        Produces the lifespan curve of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of lifespan curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        Birth = self.Birth
        Death = self.Death
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
        """
        Produces the midlife curve of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of midlife curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        Birth = self.Birth
        Death = self.Death
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
        """
        Produces the multiplicative life curve of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of multiplicative life curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        Birth = self.Birth
        Death = self.Death
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
        """
        Produces the stabilized(normalized) life curve of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of stabilized(normalized) life curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        curve = self.lifecurve(meshstart,meshstop,num_in_mesh)/self.totallife()
        return curve
    def stabilizedmidlifecurve(self, meshstart, meshstop, num_in_mesh):
        """
        Produces the stabilized(normalized) midlife curve of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of stabilized(normalized) midlife curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        curve = self.midlifecurve(meshstart,meshstop,num_in_mesh)/self.totalmidlife()
        return curve
    def stabilizedmultilifecurve(self, meshstart, meshstop, num_in_mesh):
        """
        Produces the stabilized(normalized) multiplicative life curve of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of stabilized(normalized) multiplicative life curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        curve = self.multilifecurve(meshstart,meshstop,num_in_mesh)/self.totalmultilife()
        return curve
    def lifeentropycurve(self, meshstart, meshstop, num_in_mesh):
        """
        Produces the life entropy curve (aka entropy summary function; https://arxiv.org/abs/1803.08304) of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of life entropy curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        tmp = self.stabilizedlifecurve(meshstart,meshstop,num_in_mesh)
        curve = -1*tmp*np.log(tmp)
        curve[np.isnan(curve)]=0
        return curve
    def multilifeentropycurve(self, meshstart, meshstop, num_in_mesh):
        """
        Produces the multiplicative life entropy curve of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of multiplicative life entropy curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        tmp = self.stabilizedmidlifecurve(meshstart,meshstop,num_in_mesh)
        curve[np.isnan(curve)]=0
        return curve
    def midlifeentropycurve(self, meshstart, meshstop, num_in_mesh):
        """
        Produces the midlife entropy curve of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of midlife entropy curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        tmp = self.stabilizedmultilifecurve(meshstart,meshstop,num_in_mesh)
        curve[np.isnan(curve)]=0
        return curve
    def custom_curve_at_t(self,fun,stat,t):
        Birth = self.Birth
        Death = self.Death
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
        """
        returns the sum of the lifespans of a diagram
        """
        Birth = self.Birth
        Death = self.Death
        return np.sum(Death-Birth)
    def totalmidlife(self):
        """
        returns the sum of midlifes in the diagram
        """
        Birth = self.Birth
        Death = self.Death
        return np.sum(Death+Birth)/2
    def totalmultilife(self):
        """
        returns the sum of multiplicative lifespans in the diagram
        """
        Birth = self.Birth
        Death = self.Death
        return np.sum(Death/Birth)
    def entropy(self):
        """
        returns the persistent entropy (https://arxiv.org/abs/1512.07613) of the diagram
        """
        Birth = self.Birth
        Death = self.Death
        return -np.sum((Death-Birth)/self.totallife()*np.log((Death-Birth)/self.totallife()))
    def gaussian_life(self, meshstart, meshstop, num_in_mesh, spread = 1):
        """
        Produces the gaussian life curve of the diagram
        
        Parameters:
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        spread: The width of the gaussian
        """
        Birth = self.Birth
        Death = self.Death
        T = np.linspace(meshstart, meshstop, num_in_mesh)
        Life = (Death - Birth)/np.sum(Death-Birth)
        B = np.matmul(np.array(Birth).reshape(-1,1),np.ones([1,num_in_mesh]))
        De = np.matmul(np.array(Death).reshape(-1,1),np.ones([1,num_in_mesh]))
        L = np.matmul(np.array(Life).reshape(-1,1),np.ones([1,num_in_mesh]))
        return np.sum(L*stats.norm.cdf((T-B)/spread)*(1-stats.norm.cdf((T-De)/spread)), axis=0)
    def gaussian_life_derivative(self, meshstart, meshstop, num_in_mesh, spread = 1):
        #Produces the gaussian life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        T = np.linspace(meshstart, meshstop, num_in_mesh)
        Life = (Death - Birth)/np.sum(Death-Birth)
        B = np.matmul(np.array(Birth).reshape(-1,1),np.ones([1,num_in_mesh]))
        De = np.matmul(np.array(Death).reshape(-1,1),np.ones([1,num_in_mesh]))
        L = np.matmul(np.array(Life).reshape(-1,1),np.ones([1,num_in_mesh]))
        return np.sum(L*(stats.norm.pdf((T-B)/spread) - stats.norm.pdf((T-B)/spread)*(stats.norm.cdf((T-De)/spread))-stats.norm.pdf((T-De)/spread)*(stats.norm.cdf((T-B)/spread))), axis=0)
    def gaussian_Betti(self, meshstart, meshstop, num_in_mesh,spread=1):
        #Produces the gaussian life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        T = np.linspace(meshstart, meshstop, num_in_mesh)
        Life = np.ones(Birth.shape[0])/Birth.shape[0]
        B = np.matmul(np.array(Birth).reshape(-1,1),np.ones([1,num_in_mesh]))
        De = np.matmul(np.array(Death).reshape(-1,1),np.ones([1,num_in_mesh]))
        L = np.matmul(np.array(Life).reshape(-1,1),np.ones([1,num_in_mesh]))
        return np.sum(L*stats.norm.cdf((T-B)/spread)*(1-stats.norm.cdf((T-De)/spread)), axis=0)
    def gaussian_midlife(self, meshstart, meshstop, num_in_mesh,spread=1):
        #Produces the gaussian life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        T = np.linspace(meshstart, meshstop, num_in_mesh)
        Life = (Birth+Death)/np.max(Birth+Death)
        B = np.matmul(np.array(Birth).reshape(-1,1),np.ones([1,num_in_mesh]))
        De = np.matmul(np.array(Death).reshape(-1,1),np.ones([1,num_in_mesh]))
        L = np.matmul(np.array(Life).reshape(-1,1),np.ones([1,num_in_mesh]))
        return np.sum(L*stats.norm.cdf((T-B)/spread)*(1-stats.norm.cdf((T-De)/spread)), axis=0)
    def gaussian_midlife_derivative(self, meshstart, meshstop, num_in_mesh, spread = 1):
        #Produces the gaussian life curve of the diagram
        #@param meshstart: The lowest value at which to begin the curve
        #@param meshstop: the highest value at which to stop the curve
        #@param num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        Birth = self.Birth
        Death = self.Death
        T = np.linspace(meshstart, meshstop, num_in_mesh)
        Life = (Death + Birth)/np.max(Death+Birth)
        B = np.matmul(np.array(Birth).reshape(-1,1),np.ones([1,num_in_mesh]))
        De = np.matmul(np.array(Death).reshape(-1,1),np.ones([1,num_in_mesh]))
        L = np.matmul(np.array(Life).reshape(-1,1),np.ones([1,num_in_mesh]))
        return np.sum(L*(stats.norm.pdf((T-B)/spread) - stats.norm.pdf((T-B)/spread)*(stats.norm.cdf((T-De)/spread))-stats.norm.pdf((T-De)/spread)*(stats.norm.cdf((T-B)/spread))), axis=0)