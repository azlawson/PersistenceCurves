"""
Created on Mon Nov  5 11:26:12 2018

@author: Austin Lawson (azlawson@uncg.edu)
"""

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
class Diagram:
    def __init__(self, Dgm, globalmaxdeath = None, infinitedeath=float("inf"), inf_policy="keep"):
        """
        Transforms a diagram (n by 2 NumPy array or Pandas DataFrame) into the Diagram class.
        
        Parameters
        Dgm: A persistence diagram, i.e. an n by 2 array or pandas dataframe
        
        globalmaxdeath: The maximum possible death value for the persistence process, used to replace infinite death value.
        
        infinitedeath: The value that represents an infinite death value, e.g. for Perseus that value is -1, and for gudhi the value is float("inf") which is default for this package
        
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
        self.diagram = np.stack([self.Birth, self.Death], axis = 1)
        self.shape = self.diagram.shape

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
        Birth = self.Birth.reshape([self.shape[0],1])
        Death = self.Death.reshape([self.shape[0],1])
        FUN = np.ones([self.shape[0],1])
        T = np.linspace(meshstart,meshstop,num_in_mesh)*np.ones([self.shape[0],num_in_mesh])
        curve=np.where(((T>=Birth) & (T<Death)),FUN,0).sum(axis=0)
        return curve

    def landscape(self, k, meshstart, meshstop, num_in_mesh):
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
        Birth = self.Birth.reshape([self.shape[0],1]).astype(np.float32)
        Death = self.Death.reshape([self.shape[0],1]).astype(np.float32)
        T = np.linspace(meshstart,meshstop,num_in_mesh)*np.ones([self.shape[0],num_in_mesh]).astype(np.float32)
        tmpB = (T-Birth).astype(np.float32)
        tmpD = (Death-T).astype(np.float32)
        tri=np.minimum(tmpB,tmpD).astype(np.float32)
        land=np.sort(tri,axis=0,)[self.shape[0]-k-1]
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
        Birth = self.Birth.reshape([self.shape[0],1])
        Death = self.Death.reshape([self.shape[0],1])
        FUN = Death-Birth
        T = np.linspace(meshstart,meshstop,num_in_mesh)*np.ones([self.shape[0],num_in_mesh])
        curve=np.where(((T>=Birth) & (T<Death)),FUN,0).sum(axis=0)
        return curve
    def deathcurve(self, meshstart, meshstop, num_in_mesh):
        """
        Produces the death curve of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of lifespan curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        Birth = self.Birth.reshape([self.shape[0],1])
        Death = self.Death.reshape([self.shape[0],1])
        FUN = Death
        T = np.linspace(meshstart,meshstop,num_in_mesh)*np.ones([self.shape[0],num_in_mesh])
        curve=np.where(((T>=Birth) & (T<Death)),FUN,0).sum(axis=0)
        return curve
    def birthcurve(self, meshstart, meshstop, num_in_mesh):
        """
        Produces the birth curve of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of lifespan curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        Birth = self.Birth.reshape([self.shape[0],1])
        Death = self.Death.reshape([self.shape[0],1])
        FUN = Birth
        T = np.linspace(meshstart,meshstop,num_in_mesh)*np.ones([self.shape[0],num_in_mesh])
        curve=np.where(((T>=Birth) & (T<Death)),FUN,0).sum(axis=0)
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
        Birth = self.Birth.reshape([self.shape[0],1])
        Death = self.Death.reshape([self.shape[0],1])
        FUN = 0.5*(Birth+Death)
        T = np.linspace(meshstart,meshstop,num_in_mesh)*np.ones([self.shape[0],num_in_mesh])
        curve=np.where(((T>=Birth) & (T<Death)),FUN,0).sum(axis=0)
        return curve
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
        Birth = self.Birth.reshape([self.shape[0],1])
        Death = self.Death.reshape([self.shape[0],1])
        FUN = Death/Birth
        FUN[FUN == float("inf")] = Death.max() 
        T = np.linspace(meshstart,meshstop,num_in_mesh)*np.ones([self.shape[0],num_in_mesh])
        curve=np.where(((T>=Birth) & (T<Death)),FUN,0).sum(axis=0)
        return curve
    def normalizedBetticurve(self, meshstart, meshstop, num_in_mesh):
        """
        Produces the normalized(normalized) Betti curve of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of normalized(normalized) life curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        curve = self.Betticurve(meshstart,meshstop,num_in_mesh)/self.totallife()
        return curve
    def normalizedlifecurve(self, meshstart, meshstop, num_in_mesh):
        """
        Produces the normalized(normalized) life curve of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of normalized(normalized) life curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        curve = self.lifecurve(meshstart,meshstop,num_in_mesh)/self.totallife()
        return curve
    def normalizedmidlifecurve(self, meshstart, meshstop, num_in_mesh):
        """
        Produces the normalized(normalized) midlife curve of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of normalized(normalized) midlife curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        curve = self.midlifecurve(meshstart,meshstop,num_in_mesh)/self.totalmidlife()
        return curve
    def normalizedmultilifecurve(self, meshstart, meshstop, num_in_mesh):
        """
        Produces the normalized(normalized) multiplicative life curve of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of normalized(normalized) multiplicative life curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        curve = self.multilifecurve(meshstart,meshstop,num_in_mesh)/self.totalmultilife()
        return curve
    def Bettientropycurve(self, meshstart, meshstop, num_in_mesh):
        """
        Produces the life entropy curve (aka entropy summary function; https://arxiv.org/abs/1803.08304) of the diagram
        
        Parameters:
        
        meshstart: The lowest value at which to begin the curve
        meshstop: the highest value at which to stop the curve
        num_in_mesh: The number of evenly spaced points between meshstart and meshstop at which to compute the curve values
        
        Output:
        num_in_mesh dimensional vector of life entropy curve values computed at num_in_mesh evenly spaced points starting at meshstart and ending at meshstop
        """
        tmp = self.normalizedBetticurve(meshstart,meshstop,num_in_mesh)
        tmp[tmp==0] = 1
        curve = -1*tmp*np.log(tmp)
        curve[np.isnan(curve)] = 0
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
        tmp = self.normalizedlifecurve(meshstart,meshstop,num_in_mesh)
        tmp[tmp==0] = 1
        curve = -1*tmp*np.log(tmp)
        curve[np.isnan(curve)] = 0
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
        tmp = self.normalizedmidlifecurve(meshstart,meshstop,num_in_mesh)
        tmp[tmp==0] = 1
        curve = -1*tmp*np.log(tmp)
        curve[np.isnan(curve)] = 0
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
        tmp = self.normalizedmultilifecurve(meshstart,meshstop,num_in_mesh)
        tmp[tmp==0] = 1
        curve = -1*tmp*np.log(tmp)
        curve[np.isnan(curve)] = 0
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
        will replace division by 0 with max death
        """
        Birth = self.Birth
        Death = self.Death
        FUN = Death/Birth
        FUN[FUN==float("inf")]=Death.max()
        out = np.sum(FUN)
        return out
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
        B = np.matmul(np.array(Birth).reshape(-1,1),np.ones([1,num_in_mesh]))
        De = np.matmul(np.array(Death).reshape(-1,1),np.ones([1,num_in_mesh]))
        return np.sum(stats.norm.cdf((T-B)/spread)*(1-stats.norm.cdf((T-De)/spread)), axis=0)
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
