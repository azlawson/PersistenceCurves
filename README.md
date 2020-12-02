# PersistenceCurves
A python package for computing [Persistence Curves](https://arxiv.org/abs/1904.07768)

## Brief Introduction
Computational Topology is a field of mathematics concerned with examining and utilizing the shape of a dataset to discern the shape of the underlying space. The main tool of this field is Persistent Homology. Using this tool on a dataset yields a visual summary called a Persistence Diagram, which are multisets of points. We can define a metric on these diagrams called the bottleneck distance. These diagrams are stable with respect to this distance in the sense that a small change in the dataset leads to a small change in the diagrams. The set of all persistence diagrams with this bottleneck distance forms a metric space. However, performing machine learning tasks with this space is difficult. For this reason, much work has been done towards creating useful, clever vectorizations of these diagrams that are compatible with machine learning algorithms. This package is named after the [Persistence Curve framework](https://arxiv.org/abs/1904.07768), which is a generalization of many such vectorizations.

##  The Diagram Class

The sole class of this package is ``Diagram``. This package assumes the user has already calculated the persistence diagram(s) of interest. Diagrams are collection of ordered pairs (b,d) (birth and death respectively) where d>b and d can take the value infinity. Essentially, a diagram is an array or data frame of shape (n,2). Suppose Dgm is a diagram. The code below transforms D to the Diagram class. 

```python
import persistencecurves as pc
D = pc.Diagram(Dgm =Dgm, globalmaxdeath = None, infinitedeath=None, inf_policy="keep")
```

### Global Values
This class has various global variables such as 

``diagram`` = Dgm

``Birth`` = All the Birth values of Dgm

``Death`` = All the Death values of Dgm

``globalmaxdeath`` = If one is considering multiple samples of a space where a largest possible death value, that value should be input here. For example, images have a global max death of 255.

``infinitedeath`` = The value signifying an infinite death. For most softwares, this value is -1. Thus if left unset, any negative death value is assumed to be infinite.

``shape`` = a tuple (n,2) where n is the number of points in the diagram

 ### The Built-in Curves
 The bulk of this package is the curve calculation. The Diagram class comes with some pre-set Persistence Curves as well as the option for a custom curve. All of the preset curves, except landscape, have as inputs **meshstart**, **meshstop**, and **num_in_mesh**. These numbers follow the same concept as NumPy's linspace where **meshstart** is the leftmost endpoint, **meshstop** is the right endpoint and **num_in_mesh** is the number of points in the mesh. The method returns a vector of length **num_in_mesh**. For a Diagram arising from an image we can (and do) use:

```python
D.Betticurve(meshstart=0,meshstop=255,num_in_mesh=256)
```

For landscapes, we have an additional argument **k** for the level. 

```python
D.landscape(k=1,meshstart=0,meshstop=255,num_in_mesh=256)
```

### The Custom Curve
A Persistence Curve has two main components. An inner function **fun** and a **statistic** (recall that a statistic is just a function of a multi-set). As such the ``customcurve`` method requires an input of both. 

The inner function **fun** should be a function of three variables b,d,and t that returns a single value. The statistic can be any function applicable to an entire array. Some examples follow:

The code below returns exactly the Betti curve.

```python
def fun(b,d,t):
    return 1
D.custom_curve(fun = fun, statistic = np.sum, meshstart=0, meshstop=255, num_in_mesh=256)
```

