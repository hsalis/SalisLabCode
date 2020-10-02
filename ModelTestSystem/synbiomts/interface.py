"""
Model interface to import models, write model wrappers, collect model output
and run statistics specific to model type

Copyright 2017 Alexander C. Reis, Howard M. Salis, all rights reserved.
  
"""

import numbers
import inspect
from functools import partial

'''
class Model(obj):
    def __init__(self):
        pass
'''

class Container(dict):
    __name__ = "Models"
    form = {}

    def __init__(self):
        pass

    def add(self, model, *args, **kargs):
        '''Register a model with the model Container. You may provide default
        arguments that will be passed automatically when calling the registered
        model. Fixed arguments can then be overriden at model call time.
        Inputs:
        model (object) = The model to which refer the name.
        argument       = One or more argument (and keyword argument) to pass
                         automatically to the registered model when called,
                         optional.

        The following code block is an example of how the toolbox is used. ::
            >>> def RBS_Calculator_v2_0(sequence, organism, temp):
            ...     # run model code here
            ...     print seq, organism, temperature
            ...
            >>> models = Container()
            >>> models.add(RBS_Calculator_v2_0,'ACTAGC',temp=37.0)
            >>> models.RBSCalc2('Escherichia coli')
            'ACTAGC' 'Escherichia coli' 37.0

        The registered model will be given the attributes :attr:`__name__`
        set to the name and :attr:`__doc__` set to the original model's
        documentation. The :attr:`__dict__` attribute will also be updated
        with the original model's instance dictionary, if any. '''
        
        pmodel = partial(model, *args, **kargs)
        pmodel.__doc__ = model.__doc__
        name = model.__name__
        ArgSpec = inspect.getargspec(model)
        pmodel.variables = ArgSpec[0]

        if hasattr(model, "__dict__") and not isinstance(model, type):
            # Some functions don't have a dictionary, in these cases
            # simply don't copy it. Moreover, if the model is actually
            # a class, we do not want to copy the dictionary.
            pmodel.__dict__.update(model.__dict__.copy())

        setattr(self, name, pmodel) # not currently used
        self[name] = pmodel
        self[name].set = False
        self.available = sorted(self.keys())

    def remove(self, model):
        '''Unregister model from the model Container.
        model (string or function) = Can be a string or the function itself
        '''
        name = model if isinstance(model,str) else model.__name__
        assert name in self.available, "{} is not registered with the Container.".format(name)
        delattr(self, name)
        self.pop(name)
        self.available = sorted(self.keys())

    def setform(self,modelNames,x,y,std,xScale='linear',yScale='linear',a1=None):

        options = ['linear','ln','log10']
        for model in modelNames:
            assert model in self.available, "Model {} not in container.".format(model)
        assert isinstance(x,str), "x should be a string, is {}.".format(type(x))
        assert isinstance(y,str), "y should be a string, is {}.".format(type(y))
        assert isinstance(std,str), "std should be a string, is {}.".format(type(std))
        assert xScale in options, "xScale should be one of: {}".format(options)
        assert yScale in options, "yScale should be one of: {}".format(options)
        if not a1 is None:
            assert isinstance(a1,(int,float)), "a1, {}, is the slope and should be a number.".format(a1)

        for name in modelNames:
            self[name].set = True
            self[name].x = x
            self[name].y = y
            self[name].std = std
            self[name].xScale = xScale
            self[name].yScale = yScale
            self[name].a1 = float(a1) if isinstance(a1,(int,float)) else None
            # self.form[name] = {'x': x, 'y': y, 'yScale': yScale,  'a1': float(a1)}

    def changeName(self,oldName,newName):
        self[newName] = self.pop(oldName)
        self.available = sorted(self.keys())

    def __add__(self,another):
        assert another.__name__ == "Models", "Must add two Models to combine."
        for name,pmodel in another.iteritems():
            self.add(name,pmodel)
        return self

    def __sub__(self,another):
        assert another.__name__ == "Models", "Must subtract two Models to remove."
        for name,pmodel in another.iteritems():
            self.remove(name,pmodel)
        return self


if __name__ == "__main__":
    pass