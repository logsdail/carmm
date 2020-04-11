import numpy as np

"""Module for wrapping an array without being an array.

This can be desirable because we would like atoms.cell to be like an array,
but we don't like atoms.positions @ atoms.cell to also be a cell, and as far
as I/we know, it's not possible to subclass array without that kind of thing
happening.

For most methods and attributes we can just bind the name as a property:

    @property
    def reshape(self):
        return np.asarray(self).reshape

For 'in-place' operations like += we need to make sure to return self
rather than the array though:

    def __iadd__(self, other):
        np.asarray(self).__iadd__(other)
        return self

This module provides the @arraylike decorator which does these things
for all the interesting ndarray methods.
"""


inplace_methods = ['__iadd__', '__imul__', '__ipow__', '__isub__',
                   '__itruediv__', '__imatmul__']

forward_methods = ['__abs__', '__add__', '__contains__', '__eq__',
                   '__ge__', '__getitem__', '__gt__', '__hash__',
                   '__iter__', '__le__', '__len__', '__lt__',
                   '__mul__', '__ne__', '__neg__', '__pos__',
                   '__pow__', '__radd__', '__rmul__', '__rpow__',
                   '__rsub__', '__rtruediv__', '__setitem__',
                   '__sub__', '__truediv__']


if hasattr(np.ndarray, '__matmul__'):
    forward_methods += ['__matmul__', '__rmatmul__']


def forward_inplace_call(name):
    arraymeth = getattr(np.ndarray, name)
    def f(self, obj):
        a = self.__array__()
        arraymeth(a, obj)
        return self
    # use update_wrapper()?
    f.__name__ = name
    f.__qualname__ = name
    return f


def wrap_array_attribute(name):
    wrappee = getattr(np.ndarray, name)
    if wrappee is None:  # For example, __hash__ is None
        assert name == '__hash__'
        return None
    def attr(self):
        array = np.asarray(self)
        return getattr(array, name)
    attr.__name__ = wrappee.__name__
    attr.__qualname__ == wrappee.__qualname__
    return property(attr)


def arraylike(cls):
    """Decorator for being like an array without being an array.

    Poke decorators onto cls so that getting an attribute
    really gets that attribute from the wrapped ndarray.

    Exceptions are made for in-place methods like +=, *=, etc.
    These must return self since otherwise myobj += 1 would
    magically turn into an ndarray."""
    for name in inplace_methods:
        if hasattr(np.ndarray, name) and not hasattr(cls, name):
            meth = forward_inplace_call(name)
        setattr(cls, name, meth)

    allnames = [name for name in dir(np.ndarray) if not name.startswith('_')]
    for name in forward_methods + allnames:
        if hasattr(cls, name) and not name.startswith('_'):
            continue  # Was overridden -- or there's a conflict.

        prop = wrap_array_attribute(name)
        setattr(cls, name, prop)
    return cls
