


import numpy as np
import matplotlib.pyplot as plt

from typing import Union, List, Tuple, Dict, Any, Optional, Callable, Iterable, Sequence
from numpy.typing import ArrayLike



def make_iterable( x ):
    if isinstance( x, str ) or not isinstance( x, Iterable ):
        return [ x ]
    return x

def multiple_formatter(
        denominator=2,
        number=np.pi,
        latex='\\pi',
        ):
    def gcd(a, b):
        while b:
            a, b = b, a%b
        return a
    def _multiple_formatter(x, pos):
        den = denominator
        num = np.int(np.rint(den*x/number))
        com = gcd(num,den)
        (num,den) = (int(num/com),int(den/com))
        if den==1:
            if num==0:
                return r'$0$'
            if num==1:
                 return r'$%s$'%latex
            elif num==-1:
                return r'$-%s$'%latex
            else:
                return r'$%s%s$'%(num,latex)
        else:
            if num==1:
                return r'$\frac{%s}{%s}$'%(latex,den)
            elif num==-1:
                return r'$\frac{-%s}{%s}$'%(latex,den)
            else:
                return r'$\frac{%s%s}{%s}$'%(num,latex,den)
    return _multiple_formatter

class Multiple:
    def __init__(
            self,
            denominator=2,
            number=np.pi,
            latex='\\pi',
            ):
        self.denominator = denominator
        self.number = number
        self.latex = latex
    def locator(self):
        return plt.MultipleLocator(
            self.number / self.denominator
            )
    def formatter(self):
        return plt.FuncFormatter(
            multiple_formatter(
                self.denominator,
                self.number,
                self.latex,
                )
            )
    
def strictly_increasing( L ):
    return all(
        x < y
        for x, y in zip(
            L,
            L[1:],
        )
    )