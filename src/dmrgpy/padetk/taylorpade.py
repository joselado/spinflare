from scipy.interpolate import pade
import numpy as np

raise # this does not work

def fit(zs,us,m=10,n=10):
    """Fit to a Pade approximant"""
    us = 1./us
    e_exp = np.polyfit(zs,us,(m+n)*2)
    p, q = pade(np.flip(e_exp), m,n)
#    e_exp = np.flip(e_exp)
    print(e_exp)
    e_poly = np.poly1d(e_exp)
    def f(z): return 1/(p(z)/q(z))
    return f
