import os
import sys
import numpy as np
mainpath = os.path.dirname(os.path.realpath(__file__)) + "/../src/"
sys.path.append(mainpath) # add dmrgpy library

import qtwrap # import the library with simple wrappers
from dmrgpy import spinchain
app = qtwrap.App() # this is the main interface


def execute_script(name):
  path = os.path.dirname(os.path.realpath(__file__)) + "/../utilities/"
  os.system(path+name+" &") # execute that script


def getsc(app):
    """Get the spin chain"""
    ns = int(app.get("nsites")) # number of sites
    spins = [2 for i in range(ns)] # define the spins
    sc = spinchain.Spin_Hamiltonian(spins) # return spin chain object
    # now initialize the Hamiltonian
    fj = get_exchange_function(app) # add exchange field
    fb = get_field_function(app) # get the function that adds b field
    sc.set_exchange(fj) # add the exchange
    sc.set_fields(fb) # set the magnetic field
    sc.maxm = int(app.get("maxm"))
    sc.nsweeps = int(app.get("nsweeps"))
    return sc


def get_exchange_function(app):
    """Return the exchange function"""
    jtype = app.getbox("coupling_type") # type of coupling
    jmod = app.getbox("coupling_modulation") # type of modulation
    if jtype=="Heisenberg": m = np.identity(3) # identity matrix
    elif jtype=="Ising": m = np.zeros((3,3)) ; m[2,2] = 1.0 # Ising
    elif jtype=="XXZ": m = np.identity(3) ; m[2,2] = app.get("XXZJ") # XXZ
    else: raise # not implemented
    def fj(i,j):
        if abs(i-j)==1: return m
        else: return m*0
    return fj

def get_field_function(app):
    """Function that returns the magnetic field"""
    return lambda i: [0.0,0.0,0.0]



def get_static_correlator():
    """Return the static correlator"""
    sc = getsc(app) # get the spin chain object
    pairs = get_pairs(app,sc.ns,"type_static")# get the pairs
    out = sc.get_correlator(pairs)
    np.savetxt("STATIC_CORRELATOR.OUT",np.array([range(len(pairs)),out]).T)
    execute_script("sf-correlator STATIC_CORRELATOR.OUT")


def get_pairs(app,n,name):
    """Return the desired pairs"""
    m = app.getbox(name)
    if m=="From edge":
      pairs = [(0,i) for i in range(n)]
    elif m=="From bulk":
      pairs = [(n//2,i) for i in range(n)]
    elif m=="First neighbor right":
      pairs = [(i,i+1) for i in range(n-1)]
    elif m=="First neighbor left":
      pairs = [(i+1,i) for i in range(n-1)]
    else: raise
    return pairs





signals = dict()
signals["get_static_correlator"] = get_static_correlator

app.connect_clicks(signals)
app.run()
