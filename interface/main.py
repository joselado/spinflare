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


def text2spin(app):
    """Convert a text into a list of spins"""
    text = app.gettext("spin_list") # get the text with the spins
    out = text.split("\n") # remove next line
    out2 = [] # empty list
    for o in out:
        if len(o)>0: out2.append(o) # store
    t2s = {"1/2":2,"1":3,"3/2":4,"2":5,"3/2":6}
    spins = [t2s[t] for t in out2] # spins using a dictionary
    return spins


def text2list(app,name):
    """Convert a text into a list of spins"""
    text = app.gettext(name) # convert text into a list
    out = text.split("\n") # remove next line
    out2 = [] # empty list
    for o in out:
        if len(o)>0: 
            o = o.split() # split list
            o2 = [float(io) for io in o]
            out2.append(o2) # store
    return out2 # return list





def spin2text(app,spins):
    """Convert a list of spins into a text"""
    text = "" # empty list
    s2t = {2:"1/2",3:"1",4:"3/2",5:"2",6:"5/2"}
    for s in spins:
        text += s2t[s]+"\n"
    app.settext("spin_list",text) # set that text






def initialize_spins():
    """Initialize spins accordinf to a scheme"""
    ns = int(app.get("nsites")) # number of sites
    name = app.getbox("spin_combination")
    if name=="1/2":  spins = [2 for i in range(ns)]
    elif name=="1":  spins = [3 for i in range(ns)]
    elif name=="3/2":  spins = [4 for i in range(ns)]
    elif name=="2":  spins = [5 for i in range(ns)]
    elif name=="5/2":  spins = [6 for i in range(ns)]
    else: raise
    spin2text(app,spins) # write spins in the table
    initialize_exchange() # initialize the exchange constants
    initialize_field() # initialize the exchange constants


def initialize_exchange():
    """Initialize the table of exchange constants"""
    spins = text2spin(app) # get the spins
    ctype = app.getbox("coupling_type") # type of coupling
    def getlist(jc):
        text = "" # initialize
        for i in range(len(spins)-1): # loop
            dd = app.get("dimerization_exchange")
            jj = jc*(1.0 + dd*(-1)**i)
            text += str(i)+"  "+str(i+1)+"   "+str(jj)+"\n" # store
        return text
    app.settext("jxx_table",getlist(app.get("Jxx"))) # set that text
    app.settext("jyy_table",getlist(app.get("Jyy"))) # set that text
    app.settext("jzz_table",getlist(app.get("Jzz"))) # set that text




def initialize_field():
    """Initialize the magnetic field"""
    spins = text2spin(app) # get the spins
    ns = len(spins) # number of spins
    def getlist(b):
        text = "" # initialize
        if b==0.0: return ""
        name = app.getbox("location_field")
        if name=="Everywhere":
          for i in range(len(spins)): # loop
            if app.getbox("modulation_field")=="Ferromagnetic": bb = b
            else: bb = b*(-1)**i
            text += str(i)+"  "+str(bb)+"\n" # store
          return text
        elif name=="First": return "0  "+str(b)+"\n"
        elif name=="Last": return str(ns-1)+"  "+str(b)+"\n"
        else: raise
    app.settext("bx_table",getlist(app.get("Bx"))) # set that text
    app.settext("by_table",getlist(app.get("By"))) # set that text
    app.settext("bz_table",getlist(app.get("Bz"))) # set that text







def get_exchange_function(app):
    """Return the function that creates the exchange constants"""
    spins = text2spin(app) # get the spins
    ns = len(spins) # number of spins
    # create matrices
    js = [[np.zeros((3,3)) for i in range(ns)] for j in range(ns)]
    # loop over the three coponents
    for l in text2list(app,"jxx_table"): # loop
        i,j,c = int(l[0]),int(l[1]),float(l[2])
        js[i][j][0,0] = c # store
        js[j][i][0,0] = c # store
    for l in text2list(app,"jyy_table"): # loop
        i,j,c = int(l[0]),int(l[1]),float(l[2])
        js[i][j][1,1] = c # store
        js[j][i][1,1] = c # store
    for l in text2list(app,"jzz_table"): # loop
        i,j,c = int(l[0]),int(l[1]),float(l[2])
        js[i][j][2,2] = c # store
        js[j][i][2,2] = c # store
    return lambda i,j: js[i][j] # return function


def get_field_function(app):
    """Return the function that creates the magnetic field"""
    spins = text2spin(app) # get the spins
    ns = len(spins) # number of spins
    # create matrices
    bs = [np.zeros(3) for i in range(ns)]
    # loop over the three coponents
    for l in text2list(app,"bx_table"): # loop
        i,c = int(l[0]),float(l[1])
        bs[i][0] = c # store
    for l in text2list(app,"by_table"): # loop
        i,c = int(l[0]),float(l[1])
        bs[i][1] = c # store
    for l in text2list(app,"bz_table"): # loop
        i,c = int(l[0]),float(l[1])
        bs[i][2] = c # store
    return lambda i: bs[i] # return function






def get_spins():
    """
    Return the spins
    """
    spins = text2spin(app) # get the spins
    return spins




def getsc(app):
    """Get the spin chain"""
    spins = get_spins() # return the spins
    sc = spinchain.Spin_Hamiltonian(spins) # return spin chain object
    # now initialize the Hamiltonian
    fj = get_exchange_function(app) # add exchange field
    fb = get_field_function(app) # get the function that adds b field
    sc.set_exchange(fj) # add the exchange
    sc.set_fields(fb) # set the magnetic field
    sc.maxm = int(app.get("maxm"))
    sc.nsweeps = int(app.get("nsweeps"))
    return sc





def get_static_correlator():
    """Return the static correlator"""
    sc = getsc(app) # get the spin chain object
    pairs = get_pairs(app,sc.ns,"static_arrangement") # get the pairs
    cname = app.getbox("static_type") # actual operator for the correlator
    out = sc.get_correlator(pairs,name=cname) # get the correlator
    np.savetxt("STATIC_CORRELATOR.OUT",np.array([range(len(pairs)),out.real,out.imag]).T)
    execute_script("sf-correlator STATIC_CORRELATOR.OUT")





def get_dynamical_correlator_single():
    """Dynamic correlator in a single site"""
    sc = getsc(app) # get the spin chain object
    cname = app.getbox("dynamic_type_single") # operator
    delta = app.get("smearing_dynamic") # delta
    ii = int(app.get("dynamic_site_i_single")) # delta
    (es,ds) = sc.get_dynamical_correlator(delta=delta,name=cname,i=ii,j=ii)
    np.savetxt("STATIC_CORRELATOR.OUT",
            np.array([es,ds.real,ds.imag]).T)
    
    execute_script("sf-dynamical_correlator DYNAMICAL_CORRELATOR.OUT")


def get_magnetization():
    sc = getsc(app) # get the spin chain object
    mx,my,mz = sc.get_magnetization() # get the magnetization
    inds = np.array(range(sc.ns))
    name = app.getbox("magnetization_type")
    if name=="X": mi = mx
    elif name=="Y": mi = my
    elif name=="Z": mi = mz
    else: raise
    np.savetxt("MAGNETIZATION.OUT",np.array([inds,mi]).T)
    execute_script("sf-magnetization MAGNETIZATION.OUT")

    
def initialize_interface():
    """Initialize interface"""
    initialize_spins() # initialize the spins






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
signals["initialize_spins"] = initialize_spins
signals["get_dynamical_correlator_single"] = get_dynamical_correlator_single
signals["get_magnetization"] = get_magnetization

app.connect_clicks(signals)


# now initialize the interface
initialize_interface()


app.run()
