import os
import sys
import numpy as np
mainpath = os.path.dirname(os.path.realpath(__file__)) + "/../src/"
sys.path.append(mainpath) # add dmrgpy library

import qtwrap # import the library with simple wrappers
from dmrgpy import spinchain
app = qtwrap.App() # this is the main interface



def get_dimension(spins):
    """
    Return the dimension of the Hilbert space of the spin chain,
    return None if the dimension is bigger than 10000
    """
    out = 1
    for s in spins: 
        out *= s
        if out>10000: return None
    return out



def get_mode():
    """Decide if DMRG should be used according to the size of the matrix"""
    if app.getbox("preferred_mode")=="DMRG": return "DMRG"
    else: # try to use ED
      spins = get_spins() # get the spins
      out = 1 # initialize
      for s in spins:
          out *= s 
          if out>1000: return "DMRG"
    print("ED mode can be used")
    return "ED"



# set the logo
def set_logo():
  path = os.path.dirname(os.path.realpath(__file__))+"/../"
  logopath = path+"/logo/spinflare.png"
  app.set_image("spinflare_logo",logopath)
  app.set_image("exchange_image",path+"/logo/exchange_image.png")
  app.set_image("field_image",path+"/logo/field_image.png")
  app.set_image("magnetization_image",path+"/logo/magnetization_image.png")
  app.set_image("dynamical_image",path+"/logo/dynamical_image.png")

set_logo()

def execute_script(name):
  path = os.path.dirname(os.path.realpath(__file__)) + "/../utilities/"
  os.system("python "+path+name+" &") # execute that script


def text2spin(app):
    """Convert a text into a list of spins"""
    text = app.gettext("spin_list") # get the text with the spins
    out = text.split("\n") # remove next line
    out2 = [] # empty list
    for o in out:
        if len(o)>0: out2.append(o) # store
    t2s = {"1/2":2,"1":3,"3/2":4,"2":5,"5/2":6}
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




def set_data(name):
    """
    Set the latest data in the window
    """
    m = np.genfromtxt(name) # get the data
    np.savetxt("DATA.OUT",m,fmt='%.5f') # save the data
    text = open("DATA.OUT").read() # read the text
    text = "# "+name+"\n" + text.replace(" ","   ") # name of the file
    app.settext("data_box",text) # set that text




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
    elif name=="1/2 & 1":  spins = [2+i%2 for i in range(ns)]
    elif name=="1 & 3/2":  spins = [3+i%2 for i in range(ns)]
    else: raise
    spin2text(app,spins) # write spins in the table
    initialize_exchange() # initialize the exchange constants
    initialize_field() # initialize the exchange constants


def initialize_exchange():
    """Initialize the table of exchange constants"""
    spins = text2spin(app) # get the spins
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
    ns = len(spins)
    h = 0 # initialize
    # add exchange
    for i in range(ns): # loop
      for j in range(ns): # loop
          cij = fj(i,j) # get matrix
          Si = [sc.Sx[i],sc.Sy[i],sc.Sz[i]]
          Sj = [sc.Sx[j],sc.Sy[j],sc.Sz[j]]
          for ii in range(3):
            for jj in range(3):
                if cij[ii,jj]!=0.0: 
                  h = h + cij[ii,jj]*Si[ii]*Sj[jj]
    # add zeeman
    for i in range(ns): # loop
        Si = [sc.Sx[i],sc.Sy[i],sc.Sz[i]]
        bi = fb(i)
        for ii in range(3):
            if bi[ii]!=0.0:
                h = h + bi[ii]*Si[ii]
    sc.set_hamiltonian(h)
    sc.maxm = int(app.get("maxm"))
    sc.kpmmaxm = int(app.get("maxm"))
    sc.nsweeps = int(app.get("nsweeps"))
    return sc





def get_static_correlator():
    """Return the static correlator"""
    sc = getsc(app) # get the spin chain object
    pairs = get_pairs(app,sc.ns,"static_arrangement") # get the pairs
    cname = app.getbox("static_type") # actual operator for the correlator
    # get the correlator
    if cname=="XX": ops = [sc.Sx[p[0]]*sc.Sx[p[1]] for p in pairs]
    elif cname=="YY": ops = [sc.Sy[p[0]]*sc.Sy[p[1]] for p in pairs]
    elif cname=="ZZ": ops = [sc.Sz[p[0]]*sc.Sz[p[1]] for p in pairs]
    elif cname=="SS": 
        ops = [sc.Sz[p[0]]*sc.Sz[p[1]] + sc.Sx[p[0]]*sc.Sx[p[1]] + sc.Sy[p[0]]*sc.Sy[p[1]] for p in pairs]
    else: raise
    out = np.array([sc.vev(o) for o in ops])
#    out = sc.get_correlator(pairs,name=cname) 
    np.savetxt("STATIC_CORRELATOR.OUT",np.array([range(len(pairs)),out.real,out.imag]).T)
    set_data("STATIC_CORRELATOR.OUT")
    execute_script("sf-correlator --input STATIC_CORRELATOR.OUT   --ylabel "+cname)





def get_dynamical_correlator_single():
    """Dynamic correlator in a single site"""
    sc = getsc(app) # get the spin chain object
    cname = app.getbox("dynamic_type_single") # operator
    delta = app.get("smearing_dynamic") # delta
    ii = int(app.get("dynamic_site_i_single")) # delta
    if cname=="XX": name = (sc.Sx[ii],sc.Sx[ii])
    elif cname=="YY": name = (sc.Sy[ii],sc.Sy[ii])
    elif cname=="ZZ": name = (sc.Sz[ii],sc.Sz[ii])
    else: raise
    (es,ds) = sc.get_dynamical_correlator(delta=delta,name=name)
    np.savetxt("DYNAMICAL_CORRELATOR.OUT",
            np.array([es,ds.real,ds.imag]).T)
    
    set_data("DYNAMICAL_CORRELATOR.OUT")
    execute_script("sf-dynamical_correlator --input DYNAMICAL_CORRELATOR.OUT  --ylabel "+cname)


def is_parallel():
    if app.getbox("parallelization")=="Yes": return True
    else: return False

def multievaluate(f,ps):
    """Evaluate a function for multiple parameters"""
    if is_parallel():
        from dmrgpy import parallel
        parallel.cores = parallel.maxcpu # set to the maximum value
        fs = [lambda: f(p) for p in ps] # define all the functions
        return parallel.multicall(fs)
    else: return [f(p) for p in ps] # return all the results


def get_dynamical_correlator_map():
    """Dynamic correlator in a single site"""
    sc = getsc(app) # get the spin chain object
    cname = app.getbox("dynamic_type_map") # operator
    delta = app.get("smearing_dynamic") # delta
    pairs = get_pairs(app,sc.ns,"dynamic_arrangement") # get the pairs
    def getAB(p):
        ii,jj = p[0],p[1]
        if cname=="XX": return (sc.Sx[ii],sc.Sx[jj])
        elif cname=="YY": return (sc.Sy[ii],sc.Sy[jj])
        elif cname=="ZZ": return (sc.Sz[ii],sc.Sz[jj])
        else: raise
    fo = open("DYNAMICAL_CORRELATOR_MAP.OUT","w")
    ip = 0
    #### Parallel execution ###
    if app.getbox("parallelization")=="Yes": # parallel calculation
        def fm(p): # function to call
            name = getAB(p)
            sc0 = sc.clone() # clone directory
            def f(): # function to call
                print("Computing DYNCORR for",p)
                out = sc0.get_dynamical_correlator(delta=delta,name=name)
                sc0.clean()
                return out
            return f # return dummy function
        from dmrgpy import parallel
        parallel.cores = parallel.maxcpu # set to the maximum value
        outs = parallel.multicall([fm(p) for p in pairs])
    #### Serial execution ####
    elif app.getbox("parallelization")=="No":
      sc.get_gs() # compute ground state
      def f(p): # function to call
        print("Computing DYNCORR for",p)
        name = getAB(p)
        return sc.get_dynamical_correlator(delta=delta,name=name)
      outs = [f(p) for p in pairs] # compute all the outputs
    else: raise # something wrong
    for ip in range(len(outs)): # loop over outputs
      (es,ds) = outs[ip] # get those 
      for (ei,di) in zip(es,ds):
          fo.write(str(ip)+"   ")
          fo.write(str(ei)+"   ")
          fo.write(str(di.real)+"   ")
          fo.write(str(di.imag)+"\n")
      fo.flush()
    fo.close()
    set_data("DYNAMICAL_CORRELATOR_MAP.OUT")
    execute_script("sf-dynamical_correlator-map --input DYNAMICAL_CORRELATOR_MAP.OUT  --ylabel "+cname)







def get_magnetization():
    sc = getsc(app) # get the spin chain object
#    mx,my,mz = sc.get_magnetization(mode=get_mode()) # get the magnetization
    inds = np.array(range(sc.ns))
    name = app.getbox("magnetization_type")
    if name=="X": ops = sc.Sx
    elif name=="Y": ops = sc.Sy 
    elif name=="Z": ops = sc.Sz
    else: raise
    mi = np.array([sc.vev(o) for o in ops]).real
    np.savetxt("MAGNETIZATION.OUT",np.array([inds,mi]).T)
    set_data("MAGNETIZATION.OUT")
    execute_script("sf-magnetization --input MAGNETIZATION.OUT --ylabel "+name)


def get_gs_convergence():
    """Get the ground state converge with bond dimension"""
    maxm_min = int(app.get("maxm_min"))
    maxm_max = int(app.get("maxm_max"))
    maxm_step = int(app.get("maxm_step"))
    ms = np.array(range(maxm_min,maxm_max,maxm_step))
    es = []
    for m in ms:
        sc = getsc(app) # get the spin chain object
        sc.maxm = int(m) # overwrite bond dimension
        es.append(sc.gs_energy()) # compute ground state energy
    np.savetxt("ENERGY_VS_MAXM.OUT",np.array([ms,es]).T)
    set_data("ENERGY_VS_MAXM.OUT")
    execute_script("sf-gs-energy-maxm ENERGY_VS_MAXM.OUT")



def get_excited_states():
    """Get the ground state converge with bond dimension"""
    sc = getsc(app) # get the spin chain object
    n = int(app.get("num_excited"))
    es = sc.get_excited(n=n,mode=get_mode()) # get excited states
    es = np.sort(es) # sort energies
    ndim = get_dimension(get_spins()) # dimension of the system
    if ndim is not None: # if not too big dimension
        if ndim<=n: es = es[0:ndim] # take only this energies
    rmgs = app.getbox("substract_gs_excited")=="Yes"
    args = "" # arguments for the plotting
    if rmgs: 
        es -= es[0] # remove ground state energy
        args += " --ylabel GS" 
    np.savetxt("EXCITED.OUT",np.array(es))
    set_data("EXCITED.OUT")
    execute_script("sf-excited-states --input EXCITED.OUT"+args)




    
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
    elif m=="Onsite":
      pairs = [(i,i) for i in range(n)]
    else: raise
    return pairs





signals = dict()
signals["get_static_correlator"] = get_static_correlator
signals["initialize_spins"] = initialize_spins
signals["get_dynamical_correlator_single"] = get_dynamical_correlator_single
signals["get_dynamical_correlator_map"] = get_dynamical_correlator_map
signals["get_magnetization"] = get_magnetization
signals["get_gs_convergence"] = get_gs_convergence
signals["get_excited_states"] = get_excited_states

app.connect_clicks(signals)


# now initialize the interface
initialize_interface()


app.run()
