from scipy.sparse import issparse
from scipy.sparse import csc_matrix as csc
import scipy.linalg as dlg
import scipy.sparse.linalg as slg
import numpy as np


maxsize = 3000



def braket_wAw(w,A,wi=None):
  """
  Compute the braket of a wavefunction
  """
  if wi is None: wi = w
  if issparse(A): # sparse matrices
    return (np.conjugate(wi)@A@w) # modern way
  else: # matrices and arrays
      if len(w.shape)==1: return (np.conjugate(wi)@A@w) # modern way
      else: return (np.conjugate(wi)@A@w)[0,0] # modern way






def braket_ww(w,wi):
  """
  Compute the braket of two wavefunctions
  """
  w = matrix2vector(w) # convert to vector
  wi = matrix2vector(wi) # convert to vector
  return (np.conjugate(w)@wi) # modern way




def disentangle_manifold(wfs,A):
  """
  Disentangles the wavefunctions of a degenerate manifold
  by expressing them in terms of eigenvalues of an input operator
  """
  ma = get_representation(wfs,A) # get the matrix form of the operator
  wfsout = [] # empty list
  evals,evecs = dlg.eigh(ma) # diagonalize

  evecs = evecs.transpose() # transpose eigenvectors
  for v in evecs: # loop over eigenvectors
    wf = wfs[0]*0.0j
    for (i,iv) in zip(range(len(v)),v): # loop over components
      wf += iv*wfs[i] # add contribution
    wfsout.append(wf.copy()) # store wavefunction
  return wfsout



def get_representation(wfs,A):
  """
  Gets the matrix representation of a certain operator
  """
  n = len(wfs) # number of eigenfunctions
  ma = np.zeros((n,n),dtype=np.complex) # representation of A
  sa = csc(A) # sparse matrix
  for i in range(n):
    vi = csc(np.conjugate(wfs[i])) # first wavefunction
    for j in range(n):
      vj = csc(wfs[j]).transpose() # second wavefunction
      data = (vi@sa@vj).todense()[0,0]
      ma[i,j] = data
  return ma





## routines for diagonalization ##

error = 1e-7



accelerate = False

def eigh(m):
    """Wrapper for linalg"""
    m = todense(m)
    return dlg.eigh(m)

def eigvalsh(m):
    """Wrapper for linalg"""
    m = todense(m)
    return dlg.eigvalsh(m)


def matrix2vector(v):
    """Transform a matrix into a vector"""
    if issparse(v): # sparse matrix
      v = v.todense() # convert to conventional matrix
    v = np.array(v) # convert to array
    if len(v.shape)==1: return v
    else: return v.reshape(v.shape[0]*v.shape[1])


def ground_state(h,nmax=maxsize):
  """Get a ground state"""
  info = False
  if h.shape[0]>nmax:
    if info: print("Calling ARPACK")
    eig,eigvec = slg.eigsh(h,k=10,which="SA",maxiter=100000)
    eig = np.sort(eig)
  else:
    if info: print("Full diagonalization")
    eig,eigvec = dlg.eigh(todense(h))
  return eig[0],eigvec.transpose()[0]


def todense(m):
    """Turn a matrix dense"""
    if issparse(m):
        if m.shape[0]>maxsize: raise
        else: return m.todense()
    else: return m


def lowest_eigenvalues(h,n=10,nmax=maxsize):
  """Get a ground state"""
  info = False
  if h.shape[0]>nmax: # for sparse use arpack
      eig,vs = lowest_states(h,n=n)
  else:
    if info: print("Full diagonalization")
    if ishermitian(h):
      eig = dlg.eigvalsh(h.todense())
    else:
      eig = dlg.eigvals(h.todense())
      eig = [y for (x,y) in sorted(zip(eig.real,eig))]
  return eig[0:n]


def lowest_states(h,n=10,nmax=maxsize):
  """Get a ground state"""
  info = False
  if h.shape[0]>nmax:
    if info: print("Calling ARPACK")
    if ishermitian(h): # Hermitian matrix
      eig,eigvec = slg.eigsh(h,k=n,which="SA",maxiter=100000)
      eig,eigvec = sorteigen(eig,eigvec.T)
      return (eig,eigvec)
    else: 
      eig,eigvec = slg.eigs(h,k=n,which="SR",maxiter=100000)
      eig,eigvec = sorteigen(eig,eigvec.T)
      return (eig,eigvec)
  else:
    if info: print("Full diagonalization")
    if ishermitian(h): # Hermitian matrix
      eig,vs = dlg.eigh(h.todense())
      return eig[0:n],vs.T[0:n] 
    else: # non Hermitian matrix
      eig,vs = dlg.eig(h.todense())
      eig,vs = sorteigen(eig,vs.T)
      return eig[0:n],vs[0:n]

lowest_eigenvectors = lowest_states

def sorteigen(eig,vs):
    """Return sorted eigenvalues and eigenvectors"""
    vs = [y for (x,y) in sorted(zip(eig.real,vs),key=lambda x: x[0])]
    eig = [y for (x,y) in sorted(zip(eig.real,eig),key=lambda x: x[0])]
    return np.array(eig),vs



def ishermitian(m):
    """Check if a matrix is Hermitian"""
    d = m - np.conjugate(m.T)
    if np.max(np.abs(d))>1e-6: return False
    return True

def expm(m):
    """Compute exponential"""
    m = todense(m)
    return dlg.expm(m)


def inv(m):
    """Inverse"""
    m = todense(m)
    return dlg.inv(m)


#
#def lowest_eigenvectors(h,n=10,nmax=maxsize):
#  """Get a ground state"""
#  info = False
#  if h.shape[0]>nmax:
#    if info: print("Calling ARPACK")
#    eig,eigvec = slg.eigsh(h,k=n,which="SA",maxiter=100000)
##    eigvec = [v for (e,v) in zip(eig,eigvec.T)]
#  else:
#    if info: print("Full diagonalization")
#    eig,eigvec = dlg.eigh(h.todense())
#    eigvec = eigvec.T # transpose
##  print(sorted(eig))
##  eigevec = [v for (e,v) in sorted(zip(eig,eigvec))]
#  return eigvec[0:n] # return eigenvectors
#


def expm(m):
    m = todense(m)
    return dlg.expm(m) # exponential matrix

#def expm(m):
#    m = todense(m)
#    es,vs = dlg.eig(m)
#    d = np.zeros(m.shape,dtype=np.complex)
#    for i in range(len(es)): d[i,i] = np.exp(es[i])
#    R = vs.T
#    Rh = np.conjugate(R.T)
#    U = Rh@d@R
#    return U



def ismatrix(m):
    return type(m)==np.ndarray or issparse(m) or type(m)==np.matrix


