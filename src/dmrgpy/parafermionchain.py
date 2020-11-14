from .manybodychain import Many_Body_Chain
import numpy as np


class Parafermionic_Chain(Many_Body_Chain):
    """Class for a parafermionic chain"""
    def __init__(self,n,Z=3):
        self.Z = Z # type of parafermion
        self.N = [self.get_operator("N",i) for i in range(n)]
        self.Sig = [self.get_operator("Sig",i) for i in range(n)]
        self.Sigd = [self.get_operator("SigDag",i) for i in range(n)]
        self.Tau = [self.get_operator("Tau",i) for i in range(n)]
        self.Taud = [self.get_operator("TauDag",i) for i in range(n)]
        self.Id = self.get_operator("Id",1)
        if Z==3: Many_Body_Chain.__init__(self,[-2 for i in range(n)])
        elif Z==4: Many_Body_Chain.__init__(self,[-3 for i in range(n)])
        elif Z==2: 
             print("Not with DMRG!")
             Many_Body_Chain.__init__(self,[-1 for i in range(n)])
        else: raise
        self.use_ampo_hamiltonian = True # use ampo
        self.Chi = []
        self.Dis = []
        self.Psi = []
        for i in range(n): 
            t = 1
            for j in range(i): t = t*self.Tau[j]
            self.Chi.append(t*self.Sig[i])
            self.Psi.append(t*self.Sig[i]*self.Tau[i])
        for i in range(n): 
            t = 1
            for j in range(i+1): t = t*self.Tau[j]
            self.Dis.append(t)
        self.Psid = [o.get_dagger() for o in self.Psi]
        self.Chid = [o.get_dagger() for o in self.Chi]
    def get_ED_obj(self):
        from .pyparafermion import parafermion
        obj = parafermion.Parafermion_Chain(self)
        return obj
    def get_dynamical_correlator(self,mode="DMRG",**kwargs):
        if mode=="DMRG":
            return super().get_dynamical_correlator_MB(**kwargs)
        elif mode=="ED":
            return self.get_ED_obj().get_dynamical_correlator(**kwargs)
    def test(self,**kwargs):
        return test_commutation(self,**kwargs)








def test_commutation(self):
    """Perform a test of the commutation relations"""
    Chi = self.Chi
    Chid = self.Chid
    Tau = self.Tau 
    Taud = self.Taud
    Psi = self.Psi
    Psid = self.Psid
    n = len(Chi) # number of sites
    ntries = 3 # number of tries
    omega = np.exp(1j*2*np.pi/self.Z)
    for ii in range(ntries):
        i = np.random.randint(n)
        j = np.random.randint(n)
        if i>=j: continue
        d = Chi[i]*Chi[j] - omega*Chi[j]*Chi[i]
        if not self.is_zero_operator(d):
            print("Chi,Chi failed",i,j)
            raise
        d = Psi[i]*Psi[j] - omega*Psi[j]*Psi[i]
        if not self.is_zero_operator(d):
            print("Psi,Psi failed",i,j)
            raise
        d = Chi[i]*Psi[j] - omega*Psi[j]*Chi[i]
        if not self.is_zero_operator(d):
            print("Psi,Chi failed",i,j)
            raise
    print("Commutation test passed")

