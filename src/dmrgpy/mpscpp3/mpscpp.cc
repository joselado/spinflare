#include "itensor/all.h"
#include "extra/all.h"
#include <iostream>
#include <fstream>
#include <iomanip>



using namespace itensor;
using namespace std;



#include"utilities.h" // read the different tasks
#include"check_task.h" // read the different tasks
#include"mpsalgebra.h" // functions to deal with MPS
#include"get_sweeps.h" // get the sweep info
#include"bandwidth.h"  // return the bandwidth of the hamiltonian
#include"get_gap.h" // compute the gap
#include"get_sites.h" // get the sites from a file
#include"get_ampo_operator.h" // get an arbitrary AMPO operator
#include"operators.h" // read the different tasks
#include"get_hamiltonian.h" // get the hoppings (in case there are)
#include"read_wf.h" // this does not work yet
#include"get_gs.h" // compute ground state energy and wavefunction
#include"get_excited.h" // compute excited states
#include"get_dos.h" // compute the DOS
#include"get_correlator.h" // compute correaltors betwee sites
#include"get_entropy.h" // compute entanglement entropy
#include"kpm.h" // KPM routines
#include"compute_overlap.h" // Compute overlap
#include"cvm_dynamical_correlator.h" // CVM dynamical correlator
#include"time_evolution.h" // Time evolution
#include"dynamical_correlator_excited.h" // dynamical correlator with exited
#include"vev.h" // Reduced density matrix


int 
main()
    {
    system("touch ERROR") ; // create error file


    // read the number of sites
    ifstream sfile; // file to read
    auto sites = get_sites(); // Get the different sites
    auto H = get_hamiltonian(sites) ; // get the Hamiltonian
//    test_hopping(H,sites); // test the hoppings
    auto sweeps = get_sweeps(); // get the DMRG sweeps

    if (check_task("GS")) {
      get_gs(sites,H) ; // get ground state wavefunction
    } ;
    if (check_task("correlator")) get_correlator() ; // write correlators 
    if (check_task("gap")) get_gap(H,sites,sweeps); // calculate the gap 
//    if (check_task("entropy")) get_entropy(psi,2); 
    if (check_task("excited")) {
      int nexcited = get_int_value("nexcited") ; // number of excited states
      auto wfs = get_excited(H,sites,sweeps,nexcited); // compute states 
    } ;
    if (check_task("dos")) get_dos(H,sites) ; // get the DOS
    if (check_task("dynamical_correlator"))  { // dynamical correlation
       get_moments_dynamical_correlator(sites,H);
    } ;
    if (check_task("dos"))  { // global DOS
       get_moments_dos(sites,H);
    } ;
    if (check_task("cvm"))  cvm_dynamical_correlator(sites) ; // CVM
// overlap task
    if (check_task("overlap"))  compute_overlap() ; // compute overlap
    if (check_task("time_evolution"))  quench(sites) ; // time evolution
//    if (check_task("density_matrix"))  reduced_dm() ; // DM
    if (check_task("vev"))  vev() ; // Vacuum expectation value
    if (check_task("dynamical_correlator_excited"))  
	    dynamical_correlator_excited(); // DM
    system("rm -f ERROR") ; // remove error file
    return 0;
    }
