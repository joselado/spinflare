// scale the Hamiltonian so it lies between -1 and 1
static auto bandwidth=[](auto sites, auto H) {
    auto sweeps = get_sweeps(); // get the sweeps
    // get minimum energy
    auto wf0 = randomMPS(sites) ;
    auto [emin,wfmin] = dmrg(H,wf0,sweeps,{"Quiet=",true}); 
    // get maximum energy
    auto wf1 = randomMPS(sites) ;
    auto [emax,wfmax] = dmrg(-1*H,wf1,sweeps,{"Quiet=",true}); 
    return -emax-emin ; // return the bandwidth
}
;


