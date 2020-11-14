
// Get the ground state of this Hamiltonian, and write
// the wavefunction into a file

static auto get_gs=[](auto sites, auto H) {
    // read the GS from a file
    auto psi0 = randomMPS(sites);
    if (get_bool("gs_from_file"))  {
	    psi0 = read_wf(get_str("starting_file_gs")) ;
//  check if return this wavefunction 
           if (get_bool("skip_dmrg_gs")) return psi0 ;
    };
    auto sweeps = get_sweeps(); // get the DMRG sweeps
    auto [energy,psi] = dmrg(H,psi0,sweeps); // ground state energy
    ofstream myfile; // create object
    writeToFile("psi_GS.mps",psi); // write the GS wavefunction
    writeToFile("sites.sites",sites); // write the sites
    myfile.open("GS_ENERGY.OUT"); // open file
    myfile << std::setprecision(20) << energy << endl; // write file
    psi = psi.normalize(); // normalize wavefunction
//  } ;
  return psi ; // return the ground state
}
;
