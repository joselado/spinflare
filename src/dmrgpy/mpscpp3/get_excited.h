// function to get several excited states

float get_energy_fluctuation(auto psi1, auto H) {
    return 0.0; // this has to be fixes
 //   psi1.orthogonalize(); // normalize the wavefunction
    auto psi2 = applyMPO(H,psi1,{"MaxDim",get_int_value("maxm"),
                 "Cutoff",get_float_value("cutoff")}) ;
    float de;
    de = overlap(psi1,psi2);
    de = overlap(psi2,psi2) - de*de; // energy fluctuation
    return de; // return energy fluctuation
};



static auto get_excited=[](auto H, auto sites, auto sweeps, int nexcited) {
  ofstream myfile;
  myfile.open("EXCITED.OUT"); // open file
  auto psi0 = randomMPS(sites); // first wavefunction
  auto [en0,psi] = dmrg(H,psi0,sweeps,{"Quiet=",true});
  auto de = get_energy_fluctuation(psi,H);// fluctuation
  myfile << std::setprecision(20) << en0 << "  " << de << endl;
  int i; 
  auto wfs = std::vector<MPS>(nexcited);
  for (i=0;i<nexcited;i++) wfs.at(i) = psi; // initialize 
  // lagrange multiplier
  float weight = bandwidth(sites,H)*get_float_value("scale_lagrange_excited"); 
  int numw=1; // number of wavefunctions found
  auto psi1 = randomMPS(sites) ; // new random wavefunction
  for (i=1;i<nexcited;i++)  { 
    // now compute a new excited state
    // new energy
    auto [en2,psi2] = dmrg(H,wfs,psi1,sweeps,{"Quiet=",true,"Weight=",weight}); 
    wfs.at(i) = psi2 ; // store this wavefunction
    de = get_energy_fluctuation(psi2,H);
    if (de<1e-2) { // if the fluctuation is small enough
      wfs.at(numw) = psi2 ; // store this wavefunction
      psi1 = randomMPS(sites) ; // new random wavefunction
      myfile << std::setprecision(8) << en2 << "  " << de << endl; // write
      numw += 1; // increase counter
    };

  } ;
  auto wfsout = std::vector<MPS>(numw); // wavefunctions found
  for (i=0;i<numw;i++)  wfsout.at(i) = wfs.at(i); // store
  return wfsout ;
}
;

