// perform a time evolution

static auto quench=[](auto sites) {
  auto name1 = get_str("tevol_operator_i");  // first operator
  auto name2 = get_str("tevol_operator_j");  // second operator
  auto i1 = get_int_value("tevol_site_i");  // second operator
  auto i2 = get_int_value("tevol_site_j");  // second operator
  // now get the operators
  auto A1 = get_operator(sites,i1,name1);
  auto A2 = get_operator(sites,i2,name2);
  auto H = get_hamiltonian(sites) ; // get the ampo for the Hamiltonian
  auto psi = get_gs(sites,H) ; // get the ground state
  int maxm = get_int_value("maxm") ; // bond dimension
  auto cutoff = get_float_value("cutoff") ; // cutoff
  // apply the first operator
  // now do the time evolution
  auto nt = get_int_value("tevol_nt"); // number of time steps
  auto dt = get_float_value("tevol_dt"); // delta of time
  int it;
  auto args = Args("Cutoff",cutoff,"MaxDim",maxm);
  auto args2 = Args("Cutoff",cutoff,"MaxDim",maxm,"Method","Fit");
  // compute ground state energy
  auto EGS = overlap(psi,H,psi)/overlap(psi,psi); // ground state enrgy
  // get the AutoMPO
  auto ampo = get_ampo(sites) ; // get the ampo for the Hamiltonian
  // shift by the ground state energy
  ampo += -EGS,"Id", 1; // minus ground state energy
  auto expH = toExpH(ampo,dt*Cplx_i); // get the exponential of the H
  ofstream fileevol; // file for the evolution
  fileevol.open("TIME_EVOLUTION.OUT"); // time evolution
  auto psi1 = applyMPO(A1,psi,args) ;
  auto psi2 = applyMPO(A2,psi,args) ;
//  normalize(psi1); // normalize
//  normalize(psi2); // normalize
  auto norm0 = sqrt(innerC(psi1,psi1)) ;
  auto fittd = get_bool("tevol_fit_td") ; // use fitting method
  for (it=0;it<nt;it++) { // loop
	      if (fittd) applyMPO(expH,psi1,args2) ; // evolve
	      if (not fittd) psi1 = applyMPO(expH,psi1,args); // evolve
              normalize(psi1); // normalize
	      psi1 *= norm0 ; // restore initial norm
	      auto z = innerC(psi2,psi1) ; // overlap
//	      auto z = innerC(psi,psi1) ; // overlap
	      // write in a file
	      fileevol << std::setprecision(8) << real(z) << "  "
                       << std::setprecision(8)<< imag(z) << endl;
  } ;
  fileevol.close(); // close file
};

