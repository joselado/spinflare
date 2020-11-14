



// calculate all the correlators
int get_correlator()   {
  auto sites = get_sites();
  auto H = get_hamiltonian(sites) ;
  auto psi = get_gs(sites,H) ;
  ifstream cfile; // declare
  ofstream ofile; // declare
  cfile.open("correlators.in");  // open file
  ofile.open("CORRELATORS.OUT");  // open file
  int nc; 
  cfile >> nc; // number of correlators
  float c; // declare float
  int i,j ;
  int maxm = get_int_value("maxm") ; // bond dimension for KPM
  float cutoff = get_float_value("cutoff") ; // cutoff for KPM
  for (int ic=0;ic<nc;++ic) { // loop over correlators
    cfile >> i >> j; // index of correlators
    // get the two operators
    auto opi = get_operator(sites,i,get_str("correlator_operator_i")) ;
    auto opj = get_operator(sites,j,get_str("correlator_operator_j")) ;
    auto psii = applyMPO(opi,psi,{"MaxDim",maxm,"Cutoff",cutoff}) ;
    // in case the Hamiltonian is applied in the middle
    if (get_bool("correlator_apply_hamiltonian")) {
    psii = applyMPO(H,psii,{"MaxDim",maxm,"Cutoff",cutoff});
    }
    auto psij = applyMPO(opj,psi,{"MaxDim",maxm,"Cutoff",cutoff}) ;
    auto c = innerC(prime(psii),psij);
    ofile << ic << "   "  ; // write index 
    ofile << std::setprecision(20) << real(c) << "  " << imag(c) << endl ;
  };
  ofile.close() ;
  cfile.close() ;
  return 0 ;
}
