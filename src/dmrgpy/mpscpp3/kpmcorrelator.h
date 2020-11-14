// compute the KPM moments for matrix m and vectors vi and vj
// and a shift in the energy
static auto moments_vi_vj_shift=[](auto m, auto vi, auto vj, int n, auto shift) {
  ofstream myfile;
  myfile.open("KPM_MOMENTS.OUT"); // open file
  int kpmmaxm = get_int_value("kpmmaxm") ; // bond dimension for KPM
  auto v = vi*1.0 ; // initialize
  auto am = vi*1.0 ; // initialize
  auto a = applyMPO(m,v,{"MaxDim",kpmmaxm,"Cutoff",1E-7}) ; // initialize
  a = sum(a,shift*v,{"MaxDim",kpmmaxm,"Cutoff",1E-7}) ; // shift
  auto ap = a*1.0 ; // initialize
  auto bk = innerC(vj,v) ; // overlap
  auto bk1 = innerC(vj,a) ; // overlap
  myfile << std::setprecision(20) << real(bk) << endl;
  myfile << std::setprecision(20) << real(bk1) << endl;
  int i ;
  for(i=0;i<n;i++) {
    ap = applyMPO(m,a,{"MaxDim",kpmmaxm,"Cutoff",1E-7}) ; // apply
    ap = 2.0*sum(ap,shift*a,{"MaxDim",kpmmaxm,"Cutoff",1E-7}) ; // shift
    ap = sum(ap,-1.0*am,{"MaxDim",kpmmaxm,"Cutoff",1E-7}) ; // recursion relation
    bk = innerC(vj,ap) ; // compute term 
    myfile << std::setprecision(20) << real(bk) << endl;
    am = a*1.0; // next iteration
    a = ap*1.0; // next iteration
  } ;
  return 0 ;
} ;






static auto get_moments_dynamical_correlator=[](auto sites, auto H)
{
  auto n = get_int_value("nkpm");
  auto delta = get_float_value("kpm_delta");
  // define on the fly the number of polynomials
  n = int(round(bandwidth(sites,H)/delta))*get_int_value("kpm_n_scale");
  ofstream myfile;
  myfile.open("KPM_NUM_POLYNOMIALS.OUT");
  myfile << std::setprecision(8) << n << endl;
  myfile.close(); // close file
  auto psi = get_gs(sites,H) ; // get the ground state
  auto m = scale_hamiltonian(sites,H) ; // scale this Hamiltonian
  auto m1 = Iden(sites) ; // identity
  auto m2 = Iden(sites) ; // identity
  // now read the operators to use
  // First the operator i
  if (get_bool("kpm_multioperator_i")) {
    m1 = get_multioperator("kpm_multioperator_i",sites);
  } ;
  if (not get_bool("kpm_multioperator_i")) {
    m1 = get_operator(sites,get_int_value("site_i_kpm"),
		    get_str("kpm_operator_i")); // first operator
  } ;
  // afterwards the operator j
  if (get_bool("kpm_multioperator_j")) {
    m2 = get_multioperator("kpm_multioperator_j",sites);
  } ;
  if (not get_bool("kpm_multioperator_j")) {
    m2 = get_operator(sites,get_int_value("site_j_kpm"),
		    get_str("kpm_operator_j")); // first operator
  } ;
  ///////////////////////////////////
  ////// once the operators have been read, continue
  ///////////////////////////////////
  int kpmmaxm = get_int_value("kpmmaxm") ; // bond dimension for KPM
  auto kpmcutoff = get_float_value("kpm_cutoff") ; // bond dimension for KPM
  auto psi1 = applyMPO(m1,psi,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff}) ;
  auto psi2 = applyMPO(m2,psi,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff}) ;
  moments_vi_vj(m,psi1,psi2,n) ;  //compute the KPM moments
  return 0 ;
} ;













