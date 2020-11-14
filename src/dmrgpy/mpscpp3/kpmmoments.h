// compute the KPM moments for matrix m and vectors vi and vj
static auto moments_vi_vj_full=[](auto m, auto vi, auto vj, int n) {
  // technique use to apply the mpo
//  auto fitmpo = get_bool("fitmpo_kpm") ;
//  fitmpo = false ; // this does not work ok
  ofstream myfile; // file for the moments
  ofstream entropyfile; // file for the entropies
  myfile.open("KPM_MOMENTS.OUT"); // open file
  entropyfile.open("KPM_ENTROPY.OUT"); // open file
  int kpmmaxm = get_int_value("kpmmaxm") ; // bond dimension for KPM
  auto kpmcutoff = get_float_value("kpm_cutoff") ; // bond dimension for KPM
  auto v = vi*1.0 ; // initialize
  auto am = vi*1.0 ; // initialize
  auto a = vi*0.0 ; // initialize
  a = applyMPO(m,v,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff}) ;
  auto ap = a*1.0 ; // initialize
  auto bk = innerC(vj,v) ; // overlap
  auto bk1 = innerC(vj,a) ; // overlap
  int cindex = length(v)/2 ; // central site
  myfile << std::setprecision(20) << real(bk) << "  "
                       << std::setprecision(20)<< imag(bk) << endl;
  myfile << std::setprecision(20) << real(bk1) << "  "
                       << std::setprecision(20)<< imag(bk1) << endl;
  entropyfile << entropy(v,cindex) << endl ;
  entropyfile << entropy(a,cindex) << endl ;
  int i ;
  for(i=0;i<n;i++) {
    ap = applyMPO(m,a,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff}) ;
    ap = sum(2.0*ap,-1.0*am,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff}) ; // recursion relation
    bk = innerC(vj,ap) ; // compute term 
    myfile << std::setprecision(20) << real(bk) << "  "
                       << std::setprecision(20)<< imag(bk) << endl;
    entropyfile << entropy(ap,cindex) << endl ;
    am = a*1.0; // next iteration
    a = ap*1.0; // next iteration
  } ;
  entropyfile.close();
  myfile.close();
  return 0 ;
} ;




// this is a modified technique that can be used when the
// two KPM vectors are the same
static auto moments_vi_accelerated=[](auto m, auto vi, int n) {
  // technique use to apply the mpo
  ofstream myfile; // file for the moments
  myfile.open("KPM_MOMENTS.OUT"); // open file
  int kpmmaxm = get_int_value("kpmmaxm") ; // bond dimension for KPM
  auto kpmcutoff = get_float_value("kpm_cutoff") ; // bond dimension for KPM
  auto v = vi*1.0 ; // initialize
  auto am = vi*1.0 ; // initialize
  auto a = vi*0.0 ; // initialize
  a = applyMPO(m,vi,Args("MaxDim",kpmmaxm,"Cutoff",kpmcutoff)) ;
  auto ap = a*1.0 ; // initialize
  auto bk = innerC(vi,vi) ; // overlap
  auto bk1 = innerC(vi,a) ; // overlap
  auto mu0 = bk ; // save the zeroth
  auto mu1 = bk1 ; // save the first
  myfile << std::setprecision(20) << real(bk) << "  "
                       << std::setprecision(20)<< imag(bk) << endl;
  myfile << std::setprecision(20) << real(bk1) << "  "
                       << std::setprecision(20)<< imag(bk1) << endl;
  int i ;
  for(i=0;i<n/2;i++) {
    ap = applyMPO(m,a,Args("MaxDim",kpmmaxm,"Cutoff",kpmcutoff)) ;
    ap = sum(2.0*ap,-1.0*am,Args("MaxDim",kpmmaxm,"Cutoff",kpmcutoff)) ; // recursion relation
    bk = innerC(a,a) ; // compute overlap term 
    bk1 = innerC(a,ap) ; // compute overlap term 
    bk = 2*bk - mu0; // correction due to the trick
    bk1 = 2*bk1 - mu1; // correction due to the trick
    myfile << std::setprecision(20) << real(bk) << "  "
                       << std::setprecision(20)<< imag(bk) << endl;
    myfile << std::setprecision(20) << real(bk1) << "  "
                       << std::setprecision(20)<< imag(bk1) << endl;
    am = a*1.0; // next iteration
    a = ap*1.0; // next iteration
  } ;
  myfile.close();
  return 0 ;
} ;



// select which method to use
static auto moments_vi_vj=[](auto m, auto vi, auto vj, int n) {
	// check if the default should be used
	auto kpmaccelerate = get_bool("kpm_accelerate");
	if (kpmaccelerate) {
		// if the two vectors are the same, use special function
		bool aresame = same_mps(vi,vj);
		if (same_mps(vi,vj)) {
		  cout << "KPM Accelerated mode" << endl;
	          moments_vi_accelerated(m,vi,n) ;
		return 0; // return
		}
	};
	// use the conventional method if the conditions were not met
	moments_vi_vj_full(m,vi,vj,n) ;
	return 0;
};




