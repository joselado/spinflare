

MPS conjMPS(MPS psi){
    auto tmp = psi*1.;
    for(int i = 0; i < tmp.N(); ++i){
        tmp.Aref(i+1).conj();
    }
    return tmp;
}


// CVM solver
MPS bicstab(MPO A, MPS b, double tol, int max_it, Args const& args){

    MPS x = b;
    MPS r_old = sum(b, -1 * exactApplyMPO(A, x, args));
    MPS r_new;
    MPS r_ = r_old;
    MPS p = r_old;
    MPS s;
    MPS Ap;
    MPS As;
    std::complex<double> alpha;
    std::complex<double> beta;
    std::complex<double> w;
    double res;
    int k = 0;

    while(k < max_it){

        Ap = exactApplyMPO(A, p, args);
        alpha = overlapC(conjMPS(r_old), r_) / overlapC(conjMPS(Ap), r_);
        s = sum(r_old, -alpha * Ap, args);
        As = exactApplyMPO(A, s, args);
        w = overlapC(conjMPS(As), s) / overlapC(conjMPS(As), As);
        x = sum(x, sum(alpha * p, w * s, args), args);
        r_new = sum(s, -w * As, args);
        res = sqrt(abs(overlapC(conjMPS(r_new), r_new).real()));

        if(res <= tol){
            std::cout << "Residue = " << res << std::endl;
            break;
        }

        beta = (alpha / w) * overlapC(conjMPS(r_new), r_) / overlapC(conjMPS(r_old), r_);
        p = sum(r_new, beta * sum(p, -w * Ap, args), args);
        r_old = r_new;
        k++;

    };

    return x;
}

// main CVM function
static auto spectral_function=[](MPS psi, MPO H, MPO S1, MPO S2, double omega,
		double eta, double energy, double tol, int max_it,
		int maxm, double cut, auto sites) {

    auto args = Args({"Maxm", maxm, "Cutoff", cut});
    const std::complex<double> z(omega + energy, eta);
    auto A = sum(z * Iden(sites), -1. * H, args);
    auto b =  exactApplyMPO(S2, psi, args);
    auto x = bicstab(A, b, tol, max_it, args);
    std::complex<double> G = overlapC(psi, S1, x);

    return -G.imag() / M_PI;
};



// main CVM function
static auto apply_inverse=[]() {
    int maxm = get_int_value("maxm") ; // bond dimension
    int max_it = get_int_value("cvm_nit") ; // number of iterations
    auto cutoff = get_float_value("cutoff") ; // cutoff of DMRG
    auto tol = get_float_value("cvm_tol") ; // GS energy
    auto args = Args("Cutoff",cutoff,"Maxm",maxm); // MPS arguments
    auto A = get_mpo_operator("apply_inverse_multioperator.in");
    auto wf = read_wf("apply_inverse_wf0.mps") ; 
    auto x = bicstab(A, wf, tol, max_it, args);
    writeToFile("apply_inverse_wf1.mps",x);
    return 0;
};

