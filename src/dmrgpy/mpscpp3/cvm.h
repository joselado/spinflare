

MPS conjMPS(MPS psi){
    return dag(psi);
}


// CVM solver
MPS bicstab(MPO A, MPS b, double tol, int max_it, Args const& args){

    MPS x = b;
    MPS r_old = sum(b, -1 * applyMPO(A, x, args));
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

        Ap = applyMPO(A, p, args);
        alpha = innerC(prime(conjMPS(r_old)), r_) / innerC(prime(conjMPS(Ap)), r_);
        s = sum(r_old, -alpha * Ap, args);
        As = applyMPO(A, s, args);
        w = innerC(prime(conjMPS(As)), s) / innerC(prime(conjMPS(As)), As);
        x = sum(x, sum(alpha * p, w * s, args), args);
        r_new = sum(s, -w * As, args);
        res = sqrt(abs(innerC(prime(conjMPS(r_new)), r_new).real()));

        if(res <= tol){
            std::cout << "Residue = " << res << std::endl;
            break;
        }

        beta = (alpha / w) * innerC(prime(conjMPS(r_new)), r_) / innerC(prime(conjMPS(r_old)), r_);
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

    auto args = Args({"MaxDim", maxm, "Cutoff", cut});
    const std::complex<double> z(omega + energy, eta);
    auto A = sum(z * Iden(sites), -1. * H, args);
    auto b =  applyMPO(S2, psi, args);
    auto x = bicstab(A, b, tol, max_it, args);
    std::complex<double> G = innerC(prime(psi), S1, x);

    return -G.imag() / M_PI;
};
