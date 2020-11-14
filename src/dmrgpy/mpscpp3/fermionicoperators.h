

// return a spinless creation/annhilation operator
// with Jordan-Wigner strings
static auto fermionic_operator_spinless= [](auto sites, int i, auto name) {
        auto ampo = AutoMPO(sites);  // create MPO
        if (compare_string(name,"C" )) {
        ampo += 1.0,"A",i+1; // bosonic one
        }
        if (compare_string(name,"Cdag" )) {
        ampo += 1.0,"Adag",i+1; // bosonic one
        }
        int maxm = get_int_value("maxm") ;
        auto cutoff = get_float_value("cutoff") ;
        auto m0 = toMPO(ampo) ; // create MPO
        auto m = toMPO(ampo) ; // create MPO
        for(int j=0;j<i;j++) {
          auto ampoi = AutoMPO(sites); // temporal one
          ampoi += 1.0,"F",j+1; // string operator
          auto mi = toMPO(ampoi) ; // create MPO
          m = nmultMPO(prime(mi),m0,{"MaxDim",maxm,"Cutoff",cutoff}) ; // multiply MPO
          m0 = m ; // reassign for next iteration
        }
        return m; // return MPO
};


