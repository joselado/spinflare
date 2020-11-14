#include"fermionicoperators.h" // library with fermionic operators



static auto Iden=[](auto sites) {
    auto ampo = AutoMPO(sites);
    ampo += "Id", 1;
    return toMPO(ampo);
};




static auto get_operator= [](auto sites, int i, auto name) {
	if (compare_string(name,"Id" )) return Iden(sites) ; // return identity
	auto ampo = AutoMPO(sites);
	if (site_type(i)>1)  // spin site
		ampo += 1.0,name,i+1 ;
	else if (site_type(i)==1) { // spinful fermionic site
        	if (compare_string(name,"Sx" )) {
        		ampo += 0.5,"Cdagdn",i+1,"Cup",i+1;
        		ampo += 0.5,"Cdagup",i+1,"Cdn",i+1;
        	} ;
		if (compare_string(name,"Sy" )) {
        		ampo += -0.5*1i,"Cdagdn",i+1,"Cup",i+1;
        		ampo += 0.5*1i,"Cdagup",i+1,"Cdn",i+1;
        	} ;
		if (compare_string(name,"Sz" )) {
        		ampo += 0.5,"Cdagup",i+1,"Cup",i+1;
        		ampo += -0.5,"Cdagdn",i+1,"Cdn",i+1;
        	} ;
		if (compare_string(name,"Cdag" )) {
        		ampo += 1.0,"Cdagup",i+1;
        		ampo += 1.0,"Cdagdn",i+1;
        	} ;
		if (compare_string(name,"C" )) {
        		ampo += 1.0,"Cup",i+1;
        		ampo += 1.0,"Cdn",i+1;
        	} ;
		if (compare_string(name,"Cup" )) {
        		ampo += 1.0,"Cup",i+1;
        	} ;
		if (compare_string(name,"Cdn" )) {
        		ampo += 1.0,"Cdn",i+1;
        	} ;
		if (compare_string(name,"Cdagup" )) {
        		ampo += 1.0,"Cdagup",i+1;
        	} ;
		if (compare_string(name,"Cdagdn" )) {
        		ampo += 1.0,"Cdagdn",i+1;
        	} ;
		// superconducting terms
		if (compare_string(name,"delta" )) {
        		ampo += 1.0,"Cdn",i+1,"Cup",i+1;
        	} ;
		if (compare_string(name,"deltad" )) {
        		ampo += 1.0,"Cdagdn",i+1,"Cdagup",i+1;
        	} ;
		// density terms
		if (compare_string(name,"density" )) {
        		ampo += 1.0,"Cdagdn",i+1,"Cdn",i+1;
        		ampo += 1.0,"Cdagup",i+1,"Cup",i+1;
        	} ;
		if (compare_string(name,"density_up" )) {
        		ampo += 1.0,"Cdagup",i+1,"Cup",i+1;
        	} ;
		if (compare_string(name,"density_dn" )) {
        		ampo += 1.0,"Cdagdn",i+1,"Cdn",i+1;
        	} ;
	}
	else if (site_type(i)==0) { // spinless fermionic site
		cout << "Requesting site " << i << endl;
		if ((compare_string(name,"C" )) 
				or (compare_string(name,"Cdag" ))) {
			return fermionic_operator_spinless(sites,i,name) ;
        	} ;
		// density terms
		if (compare_string(name,"density" )) {
        		ampo += 1.0,"Cdag",i+1,"C",i+1;
        	} ;
        };
        auto m = toMPO(ampo) ;	
return m ;
}
;


static auto get_occupation_operator= [](auto sites, int i) {
	auto ampo = AutoMPO(sites);
	if (site_type(i)==1)  // fermionic site
        	ampo += 1.0,"Ntot",i+1;
	if (site_type(i)==0)  // fermionic site
        	ampo += 1.0,"N",i+1;
        auto m = toMPO(ampo) ;	
return m ;
}
;



static auto get_occupation2_operator= [](auto sites, int i) {
	auto ampo = AutoMPO(sites);
	if (site_type(i)==1)  // fermionic site
        	ampo += 1.0,"Ntot",i+1,"Ntot",i+1;
	if (site_type(i)==0)  // fermionic site
        	ampo += 1.0,"N",i+1,"N",i+1;
        auto m = toMPO(ampo) ;	
return m ;
}
;




static auto get_delta_operator= [](auto sites, int i) {
	auto ampo = AutoMPO(sites);
	if (site_type(i)==1) { // fermionic site
        	ampo += 1.0,"Cdn",i+1,"Cup",i+1;
        	ampo += 1.0,"Cdagup",i+1,"Cdagdn",i+1;}
	if (site_type(i)==0)  // fermionic site
        	ampo += 1.0,"C",i+1,"C",i+1;
        auto m = toMPO(ampo) ;	
return m ;
}
;

auto get_up_creation_operator= [](auto sites, int i) {
	auto ampo = AutoMPO(sites);
	if (site_type(i)==1)  // fermionic site
        	ampo += 1.0,"Cdagup",i+1;
        auto m = toMPO(ampo) ;	
return m ;
}
;


static auto get_hopping_operator= [](auto sites, int i, int j) {
	auto ampo = AutoMPO(sites);
	if (site_type(i)==1) { // spinful fermionic site
        	ampo += 1.0,"Cdagup",i+1,"Cup",j+1;
        	ampo += 1.0,"Cdagdn",i+1,"Cdn",j+1;}
	if (site_type(i)==0)  // spinless fermionic site
        	ampo += 1.0,"Cdag",i+1,"C",j+1;
        auto m = toMPO(ampo) ;	
return m ;
}
;

static auto get_sidotsj_operator= [](auto sites, int i, int j) {
	auto ampo = AutoMPO(sites);
	if (site_type(i)>1) { // spin site
        	ampo += 1.0,"Sx",i+1,"Sx",j+1;
        	ampo += 1.0,"Sy",i+1,"Sy",j+1;
        	ampo += 1.0,"Sz",i+1,"Sz",j+1;
	}
        auto m = toMPO(ampo) ;	
return m ;
}
;



static auto add_spin_operator= [](auto ampo, auto sites, float v, int i, auto name) {
	if (site_type(i)!=1)  // spin site
		ampo += v,name,i+1 ;
	else if (site_type(i)==1) { // fermionic site
        	if (compare_string(name,"Sx" )) {
        		ampo += v*0.5,"Cdagdn",i+1,"Cup",i+1;
        		ampo += v*0.5,"Cdagup",i+1,"Cdn",i+1;
        	} ;
        	if (compare_string(name,"Sy" )) {
        		ampo += -v*0.5*1i,"Cdagdn",i+1,"Cup",i+1;
        		ampo += v*0.5*1i,"Cdagup",i+1,"Cdn",i+1;
        	} ;
        	if (compare_string(name,"Sz" )) {
        		ampo += v*0.5,"Nup",i+1;
        		ampo += -v*0.5,"Ndn",i+1;
        	} ;
        }
	return ampo ;
}
;



#include"multioperator.h" // library to read a multioperator

