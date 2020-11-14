class MBchain : public SiteSet
    {
    public:


    MBchain(Args const& args = Args::global());


    void
    read(std::istream& s);

    };


void inline MBchain::
read(std::istream& s)
    {
    int N = itensor::read<int>(s);
    if(N > 0)
        {
        auto store = SiteStore(N);
        ifstream sfile; // file to read
        sfile.open("sites.in"); // file with the coupling
	int N2 = 0;
        sfile >> N2; // read the number of sites from file
        auto sites = SiteStore(N); // get an empty list of sites
	int nm = 0; // initialize
        for(int j = 1; j <= N; ++j)
            {
            auto I = Index{}; // generate an identity for this index
            sfile >> nm ; // read the label of the site from sites.in
            I.read(s); // read the site from sites.sites
	    cout << "Reading site from sites.sites "<< nm << endl ;
	    cout << "Read from sites.sites "<< I << endl ;
            if (nm==0) store.set(j,FermionSite(I)); // use fermions
            else if(nm == (-1)) store.set(j,Z3Site(I)); // Z3
            else if(nm == 1) store.set(j,BosonSite(I)); // boson
            else if(nm == 2) store.set(j,SpinHalfSite(I));
	    else if(nm == 3) store.set(j,SpinOneSite(I));
            else if(nm == 4) store.set(j,SpinThreeHalfSite(I));
            else if(nm == 5) store.set(j,SpinTwoSite(I));
            else if(nm == 6) store.set(j,SpinFiveHalfSite(I));
            else Error(format("MBchain cannot read index of size %d",nm));
            }
	sfile.close() ;
        init(std::move(store));
        }
    }

//} //namespace itensor







inline MBchain::
MBchain(Args const& args)
    {
    ifstream sfile; // file to read
    sfile.open("sites.in"); // file with the coupling
    int N, nm;
    sfile >> N; // read the number of sites and number of projections
    auto sites = SiteStore(N); // get an empty list of sites
    for (int i=1;i<=N;i++)  {
      sfile >> nm ; // read the name of that site
      cout << "Reading  " << nm << endl; // record this
      if (nm==0) sites.set(i,FermionSite({"SiteNumber",i,"ConserveQNs",false,"ConserveNf",false})); // use spinless
      else if (nm==(-1)) sites.set(i,Z3Site({"SiteNumber",i,"ConserveQNs",false})); // use Z3
      else if (nm==1) sites.set(i,BosonSite({"SiteNumber",i,"ConserveQNs",false})); // use spinful
      else if (nm==2) sites.set(i,SpinHalfSite({"SiteNumber",i,"ConserveQNs",false})); // use spin=1/2
      else if (nm==3) sites.set(i,SpinOneSite({"SiteNumber",i,"ConserveQNs",false})); // use spin=1
      else if (nm==4) sites.set(i,SpinThreeHalfSite({"SiteNumber",i,"ConserveQNs",false})); // use spin=3/2
      else if (nm==5) sites.set(i,SpinTwoSite({"SiteNumber",i,"ConserveQNs",false})); // use spin=2
      else if (nm==6) sites.set(i,SpinFiveHalfSite({"SiteNumber",i,"ConserveQNs",false})); // use spin=5/2
      else Error(format("MBchain cannot read index of size "));
    } ;
    sfile.close(); // close file

    SiteSet::init(std::move(sites));
    }




auto generate_sites() { // function to generate the sites
    auto sites = MBchain() ; // read from file
    return sites ;
}



auto get_sites() { // function to get the sites
    auto sites = generate_sites() ;  // generate the sites
    // this was for testing
//    auto sites = BasicSiteSet<FermionSite>(6,{"ConserveQNs=",false,"ConserveNf",false}) ;
    // overwrite the sites
    if (check_task("gs_from_file")) readFromFile("sites.sites",sites);
    if (check_task("sites_from_file")) readFromFile("sites.sites",sites);
    cout << "Number of sites " << length(sites) << endl ;
    return sites ;
}



int site_type(int index) {
    ifstream sfile; // file to read
    sfile.open("sites.in"); // file with the sites
    int N, nm, out=-1;
    sfile >> N; // read the number of sites and number of projections
    for (int i=1;i<=N;i++)  {
      sfile >> nm ; // read this spin
      if (i-1==index) out = nm ; }
    sfile.close() ;
    cout << index << "  " << out << endl ;
    return out ;
}





