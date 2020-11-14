// function to get a gap

static auto get_gap = [](auto H, auto sites, auto sweeps) {
  auto psi0 = randomMPS(sites);
  auto [en0,psi3] = dmrg(H,psi0,sweeps,{"Quiet=",true});
  auto wfs = std::vector<MPS>(1);
  wfs.at(0) = psi0;
  auto psi1 = randomMPS(sites);
  auto [en1,psi2] = dmrg(H,wfs,psi1,sweeps,{"Quiet=",true,"Weight=",20.0});
  ofstream myfile; // create object
  myfile.open("GAP.OUT"); // open file
  myfile << en1-en0 << endl; // write file
  return en1-en0 ;
};
