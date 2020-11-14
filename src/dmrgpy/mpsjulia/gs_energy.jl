using ITensors
include("common.jl")

let
  sites = get_sites()
  ampo = read_operator("hamiltonian.in")
  H = MPO(ampo,sites)
  psi0 = randomMPS(sites)
  sweeps = Sweeps(get_int("nsweeps"))
  maxdim!(sweeps, get_int("maxm"))
  cutoff!(sweeps, get_float("cutoff"))
  @show sweeps
  energy, psi = dmrg(H,psi0, sweeps)
  println("Final energy = $energy")
  write_in_file("GS_ENERGY.OUT", energy)
end
