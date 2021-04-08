# SpinFlare

## User interface to compute many body spin systems
![Alt text](images/interface.png?raw=true "Main window of SpinFlare")

## Capabilities
Local magnetization

Static correlation functions

Dynamical correlator function

Excited states

## How to install

## Linux
Execute the script "install.sh", afterwards you just have to write in the
terminal
"spinflare"

The compilation requires a modern C++ compiler and lapack libraries

The library relies on several Python libraries, a simple way of getting
all the dependencies is by installing Anaconda Python.


## Windows ##

For using this program in Windows, the easiest solution is to create a virtual
machine using [Virtual Box](https://www.virtualbox.org/), installing
a version of [Ubuntu](https://releases.ubuntu.com/20.04/)
in that virtual machine, and following the previous
instructions.


## Dependencies
At least g++ version 6
### How to install all the C++ dependencies in Ubuntu
sudo add-apt-repository ppa:ubuntu-toolchain-r/test

sudo apt-get update

sudo apt-get install g++-6

sudo apt-get install liblapack-dev

## Example
### Dynamical correlator of a topological spin chain (S=1)
The example below shows the calculations of the dynamical structure factor of a Heisenberg chain with S=1 and 48 sites. This computation shows the emergence of the bulk Haldane gap and the gapless fractionalizated edge modes.
![Alt text](images/dyncorr_Haldane.png?raw=true "Dynamical correlator of the Haldane Heisenberg model")

### Local magnetization of an S=1/2 quantum chain with an edge field
The example below shows the local magnetization for a Heisenberg chain of S=1/2, with a local magnetic field in the x direction applied on the first site. Due to the quantum disordered ground state, it is observed a decay of the magnetization as we go further from the edge site. The decay is a power-law due to the gapless nature of the system.
![Alt text](images/magnetization_s12.png?raw=true "Local magnetization of a S=1/2 Heisenberg chain with a local field")

### Non-local correlator of an S=1 chain (with S=1/2 on the edge)
The example below shows the non-local correlator (in log scale) between the edge and the different sites for a Heisenberg chain of S=1 (with S=1/2 in the edge to lift the topological modes). Due to the many-body gap of the system, the correlator decays exponentially with distance as seen in the plot.
![Alt text](images/correlator_s1.png?raw=true "Non-local correlator in log scale")

