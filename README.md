# SpinFlare

## User interface to compute many body spin systems
![Alt text](images/interface.png?raw=true "Main window of SpinFlare")

## Capabilities
Local magnetization
Static correlators
Dynamical correlators

## How to install
Execute the script "install.sh", afterwards you just have to write in the
terminal
"spinflare"

The compilation requires a modern C++ compiler and lapack libraries

The library relies on several Python libraries, a simple way of getting
all the dependencies is by installing Anaconda Python.

## Dependencies
At least g++ version 6
### How to install all the C++ dependencies in Ubuntu
sudo add-apt-repository ppa:ubuntu-toolchain-r/test

sudo apt-get update

sudo apt-get install g++-6

sudo apt-get install liblapack-dev

## Example
The example below shows the calculations of the dynamical structure factor of a Heisenberg chain with S=1 and 48 sites. This computations shows The emergence of the bulk Haldane gap and the gapless fractionalizated edge modes.
![Alt text](images/dyncorr_Haldane.png?raw=true "Dynamical correlator of the Haldane Heisenberg model")

