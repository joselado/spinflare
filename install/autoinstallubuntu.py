import os

def install(cpp="g++-10"):
  print("This script will install all the depencies for Ubuntu")
  print("### It will ask for your sudo password ###")
  os.system("sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y")
  os.system("sudo apt-get update -y")
  os.system("sudo sudo apt-get install "+cpp+" -y")
  os.system("sudo apt-get install liblapack-dev -y")


if  __name__=="__main__":
    install() # execute installation
