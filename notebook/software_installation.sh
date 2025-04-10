# Document software installtion in the stats servers at UW-Madison 

#-------------- Install julia --------------# 
curl -fsSL https://install.julialang.org | sh # choose to proceed
. /u/b/i/bingl/.bashrc # Depending on the path 
. /u/b/i/bingl/.bash_profile 
# Install packages used in Julia 
julia
import Pkg
Pkg.add("PhyloPlots")
Pkg.add("PhyloNetworks")
Pkg.add("ArgParse") 
Pkg.add("StatsBase") 
Pkg.add("StableRNGs")

#-------------- Download simphy --------------# 
cd ~/private/software/ 
wget https://github.com/adamallo/SimPhy/releases/download/v1.0.2/SimPhy_1.0.2.tar.gz 
# simphy 1.0.2 is the most updated version by Feb 5 2025 
tar -xvf SimPhy_1.0.2.tar.gz
chmod +x SimPhy_1.0.2/bin/simphy_lnx64
cd ~/PATH/TO/simulation-reptiles # go back to simulation-reptiles folder 
ln -s ~/private/software/SimPhy_1.0.2/bin/simphy_lnx64 executables/simphy 
./executables/simphy -h # test 

#-------------- Download Seq-Gen--------------# 
cd ~/private/software 
wget https://github.com/rambaut/Seq-Gen/archive/refs/tags/v1.3.5.tar.gz 
# Seq-Gen v.1.3.5 is the most updated version by Feb 6 2025 
tar -xvf v1.3.5.tar.gz  
cd Seq-Gen-1.3.5/source 
make 
cd ~/PATH/TO/simulation-reptiles 
ln -s ~/private/software/Seq-Gen-1.3.5/source/seq-gen executables/seq-gen 
executables/seq-gen -h 

#-------------- Download IQtree --------------# 
cd ~/private/software/
wget https://github.com/iqtree/iqtree2/releases/download/v2.4.0/iqtree-2.4.0-Linux-intel.tar.gz 
tar -xvf iqtree-2.4.0-Linux-intel.tar.gz # Latest version on Feb 7 2025 
# Interestingly, iqtree 2.4.0 was released on Feb 7 2025. Based on the new manual, the executable is located in bin/ and could be directly used. 
cd iqtree-2.4.0-Linux-intel 
bin/iqtree2 -s example.phy # run an example to test installation --> work okay 
cd ~/PATH/TO/simulation-reptiles 
ln -s ~/private/software/iqtree-2.4.0-Linux-intel/bin/iqtree2 executables/iqtree2 
./executables/iqtree2 -h # test 

#-------------- Download ASTER/Astral --------------#
cd ~/private/software/
wget https://github.com/chaoszhang/ASTER/archive/refs/heads/Linux.zip # Version: v1.20.3.6 -- Newest version on Feb 14 2025 
unzip Linux.zip  
cd ASTER-Linux/ && make 
cd ~/PATH/TO/simulation-reptiles  
ln -s ~/private/software/ASTER-Linux/bin/astral-pro3 executables/astral-pro3 

#-------------- Make sure to have Python 3 --------------# 
# Here I installed python 3.12 
cd ~/private/software/
wget https://www.python.org/ftp/python/3.12.0/Python-3.12.0.tgz 
cd Python-3.12 
./configure --prefix=~/private/software/python3.12
make -j$(nproc) 
echo 'export PATH=/u/b/i/bingl/private/software/Python-3.12.0:$PATH' >> ~/.bashrc # where the python binary is. 
source ~/.bashrc
python --version # double check if it is python 3.12. 

# Install python dependencies: 
python3.12 -m ensurepip --upgrade
python3.12 -m pip install --upgrade pip
python3.12 -m pip install biopython
python -c  "import Bio.SeqIO, os, sys, collections" # double check if modules are installed successfully 
 

#--------------Download snp-sites --------------# 
# snp-sites calls SNPs from molecular sequences from Seq-Gen 
wget https://github.com/sanger-pathogens/snp-sites/archive/refs/tags/v2.5.1.tar.gz # version 2.5.1 newest version on March 7 2025 
tar -xvf v2.5.1.tar.gz  
cd snp-sites-2.5.1/
autoreconf -i -f
./configure
make
cd PATH/TO/MY/REPO/simulation-reptiles/
ln -s ~/private/software/snp-sites-2.5.1/src/snp-sites executables/snp-sites 
executables/snp-sites --help # Check if installed successfully  

#--------------Install required R packages-----------# 
R # open the interactivate mode 
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
  devtools::install_github("uqrmaie1/admixtools")
} 
install.packages("dplyr") 
install.packages("optparse") 













