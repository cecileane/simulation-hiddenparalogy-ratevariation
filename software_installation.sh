# Document software installtion in Franklin00 

# Download simphy 
cd ~/private/software/ 
wget https://github.com/adamallo/SimPhy/releases/download/v1.0.2/SimPhy_1.0.2.tar.gz 
# simphy 1.0.2 is the most updated version by Feb 5 2025 
tar -xvf SimPhy_1.0.2.tar.gz
chmod +x SimPhy_1.0.2/bin/simphy_lnx64
cd ~/PATH/TO/simulation-reptiles # go back to simulation-reptiles folder 
ln -s ~/private/software/SimPhy_1.0.2/bin/simphy_lnx64 executables/simphy 
./executables/simphy -h # test 

# Download Seq-Gen 
cd ~/private/software 
wget https://github.com/rambaut/Seq-Gen/archive/refs/tags/v1.3.5.tar.gz 
# Seq-Gen v.1.3.5 is the most updated version by Feb 6 2025 
tar -xvf v1.3.5.tar.gz  
cd Seq-Gen-1.3.5/source 
make 
cd ~/PATH/TO/simulation-reptiles 
ln -s ~/private/software/Seq-Gen-1.3.5/source/seq-gen executables/seq-gen 
executables/seq-gen -h 

# Download IQtree
cd ~/private/software/
wget https://github.com/iqtree/iqtree2/archive/refs/tags/v2.3.6.tar.gz 
tar -xvf v2.3.6.tar.gz 
mkdir build && cd build 


