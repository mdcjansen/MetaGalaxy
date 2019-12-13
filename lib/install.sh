#!/bin/bash
git clone https://github.com/rrwick/Filtlong
git clone https://github.com/fenderglass/Flye
git clone https://github.com/samtools/htslib
git clone https://github.com/samtools/samtools
git clone https://github.com/DerrickWood/kraken2
git clone https://github.com/marbl/Krona
git clone https://github.com/lh3/minimap2
git clone https://github.com/matsen/pplacer
git clone https://github.com/hyattpd/Prodigal
git clone https://github.com/ablab/quast
wget https://github.com/bbuchfink/diamond/releases/download/v0.9.29/diamond-linux64.tar.gz
wget https://bitbucket.org/berkeleylab/metabat/get/master.tar.gz
wget https://github.com/rrwick/Bandage/releases/download/v0.8.1/Bandage_Ubuntu_static_v0_8_1.zip
unzip Bandage_Ubuntu_static_v0_8_1.zip
tar xzf diamond-linux64.tar.gz
mv diamond* ./diamond
tar xzvf master.tar.gz
mv berkeleylab-metabat* ./metabat
mv quast ./Quast
cd Filtlong
make -j
cd ../Flye
make
cd ../htslib
autoheader
autoconf
./configure
make
make install --prefix=$(pwd)
cd ../samtools
autoheader
autoconf
./configure
make
make install --prefix=$(pwd)
cd ../kraken2
./install_kraken2.sh ./
cd ../Krona/KronaTools
./install.pl --prefix ./
./updateTaxonomy.sh
./undateAccessions.sh
cd ../../minimap2 && make
cd ../pplacer/scripts
python setup.py install
cd ../../Prodigal
make
cd ../Quast
./setup.py install_full
cd ../metabat
mkdir build ; cd build && cmake -DCMAKE_INSTALL_PREFIX=$(pwd) .. && make && make install && mv bin ../ && rm -rf ../build
