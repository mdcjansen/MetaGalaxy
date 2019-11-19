#!/bin/bash

cd ~/MetaGalaxy/lib/Bandage
export QT_SELECT=5
qmake
make
cd ~/MetaGalaxy/lib/Prodigal
make install
cd ~/MetaGalaxy/lib/Filtlong
make -j
cd ~/MetaGalaxy/lib/Flye
make
cd ~/MetaGalaxy/lib/htslib
autoheader && autoconf && ./configure && make && make install
cd ~/MetaGalaxy/lib/kraken2
./isntall_kraken2.sh cd ~/MetaGalaxy/lib/kraken2
cd ~/MetaGalaxy/lib/Krona
./install.pl --prefix cd ~/MetaGalaxy/lib/Krona && sh updateTaxonomy.sh && updateAccenssions.sh
mkdir ~/MetaGalaxy/lib/metabat/build
cd ~/MetaGalaxy/lib/metabat/build
cmake -DCMAKE_INSTALL_PREFIX=$~/MetaGalaxy/lib/metabat .. && make && make install
cd ~/MetaGalaxy/lib/minimap2
make
cd ~/MetaGalaxy/lib/pplacer/scripts
python setup.py install
cd ~/MetaGalaxy/lib/Quast
python setup.py install_full
cd ~/MetaGalaxy/lib/samtools
autoheader && autoconf &&./configure && make && make install
