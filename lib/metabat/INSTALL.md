#Install instructions for supported Operating Systems

# Recent Linux Distributions with packages for MetaBAT pre-requisites:

## Docker:
-----------------
```
git clone https://bitbucket.org/berkeleylab/metabat.git
cd metabat
docker build --tag metabat .
docker run metabat runMetaBat.sh ...

```

### Prerequisites for Linux Ubuntu 16.04
----------------
```

# install boost and a build environment
sudo apt-get update 
sudo apt-get install -y build-essential libboost-all-dev git cmake curl libncurses5-dev zlib1g-dev

mkdir build ; cd build && cmake -DCMAKE_INSTALL_PREFIX=$HOME/metabat .. && make && make install
```

### Prerequisties for Linux Fedora 20
--------------------

```

#### install g++, boost and other build dependencies
sudo yum install gcc-c++ boost.x86_64 boost-devel.x86_64 zlib-devel.x86_64 libstdc++-static cmake

mkdir build ; cd build && cmake -DCMAKE_INSTALL_PREFIX=$HOME/metabat .. && make && make install
```

### Prerequisties for MacOS X (10.14 Mojave) : (  using Homebrew http://brew.sh/ )
---------------------

```

# First install Xcode from the App Store (version 10.2) 
# Second install Homebrew 
# Third install llvm with openmpi and boost and cmake
brew tap homebrew/versions
brew install llvm libomp boost cmake
brew link libomp

# use the latest llvm compiler and flags
export CPPFLAGS="-I/usr/local/opt/llvm/include"
export LDFLAGS="-L/usr/local/opt/llvm/lib"
export CC=/usr/local/opt/llvm/bin/clang
export CXX=/usr/local/opt/llvm/bin/clang++

mkdir build ; cd build && cmake -DCMAKE_INSTALL_PREFIX=$HOME/metabat .. && make && make install
```


### Older distributions must build and install:
```
gcc/g++ >= 4.9 or intel >= 18.0.1 or llvm >= 8.0
boost >= 1.53
cmake >= 3.8.2
make >= 4.1
```

# Build and install MetaBAT

```
git clone https://bitbucket.org/berkeleylab/metabat.git
cd metabat
mkdir build 
cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/metabat ..
make
make install
cd ..
rm -rf build

```


