# Differential di-jet imbalance paper 
### Nick Elsey

## Compiling

The project is built using cmake - in-source builds are not supported. Building should be as easy as 
```
cd /path/to/source
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/desired/install/location -DBUILD_32BIT=on -DBUILD_TEST=ON ..
make install
```
The configure and build will take a decent amount of time; the project will fetch and build all of its dependencies, including fastjet 3.3, which is a requirement but not present on RCF. Once the project is installed into the desired location, there is one final step.
```
cd /path/to/install
setenv LD_LIBRARY_PATH `pwd`/lib:${LD_LIBRARY_PATH} (if using csh)
export LD_LIBRARY_PATH=`pwd`/lib:${LD_LIBRARY_PATH} (if using bash)
```

