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
The configure and build will take a decent amount of time; the project will fetch and build all of its dependencies, including fastjet 3.3, which is a requirement but not present on RCF. Once the project is installed into the desired location, there is one final step to put the libdijetcore.so library into the linker search path
```
cd /path/to/install
setenv LD_LIBRARY_PATH `pwd`/lib:${LD_LIBRARY_PATH} (if using csh)
export LD_LIBRARY_PATH=`pwd`/lib:${LD_LIBRARY_PATH} (if using bash)
```

In the install directory, all the analyses are standalone binaries, in the bin directory. the lib directory holds libdijectore, a library of analysis tools. The submit directory holds submit scripts for submitting batches of jobs, and  resources holds needed information such as efficiency curves, dataset lists, bad tower lists, etc. 

## Data

The di-jet imbalance paper uses TStarJetPico trees, not standard MuDst or PicoDst trees. I have put the relevant data on RCF at 
```
/gpfs01/star/pwg/nelsey/papers/psn0729/data
```
There are three relevant directories - TStarJetPico_AuAuHT,  TStarJetPico_AuAuMB, TStarJetPico_pp, which contain the AuAu high-tower data used for the AuAu triggered sample, AuAu min-bias used for pp embedding and for off-axis background calculations, and the pp high-tower data for the pp triggered sample.

I also copied over fully processed trees, because the analysis can be quite lengthy, especially for the p+p and the systematic uncertainty estimations. This can be found at the same path under the dijetworker_output directory, which contains the two folders y6 and y7. y6 contains the embedded p+p di-jet trees, and y7 contains the y7 di-jet trees.

# Running the di-jet finder

The jet finding is done by the two binaries
```
bin/dijet_imbalance/differential_aj_auau_y7
bin/dijet_imbalance/differential_aj_pp_y6_y7_embed
```

for Au+Au and p+p, respectively. To see the command line options for each binary, use the command line flag `--helpshort`.

Below are examples of a command one would use to perform the analysis to reproduce the di-jet trees used in this analysis.

Au+Au
```
./bin/dijet_imbalance/differential_aj_auau_y7 --name=auau_y7_ --id=0 --input=/gpfs01/star/pwg/nelsey/papers/psn0729/data/TStarJetPico_AuAuHT/Clean809.root --outputDir= results/auau/grid_16_8 --leadJetPt=16 --subJetPt=8 --leadR=0.2,0.25,0.3,0.35,0.4 --subR=0.2,0.25,0.3,0.35,0.4 --leadMatchR=0.2,0.25,0.3,0.35,0.4 --subMatchR=0.2,0.25,0.3,0.35,0.4 --leadConstPt=1.0,1.5,2.0,2.5,3.0 --subConstPt=1.0,1.5,2.0,2.5,3.0 --towList=resources/bad_tower_lists/y7_y6_bad_tower_strict.txt --offAxisInput=resources/data_lists/y7_mb_file_list_rcf.txt --logtostderr --forceMatchJetResolutionEquality=true
```

p+p
```
./bin/dijet_imbalance/differential_aj_pp_y6_y7_embed --name=pp_y6_ --id=0 --input=resources/data_lists/pp_ht_rcf/pp_00.list --outputDir= results/pp/grid_16_8 --leadJetPt=16 --subJetPt=8 --leadR=0.2,0.25,0.3,0.35,0.4 --subR=0.2,0.25,0.3,0.35,0.4 --leadMatchR=0.2,0.25,0.3,0.35,0.4 --subMatchR=0.2,0.25,0.3,0.35,0.4 --leadConstPt=1.0,1.5,2.0,2.5,3.0 --subConstPt=1.0,1.5,2.0,2.5,3.0 --towList=resources/bad_tower_lists/y7_y6_bad_tower_strict.txt --embedInput=resources/data_lists/y7_mb_file_list_rcf.txt --efficiencyFile=resources/efficiencies/y7_effic.root --logtostderr --forceMatchJetResolutionEquality=true
```

These two examples are for the grid measurements (hard-core and matched, as well as off-axis background estimation). To perform measurements for the trees that produce the radial scan plot, simply change the following options to the specified values (along with appropriate changes to the output directory, name, etc):
```
--leadR=0.2
--subR=0.2
--forceMatchJetResolutionEquality=false
```

As you can see from the `--id` tag and that each of the analysis routines only takes a single input, these must be run multiple times, once for each input file, and in the case of the p+p, all input files must be run 5 times, once for the nominal case and four times for the systematic variations. Instead of submitting by hand to a scheduler, there are python submit scripts in the submit directory of the install location. These will iterate over the given input files, give each input a
unique id, and submit the job. Currently these are designed to work with PBS scheduling software, but should be easy to rewrite to work with condor. The two submit files that are necessary are:
```
submit/pbs_submit_dijet_imbalance_auau_y7.py
submit/pbs_submit_dijet_imbalance_pp_y6.py
```

The submit scripts use command line arguments that are similar to the command line options in the actual binaries, and can generally be translated one-to-one. The only major difference is that instead of using a command like `--input=`, you can use a wildcard for your input as such:
```
python submit/pbs_submit_dijet_imbalance_pp_y6.py resources/data_lists/pp_ht_rcf/*
```

Finally, the plot making binary expects the Au+Au data to be hadded into a single root file, and the  p+p data to be in a certain format. If you go into the output directory you specified for the `submit/pbs_submit_dijet_imbalance_pp_y6.py`, you should find five directories labeled 
```
tow_0_trk_0
tow_1_trk_0
tow_-1_trk_0
tow_0_trk_1
tow_0_trk_-1
```

tow_0_trk_0 is the nominal run (the datapoints used in the analysis), while the other four are our systematic variations. The plot-making binary expects five files to exist in a single directory:
```
tow_0_trk_0.root
tow_1_trk_0.root
tow_-1_trk_0.root
tow_0_trk_1.root
tow_0_trk_-1.root
```

Where tow_0_trk_0.root is the hadd'ed contents of the tow_0_trk_0 subdirectory of the p+p results.

## Producing plots

To produce all four grid plots from the paper, the submit scripts must be run for both the grid and the radial scan trees. Once all the data has been hadded (so you have two Au+Au files, one for the grid and one for the radial scan, and two p+p directories each with the required five root files), then producing the plots is relatively easy, with the following two commands.
```
./bin/dijet_imbalance/print/print_aj_auau_scan_paper --auau=/path/to/auau/grid/root/file --ppDir=/path/to/pp/grid/directory --setScanning=false --outputDir=results/grid_16_8_plots
./bin/dijet_imbalance/print/print_aj_auau_scan_paper --auau=/path/to/auau/scanning/root/file --ppDir=/path/to/pp/scanning/directory --setScanning=false --outputDir=results/scanning_16_8_plots
```

In the resulting directories grid_print and test_results, which hold the final plots and the ks test values respectively. The other directories all have individual plots for each di-jet definition.
