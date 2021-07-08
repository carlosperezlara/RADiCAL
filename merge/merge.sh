!/bin/sh
run=${1}
echo $run
fileDRS="/Volumes/uva/testbeam_2021_06_data/DRS_data/Run_${run}.dat"
fileOut="/Volumes/uva/testbeam_2021_06_data/merged/Run_${run}.root"
echo $fileDRS
echo $fileOut
cp $fileDRS ./drs.dat
./maketree --inputFileName=drs.dat --outputFileName=out.root --nEvents=1000000
rm drs.dat
mv ./out.root ${fileOut}
