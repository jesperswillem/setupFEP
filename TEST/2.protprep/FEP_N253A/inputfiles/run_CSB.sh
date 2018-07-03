#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#              d-hh:mm:ss
#SBATCH --time=0-03:00:00

export TMP=/tmp
export TEMP=/tmp
export TMPDIR=/tmp
## Load modules for qdynp
module load openmpi-x86_64


## define qdynp location
qdyn=/home/apps/q-5.06/qdynp
fepfiles=(FEP1.fep FEP2.fep FEP3.fep FEP4.fep FEP5.fep)
temperature=298
run=10
finalMDrestart=md_0000_1000.re

workdir=/home/jespers/adenosine/1.A1-A2A_selectivity/A1/5.FEP/holo
inputfiles=/home/jespers/adenosine/1.A1-A2A_selectivity/A1/5.FEP/holo/inputfiles
length=${#fepfiles[@]}
length=$((length-1))
for index in $(seq 0 $length);do
fepfile=${fepfiles[$index]}
fepdir=$workdir/FEP$((index+1))
mkdir -p $fepdir
cd $fepdir
tempdir=$fepdir/$temperature
mkdir -p $tempdir
cd $tempdir

rundir=$tempdir/$run
mkdir -p $rundir
cd $rundir

cp $inputfiles/md*.inp .
cp $inputfiles/*.top .
cp $inputfiles/$fepfile .
cp $inputfiles/run_0500-1000.sh .
cp $inputfiles/run_0500-0000.sh .

if [ $index -lt 1 ]; then
cp $inputfiles/eq*.inp .
sed -i s/SEED_VAR/"$[1 + $[RANDOM % 9999]]"/ eq1.inp
else
lastfep=FEP$index
cp $workdir/$lastfep/$temperature/$run/$finalMDrestart $rundir/eq5.re
fi

sed -i s/T_VAR/"$temperature"/ *.inp
sed -i s/FEP_VAR/"$fepfile"/ *.inp
if [ $index -lt 1 ]; then
#time mpirun -np 16 $qdyn eq1.inp > eq1.log
#EQ_FILES
time mpirun -np 16 $qdyn eq1.inp> eq1.log
time mpirun -np 16 $qdyn eq2.inp> eq2.log
time mpirun -np 16 $qdyn eq3.inp> eq3.log
time mpirun -np 16 $qdyn eq4.inp> eq4.log
time mpirun -np 16 $qdyn eq5.inp> eq5.log
fi
#RUN_FILES
time mpirun -np 16 $qdyn md_1000_0000.inp> md_1000_0000.log
time mpirun -np 16 $qdyn md_0950_0050.inp> md_0950_0050.log
time mpirun -np 16 $qdyn md_0900_0100.inp> md_0900_0100.log
time mpirun -np 16 $qdyn md_0850_0150.inp> md_0850_0150.log
time mpirun -np 16 $qdyn md_0800_0200.inp> md_0800_0200.log
time mpirun -np 16 $qdyn md_0750_0250.inp> md_0750_0250.log
time mpirun -np 16 $qdyn md_0700_0300.inp> md_0700_0300.log
time mpirun -np 16 $qdyn md_0650_0350.inp> md_0650_0350.log
time mpirun -np 16 $qdyn md_0600_0400.inp> md_0600_0400.log
time mpirun -np 16 $qdyn md_0550_0450.inp> md_0550_0450.log
time mpirun -np 16 $qdyn md_0500_0500.inp> md_0500_0500.log
time mpirun -np 16 $qdyn md_0450_0550.inp> md_0450_0550.log
time mpirun -np 16 $qdyn md_0400_0600.inp> md_0400_0600.log
time mpirun -np 16 $qdyn md_0350_0650.inp> md_0350_0650.log
time mpirun -np 16 $qdyn md_0300_0700.inp> md_0300_0700.log
time mpirun -np 16 $qdyn md_0250_0750.inp> md_0250_0750.log
time mpirun -np 16 $qdyn md_0200_0800.inp> md_0200_0800.log
time mpirun -np 16 $qdyn md_0150_0850.inp> md_0150_0850.log
time mpirun -np 16 $qdyn md_0100_0900.inp> md_0100_0900.log
time mpirun -np 16 $qdyn md_0050_0950.inp> md_0050_0950.log
time mpirun -np 16 $qdyn md_0000_1000.inp> md_0000_1000.log
done
