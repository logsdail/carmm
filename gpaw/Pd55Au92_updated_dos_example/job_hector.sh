#!/bin/bash --login
#PBS -N gpaw
#PBS -A e05-react-sok
#PBS -l mppwidth=1024
#PBS -l mppnppn=32
#PBS -l walltime=03:00:00

module del PrgEnv-pgi
module del cce
module add PrgEnv-gnu/4.0.30
module swap xt-libsci/12.0.00 xt-libsci/11.0.06
module add fftw/3.3.0.0
module use /work/y07/y07/nag/packages/gccmod/modules
module add xe-cblas/2003
module add xe-expat/2.0.1
module add xe-zlib/1.2.6
module add xe-libpng/1.5.10
module add xe-freetype/2.4.9
module add xe-ghostscript/9.05
module add xe-nose/1.1.2
module add xe-numpy/1.6.1
module add xe-matplotlib/1.1.0
module add xe-ase/3.6.0
module add xe-gpaw/0.9.0.9624

export LD_LIBRARY_PATH=/opt/xt-libsci/11.0.06/gnu/46/interlagos/lib/:$LD_LIBRARY_PATH
export PATH=.:$PATH
cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=1

export MYJOB=`qstat -f $PBS_JOBID | grep Job_Name | awk '{print $3}'`
export NTASK=`qstat -f $PBS_JOBID | grep mppnppn  | awk '{print $3}'`
export NPROC=`qstat -f $PBS_JOBID | grep mppwidth | awk '{print $3}'`

aprun -n $NPROC -N $NTASK gpaw-python cube_and_dos.py 
