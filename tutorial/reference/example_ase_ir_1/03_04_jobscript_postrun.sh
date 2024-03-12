#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q xs2.q
#$ -pe x16 16
#$ -j y
#$ -N calc_03_04_postrun

module load qe/7.2-mpi
mymamba
mamba activate mtkenv

export I_MPI_PIN=1
export OMP_NUM_THREADS=1
export I_MPI_FABRICS=shm:ofi

cat $PE_HOSTFILE | awk '{ print $1":"$2/ENVIRON["OMP_NUM_THREADS"] }' > hostfile.$JOB_ID

echo "========= Job started  at `date` =========="

python 03_dipole_collector.py
python 04_postprocess.py

echo "========= Job finished at `date` =========="


# cleanup
rm -rf hostfile.$JOB_ID

