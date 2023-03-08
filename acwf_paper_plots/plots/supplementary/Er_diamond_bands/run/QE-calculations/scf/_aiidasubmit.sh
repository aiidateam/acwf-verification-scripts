#!/bin/bash
exec > _scheduler-stdout.txt
exec 2> _scheduler-stderr.txt


eval "$(conda shell.posix hook)"
conda activate /home/jovyan/.conda/envs/quantum-espresso-7.0
export OMP_NUM_THREADS=1

'mpirun' '-np' '8' '/home/jovyan/.conda/envs/quantum-espresso-7.0/bin/pw.x' '-in' 'aiida.in'  > 'aiida.out' 
