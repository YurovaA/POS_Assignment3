#!/bin/bash 
#@ wall_clock_limit = 00:$(mins):00
#@ job_name = pos-cannon-mpi-ibm
#@ job_type = Parallel
#@ output = cannon_64_$(jobid).out
#@ error = cannon_64_$(jobid).out
#@ class = test
#@ node = 4
#@ total_tasks = 64
#@ node_usage = not_shared
#@ energy_policy_tag = cannon
#@ minimize_time_to_solution = yes
#@ notification = never
#@ island_count = 1
#@ queue


. /etc/profile
. /etc/profile.d/modules.sh
module load python/2.7.5

date
python execute_n.py $(times) "mpiexec -n 64 ./$(exe) ../cannon_matrices/$(size)x$(size)-1.in ../cannon_matrices/$(size)x$(size)-2.in" $(mins_stop)
date
