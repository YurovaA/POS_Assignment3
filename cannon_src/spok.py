import os
import os.path
import subprocess;
from subprocess import call
from time import sleep
from random import randint
import re
from math import sqrt

sigmas = 1
runs = 20

def means_devs_ci(dct_of_lst):
	print('FUNC ' * 10)
	dev_dct = {}
	mean_dct = {}
	ci_dct = {}
	for key,lst in dct_of_lst.items():
		lst_mean = sum(lst) / len(lst)
		lst_sq = [x**2 for x in lst]
		lst_sq_mean = sum(lst_sq) / len(lst_sq)
		print('{0}: {1}'.format(key, lst_sq_mean - lst_mean**2))
		dev_dct[key] = sqrt(max(lst_sq_mean - lst_mean**2, 0))
		mean_dct[key] = lst_mean
		# percentage of 95% confidence interval
		ci_dct[key] = (sigmas * dev_dct[key] / sqrt(len(lst))) / lst_mean * 100
			
	return mean_dct, dev_dct, ci_dct




print("Start Spok")
#def_batch_name = 'll-64-ibm.sh'
#def_batch_name = 'batch_optim.sh'
def_batch_name = 'batch_small.sh'

comp_results_dct = {}
mpi_results_dct = {}
exit = False
counter = 0


while not exit:
	rnd = randint(0, 100100)
	new_batch_name = def_batch_name.replace('.sh', '_{0}.sh'.format(rnd))
	call(['cp', def_batch_name, new_batch_name])
	# Here we have created the copy of a batch file. Let's chage the name of the output file
	b_file = open(new_batch_name, 'r')
	cont = b_file.read()
	cont = cont.replace('$(jobid)', '{0}'.format(rnd))

	b_file.close()
	# Rewriting the content
	b_file = open(new_batch_name, 'w')
	b_file.write(cont)
	b_file.close()


	#submit = subprocess.check_output("llsubmit " + new_batch_name + " > dummy ", shell=True)
	print(60 * 'o')
	call(['llsubmit', new_batch_name])	
	print(60 * 'o')

	running = True
	# Waiting for output
	while running:
		sleep(1)
		out = subprocess.check_output("llq -u h039y11", shell=True)
		if (out.startswith('llq: The')):
			running = False
		
	# Parsing output
	output_file_name = 'cannon_64_{0}.out'.format(rnd)
	while not os.path.isfile(output_file_name):
		sleep(1)
	with open(output_file_name, 'r') as output_file:
		# Variable containing the current matrix size
		cur_matrix_size = 0
		for line in output_file:
			if line.startswith('('):
				cur_matrix_size = int(re.findall(r'[1-9][0-9]*', line)[0])
				if not cur_matrix_size in comp_results_dct:
					comp_results_dct[cur_matrix_size] = []
					mpi_results_dct[cur_matrix_size] = []
			elif line.startswith('Computation time:'):
				comp_time = float(re.findall('\d+\.\d+', line)[0])
				comp_results_dct[cur_matrix_size].append(comp_time)
			elif line.startswith('MPI time:'):
				mpi_time = float(re.findall('\d+\.\d+', line)[0])
				mpi_results_dct[cur_matrix_size].append(mpi_time)
	#print(comp_results_dct)
	#print(mpi_results_dct)
	
	counter += 1
	# compute intervals
	if counter > 2:
		print(60 * '*')
		print(counter)
		c_means, c_devs, c_ci = means_devs_ci(comp_results_dct)
		m_means, m_devs, m_ci = means_devs_ci(mpi_results_dct)		
		print('Max confidence interval (%) in computation time is {0}'.format(max(c_ci.values() ) ) )
		print('Max confidence interval (%) in mpi time is {0}'.format(max(m_ci.values() ) ) )
		#print(m_ci.values())
		print(m_means.values())
		#print(m_devs.values())
		print(60 * '*')
		
	if counter == runs:
		exit = True
	
c_means, c_devs, c_ci = means_devs_ci(comp_results_dct)
m_means, m_devs, m_ci = means_devs_ci(mpi_results_dct)	
with open('stat.csv', 'w') as stat_file:
	stat_file.write('\n')
	stat_file.write(';Computation time statistics\n;')
	for key, val in c_means.items():#:#:
		stat_file.write(str(key) + ';')
	stat_file.write('\nmeans:;')
	for key, val in c_means.items():#:#:
		stat_file.write(str(val) + ';')
	stat_file.write('\ndevs:;')
	for key, val in c_devs.items():#:#:
		stat_file.write(str(val) + ';')
	stat_file.write('\nconf.int({0}):;'.format(sigmas))
	for key, val in c_ci.items():#:#:
		stat_file.write(str(val) + ';')
		
	stat_file.write('\n\n\n')
	stat_file.write(';MPI time statistics\n;')
	for key, val in c_means.items():#:
		stat_file.write(str(key) + ';')
	stat_file.write('\nmeans:;')
	for key, val in m_means.items():#:
		stat_file.write(str(val) + ';')
	stat_file.write('\ndevs:;')
	for key, val in m_devs.items():#:
		stat_file.write(str(val) + ';')
	stat_file.write('\nconf.int({0}):;'.format(sigmas))
	for key, val in m_ci.items():#:
		stat_file.write(str(val) + ';')
	
	stat_file.write('\n\n\n')

				
print("End Spok")

