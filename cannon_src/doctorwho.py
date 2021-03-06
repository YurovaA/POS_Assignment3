import sys
import subprocess;
from subprocess import call
from time import sleep
import os
import re
from math import sqrt

# STAT_{NAME}


# times is the number of times we want to start a batch
# size is the size of the matrix we want to test
def gen_batch(path, path_new, mins, times, size, exe):
	jobid = '{0}_{1}_{2}'.format(mins, times, size)
	batch_content = ''
	with open(path, 'r') as batch_template:
		batch_content = batch_template.read()
		str_mins = str(mins)
		if len(str_mins) == 1:
			str_mins = '0' + str_mins
		batch_content = batch_content.replace("$(mins)", str_mins).replace("$(times)", str(times)).replace("$(size)", str(size)).replace("$(jobid)", str(jobid)).replace("$(exe)", exe).replace("$(mins_stop)", str(max(20, mins - 3)))
		with open(path_new, 'w') as batch_new:
			batch_new.write(batch_content)
	return 'cannon_64_' + jobid + '.out'

def sys_args_to_dct():
	ret_dct = {}
	argc = len(sys.argv)
	for i in range(1, argc):
		arg = sys.argv[i]
		comps = arg.split('=')
		if len(comps) == 2 and comps[0].startswith('-'):		
			ret_dct[comps[0][1:]] = comps[1]
	return ret_dct

def means_devs_ci(dct_of_lst):
	dev_dct = {}
	mean_dct = {}
	ci_dct = {}
	for key,lst in dct_of_lst.items():
		lst_mean = sum(lst) / len(lst)
		lst_sq = [x**2 for x in lst]
		lst_sq_mean = sum(lst_sq) / len(lst_sq)
		#print('{0}: {1}'.format(key, lst_sq_mean - lst_mean**2))
		dev_dct[key] = sqrt(max(lst_sq_mean - lst_mean**2, 0))
		mean_dct[key] = lst_mean
		# percentage of 95% confidence interval
		ci_dct[key] = (sigmas * dev_dct[key] / sqrt(len(lst))) / lst_mean * 100
			
	return mean_dct, dev_dct, ci_dct		

def ParseDcts(file_name, comp_results_dct, mpi_results_dct):
	while not os.path.isfile(file_name):
		sleep(1)
	with open(file_name, 'r') as output_file:
		# Variable containing the current matrix size
		cur_matrix_size = 0
		for line in output_file:
			if line.startswith('('):
				cur_matrix_size = int(re.findall(r'[1-9][0-9]*', line)[0])				
				if not cur_matrix_size in comp_results_dct:
					print(cur_matrix_size)
					comp_results_dct[cur_matrix_size] = []
					mpi_results_dct[cur_matrix_size] = []
			elif line.startswith('Computation time:'):
				comp_time = float(re.findall('\d+\.\d+', line)[0])
				comp_results_dct[cur_matrix_size].append(comp_time)
			elif line.startswith('MPI time:'):
				mpi_time = float(re.findall('\d+\.\d+', line)[0])
				mpi_results_dct[cur_matrix_size].append(mpi_time)
	return comp_results_dct, mpi_results_dct

def ParseDctsOfSize(file_name, comp_results_dct, mpi_results_dct, size, times):
	while not os.path.isfile(file_name):
		sleep(1)
	if not size in comp_results_dct:
		print("Key {0} created".format(size))
		comp_results_dct[size] = []
		mpi_results_dct[size] = []
	
	times_in_file = 0
	while times_in_file != times:
		times_in_file = 0
		with open(file_name, 'r') as output_file:
			for line in output_file:
				if line.startswith('MPI time:'):
					times_in_file += 1		
				
	
	with open(file_name, 'r') as output_file:
		for line in output_file:
			if line.startswith('Computation time:'):
				comp_time = float(re.findall('\d+\.\d+', line)[0])
				comp_results_dct[size].append(comp_time)
			elif line.startswith('MPI time:'):
				mpi_time = float(re.findall('\d+\.\d+', line)[0])
				mpi_results_dct[size].append(mpi_time)
	return comp_results_dct, mpi_results_dct


def AdjustList(mpi_res, mtx_size):
	if len(mpi_res[mtx_size]) <= 1:
		raise NameError("List too small")
	u = max(mpi_res[mtx_size])
	d = min(mpi_res[mtx_size])
	count_u = -1
	count_new = 0
	lst_u = []
	lst_d = []
	while count_u != count_new:
		count_u = count_new
		lst_u, lst_d = SeparateList(mpi_res[mtx_size], u, d)
		u = sum(lst_u) / len(lst_u)
		d = sum(lst_d) / len(lst_d)
		count_new = len(lst_u)
	lst_u = list(map(lambda x: x - (u - d), lst_u))
	ret_lst = lst_u + lst_d
	return ret_lst
	
def SeparateList(lst, u, d):
	lst_u = []
	lst_d = []
	for el in lst:
		if abs(el - u) <= abs(el - d):
			lst_u.append(el)
		else:
			lst_d.append(el)
	return lst_u, lst_d

def WaitOutput():
	running = True
	# Waiting for output
	while running:
		sleep(1)
		out = subprocess.check_output("llq -u h039y11", shell=True)
		if (out.startswith('llq: The')):
			running = False

precision = 0.01
maxmins = 20
mtx_size = 64
exe_file = "cannon"
sigmas = 1
merge = 0
prelim_size = 15

print("Command line arguments: ")
arg_dct = sys_args_to_dct()
print(arg_dct)
if "prec" in arg_dct:
	precision = float(arg_dct["prec"])
	print("Precision set to {0}".format(precision))	
if "maxmins" in arg_dct:
	maxmins = int(arg_dct["maxmins"])
	print("Max time set to {0} minutes".format(maxmins))
if "size" in arg_dct:
	mtx_size = int(arg_dct["size"])
	print("Max matrix size set to {0}".format(mtx_size))
if "exe" in arg_dct:
	exe_file = arg_dct["exe"]
	print("Exec file set to " + exe_file)
if "sigmas" in arg_dct:
	sigmas = float(arg_dct["sigmas"])
	print("Sigmas set to {0}".format(sigmas))
if "merge" in arg_dct:
	merge = int(arg_dct["merge"])
	print("Merge set to {0}".format(merge))
if "prelim" in arg_dct:
	prelim_size = int(arg_dct["prelim"])
	print("Preliminary size set to {0}".format(prelim_size))



# Start generating preliminary statistics	
batch_new = "batch_prelim.sh"
output_prelim = gen_batch("batch_template.sh", batch_new, maxmins, prelim_size, mtx_size, exe_file)
print(60 * 'o')
call(['rm', output_prelim])
call(['llsubmit', batch_new])	
print(60 * 'o')

WaitOutput()

comp_res = {}
mpi_res = {}
comp_res, mpi_res = ParseDctsOfSize(output_prelim, comp_res, mpi_res, mtx_size, prelim_size)
mpi_dct = mpi_res
print("SIZE: {0}".format(len(mpi_dct[mtx_size])))
if merge == 1:
	new_mpi = AdjustList(mpi_res, mtx_size)
	print(mpi_res[mtx_size])
	print(new_mpi)
	mpi_dct = { mtx_size: new_mpi }
	print("SIZ in: {0}".format(len(mpi_dct[mtx_size])))
# Now we can estimate the variation
m_means, m_devs, m_ci = means_devs_ci(mpi_dct)
sigma = m_devs[mtx_size]
mean = m_means[mtx_size]
# Estimated sample size to reach a desired precision
est_n = (sigma / (precision * mean))**2
print('Answer ' * 8)
print(est_n)
print('Answer ' * 8)

# Now we are doing the main run
batch_main = "batch_main_{0}.sh".format(mtx_size)
if prelim_size < est_n:
	main_runs = int(est_n) - prelim_size
	output_main = 'main_' + gen_batch("batch_template.sh", batch_main, maxmins, main_runs, mtx_size, exe_file)
	print(60 * 'q')
	call(['rm', output_main])
	call(['llsubmit', batch_main])	
	print(60 * 'q')
	
	WaitOutput()
	comp_res, mpi_res = ParseDctsOfSize(output_main, comp_res, mpi_res, mtx_size, main_runs)
	mpi_dct = mpi_res
	print("SIZE main: {0}".format(len(mpi_dct[mtx_size])))
	if merge == 1:
		new_mpi = AdjustList(mpi_res, mtx_size)
		print(mpi_res[mtx_size])
		print(new_mpi)
		mpi_dct = { mtx_size: new_mpi }
		print("SIZE main in: {0}".format(len(mpi_dct[mtx_size])))

c_means, c_devs, c_ci = means_devs_ci(comp_res)
print("SIZE: {0}".format(len(mpi_dct[mtx_size])))
m_means, m_devs, m_ci = means_devs_ci(mpi_dct)	
with open('stat_{0}.csv'.format(mtx_size), 'w') as stat_file:
	stat_file.write('\n')
	stat_file.write(';Runs:;{0};;Size:;{1}'.format(int(est_n), mtx_size))
	stat_file.write('\n;Computation time statistics\n;')
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

