import sys
import argparse
import os
import re
import csv


__author__ = 'anna'
parser=argparse.ArgumentParser(
     description=''' argument1 - folder with stat files, argument2 - name of the output''')
parser.add_argument('folder', type=str, help='folder with stat files')
parser.add_argument('out_name', type=str, help = 'the name of the output file')

args = parser.parse_args()

folder = args.folder
output = args.out_name
ct_means = {}
ct_devs = {}
ct_confint = {}
mpi_means = {}
mpi_devs = {}
mpi_confint = {}
            		
results = {}
for root, dirs, filenames in os.walk(folder):
    for f in filenames:
    	if (re.match('^stat', f)):
    		print(f)
        	size = f.split('_')[1].split('.')[0]
        	with open(folder + f, 'r') as thefile:
            		lines = [x.strip().upper() for x in thefile.readlines()]
            		#print(lines)
            		#//current_line = lines.index(student.upper()) + 1
            		current_line = 0
            		while current_line < len(lines) and not re.match('^;MPI', lines[current_line]):
                		if re.match('^MEANS', lines[current_line]):
                	    		ct_means[size] = lines[current_line].split(';')[1]
                		if re.match('^DEVS', lines[current_line]):
                	    		ct_devs[size] = lines[current_line].split(';')[1]
                		if re.match('^CONF', lines[current_line]):
                	    		ct_confint[size] = lines[current_line].split(';')[1]
                		current_line += 1
            		
            		while current_line < len(lines):
                		if re.match('^MEANS', lines[current_line]):
                	    		mpi_means[size] = lines[current_line].split(';')[1]
                		if re.match('^DEVS', lines[current_line]):
                	    		mpi_devs[size] = lines[current_line].split(';')[1]
                		if re.match('^CONF', lines[current_line]):
                	    		mpi_confint[size] = lines[current_line].split(';')[1]
                		current_line += 1

keys = ct_means.keys()
file = open(output, 'w')
file.write('Computational time: \n')
file.write('key, mean, dev, conf\n')
for k in keys:
    line = str(k) + ',' + str(ct_means.get(k)) + ',' + str(ct_devs.get(k)) + ',' + str(ct_confint.get(k)) + '\n'
    print(line)
    file.write(line)

file.write('MPI time: \n')
file.write('key, mean, dev, conf\n')
for k in keys:
    line = str(k) + ',' + str(mpi_means.get(k)) + ',' + str(mpi_devs.get(k)) + ',' + str(mpi_confint.get(k)) + '\n'
    file.write(line)


