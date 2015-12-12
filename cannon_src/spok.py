import os
import subprocess;
from subprocess import call
from time import sleep

print("Start Spok")
call(['llsubmit', 'll-64-ibm.sh'])
#subprocess.call(['ls', '-l'])

running = True
while running:
	sleep(1)
	#qq = subprocess.check_output("cat Makefile | grep cannon", shell=True)
	output = subprocess.check_output("llq | grep  h039y11", shell=True)
	if (len(output) > 0):
		running = False
		
print("End Spok")

