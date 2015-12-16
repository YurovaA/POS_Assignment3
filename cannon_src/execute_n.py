import sys
import subprocess;
from subprocess import call
from datetime import datetime as dt

num = int(sys.argv[1])
command = sys.argv[2]
max_secs = int(sys.argv[3]) * 60

start = dt.now()
for i in range(num):
	#out = subprocess.check_output(command, shell=True)
	call(command.split(' '))
	if (dt.now() - start).seconds > max_secs:
		break
