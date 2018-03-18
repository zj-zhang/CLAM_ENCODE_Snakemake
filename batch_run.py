import os
import sys

for project in os.listdir('projects'):
	#if not project.startswith('PTB'):
	if not project.startswith('HNRNPC'):
		continue
	os.system("./prepare_individual_run.sh "+project)
	#break
