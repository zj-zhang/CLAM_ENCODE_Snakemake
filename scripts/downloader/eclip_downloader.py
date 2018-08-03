## read the url from stdin and call download_py
## Zijun Zhang
## 3.7.2018

import sys
import os
from download_py import *

max_try = 5

for line in sys.stdin:
	ele = line.strip().split()
	url, file_path, hash = ele[0:3]
	dirname = os.path.dirname(os.path.realpath(file_path))
	if not os.path.isdir(dirname):
		os.mkdir(dirname)
	for i in range(max_try):
		print("%s, try %i.."%(file_path, i) )
		try:
			download_with_resume(url, file_path, hash)
		except:
			print("failed. resume")
		finally:
			break