#expects: script.py burnin
#input: beast .trees files in cwd
#output: filename_resampled_clean.trees in cwd


import os
import glob
import sys

burnin = int(sys.argv[1])

treeslist = glob.glob('*.trees')


for tree in treeslist:
	filestem = tree.split('/')[-1].split('.')[0]
	cmd = 'logcombiner -trees -burnin %d %s ./%s_resampled.trees'%(burnin, tree, filestem)
	os.system(cmd)
