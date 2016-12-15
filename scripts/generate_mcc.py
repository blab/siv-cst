#input: beast .trees files in cwd
#does: returns .mcc for each file; burnin specified as number of states as cmd line argument.
#expects: cwd$ python script.py N hosttemplate
# where N is the integer state of burnin to remove, and hosttemplate is optional; if included, will append the figtree template
#requires: os, glob, sys

import sys
import os
import glob

cwd = os.getcwd()

treeslist = glob.glob('*.trees')

print 'found %d tree files'%len(treeslist)

figtree_template = '/Users/Sidney/Dropbox/siv/scripts/beastAnalysis/figtreeTemplates/host_figtree_template.txt'

for tree in treeslist:
	analysis, segment_id, date = tree.split('/')[-1].split('.')[0].split('_')
	segment_id = str(int(segment_id)-1)
	filestem = '%s_%s_%s'%('discreteTraits', segment_id, date)

	cmd = 'treeannotator -burnin %d -heights mean %s %s.mcc'%(int(sys.argv[1]), tree, filestem)
	os.system(cmd)
	print 'Found mcc tree for %s'%tree

	if 'hosttemplate' in sys.argv:
		mcc = '%s_temp.mcc'%filestem
		mcc_with_template = '%s.mcc'%filestem
		os.system('cat %s %s > %s'%(mcc, figtree_template, mcc_with_template))
		os.system('rm %s'%mcc)

print 'done'
