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

# 	
# resampledlist = glob.glob('*_resampled.trees')
# for tree in resampledlist:
# 	treefile = open(tree, 'rw')
# 	contents = treefile.read()
# 	newcontents = contents.replace("COTED'IVOIRE", "COTEDIVOIRE").replace("M29973|SE_id250619|Cercopithecus_aethiops_aethiops|sub_given|GRV|GRI_677_gri_1_677, Grived (gr-1)|len2438|AFRSSA|ETHIOPIA||", "M29973|SE_id250619|Cercopithecus_aethiops_aethiops|sub_given|GRV|GRI_677_gri_1_677_Grived_(gr-1)|len2438|AFRSSA|ETHIOPIA||")
# 	outfileName = tree.split('resampled')[0]+'_clean_resampled.trees'
# 	outfile = open(outfileName, 'w')
# 	outfile.write(newcontents)
# 	os.system('rm '+tree)