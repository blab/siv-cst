#input: beast .log file with actualRates and indicator parameters; xml with list of state classes (in order)
#output: rates and bf: csv with formatted matrices, heatmaps. primitive network diagram.
#burnin must be a sampled state; removal is inclusive

print 'loading modules'
import argparse
import elementtree.ElementTree as ET
import numpy as np
import pandas as pd
import beastmatrix as bm
import pickle

#################### Parse and reformat input #####################
print 'parsing input'

parrot = argparse.ArgumentParser(description="Analyze log files; make figures.")
parrot.add_argument('-burnin', default=1000000, type=int, help="Number of states to remove as burnin.")
parrot.add_argument('-cutoff', default=None, type=float, help="BF Cutoff for Network Diagram.")
parrot.add_argument('-logfile', type=str, help=".log file to be analyzed.")
parrot.add_argument('-xml', type=str, help=".xml to pull states and priors from.")
parrot.add_argument('-trait', default='host', type=str, help="trait to calculate transition rates for")
parrot.add_argument('-prior', default=None, type=float, help="prior expectation for trait.nonzerorates parameter")
parrot.add_argument('-make_figs', default=False, type=bool, help="generate heatmaps for AR & indicators, and network diagram?")

args = vars(parrot.parse_args())
burnin, logfile, xmlfile, cutoff, trait, prior, make_figs = args['burnin'], open(args['logfile'], 'r'), open(args['xml'], 'r'), args['cutoff'], args['trait'], args['prior'], args['make_figs']
common_names = pickle.load(open('/Users/Sidney/Dropbox/siv-manuscript/data/hosts/common_names.p', 'rb'))

log_data = pd.read_csv(logfile, comment='#', sep="\t", index_col='state')	#parse the log file as a df
logfile.close()

xml = ET.parse(xmlfile)
root = xml.getroot()												#get the list of states, in order, from the xml
posterior = root.find('mcmc').find('posterior')
n_segments = float(len(posterior.find('likelihood')))
priors = posterior.find('prior')
states = [ i for i in root.findall('generalDataType') if i.get('id').startswith(trait) ][0]
state_list = [ common_names[state.get('code') ] for state in states]
n_demes = len(state_list)
print 'Found %d demes for trait %s. Order:'%(n_demes, trait)
print state_list

if prior == None:
	for p in priors:	# find the prior expectation for state.nonzerorates.
		ID = p[0].get('idref')
		if ID =='%s.nonZeroRates'%trait:
			try:
				prior_expectation = float(p.attrib['mean'])
				print 'found prior expectation of %.1f from block with ID %s'%(prior_expectation, ID)
			except:
				prior_expectation = n_demes-1
				print 'ERROR: non-exponential prior on %s.nonZeroRates;\nusing prior expectation of %d for now, but you should rerun with a provided prior expectation for accurate BFs.'%(trait, prior_expectation)
			break
		else:
			continue
else:
	prior_expectation = prior


if cutoff == None:												#Set a cutoff for the network diagram with 3.2*number of segments
	cutoff = 3.0*n_segments
	print 'Found %d segments; using prior cutoff of %f'%(n_segments, cutoff)
	print 'Note this assumes that the length of the likelihood block == number of nonindependent genomic segments contributing to the likelihood.'

else:
	print 'Using provided cutoff of %.1f, after accounting for %d nonindependent segments.'%(cutoff, n_segments)
	cutoff=cutoff*n_segments

xmlfile.close()

################# Deal with burnin, make data structures ###############################

if burnin not in log_data.index:		#drop (inclusively) all states before the specified burnin state.
 	assert burnin in log_data.index, ('ERROR: not a valid burnin value. Valid options:', log_data.index)
else:
	burnin_indices = range(0, log_data.index.get_loc(burnin)+1)
	log_data.drop(log_data.index[burnin_indices], inplace=True)

log_data = log_data.append(pd.DataFrame(log_data.mean(), columns=['avg']).T)	#get the average posterior values for each column

state_actualRates = pd.DataFrame(dtype = 'float', index=state_list, columns=state_list)	#make empty DFs to hold the actual rates and the bf
bf = pd.DataFrame(dtype = 'float', index=state_list, columns=state_list)

#############  Fill data structures ###################
print 'filling data structures'

state_actualRates_series = []
bf_series = []

for column, series in log_data.iteritems():												#pull the average posterior values from appropriate columns in order
	if 'actualRates' in column and trait in column:															#we will then use this to fill our matrices.
		state_actualRates_series.append(series['avg'])
	elif 'indicator' in column and trait in column:
		bf_series.append(bm.find_bf(series['avg'], prior_expectation, n_demes))		#convert indicators to bayes factors while we're at it
	else:
		continue

state_actualRates = bm.fill_rates_matrix(state_actualRates, pd.Series(state_actualRates_series))				#fill our matrices
bf = bm.fill_rates_matrix(bf, pd.Series(bf_series))

file_stem = args['logfile'].split('/')[-1].split('.')[0]

state_actualRates.to_csv(file_stem+'_actualRates.csv')									#and write them to file
bf.to_csv(file_stem+'_bf.csv')

######################### Make some pretty heatmaps ################################################
if not make_figs:
	import sys
	sys.exit()
else:
	import networkx as nx
	import matplotlib.pyplot as plt
	import seaborn as sns

bfname = "%s_bf.png"%file_stem
actualratesname = "%s_AR.png"%file_stem
mapname = "%s_transMap.png"%file_stem

print 'making heatmaps'

def make_heatmap(matrix, title, outfilename):		#make simple heatmaps; save to file.
    plt.figure(figsize=(13,9))
    heatmap = sns.heatmap(matrix, cmap = plt.cm.GnBu, square = True, linecolor = 'black', robust=True, linewidths = 0.5)
    heatmap.axes.set_title(title)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    plt.tight_layout(pad=1)
    plt.savefig(outfilename)
    plt.clf()

make_heatmap(state_actualRates, 'Actual Rates of %s Transition'%trait, actualratesname)
make_heatmap(bf, 'Bayes Factors: %s Transitions'%trait, bfname)

######################### Make a pretty [ugly] network diagram ########################################

print 'making network map'

n_demes = len(state_actualRates.columns)
nparams = len(state_actualRates_series)

if nparams > ((n_demes**2)-n_demes)/2:
	G = nx.DiGraph()
	directed = True
else:
	G = nx.Graph()
	directed = False

															#initialize a networkx graph
G.add_nodes_from(state_actualRates.columns)					#each state is a node

for state, rates in state_actualRates.iteritems():
	for state2 in state_actualRates.index:
		BF = bf[state][state2]								#pull rate and bf values from the two dataframes fro each pair of states.
		rate = float(state_actualRates[state][state2])
		if BF >= cutoff and rate > 0:
			G.add_edge(state2,state,{'rate': '%.3f'%rate, 'support': BF})
		else:
			continue


#networkx apparently can't access its native node names when graphing.
node_labels = {}
for node,d in G.nodes(data=True):
	node_labels[node] = node

#nor can it access edge attributes as weights. so, pull them manually.
bf = [float(G[u][v]['support']) for u,v in G.edges()]
rates = [float(G[u][v]['rate'])*1.5 for u,v in G.edges()]

#make network; show, rather than saving to file.
plt.figure(figsize=(15,10))
plt.figtext(0.1,0.1,'BF >= %.2f'%cutoff)
pos = nx.spring_layout(G, iterations=20, k=2)	#set layout
nx.draw(G, pos, edges=G.edges(), width=rates, edge_color='dimgray', alpha=0.6, node_size = 2000)
nx.draw_networkx_labels(G, pos, node_labels=node_labels, fontsize=28, font_weight='bold')
plt.bbox_inches="tight"
plt.savefig(mapname)
