#input: script.py master.xml (.xml file you want to pull parameters from); additional .xml files in the cwd
#output: combined .xml file with independent <blocks> for taxa, empiricaltrees, ancestraltreelikelihood
#cwd: should contain masterfile.xml, all other .xml input files. Will check and warn for corresponding .trees files, but not required to construct xml. BEAST will require the xml and .trees files in the same cwd.
#requires: sys, glob, elementtree, datetime

############################## Initialize, find input ###############################
import sys
import elementtree.ElementTree as ET
from datetime import date
from glob import glob

dateTag = str(date.today())
masterfilename = sys.argv[1].split('/')[-1]

treeslist = glob('*.trees')
xmllist = glob('*.xml')

for file in xmllist:
	if file == masterfilename or file.startswith('discreteTraits'):
		xmllist.remove(file)
	else:
		pass

assert xmllist != []

if len(xmllist)!=len(treeslist):
	print 'Warning: you may be missing .trees files or have extra .xml files in your cwd.'

if xmllist[0].startswith('h') and len(xmllist[0].split('_')) >= 3:
	xmllist.sort(key=lambda x: int(x.split('_')[1]))
else:
	xmllist.sort()

########### Define our parsing functions #####################################

#takes a file and the name of the block we want; returns the block(s) as an element object
def getblock(file, name):
	input = open(file, 'r')
	tree = ET.parse(input)
	root = tree.getroot()
	return root.find(name)
	input.close()

#takes a block, the tag we want to rename, and a counter(#); returns the block with '.#' appended to the tag name.
def renameblock(block, tag, counter):
	oldname = block.get(tag)
	block.set(tag, oldname+'.'+str(counter))
	return block 

#takes the .xml filename, list of desired blocks, and the counter # we want to use to rename the tags.
#returns a dictionary of filename = {'blockA' = <blockA.#>, 'blockB' = <blockB.#>, etc... 'treefile = '#treefile.trees', 'counter' = #, 'hostlist' = ['species'], 'taxalist' = ['accession|blah|host|....']}
def parsexml(filename, blocklist, counter):
	xmldict = {}

	for block in blocklist:
	#first retrieve the block with getblock; then rename its 'id' tag with renameblock
		oldblock = getblock(file, block)
		newblock = renameblock(oldblock, 'id', counter)
	#add it to our dictionary for this file
		xmldict[block] = newblock
		
	#also make an entry for the name of the tree file and for the file id number (counter)
	xmldict['treefile'] = filename.split('.')[0]+'_resampled.trees'
	xmldict['counter'] = str(counter)
	#figure out which hosts and accession numbers were represented in this file.
	hostlist = []
	taxalist = []
	alignblock = getblock(file, 'alignment')
	
	for child in alignblock:					#parse host species
		description = child[0].get('idref')
		if description == 'pSIVgml':
			continue
		else:
			host = description.split('|')[2]
			if host in hostlist:
				continue
			else:
				hostlist.append(host)
				
	for child in alignblock:					#parse taxa
		description = child[0].get('idref')
		if description in taxalist:
			continue
		else:		
			taxalist.append(description)
	xmldict['hostlist'] = hostlist
	xmldict['taxalist'] = taxalist				#we'll deal with the pSIVgml outgroup and any taxa without a host tag later.
	return xmldict



##############Parse the input, non-master files for the independent blocks ###############################

comboList = []									#make a list to hold the dictionaries for each of the non-master files.
independentBlocks = ['taxa']					#which blocks do we want to snag from each non-master file? 
counter = 0										#how many files have we already seen? use this to keep each file separate; assigned in alphabetical order by file name.
for file in xmllist:							#parse each non-master .xml file, creating a dictionary of filename = {blockA = <blockA.#>, blockB = <block2.#>, etc.}
	filedict = parsexml(file, independentBlocks, counter)
	counter += 1								#add it to our list of file dictionaries to combine, and up our counter to mark that we've seen it.
	comboList.append(filedict)

############# Make blocks for the final xml ######################################################################		
#We now have a master template and a data structure for the files we want to combine. 
# combolist = [{'blockA' = <blockA.1>, 'blockB' = <blockB.1>, 'treefile' = 'file1.trees', 'counter' = 1 ...}, {'blockA' = <blockA.2>, ...}, ... ]

for dict in comboList:
#make the empiricalTreeDistributionModel block and taxa subblock
	treeModID = 'treeModel.'+ dict['counter']
	treeModel = ET.Element('empiricalTreeDistributionModel')
	treeModel.set('id', treeModID)
	treeModel.set('fileName', dict['treefile'])
	addTaxa = ET.SubElement(treeModel, 'taxa')
	addTaxa.set('idref', 'taxa')
	dict['treeModel'] = treeModel
#make the accompanying current tree statistic block and treemodel subblock
	currentTree = ET.Element('statistic')
	currentTree.set('id', 'treeModel.'+dict['counter']+'.currentTree')
	currentTree.set('name', 'Current Tree')
	addTrees = ET.SubElement(currentTree, 'empiricalTreeDistributionModel')
	addTrees.set('idref', treeModID)
	dict['currentTree'] = currentTree
	
#make the ancestralTreeLikelihood block
	ancesTree = ET.Element('ancestralTreeLikelihood')
	ancesTree.set('id', 'host.treeLikelihood.'+dict['counter'])
	ancesTree.set('stateTagName', 'host.states')
	attribpatt = ET.SubElement(ancesTree, 'attributePatterns')
	attribpatt.set('idref', 'host.pattern')
	addTrees2 = ET.SubElement(ancesTree, 'treeModel')
	addTrees2.set('idref', treeModID)
	addSites = ET.SubElement(ancesTree, 'siteModel')
	addSites.set('idref', 'host.siteModel')
	addSubModel = ET.SubElement(ancesTree, 'generalSubstitutionModel')
	addSubModel.set('idref', 'host.model')
	dict['ancesTree'] = ancesTree

#make the logTree block
	logtree = ET.Element('logTree')
	logtree.set('fileName', 'discreteTraits_%s_%s.trees'%(dict['counter'], dateTag))	#unique output .trees files, tagged with date.
	logtree.set('id', 'treeFileLog.%s'%dict['counter'])
	logtree.set('logEvery', '5000')
	logtree.set('nexusFormat', 'true')
	logtree.set('sortTranslationTable', 'true')
	addTreeModel1 = ET.SubElement(logtree, 'treeModel')
	addTreeModel1.set('idref', treeModID)
	addpost = ET.SubElement(logtree, 'posterior')
	addpost.set('idref', 'posterior')
	addtraits = ET.SubElement(logtree, 'trait')
	addtraits.set('name', 'host.states')
	addtraits.set('tag', 'host')
	addancestree = ET.SubElement(addtraits, 'ancestralTreeLikelihood')
	addancestree.set('idref', 'host.treeLikelihood.'+dict['counter'])
	dict['logtree'] = logtree
	
#make the empiricalTreeDistribution / treemodel operator
	treeModelOp = ET.Element('empiricalTreeDistributionOperator')
	treeModelOp.set('weight', '1')
	addTreeModel = ET.SubElement(treeModelOp, 'empiricalTreeDistributionModel')
	addTreeModel.set('idref', treeModID)
	dict['treeModelOp'] = treeModelOp
	
#make the mcmc ancestral tree likelihood block
	likelihood = ET.Element('likelihood')
	likelihood.set('id', 'likelihood.'+dict['counter'])
	addAncesTree = ET.SubElement(likelihood, 'ancestralTreeLikelihood')
	addAncesTree.set('idref', 'host.treeLikelihood.'+dict['counter'])
	dict['likelihood'] = likelihood 
	
#compile a list of all the taxa and hosts we've seen 
masterhostlist = []
for dict in comboList:	
	for host in dict['hostlist']:
		if host in masterhostlist:
			continue
		else:
			masterhostlist.append(host)
			
mastertaxalist = []
taxacount = 0
for dict in comboList:
	for taxa in dict['taxalist']:
		if taxa in mastertaxalist:
			continue
		else:
			mastertaxalist.append(taxa)
			taxacount += 1

#make the shared host classes block
hostClasses = ET.Element('generalDataType')
hostClasses.set('id', 'host.dataType')

n=0 														#number of taxa without a specified host
for host in masterhostlist:
	if host != ' ' and host != '':
		newhost = ET.SubElement(hostClasses, 'state')		#add a state code for each new host species we encounter
		newhost.set('code', host)
	else:
		n+=1

#find the number of unique host species and calculate some dimension parameters.
hostCount = len(masterhostlist) - n	
rateDimension = (hostCount*(hostCount-1))

#check we made the right number of host state codes
assert len(hostClasses) == hostCount	

#make the shared taxa block.
taxaShared = ET.Element('taxa')
taxaShared.set('id', 'taxa')
for taxa in mastertaxalist:
	if (taxa != 'pSIVgml') and (taxa.split('|')[2] != '' and (taxa.split('|')[2] != ' ')):
		newTaxa = ET.SubElement(taxaShared, 'taxon')
		newTaxaName = taxa.replace("COTED'IVOIRE", "COTEDIVOIRE").replace("M29973|SE_id250619|Cercopithecus_aethiops_aethiops|sub_given|GRV|GRI_677_gri_1_677, Grived (gr-1)|len2438|AFRSSA|ETHIOPIA||", "M29973|SE_id250619|Cercopithecus_aethiops_aethiops|sub_given|GRV|GRI_677_gri_1_677_Grived_(gr-1)|len2438|AFRSSA|ETHIOPIA||")
		newTaxa.set('id', newTaxaName)
		addHost = ET.SubElement(newTaxa, 'attr')
		addHost.set('name', 'host')
		addHost.text = taxa.split('|')[2]						#set each taxa's host state
	else:
		newTaxa = ET.SubElement(taxaShared, 'taxon')
		newTaxa.set('id', taxa)
		addHost = ET.SubElement(newTaxa, 'attr')
		addHost.set('name', 'host')
		addHost.text = '?'										#if no host state given, assign '?' as state code.



############################### we've got everything we need to start plugging things in. #######################################
############################### finally, parse the mastertemplate so we have somewhere to put things. ###########################

masterfile = open(sys.argv[1], 'r')
master = ET.parse(masterfile)
masterfile.close()
root = master.getroot()

#isolate some specific blocks we'll need later
operatorElement = root.find('operators')
mcmcElement = root.find('mcmc')
posteriorElement = mcmcElement.find('posterior')
mcmcLikelihoodElement = posteriorElement.find('likelihood')
fileLogElement = mcmcElement.findall('log')[1]

#update our .log and .ops file names so they're tagged with the date.
fileLogElement.set('fileName', "discreteTraits_%s.log"%dateTag) 
mcmcElement.set('operatorAnalysis', "discreteTraits_%s.ops"%dateTag)

#make sure we've given the right dimension parameters for the matrix, rates, and rateindicators
matrixdimensions = root.find('generalSubstitutionModel')[1][0][1][0]						#host.frequencies
matrixdimensions.set('dimension', str(hostCount))
ratedim = root.find('generalSubstitutionModel')[2][0]										#host.rates
ratedim.set('dimension', str(rateDimension))
indicatordim = root.find('generalSubstitutionModel')[3][0]									#host.indicators
indicatordim.set('dimension', str(rateDimension))

#keep track of where we need to insert things in an iterative way
#takes element object, name of the block that you want to insert things *in front of*
#returns appropriate index to insert things just before the indicated block
def index(element, blockname):
	return list(element).index(element.find(blockname)) 

#insert our newly constructed shared blocks (hostClasses and taxa) at the root
root.insert(0, taxaShared)
root.insert(0, hostClasses)

#insert blocks in the appropriate locations
for dict in comboList:
	root.insert(index(root, 'generalSubstitutionModel'), dict['treeModel'])					#empiricalTreesDistribution
	root.insert(index(root, 'generalSubstitutionModel'), dict['currentTree'])				#currentTree statistic
	
	root.insert(index(root, 'operators'), dict['ancesTree'])								#ancestralTreeLikelihood operator
	operatorElement.insert(index(operatorElement, 'scaleOperator'), dict['treeModelOp'])	#empiricalTreeDistribution/treeModel operator
	mcmcLikelihoodElement.append(dict['likelihood'][0])										#ancestralTreeLikelihood - mcmc
	fileLogElement.insert(index(fileLogElement, 'posterior'), dict['likelihood'][0])		#ancestralTreeLikelihood - file log
	mcmcElement.insert(index(mcmcElement, 'log'), dict['logtree'])							#logTree - mcmc
	
#write it to a new clean file. deal with the dumb pretty printing issue (why, ET, why?)
def indent(elem, level=0):
    i = "\n" + level*"    "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

outname = 'discreteTraits_'+dateTag+'.xml'
indent(root)
output = open(outname, 'w')
master.write(output)
output.close()
