######################		Initialize		############################

import elementtree.ElementTree as ET
from Bio import SeqIO
from datetime import date
import sys
from glob import glob
import copy

dateTag = str(date.today())					#	Keep track of when we created this xml

alignfiles = glob('*.fasta')
alignfiles.sort(key = lambda x: int(x.split('_')[-2]))

alignments = { file:SeqIO.to_dict(SeqIO.parse(file, 'fasta'), key_function = lambda rec : str(rec.description)) for file in alignfiles }

######################		Utilities		############################

def index(element, blockname):				#	Returns index of blockname in an element. Used below for appropriate insertion indices. 
	return list(element).index(element.find(blockname)) 

def set_attribs(element, listOfTuples):		#	Given [('add', 'attributes'), ...] & <xml_block />, return <xml_block add="attributes"/>
	for attrib in listOfTuples:
		element.set(attrib[0], attrib[1])
	return element

def indent(elem, level=0):					#	Make the output legible; thanks, stackoverflow.
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


###########		Make taxa and alignment blocks	###############
def make_alignment(seqobjectdict):	#	[seqobjects], N --> <alignment_block, id='N'> <taxon and sequence for each seqobject /> </alignment_block>
	alignblock = ET.Element('alignment')
	alignblock.set('id', 'alignment')
	alignblock.set('dataType', 'nucleotide')

	for header, seqobject in seqobjectdict.items():
		seqblock = ET.SubElement(alignblock, 'sequence')
		taxon = ET.SubElement(seqblock, 'taxon')
		taxon.set('idref', header)
		taxon.tail= str(seqobject.seq)
	return alignblock


def make_taxa(seqobjectdict):
	taxablock = ET.Element('taxa')
	taxablock.set('id', 'taxa')
	
	for header, seqobject in seqobjectdict.items():
		taxon = ET.SubElement(taxablock, 'taxon')
		taxon.set('id', header)
	
	return taxablock
	
##########	Make xmls	#######################
	
for segment, sequences in alignments.items():
	filename = '%s_%d_empTrees'%(segment.split('_')[0], alignfiles.index(segment))
	
	masterfile = open(sys.argv[1], 'r')
	master = ET.parse(masterfile)	#	Parse the template
	masterfile.close()
	root = master.getroot()
	
	taxa = make_taxa(sequences)
	align = make_alignment(sequences)
	
	root.insert(0, align)
	root.insert(0, taxa)
	
	opslog = root.find('mcmc')
	opslog.set('operatorAnalysis', '%s.txt'%filename)
	
	logfile = opslog.find('log')
	logfile.set('fileName', '%s.log'%filename)
	
	treelog = opslog.find('logTree')
	treelog.set('fileName', '%s.trees'%filename)
	
	indent(root)
	
	outfile = open('%s.xml'%filename, 'w')
	master.write(outfile)
	outfile.close()
		
		
		
		
		
		
		
		
