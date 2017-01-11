import pandas as pd
import xml.etree.ElementTree as ET

def find_bf(indicator_avg, prior_expectation, n_demes, len_chain=10000000):
	priorProbabilityNumerator = prior_expectation
	priorProbabilityDenominator = (n_demes**2 - n_demes)			#total number of parameters estimated that could be on
	priorProbability = float(priorProbabilityNumerator)/float(priorProbabilityDenominator)
	priorOdds = float(priorProbability) / float(1-priorProbability)
	posteriorOdds=(((indicator_avg-(1/float(len_chain)))/float((1-(indicator_avg-(1/float(len_chain)))))))	#calculate posterior odds from the average bf values, adjust slightly to avoid values of 1
	bf = float(posteriorOdds) / float(priorOdds)
	return bf


def fill_counts_matrix(df, series):
	'''
	Given an empty data frame and series of values, iterates per BEAST 1.x jump count matrices (and all BEAST 2.x matrices)
	[horizontally, skip diagonal]
	to fill and return the dataframe/matrix.
	'''
	assert len(df.columns) == len(df.index)			#did we give it an appropriately sized df and number of parameters?
	assert len(series) == len(df.columns)**2
	n_hosts = len(df.columns)
	for index, row in df.iterrows():
		row[:n_hosts] = series[:n_hosts]
		series = series[n_hosts:]
	assert len(series) == 0							#we should have used all the values now; check this.
	return df

def fill_rates_matrix(df, series):
	'''
	Given an empty data frame and series of values, iterates per BEAST host rates (and GLM predictor value) matrices
	[horizontally right of the diagonal; vertically left of the diagonal]
	to fill and return the dataframe/matrix. Accommodates both (a)symmetric models.
	'''

	#did we give it an appropriately sized df?
	assert len(df.columns) == len(df.index)
	n_hosts = len(df.columns)
	n_params = len(series)
	#do we have a valid number of parameters?
	if (n_params == ((n_hosts*n_hosts) - n_hosts) / 2):		#if we only have half a matrix, just use the same data to fill both halves.
		series = pd.concat([series, series])
	elif (n_params == ((n_hosts*n_hosts) - n_hosts)):
		pass
	else:
		print 'ERROR: invalid number of parameters'
		assert (n_params == ((n_hosts*n_hosts) - n_hosts)) or (n_params == ((n_hosts*n_hosts) - n_hosts) / 2), (n_hosts, n_params)

	skip_to = 1									#the [0,0] index of the df is on the diagonal; skip it, and fill to the end fo the row.
	for index, row in df.iterrows():
		if skip_to == len(row):					#unless we've run off the end of the matrix....
			break
		else:
			fill_length = len(row) - skip_to	#we need to fill from the skip_to starting point to the end of the row.
			row[skip_to:len(row)] = series[0:fill_length]	#fill those in from the series
			series = series[fill_length:]		#trim the already-used values off the series.
			skip_to += 1						#bump over one more when we start the next row.


	skip_to = 1									#the [0,0] index of the df is on the diagonal; skip it, and fill to the end fo the row.
	for index, column in df.iteritems():
		if skip_to == len(column):					#unless we've run off the end of the matrix....
			break
		else:
			fill_length = len(column) - skip_to	#we need to fill from the skip_to starting point to the end of the row.
			column[skip_to:len(column)] = series[0:fill_length]	#fill those in from the series
			series = series[fill_length:]		#trim the already-used values off the series.
			skip_to += 1						#bump over one more when we start the next row.
	assert len(series) == 0						#there's something wrong if we still have values left.
	return df

def make_parameter_block(id, value):
	'''
	return ET.Element like <parameter id='id', value='value' />
	'''
	block = ET.Element('parameter')
	block.set('id', id)
	block.set('value', value)
	return block

def get_predictor_vectors(predictors, hostorder):
	'''
	Given a list of predictors, generate vectors for the GLM xml's designMatrix block.
	N.B.: Requires predictor.csv in cwd; for 'hosts', returns pairwise vectors (all 0 except for that 1 at that From/To host cell)
	Returns ET.Element objects like [  <parameter id='predictor', value='000...1...0' /> ]
	'''
	predictor_blocks = []
	hostCount = len(hostorder)
	nparams = hostCount**2 - hostCount
	templatematrix = pd.DataFrame(columns=hostorder, index=hostorder)
	templatematrix = fill_rates_matrix(templatematrix, range(nparams))

	for predictor in predictors:
		if predictor == 'hosts':
			for From in hostorder:	#	From each host...
				for To in hostorder:#	...To each other host
					if From == To:
						continue
					else:
						vector = ['0' for k in range(nparams)]	#	Vector of nparams 0s, insert a 1 in the appropriate matrix index.
						index = templatematrix.at[From, To]
						vector[index] = '1'
						assert len(vector) == nparams			#	Double check we didn't FUBAR this
						assert vector.count('0') == nparams -1
						assert vector.count('1') == 1
						predictor_blocks.append( make_parameter_block('From_%s_To_%s_glm'%(From, To) , ' '.join(vector)))#	Add a `_glm` suffix to the id to avoid parameter id conflicts with jump counts
		else:
			data = pd.DataFrame.from_csv('%s.csv'%predictor)	#heavily relies on having the data in files named predictor.csv and formatted as rows==From_hosts, columns==To_hosts
			vector = ['0' for k in range(hostCount**2 - hostCount)]
			for From in hostorder:
				for To in hostorder:
					if From != To:							#skip the diagonal
						vector[templatematrix.at[From, To]]=(str(data.at[From, To]))	#pull data by host name
			predictor_blocks.append( make_parameter_block(predictor, ' '.join(vector)) )
	return predictor_blocks


def get_jump_vectors(hostorder):
	'''
	Similar to get_predictor_vectors, but operates per jump counts indexing.
	Returns [ ET.Elements <parameter id='From_x_To_y', value='000...1...0' />
	'''
	#make the jump vectors for the markovJumpsTreeLikelihood blocks
	hostCount = len(hostorder)
	nparams = hostCount**2 #- hostCount 	N.B.: jump count vectors do not skip the diagonal.

	templatematrix = pd.DataFrame(columns=hostorder, index=hostorder)
	templatematrix = fill_counts_matrix(templatematrix, range(nparams))

	jump_vectors = []

	for From in hostorder:							#from each host....
		for To in hostorder:						#to each other host.....
			if From==To:									#[except to itself]
				continue
			else:
				vector = ['0' for k in range(hostCount**2)]	#make a vector of 0s ('not recording value'), with length = (number of hosts)^2
				index = templatematrix.at[From, To]
				vector[index] = '1'							#mark which value we're interested in with a 1 ('recording this value'); pull this index from the templatematrix above.
				vector = ' '.join(vector)
				block = make_parameter_block('From_%s_To_%s'%(From, To), vector)
				jump_vectors.append(block)
	return jump_vectors
