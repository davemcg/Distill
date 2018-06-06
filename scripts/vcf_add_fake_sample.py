#!/usr/local/Anaconda/envs/py3.5/bin/python

import fileinput
import sys
import gzip

sample_name = 'FAUX_GNOMAD'

for line in fileinput.input():
	#line = line.decode('utf-8')
	if line[0:2] == '##':
		print(line[:-1])
	elif line[0:2] == '#C':
		print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
		print('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">')
		print('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">')
		print(line[:-1] + '\tFORMAT\t' + sample_name)
	else:
		line = line.split()
		line = line[:-1]
		if line[0] == 'chrM':
			line[0] = 'MT'
		if line[0][0:3] == 'chr':
			line[0] = line[0][3:]
		print('\t'.join(line) + '\t.\tGT:GQ:DP\t0/1:100:100')
