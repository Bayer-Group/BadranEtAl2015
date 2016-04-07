#!/usr/bin/env python
import argparse
from string import split
from Bio.Seq import Seq
import re

# TODO - export qualities

parser = argparse.ArgumentParser()
parser.add_argument('-s', dest='samfile', required=True, help='input SAM file', metavar='FILE')
parser.add_argument('-r', dest='reffile', required=True, help='reference sequence FASTA file', metavar='FILE')
parser.add_argument('-o', dest='outfile', required=True, help='output aligned FASTA file', metavar='FILE')
args = parser.parse_args()

refname = ''
refseq = ''
fr = open(args.reffile)
for line in fr:
	if line[0] == '>':
		refname += line[1:].rstrip()
	else:
		refseq += line.rstrip()
fr.close()
strpositions = range(1, len(refseq))

names = []
seqs = []

def insertgap(seq, pos):
	return seq[:pos] + '-' + seq[pos:]
def filltolength(seq, length):
	return seq + '-' * (length - len(seq))

fs = open(args.samfile)
for line in fs:
	if (line[0] != '@'):
		linea = split(line, '\t')
		FLAG = linea[1]
		POS = int(linea[3])
		SEQ = linea[9]
		if ((int(FLAG) & 0x4 != 0x4) & (int(FLAG) & 0x100 != 0x100) & (len(SEQ) + POS <= len(refseq))):
			names += [linea[0]]
			CIGAR = linea[5]
			seqstrpos = 0
			while (strpositions[seqstrpos] < POS):
				seqstrpos += 1
				SEQ = '-' + SEQ
			m = re.findall('(\d+)([MID])', CIGAR)
			for i in xrange(0, len(m)):
				if (m[i][1] == 'M'):
					endpos = strpositions[seqstrpos] + int(m[i][0])
					while ((strpositions[seqstrpos] < endpos) & (seqstrpos < len(strpositions))):
						if (refseq[seqstrpos] == '-'):
							SEQ = SEQ[:seqstrpos] + '-' + SEQ[seqstrpos:]
						seqstrpos += 1
				elif (m[i][1] == 'D'):
					endpos = strpositions[seqstrpos] + int(m[i][0])
					while ((strpositions[seqstrpos] < endpos) & (seqstrpos < len(strpositions))):
						SEQ = SEQ[:seqstrpos] + '-' + SEQ[seqstrpos:]
						seqstrpos += 1
						if (len(strpositions) < seqstrpos):
							strpositions += [strpositions[len(strpositions)-1]]
				elif (m[i][1] == 'I'):
					endpos = seqstrpos + int(m[i][0])
					while (seqstrpos < endpos):
						if (refseq[seqstrpos] != '-'):
							refseq = refseq[:seqstrpos] + '-' + refseq[seqstrpos:]
							strpositions = strpositions[:seqstrpos] + [strpositions[seqstrpos]] + strpositions[seqstrpos:]
							seqs = map(insertgap, seqs, (seqstrpos,)*len(seqs))
						seqstrpos += 1
						if (len(strpositions) < seqstrpos):
							strpositions += [strpositions[len(strpositions)-1]]
			seqs += [SEQ]
			#if (len(seqs) >= 100):
			#	break
fs.close()

seqs = map(filltolength, seqs, (len(refseq),)*len(seqs))

fo = open(args.outfile,'w')
fo.write('>' + refname + '\n' + refseq + '\n')
for i in xrange(0,len(names)):
	fo.write('>' + names[i] + '\n' + seqs[i] + '\n')
fo.close()

while (len(strpositions) < len(refseq)):
	strpositions += [strpositions[len(strpositions)-1]]

fop = open(args.outfile + '.positions', 'w')
for i in xrange(0,len(strpositions)):
	fop.write(str(i) + '\t' + str(strpositions[i]) + '\n')
fop.close()
