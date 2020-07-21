import sys
import gzip
import time

# samtools view input.bam | python bam_2_fq.py out_r1.fq out_r2.fq

OUT_R1  = sys.argv[1]
OUT_R2  = sys.argv[2]

def get_file_handle(fn, rw='r'):
	if fn[-6:].lower() == '.fq.gz' or fn[-9:].lower() == '.fastq.gz':
		return gzip.open(fn, rw)
	elif fn[-3:].lower() == '.fq' or fn[-6:].lower() == '.fastq':
		return open(fn, rw)
	else:
		print 'Input must be .fq, .fastq, .fq.gz or .fastq.gz'
		exit(1)

# return the reverse complement of a string
RC_DICT = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
def RC(s):
	return ''.join(RC_DICT[n] for n in s[::-1])

if not sys.stdin.isatty():
	input_stream = sys.stdin
else:
	print 'No input.'
	exit(1)

fo_1 = get_file_handle(OUT_R1, 'w')
fo_2 = get_file_handle(OUT_R2, 'w')

tt = time.time()
rDict = {}
unique_readnames = {}
nWritten = 0
for line in input_stream:
	if len(line) and line[0] != '#':
		splt  = line.strip().split('\t')
		rnm   = splt[0]
		unique_readnames[rnm] = True
		if rnm not in rDict:
			rDict[rnm] = [[],[]]
		flag  = int(splt[1])

		if flag&2048:	# skip supplementary alignments
			continue

		if flag&16:
			rdat = RC(splt[9])
			qdat = splt[10][::-1]
		else:
			rdat = splt[9]
			qdat = splt[10]

		if flag&1:
			if flag&64:
				rDict[rnm][0] = [rdat, qdat]
			elif flag&128:
				rDict[rnm][1] = [rdat, qdat]

		if len(rDict[rnm][0]) and len(rDict[rnm][1]):
			fo_1.write('@' + rnm + '/1\n')
			fo_1.write(rDict[rnm][0][0] + '\n')
			fo_1.write('+\n')
			fo_1.write(rDict[rnm][0][1] + '\n')
			fo_2.write('@' + rnm + '/2\n')
			fo_2.write(rDict[rnm][1][0] + '\n')
			fo_2.write('+\n')
			fo_2.write(rDict[rnm][1][1] + '\n')
			del rDict[rnm]
			nWritten += 1
			if nWritten%1000000 == 0:
				print nWritten, 'read pairs written', '('+str(int(time.time()-tt))+' sec)'
print nWritten, 'read pairs written', '('+str(int(time.time()-tt))+' sec)'

fo_2.close()
fo_1.close()
