import sys

if not sys.stdin.isatty():
	input_stream = sys.stdin
else:
	print 'No input.'
	exit(1)

pDict = {}
currentRef = -1
prevRef    = None
for line in input_stream:
	splt = line.strip().split('\t')
	if splt[0] != prevRef:
		prevRef = splt[0]
		currentRef += 1
	myKey = tuple([currentRef] + sorted([int(splt[1]), int(splt[2])]) + [splt[0]])
	if myKey not in pDict:
		pDict[myKey] = []
	pDict[myKey].append(splt[3])

for k in sorted(pDict.keys()):
	sys.stdout.write('\t'.join([k[3], str(k[1]), str(k[2]), ';'.join(list(set(pDict[k])))]) + '\n')
