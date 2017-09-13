import sys
import itertools

removalfine = sys.argv[1]

removelist = set()
with open(removalfine) as f:
    for line in f:
        line = line.strip()
        removelist.add(line)

fastqcycle = itertools.cycle([1, 2, 3, 4])
fastq_holder = []
for line in sys.stdin:
    fastq_holder.append(line.strip())
    if fastqcycle.next() == 4:
        readname = fastq_holder[0].split(' ')[0][1:] # get only text before space and don't take @
        if readname.endswith('/1') or readname.endswith('/2'):
            readname = readname[:-2]
        if readname not in removelist:
            print '\n'.join(fastq_holder)
        else:
            pass
            #removelist.remove(readname)
        fastq_holder = []
