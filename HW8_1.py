import os

samples = []

for file in os.listdir("/mnt/c/Users/SMG/Desktop/Sequencing_class/HW8/Yeast_replicates_files"):
    if file.endswith(".fastq.gz"):
        samples.append(file)

for sample in samples:
	command = "kallisto quant -i transcripts.idx -o output%s --single -l 180 -s 20 %s" % (samples.index(sample), sample)
	print "Currently running: %s" % command
	os.system(command) 
