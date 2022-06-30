#!/packages/7x/anaconda3/5.3.0/bin/python
import pysam
import numpy
from numpy.random import default_rng
import argparse

#setup argparse
parser = argparse.ArgumentParser(description="Downsamples a BAM file to a certain number of reads.") 
parser.add_argument("-t","--total-reads",dest="total_reads",help="The total number of reads within your BAM, often found with samtools flagstat")
parser.add_argument("-n","--num-to-downsample",dest="read_count",help="Number of reads you wish to have in your output BAM")
parser.add_argument("-i","--input-bam",dest="input_file",help="The undownsampled, input BAM file.")
parser.add_argument("-o","--output-bam",dest="output_file",help="Output BAM file.")
args = parser.parse_args()

total_reads = args.total_reads
read_count = args.read_count
input_file = args.input_file
output_file = args.output_file

counter=0
rng=default_rng()

lines_to_pull = rng.choice(int(total_reads),size=int(read_count),replace=False)
max_line = max(lines_to_pull)
lines_to_pull=numpy.sort(lines_to_pull)

inbam=pysam.AlignmentFile(input_file, "rb")
outbam=pysam.AlignmentFile(output_file, "wb",template=inbam)
i=0
print("Number of reads that will be in output BAM:"+str(read_count))
print("Downsample Proportion:"+str(float(read_count)/float(total_reads)))
for read in inbam:
    if((counter > max_line) or i == len(lines_to_pull)):
        break
    if(counter == lines_to_pull[i]):    
        outbam.write(read)
        i=i+1
    counter= counter+1  
