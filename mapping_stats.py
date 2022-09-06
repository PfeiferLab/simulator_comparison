import pysam
from sys import argv

input_file=argv[1]
golden_bam=argv[2]

inbam=pysam.AlignmentFile(input_file, "rb")                                     
gbam=pysam.AlignmentFile(golden_bam, "rb")

output_file=open(input_file.split(".")[0]+".mapstats.txt",'w')

num_mapped_only_to_right_chrom=0
num_mapped_perfectly=0
num_mapped_to_wrong_chrom=0
total=0.0
for read,gread in zip(inbam,gbam):
    print(str(read.query_name),str(read.reference_name),str(read.reference_start))
    if read.query_name == gread.query_name and read.is_read1 == gread.is_read1:

        if read.reference_name == gread.reference_name: 
            if read.reference_start == gread.reference_start:
                num_mapped_perfectly+=1
            else:
                num_mapped_only_to_right_chrom+=1
        else:
            num_mapped_to_wrong_chrom+=1
    total+=1
output_file.write("Mapped Perfectly:"+str(num_mapped_perfectly)+'\t'+str(num_mapped_perfectly/total)+'\n')
output_file.write("Mapped only to right chromosome:"+str(num_mapped_only_to_right_chrom)+'\t'+str(num_mapped_only_to_right_chrom/total)+'\n')
output_file.write("Mismapped:"+str(num_mapped_to_wrong_chrom)+'\t'+str(num_mapped_to_wrong_chrom/total)+'\n')
output_file.write("Total:"+str(total))
