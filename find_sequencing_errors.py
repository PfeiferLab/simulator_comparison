import pysam
from sys import argv

print_bam=argv[1]
#print_fa=argv[2]
inbam=pysam.AlignmentFile(print_bam, "rb")
infa=pysam.FastaFile("/scratch/mmilhave/BIO514/ref/GCA_000146045.2_R64_genomic.fna")

output_file_forward=open(print_bam.split(".")[0]+".errors.forward.txt",'w')
output_file_reverse=open(print_bam.split(".")[0]+".errors.reverse.txt",'w')
output_file_forward.write("NAME\tNUM_SUBS\tSUB_POS\tSUB_QUALITIES\tNUM_INSERTIONS\tINSERTION_POS\tINSERTION_LENGTHS\tINSERTION_QUALITIES\tNUM_DELETIONS\tDELETION_POS\tDELETION_LENGTHS\tIS_R1\n")
output_file_reverse.write("NAME\tNUM_SUBS\tSUB_POS\tSUB_QUALITIES\tNUM_INSERTIONS\tINSERTION_POS\tINSERTION_LENGTHS\tINSERTION_QUALITIES\tNUM_DELETIONS\tDELETION_POS\tDELETION_LENGTHS\tIS_R1\n")
total=0
hasError=False
for read in inbam:
    print(read.query_name,read.reference_name,read.reference_start)
    cigartuples=read.cigartuples
    pos = 0
    offset = 0
    seq_errors = []
    seq_qualities = []
    insertions = []
    insertion_lengths = []
    insertion_qualities = []
    deletions = []
    deletion_lengths = []
    isIndel = False
    hasIndel = False
    for tup in cigartuples:
        ref_seq=[]
        read_seq=[]
        if tup[0] == 0:
        #    continue
            print("PERFECT MATCH:"+str(tup))
        elif tup[0] == 1:
            print("INSERTION:"+str(tup))
            offset += -tup[1]
            insertions.append(pos+1)
            insertion_lengths.append(tup[1])
            insertion_qualities.append(read.query_qualities[pos:pos+tup[1]])
            isIndel = True
            hasIndel = True
        elif tup[0] == 2:
            print("DELETION:"+str(tup))
            offset += tup[1] 
            deletions.append(pos+1)
            deletion_lengths.append(tup[1])
            isIndel = True
            hasIndel = True
        else:
            #print("PROBLEM")
            print("PROBLEM:"+str(tup))
        #elif tup[0] == 3:
        #elif tup[0] == 4:
        #elif tup[0] == 5:
        #elif tup[0] == 6:
        #elif tup[0] == 7:
        #elif tup[0] == 8:
        #elif tup[0] == 9:
        if not isIndel:
            ref_seq = infa.fetch(read.reference_name,read.reference_start+pos + offset, read.reference_start+pos + tup[1] + offset)
            read_seq = read.query_sequence[pos:pos + tup[1]]
            print(ref_seq.upper())
            print(read_seq)
            print("About to check errors in region...")
            for ref_nuc, read_nuc in zip(ref_seq,read_seq):
                ref_nuc = ref_nuc.upper()
                read_nuc = read_nuc.upper()
                if ref_nuc != read_nuc:
                    print("ERROR AT POS:"+str(pos+1))
                    seq_errors.append(pos+1)
                    seq_qualities.append(read.query_qualities[pos])
                    print(seq_qualities)
                    hasError = True
                pos = pos + 1
        else:
            #print(insertion_qualities[0][0])
            print(read.query_qualities[0])
            print("Skipping Over Indel...")
            pos = pos + tup[1]
            isIndel = False

    if hasError or hasIndel:
        num_errors = len(seq_errors)
        num_insertions = len(insertions)
        num_deletions = len(deletions)
        insertion_line = ""
        insertion_lengths_line = ""
        insertion_qualities_line = ""
        deletion_line = ""
        deletion_lengths_line = ""
        seq_error_line = ""
        seq_qualities_line=""
        for item,item2 in zip(insertions,insertion_lengths):
            insertion_line = insertion_line + str(item) + ","
            insertion_lengths_line = insertion_lengths_line + str(item2) + ","
        if len(insertion_qualities) > 0:
            for item in insertion_qualities[0]:
                insertion_qualities_line = insertion_qualities_line + str(item) + ","
        for item,item2 in zip(deletions,deletion_lengths):
            deletion_line = deletion_line + str(item) + ","
            deletion_lengths_line = deletion_lengths_line + str(item2) + ","
        for item,item2 in zip(seq_errors,seq_qualities):
            seq_error_line = seq_error_line + str(item) + ","
            seq_qualities_line = seq_qualities_line + str(item2) + ","
        insertion_line = insertion_line[:-1]
        insertion_lengths_line = insertion_lengths_line[:-1]
        insertion_qualities_line = insertion_qualities_line[:-1]
        deletion_line = deletion_line[:-1]
        deletion_lengths_line = deletion_lengths_line[:-1]
        seq_error_line = seq_error_line[:-1]
        seq_qualities_line = seq_qualities_line[:-1]
        if seq_error_line == "":
            seq_error_line = "N/A"
            seq_qualities_line = "N/A"
        if insertion_line == "":
            insertion_line = "N/A"
            insertion_lengths_line = "N/A"
            insertion_qualities_line = "N/A"
        if deletion_line == "":
            deletion_line = "N/A"
            deletion_lengths_line = "N/A"
        print("NAME\tNUM_SUBS\tSUB_POS\tSUB_QUALITIES\tNUM_INSERTIONS\tINSERTION_POS\tINSERTION_LENGTHS\tINSERTION_QUALITIES\tNUM_DELETIONS\tDELETION_POS\tDELETION_LENGTHS\tIS_R1")
        print("OUTPUT_LINE:"+read.query_name+'\t'+str(num_errors)+'\t'+seq_error_line+'\t'+seq_qualities_line+'\t'+str(num_insertions)+'\t'+insertion_line+'\t'+insertion_lengths_line+'\t'+insertion_qualities_line+'\t'+str(num_deletions)+'\t'+deletion_line+'\t'+deletion_lengths_line+'\t'+str(read.is_read1))
        if read.is_reverse:
            output_file_reverse.write(read.query_name+'\t'+str(num_errors)+'\t'+seq_error_line+'\t'+seq_qualities_line+'\t'+str(num_insertions)+'\t'+insertion_line+'\t'+insertion_lengths_line+'\t'+insertion_qualities_line+'\t'+str(num_deletions)+'\t'+deletion_line+'\t'+deletion_lengths_line+'\t'+str(read.is_read1)+'\n')
        else:
            output_file_forward.write(read.query_name+'\t'+str(num_errors)+'\t'+seq_error_line+'\t'+seq_qualities_line+'\t'+str(num_insertions)+'\t'+insertion_line+'\t'+insertion_lengths_line+'\t'+insertion_qualities_line+'\t'+str(num_deletions)+'\t'+deletion_line+'\t'+deletion_lengths_line+'\t'+str(read.is_read1)+'\n')
        hasError=False
    #else:
    #    print("FLAWLESS")
    #total+=1
#    print("Total:"+str(total))
    #break
