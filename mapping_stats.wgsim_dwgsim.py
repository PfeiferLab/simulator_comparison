import pysam                                                                    
from sys import argv                                                            
                                                                                
print_file=argv[1]                                                              
                                                                                
inbam=pysam.AlignmentFile(print_file, "rb")                                     
                                                                                
output_file=open(print_file.split(".")[0]+".mapstats.txt",'w')                  
                                                                                
num_mapped_only_to_right_chrom=0                                                
num_mapped_perfectly=0                                                          
num_mapped_to_wrong_chrom=0                                                     
total=0.0                                                                       
for read in inbam:                                              
    true_coords=read.query_name.split("_")
    true_chrom=true_coords[0]
    true_pos=true_coords[1]
    print("True Coords:"+str(true_chrom)+":"+str(true_pos))
    print("Mapped Coords:"+str(read.reference_name)+":"+str(read.reference_start+1))
    if read.reference_name == true_chrom:                         
        if read.reference_start+1 == int(true_pos):
            print("Perfect Map")                   
            num_mapped_perfectly+=1                                         
        else:
            print("Mapped to right chrom")                                                               
            num_mapped_only_to_right_chrom+=1                               
    else:
        print("Mismap")                                                                   
        num_mapped_to_wrong_chrom+=1                                        
    total+=1                                                                    
output_file.write("Mapped Perfectly:"+str(num_mapped_perfectly)+'\t'+str(num_mapped_perfectly/total)+'\n')
output_file.write("Mapped only to right chromosome:"+str(num_mapped_only_to_right_chrom)+'\t'+str(num_mapped_only_to_right_chrom/total)+'\n')
output_file.write("Mismapped:"+str(num_mapped_to_wrong_chrom)+'\t'+str(num_mapped_to_wrong_chrom/total)+'\n')
output_file.write("Total:"+str(total))
