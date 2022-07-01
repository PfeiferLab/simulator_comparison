import pysam
from sys import argv

input_file = pysam.AlignmentFile(argv[1], 'rb')

tlen_dist={}
negative_count = 0
positive_count = 0
zero_count = 0

for item in input_file:
    
    tlen=item.template_length
    if tlen < 0:
        negative_count +=1
    elif tlen > 0:
        print(str(tlen))
        positive_count+=1
        if str(abs(tlen)) in tlen_dist:
            tlen_dist[str(abs(tlen))] = tlen_dist[str(abs(tlen))] + 1
        else:
            tlen_dist[str(tlen)] = 1
