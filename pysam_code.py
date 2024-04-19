# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pysam
import os
import time
import pandas as pd

start_time = time.time()
PATH_BAM= "/home/bhavika/Desktop/Nanograd/Data/subset.bam"
PATH_SAM= "/home/bhavika/Desktop/Nanograd/Data/subset_1.sam"
BAM_ALIGNMENT = pysam.AlignmentFile(PATH_BAM, "rb")
HEADER_BAM = BAM_ALIGNMENT.header.copy()
READ_INDEXED = pysam.IndexedReads(BAM_ALIGNMENT, multiple_iterators= True)
READ_INDEXED.build()

data_degraded=[]

with open(PATH_SAM) as sam:
    for line in sam:
        read_name = line.strip().split("\t")[0]
        try:
            read_object = next(READ_INDEXED.find(read_name))
        except:
            print(read_name, "not in the bam")
        read_bam_name = read_name + ".bam"
        with pysam.Samfile(read_bam_name, "wb", header=HEADER_BAM) as read_bam_file:  # write a bam containing the read
            read_bam_file.write(read_object)
        pysam.index(read_bam_name)  # index the read bam
        read_alignment = pysam.AlignmentFile(read_bam_name, "rb")
        count = 0
        start = None
        end = None

        for pileupcol in read_alignment.pileup( min_base_quality=0,
                                               max_depth=9999999):  # go through reference positions
            ref_position = pileupcol.pos
            for pileupread in pileupcol.pileups:  # go through reads at that position (only one)
                if pileupread.query_position:     # if not a deletion at this position
                    if count == 0:
                        start = ref_position
                    else:
                        end = ref_position
                    count+=1
        if start and end:
            data_degraded.append(
                    {
                            'ReadName' : read_name,
                            'Start' : start,
                            'End' : end,
                            'ReadLength' : count})    # data frame output that includes read name, start, end, read_length
        os.remove(read_bam_name)
        os.remove(read_bam_name + ".bai")

df = pd.DataFrame(data_degraded)
df.to_csv("/home/bhavika/Desktop/Nanograd/Data/pysam_output.csv", index = True, header = True)



print("done in " , start_time - time.time())

