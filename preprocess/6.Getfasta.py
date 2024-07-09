import os
     
command = "bedtools getfasta -fi /home/ycmei/data/GRCh38.p14.genome.fa -bed union.peaks.bed -fo union.peaks.pad1k.fa"
        
os.system(command)