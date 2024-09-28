import os
"""
Note: This command needs to be run on Linux, and the bedtools command is already in the environment variable; The file name should not have spaces
"""

# Set the bed folder path and output folder path
bed_folder = '/home/ycmei/model_demo/data/bed_class'
output_folder = '/home/ycmei/model_demo/data/output_bed_class'
#output_folder = '/home/ycmei/model_demo/data/output_bed_regress'

os.makedirs(output_folder, exist_ok=True)

for bed_file in sorted(os.listdir(bed_folder), key=lambda x: int(x.split('.')[0])):
    if bed_file.endswith('.bed'):
        input_bed = os.path.join(bed_folder,bed_file)
        
        output_bed = os.path.join(output_folder, f'{os.path.splitext(bed_file)[0]}.bed')
        
        command = f"bedtools intersect -wa -f 0.3 -c -a /home/ycmei/model_demo/data/whole_genome_200bp.bed -b {input_bed} > {output_bed}"
        #command = f"bedtools intersect -wa -wb -f 0.3 -a /home/ycmei/MeiProject_09_bedtools/datasets/whole_genome_200bp.bed -b {input_bed} > {output_bed}"
        os.system(command)

