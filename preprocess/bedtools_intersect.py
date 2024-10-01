import os
import sys

if len(sys.argv) < 4:
    print("Usage: python bedtools_intersect.py <task> <input_data_directory> <output_data_directory>")
    sys.exit(1)

task = sys.argv[1]
input_dir = sys.argv[2]
output_dir = sys.argv[3]

os.makedirs(output_dir, exist_ok=True)

for bed_file in sorted(os.listdir(input_dir), key=lambda x: int(x.split('.')[0])):
    if bed_file.endswith('.bed'):
        input_bed = os.path.join(input_dir,bed_file)

        if task == classification:
            output_bed = os.path.join(output_dir,"class", f'{os.path.splitext(bed_file)[0]}_classification.bed')
            command = f"bedtools intersect -wa -f 0.3 -c -a input_dir/whole_genome_200bp.bed -b {input_bed} > {output_bed}"
            
        if task == regression:
            output_bed = os.path.join(output_dir,"regress", f'{os.path.splitext(bed_file)[0]}_regression.bed')
            command = f"bedtools intersect -wa -wb -f 0.3 -a input_dir/whole_genome_200bp.bed -b {input_bed} > {output_bed}

        os.system(command)

