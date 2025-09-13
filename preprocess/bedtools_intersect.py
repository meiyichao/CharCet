#!/usr/bin/env python3
"""
Enhanced BED processing script with parameterization and logging
Note: Requires Linux environment with bedtools in PATH
"""

import os
import sys
import time
import logging
import subprocess
import argparse
from pathlib import Path

# ===== 1. Initialization and Argument Parsing =====
start_time = time.time()
parser = argparse.ArgumentParser(description='Process BED files using bedtools intersect')
parser.add_argument('-b', '--bed_folder', required=True, 
                    help='Input directory containing BED files')
parser.add_argument('-o', '--output_folder', required=True,
                    help='Output directory for processed files')
parser.add_argument('-r', '--reference_bed', required=True,
                    help='Reference BED file for genome segmentation')
parser.add_argument('-f', '--overlap_fraction', type=float, default=0.3,
                    help='Minimum overlap fraction (default: 0.3)')
parser.add_argument('-t', '--task', choices=["classification", "regression"], 
                    default="classification",
                    help='Task type: classification (count) or regression (detailed overlaps)')
args = parser.parse_args()

# Convert to absolute paths
bed_folder = os.path.abspath(args.bed_folder)
output_folder = os.path.abspath(args.output_folder)
reference_bed = os.path.abspath(args.reference_bed)

# Create output directory
os.makedirs(output_folder, exist_ok=True)

# ===== 2. Logging System Configuration =====
log_file = os.path.join(output_folder, "bed_processing.log")
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

logger.info("===== BED File Processing Started =====")
logger.info(f"Input Directory: {bed_folder}")
logger.info(f"Output Directory: {output_folder}")
logger.info(f"Reference BED File: {reference_bed}")
logger.info(f"Overlap Threshold: {args.overlap_fraction}")
logger.info(f"Task Type: {args.task}")

# ===== 3. Input Validation =====
if not os.path.exists(reference_bed):
    logger.error(f"Reference BED file not found: {reference_bed}")
    sys.exit(1)

if not os.path.exists(bed_folder) or not os.listdir(bed_folder):
    logger.error(f"Input directory is empty or does not exist: {bed_folder}")
    sys.exit(1)

# ===== 4. Collect and Sort BED Files =====
try:
    # Collect BED files and sort numerically by filename prefix
    bed_files = sorted(
        [f for f in os.listdir(bed_folder) if f.endswith('.bed')],
        key=lambda x: int(x.split('.')[0])
    )
except Exception as e:
    logger.error(f"File sorting error: {str(e)}")
    sys.exit(1)

if not bed_files:
    logger.error("No BED files found in input directory")
    sys.exit(1)

logger.info(f"Found {len(bed_files)} BED files for processing")

# ===== 5. Process Each BED File =====
processed_count = 0
error_count = 0

for i, bed_file in enumerate(bed_files, 1):
    input_path = os.path.join(bed_folder, bed_file)
    output_path = os.path.join(output_folder, f"{os.path.splitext(bed_file)[0]}.bed")
    
    logger.info(f"Processing file {i}/{len(bed_files)}: {bed_file}")
    
    # Build bedtools command
    cmd = [
        "bedtools", "intersect",
        "-wa",  # Write original feature in A for each overlap
        "-f", str(args.overlap_fraction),  # Minimum overlap fraction
        "-a", reference_bed,  # Reference genome segments
        "-b", input_path     # Input BED file
    ]
    
    # Add task-specific parameter
    if args.task == "classification":
        cmd.append("-c")  # Count overlaps
    elif args.task == "regression":
        cmd.append("-wb")  # Write detailed overlap information
    
    logger.debug(f"Executing command: {' '.join(cmd)}")
    
    try:
        # Execute command and capture output
        with open(output_path, 'w') as outfile:
            result = subprocess.run(
                cmd,
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        
        # Verify output file content
        if os.path.getsize(output_path) == 0:
            logger.warning(f"Generated empty output file: {bed_file}")
        else:
            processed_count += 1
            logger.info(f"Successfully processed: {bed_file}")
            
    except subprocess.CalledProcessError as e:
        error_count += 1
        logger.error(f"Processing failed: {bed_file}")
        logger.error(f"Error details: {e.stderr.strip()}")
    except Exception as e:
        error_count += 1
        logger.error(f"Unexpected error processing: {bed_file}")
        logger.exception(e)

# ===== 6. Processing Statistics =====
logger.info("===== Processing Summary =====")
logger.info(f"Total Files: {len(bed_files)}")
logger.info(f"Successfully Processed: {processed_count}")
logger.info(f"Failed Files: {error_count}")

# ===== 7. Resource Cleanup and Timing =====
elapsed = time.time() - start_time
hours, remainder = divmod(elapsed, 3600)
minutes, seconds = divmod(remainder, 60)
logger.info(f"Processing completed! Total duration: {int(hours)}h {int(minutes)}m {seconds:.2f}s")

# Exit with appropriate status code
sys.exit(0 if error_count == 0 else 1)