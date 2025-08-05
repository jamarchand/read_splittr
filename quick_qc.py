#!/usr/bin/env python3
"""
Quick QC Script

A simple script to quickly analyze FASTQ/FASTA files and provide basic statistics:
- Number of reads
- Average read length
- Format detection

Usage:
    python quick_qc.py input.fastq
    python quick_qc.py input.fasta
"""

import sys
from pathlib import Path
from Bio import SeqIO
import argparse


def detect_format(file_path):
    """Detect if file is FASTA or FASTQ format."""
    with open(file_path, 'r') as f:
        first_line = f.readline().strip()
        if first_line.startswith('@'):
            return 'fastq'
        elif first_line.startswith('>'):
            return 'fasta'
        else:
            raise ValueError(f"Unknown file format: {file_path}")


def analyze_file(file_path):
    """Analyze a FASTQ/FASTA file and return statistics."""
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    # Auto-detect format
    try:
        file_format = detect_format(file_path)
        print(f"Detected format: {file_format.upper()}")
    except ValueError as e:
        print(f"Error: {e}")
        return
    
    # Count reads and calculate lengths
    read_count = 0
    total_length = 0
    lengths = []
    
    try:
        for record in SeqIO.parse(file_path, file_format):
            read_count += 1
            read_length = len(record.seq)
            total_length += read_length
            lengths.append(read_length)
            
            # Print progress every 10000 reads
            if read_count % 10000 == 0:
                print(f"Processed {read_count:,} reads...")
    
    except Exception as e:
        print(f"Error reading file: {e}")
        return
    
    # Calculate statistics
    if read_count == 0:
        print("No reads found in file!")
        return
    
    avg_length = total_length / read_count
    min_length = min(lengths)
    max_length = max(lengths)
    
    # Print results
    print("\n" + "="*50)
    print("QUICK QC RESULTS")
    print("="*50)
    print(f"File: {file_path}")
    print(f"Format: {file_format.upper()}")
    print(f"Total reads: {read_count:,}")
    print(f"Total bases: {total_length:,}")
    print(f"Average length: {avg_length:.1f} bp")
    print(f"Min length: {min_length} bp")
    print(f"Max length: {max_length} bp")
    print("="*50)


def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Quick QC analysis of FASTQ/FASTA files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python quick_qc.py input.fastq
  python quick_qc.py input.fasta
  python quick_qc.py /path/to/sequencing_data.fastq
        """
    )
    
    parser.add_argument('input_file', help='Input FASTQ or FASTA file')
    parser.add_argument('--version', action='version', version='quick_qc.py 1.0')
    
    args = parser.parse_args()
    
    try:
        analyze_file(args.input_file)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 