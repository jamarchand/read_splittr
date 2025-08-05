#!/usr/bin/env python3
"""
Nanopore Concatemer Splitter

This script processes nanopore sequencing data containing concatemers by identifying
primer sequences and splitting reads before the primer matches. It handles the lower
accuracy of nanopore reads by using appropriate BLAST parameters.

Author: Assistant
Date: 2024
"""

import os
import sys
import subprocess
import tempfile
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
import logging

# =============================================================================
# TUNABLE PARAMETERS - Modify these as needed
# =============================================================================
#
# HOW TO USE THESE PARAMETERS:
# 1. Modify the values below to adjust the behavior of the script
# 2. Command line arguments will override these defaults
# 3. For BLAST parameters, see: https://www.ncbi.nlm.nih.gov/books/NBK279675/
#
# BLAST_PARAMS: Control BLAST search sensitivity and specificity
# PARSING_PARAMS: Control which hits are accepted for splitting
# SPLITTING_PARAMS: Control how reads are split into fragments
# DEBUG_PARAMS: Control debugging output and file saving
#   - verbose: True = detailed output for each split, False = summary only
#     Use False for large files (100,000+ reads) to reduce output
#
# =============================================================================

# BLAST parameters optimized for short primer sequences
BLAST_PARAMS = {
    'identity': 80,  # 80% identity for nanopore reads
    'query_coverage': 90,  # 90% query coverage
    'evalue': 1e-5,  # More permissive for short sequences
    'word_size': 4,  # Smaller word size for short primers
    'dust': 'no',  # Disable dust filtering
    'soft_masking': 'false',
    'gapopen': 5,  # Lower gap opening penalty for short sequences
    'gapextend': 2,  # Lower gap extension penalty
    'penalty': -1,  # Mismatch penalty
    'reward': 1,  # Match reward
    'max_target_seqs': 1000,  # Allow multiple hits per query
    'task': 'blastn-short'  # Optimized for short sequences
}

# BLAST output parsing thresholds
PARSING_PARAMS = {
    'identity_threshold': 80,  # Minimum identity percentage for hits
    'min_primer_overlap': 0.9,  # Minimum fraction of primer that must overlap (0.8 = 80% of primer length)
}

# Splitting parameters
SPLITTING_PARAMS = {
    'min_fragment_length': 100,  # Minimum fragment length to keep
    'split_before_primer': True,  # Split before primer (True) or after primer (False)
}

# Debug and output parameters
DEBUG_PARAMS = {
    'save_blast_outputs': False,  # Save BLAST output files for debugging
    'save_round_outputs': True,  # Save intermediate files for each round
    'use_reversed_blast': False,  # Use reversed BLAST (reads as DB, primer as query)
    'verbose': False,  # Print detailed output for each split (False = summary only)
}

# =============================================================================
# END TUNABLE PARAMETERS
# =============================================================================

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class NanoporeConcatemerSplitter:
    def __init__(self, input_file, primer_file, output_file, blast_db_dir=None, save_blast_outputs=None, use_reversed_blast=None, verbose=None):
        """
        Initialize the splitter with input files and parameters.
        
        Args:
            input_file (str): Path to input FASTA/FASTQ file
            primer_file (str): Path to primer FASTA file
            output_file (str): Path to output file
            blast_db_dir (str): Directory for BLAST databases (optional)
            save_blast_outputs (bool): Save BLAST output files for debugging (overrides DEBUG_PARAMS)
            use_reversed_blast (bool): Use reversed BLAST (overrides DEBUG_PARAMS)
            verbose (bool): Print detailed output for each split (overrides DEBUG_PARAMS)
        """
        self.input_file = Path(input_file)
        self.primer_file = Path(primer_file)
        self.output_file = Path(output_file)
        self.blast_db_dir = Path(blast_db_dir) if blast_db_dir else Path.cwd()
        
        # Use provided parameters or defaults from DEBUG_PARAMS
        self.save_blast_outputs = save_blast_outputs if save_blast_outputs is not None else DEBUG_PARAMS['save_blast_outputs']
        self.use_reversed_blast = use_reversed_blast if use_reversed_blast is not None else DEBUG_PARAMS['use_reversed_blast']
        self.verbose = verbose if verbose is not None else DEBUG_PARAMS['verbose']
        
        # Use parameters from the top of the script
        self.blast_params = BLAST_PARAMS.copy()
        self.parsing_params = PARSING_PARAMS.copy()
        self.splitting_params = SPLITTING_PARAMS.copy()
        
        # Validate input files
        self._validate_inputs()
        
    def _validate_inputs(self):
        """Validate that input files exist and are readable."""
        if not self.input_file.exists():
            raise FileNotFoundError(f"Input file not found: {self.input_file}")
        if not self.primer_file.exists():
            raise FileNotFoundError(f"Primer file not found: {self.primer_file}")
        
        # Check if input is FASTA or FASTQ
        self.input_format = self._detect_format(self.input_file)
        logger.info(f"Detected input format: {self.input_format}")
        
    def _detect_format(self, file_path):
        """Detect if file is FASTA or FASTQ format."""
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('@'):
                return 'fastq'
            elif first_line.startswith('>'):
                return 'fasta'
            else:
                raise ValueError(f"Unknown file format: {file_path}")
    
    def _create_blast_db(self, primer_file, db_name):
        """Create a BLAST database from primer sequences using simple approach."""
        # Create the database in the blast_db_dir with just the base name
        db_base_path = self.blast_db_dir / db_name
        
        logger.info(f"Creating BLAST database: {db_base_path}")
        logger.info(f"Input primer file: {primer_file}")
        logger.info(f"Database directory: {self.blast_db_dir}")
        
        # Check if primer file exists and has content
        if not primer_file.exists():
            raise FileNotFoundError(f"Primer file not found: {primer_file}")
        
        # Read and log primer file content
        with open(primer_file, 'r') as f:
            primer_content = f.read()
            logger.info(f"Primer file content:\n{primer_content}")
        
        # Use simple approach like simple_blast.py
        cmd = [
            'makeblastdb',
            '-in', str(primer_file),
            '-dbtype', 'nucl',
            '-out', str(db_base_path)
        ]
        
        logger.info(f"BLAST command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"BLAST database creation STDOUT: {result.stdout}")
            logger.info(f"BLAST database creation STDERR: {result.stderr}")
            logger.info(f"Created BLAST database: {db_base_path}")
            
            # Check if database files were created
            db_files = list(self.blast_db_dir.glob(f"{db_name}.*"))
            logger.info(f"Database files created: {[f.name for f in db_files]}")
            
            return db_base_path
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to create BLAST database: {e}")
            logger.error(f"STDOUT: {e.stdout}")
            logger.error(f"STDERR: {e.stderr}")
            raise
    
    def _run_blast(self, query_file, db_path, output_file, primer_name=None):
        """Run BLAST search with simple approach (like simple_blast.py)."""
        
        # Check if query file exists and has content
        if not query_file.exists():
            raise FileNotFoundError(f"Query file not found: {query_file}")
        
        # Check query file size and content
        query_size = query_file.stat().st_size
        logger.info(f"Query file size: {query_size} bytes")
        
        # Read first few lines of query file
        with open(query_file, 'r') as f:
            first_lines = [f.readline().strip() for _ in range(6)]
            logger.info(f"Query file first lines:\n{chr(10).join(first_lines)}")
        
        # Check if database files exist
        db_files = list(self.blast_db_dir.glob(f"{db_path.name}.*"))
        logger.info(f"Database files for {db_path.name}: {[f.name for f in db_files]}")
        
        # Build command with simple approach (like simple_blast.py)
        cmd = [
            'blastn',
            '-query', str(query_file),
            '-db', str(db_path),
            '-out', str(output_file),
            '-outfmt', '10 qseqid sseqid pident qcovs qstart qend sstart send evalue bitscore sstrand'  # Added sstrand for strand info
        ]
        
        logger.info(f"BLAST command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            logger.info(f"BLAST search STDOUT: {result.stdout}")
            logger.info(f"BLAST search STDERR: {result.stderr}")
            logger.info(f"BLAST search completed: {output_file}")
            
            # Check output file size
            if output_file.exists():
                output_size = output_file.stat().st_size
                logger.info(f"BLAST output file size: {output_size} bytes")
                
                if output_size == 0:
                    logger.warning("BLAST output file is empty - no hits found or error occurred")
                    # Read first few lines to see if there's any content
                    with open(output_file, 'r') as f:
                        content = f.read(100)
                        logger.info(f"BLAST output content (first 100 chars): '{content}'")
                else:
                    # Read first few lines of output
                    with open(output_file, 'r') as f:
                        first_lines = [f.readline().strip() for _ in range(5)]
                        logger.info(f"BLAST output first lines:\n{chr(10).join(first_lines)}")
            else:
                logger.error("BLAST output file was not created")
            
            # Save BLAST output to local directory if requested
            if self.save_blast_outputs and primer_name:
                # Create blast_outputs directory if it doesn't exist
                blast_outputs_dir = Path("blast_outputs")
                blast_outputs_dir.mkdir(exist_ok=True)
                
                # Copy BLAST output to local directory
                local_blast_output = blast_outputs_dir / f"blast_{primer_name}_output.txt"
                shutil.copy2(output_file, local_blast_output)
                logger.info(f"Saved BLAST output to: {local_blast_output}")
                
                # Also save the query file for reference
                local_query_file = blast_outputs_dir / f"blast_{primer_name}_query.fasta"
                shutil.copy2(query_file, local_query_file)
                logger.info(f"Saved query file to: {local_query_file}")
                
                # Create summary file with BLAST parameters
                summary_file = blast_outputs_dir / f"blast_{primer_name}_summary.txt"
                with open(summary_file, 'w') as f:
                    f.write(f"BLAST Search Summary for Primer: {primer_name}\n")
                    f.write("=" * 50 + "\n\n")
                    f.write("BLAST Parameters:\n")
                    f.write(f"  Output format: 10 (CSV)\n")
                    f.write("\nCommand executed:\n")
                    f.write(f"  {' '.join(cmd)}\n\n")
                    f.write("Output format: CSV (format 10)\n")
                    f.write(f"\nQuery file size: {query_size} bytes\n")
                    f.write(f"Output file size: {output_file.stat().st_size if output_file.exists() else 0} bytes\n")
                    f.write(f"\nBLAST return code: {result.returncode}\n")
                    f.write(f"BLAST STDOUT: {result.stdout}\n")
                    f.write(f"BLAST STDERR: {result.stderr}\n")
                logger.info(f"Saved BLAST summary to: {summary_file}")
                
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST search failed: {e}")
            logger.error(f"STDOUT: {e.stdout}")
            logger.error(f"STDERR: {e.stderr}")
            raise

    def _run_blast_reversed(self, reads_file, primer_file, output_file, primer_name=None):
        """Run BLAST search with reversed query-subject order (reads as DB, primer as query)."""
        
        # Check if files exist and have content
        if not reads_file.exists():
            raise FileNotFoundError(f"Reads file not found: {reads_file}")
        if not primer_file.exists():
            raise FileNotFoundError(f"Primer file not found: {primer_file}")
        
        # Check file sizes
        reads_size = reads_file.stat().st_size
        primer_size = primer_file.stat().st_size
        logger.info(f"Reads file size: {reads_size} bytes")
        logger.info(f"Primer file size: {primer_size} bytes")
        
        # Read first few lines of files
        with open(reads_file, 'r') as f:
            first_lines = [f.readline().strip() for _ in range(4)]
            logger.info(f"Reads file first lines:\n{chr(10).join(first_lines)}")
        
        with open(primer_file, 'r') as f:
            primer_content = f.read()
            logger.info(f"Primer file content:\n{primer_content}")
        
        # Create database from reads using simple approach
        reads_db_path = reads_file.parent / f"{reads_file.stem}_db"
        db_cmd = [
            'makeblastdb',
            '-in', str(reads_file),
            '-dbtype', 'nucl',
            '-out', str(reads_db_path)
        ]
        
        logger.info(f"Creating reads database: {' '.join(db_cmd)}")
        try:
            db_result = subprocess.run(db_cmd, capture_output=True, text=True)
            logger.info(f"Database creation STDOUT: {db_result.stdout}")
            logger.info(f"Database creation STDERR: {db_result.stderr}")
            if db_result.returncode != 0:
                logger.error(f"Database creation failed: {db_result.stderr}")
                raise subprocess.CalledProcessError(db_result.returncode, db_cmd)
            logger.info("✓ Reads database created successfully")
            
            # Check if database files were created
            db_files = list(reads_file.parent.glob(f"{reads_db_path.name}.*"))
            logger.info(f"Database files created: {[f.name for f in db_files]}")
            
        except Exception as e:
            logger.error(f"Database creation error: {e}")
            raise
        
        # Run BLAST search with simple approach (like simple_blast.py)
        cmd = [
            'blastn',
            '-query', str(primer_file),
            '-db', str(reads_db_path),
            '-out', str(output_file),
            '-outfmt', '10 qseqid sseqid pident qcovs qstart qend sstart send evalue bitscore sstrand'  # Added sstrand for strand info
        ]
        
        logger.info(f"BLAST command (reversed): {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            logger.info(f"BLAST search STDOUT: {result.stdout}")
            logger.info(f"BLAST search STDERR: {result.stderr}")
            logger.info(f"BLAST return code: {result.returncode}")
            logger.info(f"BLAST search completed: {output_file}")
            
            # Check output file size
            if output_file.exists():
                output_size = output_file.stat().st_size
                logger.info(f"BLAST output file size: {output_size} bytes")
                
                if output_size == 0:
                    logger.warning("BLAST output file is empty - no hits found or error occurred")
                    # Read first few lines to see if there's any content
                    with open(output_file, 'r') as f:
                        content = f.read(100)
                        logger.info(f"BLAST output content (first 100 chars): '{content}'")
                else:
                    # Read first few lines of output
                    with open(output_file, 'r') as f:
                        first_lines = [f.readline().strip() for _ in range(5)]
                        logger.info(f"BLAST output first lines:\n{chr(10).join(first_lines)}")
            else:
                logger.error("BLAST output file was not created")
            
            # Save BLAST output to local directory if requested
            if self.save_blast_outputs and primer_name:
                # Create blast_outputs directory if it doesn't exist
                blast_outputs_dir = Path("blast_outputs")
                blast_outputs_dir.mkdir(exist_ok=True)
                
                # Copy BLAST output to local directory
                local_blast_output = blast_outputs_dir / f"blast_{primer_name}_output_reversed.txt"
                shutil.copy2(output_file, local_blast_output)
                logger.info(f"Saved BLAST output to: {local_blast_output}")
                
                # Also save the files for reference
                local_primer_file = blast_outputs_dir / f"blast_{primer_name}_primer.fasta"
                shutil.copy2(primer_file, local_primer_file)
                logger.info(f"Saved primer file to: {local_primer_file}")
                
                # Create summary file with BLAST parameters
                summary_file = blast_outputs_dir / f"blast_{primer_name}_summary_reversed.txt"
                with open(summary_file, 'w') as f:
                    f.write(f"BLAST Search Summary for Primer: {primer_name} (Reversed)\n")
                    f.write("=" * 60 + "\n\n")
                    f.write("BLAST Parameters:\n")
                    f.write(f"  Output format: 10 (CSV)\n")
                    f.write("\nCommand executed:\n")
                    f.write(f"  {' '.join(cmd)}\n\n")
                    f.write("Output format: CSV (format 10)\n")
                    f.write(f"\nReads file size: {reads_size} bytes\n")
                    f.write(f"Primer file size: {primer_size} bytes\n")
                    f.write(f"Output file size: {output_file.stat().st_size if output_file.exists() else 0} bytes\n")
                    f.write("\nNote: This uses REVERSED query-subject order (primer as query, reads as database)\n")
                    f.write(f"\nBLAST return code: {result.returncode}\n")
                    f.write(f"BLAST STDOUT: {result.stdout}\n")
                    f.write(f"BLAST STDERR: {result.stderr}\n")
                logger.info(f"Saved BLAST summary to: {summary_file}")
                
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST search failed: {e}")
            logger.error(f"STDOUT: {e.stdout}")
            logger.error(f"STDERR: {e.stderr}")
            raise
    
    def _calculate_primer_overlap(self, hit_start, hit_end, primer_length):
        """Calculate the fraction of primer that overlaps with the hit."""
        # Convert to 0-based coordinates if needed
        if hit_start > hit_end:
            hit_start, hit_end = hit_end, hit_start
        
        # Calculate overlap length
        overlap_length = hit_end - hit_start + 1
        
        # Calculate fraction of primer that overlaps
        overlap_fraction = overlap_length / primer_length
        
        return overlap_fraction, overlap_length
    
    def _parse_blast_results(self, blast_output, reversed_order=False):
        """Parse BLAST results and group by query sequence."""
        hits = {}
        
        if not os.path.exists(blast_output) or os.path.getsize(blast_output) == 0:
            logger.warning(f"BLAST output file is empty or doesn't exist: {blast_output}")
            return hits
        
        logger.info(f"Parsing BLAST results from: {blast_output}")
        logger.info(f"Reversed order: {reversed_order}")
        
        # Use parameters from the top of the script
        identity_threshold = self.parsing_params['identity_threshold']
        min_primer_overlap = self.parsing_params['min_primer_overlap']
        
        # Get primer lengths for overlap calculation
        primer_lengths = {}
        for record in SeqIO.parse(self.primer_file, 'fasta'):
            primer_lengths[record.id] = len(record.seq)
        
        logger.info(f"Primer lengths: {primer_lengths}")
        
        with open(blast_output, 'r') as f:
            line_count = 0
            for line in f:
                line_count += 1
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                
                # Try to detect format - check if it's CSV (format 10) or tabular (format 6)
                if ',' in line:
                    # CSV format (format 10)
                    fields = line.split(',')
                    logger.debug(f"Detected CSV format, {len(fields)} fields")
                else:
                    # Tabular format (format 6)
                    fields = line.split('\t')
                    logger.debug(f"Detected tabular format, {len(fields)} fields")
                
                if len(fields) >= 10:  # Now expecting at least 10 fields with strand info
                    if reversed_order:
                        # In reversed order: qseqid=primer, sseqid=read
                        qseqid = fields[0]  # primer name
                        sseqid = fields[1]  # read name
                        pident = float(fields[2])
                        qcovs = float(fields[3])
                        qstart = int(fields[4])  # primer start
                        qend = int(fields[5])    # primer end
                        sstart = int(fields[6])  # read start
                        send = int(fields[7])    # read end
                        evalue = float(fields[8])
                        bitscore = float(fields[9]) if len(fields) > 9 else 0.0
                        sstrand = fields[10] if len(fields) > 10 else 'plus'  # strand info
                        
                        logger.debug(f"Parsed line {line_count}: primer={qseqid}, read={sseqid}, identity={pident}%, coverage={qcovs}%, read_pos={sstart}-{send}, strand={sstrand}")
                        
                        # Calculate primer overlap
                        if sseqid in primer_lengths:
                            primer_length = primer_lengths[sseqid]
                            overlap_fraction, overlap_length = self._calculate_primer_overlap(qstart, qend, primer_length)
                            
                            # Filter by identity and primer overlap
                            if pident >= identity_threshold and overlap_fraction >= min_primer_overlap:
                                if qseqid not in hits:  # Group by read (query)
                                    hits[qseqid] = []
                                hits[qseqid].append({
                                    'primer': sseqid,
                                    'start': read_start,  # Use read coordinates (query)
                                    'end': read_end,
                                    'identity': pident,
                                    'overlap_fraction': overlap_fraction,
                                    'overlap_length': overlap_length,
                                    'evalue': evalue,
                                    'strand': sstrand,
                                    'is_forward': sstrand == 'plus'  # True if forward strand
                                })
                                if self.verbose:
                                    logger.debug(f"Added hit for read {qseqid} at position {read_start}-{read_end} (strand: {sstrand}, overlap: {overlap_fraction:.2f})")
                            else:
                                if self.verbose:
                                    logger.debug(f"Hit filtered out: identity={pident}% (threshold={identity_threshold}%), overlap={overlap_fraction:.2f} (threshold={min_primer_overlap})")
                        else:
                            logger.warning(f"Primer {sseqid} not found in primer file")
                    else:
                        # Original order: qseqid=read, sseqid=primer
                        qseqid = fields[0]  # read name
                        sseqid = fields[1]  # primer name
                        pident = float(fields[2])
                        qcovs = float(fields[3])
                        qstart = int(fields[4])  # query start (read start) - CORRECT!
                        qend = int(fields[5])    # query end (read end) - CORRECT!
                        sstart = int(fields[6])  # subject start (primer start)
                        send = int(fields[7])    # subject end (primer end)
                        evalue = float(fields[8])
                        sstrand = fields[10] if len(fields) > 10 else 'plus'  # strand info
                        
                        # Use READ coordinates for split positions (qstart, qend)
                        read_start = qstart  # Use query coordinates for read position
                        read_end = qend
                        
                        logger.debug(f"Parsed line {line_count}: read={qseqid}, primer={sseqid}, identity={pident}%, coverage={qcovs}%, read_pos={read_start}-{read_end}, strand={sstrand}")
                        
                        # Calculate primer overlap
                        if sseqid in primer_lengths:
                            primer_length = primer_lengths[sseqid]
                            overlap_fraction, overlap_length = self._calculate_primer_overlap(qstart, qend, primer_length)
                            
                            # Filter by identity and primer overlap
                            if pident >= identity_threshold and overlap_fraction >= min_primer_overlap:
                                if qseqid not in hits:  # Group by read (query)
                                    hits[qseqid] = []
                                hits[qseqid].append({
                                    'primer': sseqid,
                                    'start': read_start,  # Use read coordinates (query)
                                    'end': read_end,
                                    'identity': pident,
                                    'overlap_fraction': overlap_fraction,
                                    'overlap_length': overlap_length,
                                    'evalue': evalue,
                                    'strand': sstrand,
                                    'is_forward': sstrand == 'plus'  # True if forward strand
                                })
                                if self.verbose:
                                    logger.debug(f"Added hit for read {qseqid} at position {read_start}-{read_end} (strand: {sstrand}, overlap: {overlap_fraction:.2f})")
                            else:
                                if self.verbose:
                                    logger.debug(f"Hit filtered out: identity={pident}% (threshold={identity_threshold}%), overlap={overlap_fraction:.2f} (threshold={min_primer_overlap})")
                        else:
                            logger.warning(f"Primer {sseqid} not found in primer file")
                else:
                    logger.warning(f"Invalid line {line_count} in BLAST output: {line.strip()}")
        
        # Sort hits by position for each read
        for read_id in hits:
            hits[read_id].sort(key=lambda x: x['start'])
            if self.verbose:
                logger.info(f"Read {read_id} has {len(hits[read_id])} hits: {[(h['start'], h['end']) for h in hits[read_id]]}")
        
        logger.info(f"Total reads with hits: {len(hits)}")
        return hits
    
    def _split_read_at_primer(self, record, primer_hits, primer_name):
        """Split a read at primer positions and return fragments."""
        fragments = []
        seq = str(record.seq)
        read_len = len(seq)
        
        logger.debug(f"Splitting read {record.id} (length: {read_len}) with {len(primer_hits)} hits")
        
        # Get quality scores if available
        qual = None
        if hasattr(record, 'letter_annotations') and 'phred_quality' in record.letter_annotations:
            qual = record.letter_annotations['phred_quality']
        
        if not primer_hits:
            if self.verbose:
                logger.warning(f"No primer hits found for read {record.id}")
            return fragments
        
        # Sort hits by position to process them in order
        sorted_hits = sorted(primer_hits, key=lambda x: x['start'])
        if self.verbose:
            logger.debug(f"Sorted hits for {record.id}: {[(h['start'], h['end']) for h in sorted_hits]}")
        
        # Start with the full sequence
        current_seq = seq
        current_qual = qual
        current_offset = 0  # Track position offset as we split
        
        for i, hit in enumerate(sorted_hits):
            # Calculate split position relative to the current sequence
            # hit['start'] and hit['end'] are now the READ coordinates (not primer coordinates)
            if hit['is_forward']:
                # Forward primer hit: split BEFORE the primer
                split_pos = hit['start'] - current_offset
                if self.verbose:
                    logger.debug(f"Forward primer hit: splitting BEFORE primer at position {split_pos} (read pos {hit['start']}, offset {current_offset})")
            else:
                # Reverse primer hit: split AFTER the primer
                split_pos = hit['end'] - current_offset + 1
                if self.verbose:
                    logger.debug(f"Reverse primer hit: splitting AFTER primer at position {split_pos} (read pos {hit['end']}, offset {current_offset})")
            
            if self.verbose:
                logger.debug(f"Processing hit {i+1}/{len(sorted_hits)}: start={hit['start']}, end={hit['end']}, split_pos={split_pos}, strand={hit['strand']}, is_forward={hit['is_forward']}")
                logger.debug(f"Current sequence length: {len(current_seq)}, current_offset: {current_offset}")
            
            # Validate split position
            if split_pos < 0 or split_pos > len(current_seq):
                if self.verbose:
                    logger.warning(f"Invalid split position {split_pos} for current sequence (length: {len(current_seq)})")
                continue
            
            # Create left fragment (before split position)
            left_seq = current_seq[:split_pos]
            if len(left_seq) >= self.splitting_params['min_fragment_length']:
                left_qual = current_qual[:split_pos] if current_qual else None
                left_record = self._create_fragment_record(
                    record, left_seq, left_qual, primer_name, 0, hit['start']
                )
                fragments.append(left_record)
                if self.verbose:
                    logger.debug(f"Created left fragment: {left_record.id} (length: {len(left_seq)})")
            else:
                if self.verbose:
                    logger.debug(f"Left fragment too short (length: {len(left_seq)}, min: {self.splitting_params['min_fragment_length']}), skipping")
            
            # Update current sequence to be the right fragment (after split position)
            current_seq = current_seq[split_pos:]
            current_qual = current_qual[split_pos:] if current_qual else None
            current_offset = hit['start'] if hit['is_forward'] else hit['end'] + 1
            
            if self.verbose:
                logger.debug(f"Updated current sequence length: {len(current_seq)}")
        
        # Add the final fragment (everything after the last primer)
        if len(current_seq) >= self.splitting_params['min_fragment_length']:
            final_record = self._create_fragment_record(
                record, current_seq, current_qual, primer_name, 1, current_offset
            )
            fragments.append(final_record)
            if self.verbose:
                logger.debug(f"Created final fragment: {final_record.id} (length: {len(current_seq)})")
        else:
            if self.verbose:
                logger.debug(f"Final fragment too short (length: {len(current_seq)}, min: {self.splitting_params['min_fragment_length']}), skipping")
        
        if self.verbose:
            logger.info(f"Successfully split read {record.id} into {len(fragments)} fragments")
        logger.debug(f"Returning {len(fragments)} fragments for read {record.id}")
        return fragments
    
    def _create_fragment_record(self, original_record, seq_str, qual, primer_name, fragment_side, split_pos):
        """Create a new SeqRecord for a fragment."""
        # Create new ID with primer and order information
        start_pos = split_pos
        end_pos = split_pos + len(seq_str)
        
        # Include primer name and fragment order in the header
        if fragment_side == 0:
            # First fragment (before primer)
            new_id = f"{original_record.id}[{start_pos}:{end_pos}][{primer_name}]"
        else:
            # Subsequent fragments (after primer)
            new_id = f"{original_record.id}[{primer_name}][{start_pos}:{end_pos}]"
        
        # Create new record (always FASTA, no quality scores)
        new_record = SeqRecord(
            seq=Seq(seq_str),
            id=new_id,
            description=""  # Empty description to avoid duplicate information
        )
        
        return new_record
    
    def _count_reads(self, file_path):
        """Count the number of reads in a FASTA/FASTQ file."""
        count = 0
        # Detect format dynamically
        file_format = self._detect_format(file_path)
        for _ in SeqIO.parse(file_path, file_format):
            count += 1
        return count
    
    def process(self):
        """Main processing function."""
        logger.info("Starting nanopore concatemer splitting process")
        
        # Create local temp directory in current working directory
        temp_dir = Path("tmp_nanopore_splitter")
        temp_dir.mkdir(exist_ok=True)
        logger.info(f"Using local temp directory: {temp_dir}")
        
        # Create round_outputs directory for debugging
        round_outputs_dir = Path("round_outputs")
        if DEBUG_PARAMS['save_round_outputs']:
            round_outputs_dir.mkdir(exist_ok=True)
            logger.info(f"Round outputs will be saved to: {round_outputs_dir}")
        else:
            logger.info("Round outputs disabled")
        
        try:
            # Update blast_db_dir to use temp directory if not specified
            if self.blast_db_dir == Path.cwd():
                self.blast_db_dir = temp_dir
            
            # Copy input file to temp directory and convert to FASTA if needed
            current_input = temp_dir / f"current_input.fasta"
            if self.input_format == 'fastq':
                # Convert FASTQ to FASTA
                with open(current_input, 'w') as f:
                    for record in SeqIO.parse(self.input_file, 'fastq'):
                        SeqIO.write(record, f, 'fasta')
                logger.info(f"Converted FASTQ input to FASTA: {current_input}")
            else:
                # Copy FASTA file
                shutil.copy2(self.input_file, current_input)
                logger.info(f"Copied FASTA input: {current_input}")
            
            # Save initial input to round_outputs for debugging
            initial_round_file = round_outputs_dir / "round_0_initial.fasta"
            shutil.copy2(current_input, initial_round_file)
            logger.info(f"Saved initial reads to: {initial_round_file}")
            
            # Count initial reads
            initial_count = self._count_reads(current_input)
            logger.info(f"Initial read count: {initial_count}")
            
            # Process each primer
            primer_records = list(SeqIO.parse(self.primer_file, 'fasta'))
            total_splits = 0
            
            for i, primer_record in enumerate(primer_records):
                primer_name = primer_record.id
                primer_seq = str(primer_record.seq)
                
                logger.info(f"Processing primer {i+1}/{len(primer_records)}: {primer_name}")
                logger.info(f"Primer sequence: {primer_seq}")
                
                # Save current input before processing this primer
                before_primer_file = round_outputs_dir / f"round_{i+1}_before_{primer_name}.fasta"
                shutil.copy2(current_input, before_primer_file)
                logger.info(f"Saved reads before {primer_name} processing to: {before_primer_file}")
                
                # Create temporary primer file
                temp_primer_file = temp_dir / f"primer_{i}.fasta"
                with open(temp_primer_file, 'w') as f:
                    f.write(f">{primer_name}\n{primer_seq}\n")
                
                # Run BLAST search
                blast_output = temp_dir / f"blast_output_{i}.txt"
                if self.use_reversed_blast:
                    # Use reversed order (reads as DB, primer as query)
                    logger.info(f"Using reversed BLAST approach for {primer_name}")
                    self._run_blast_reversed(current_input, temp_primer_file, blast_output, primer_name)
                    hits = self._parse_blast_results(blast_output, reversed_order=True)
                else:
                    # Use original order (primer as DB, reads as query)
                    logger.info(f"Using original BLAST approach for {primer_name}")
                    db_name = f"primer_db_{i}"
                    db_path = self._create_blast_db(temp_primer_file, db_name)
                    self._run_blast(current_input, db_path, blast_output, primer_name)
                    hits = self._parse_blast_results(blast_output, reversed_order=False)
                
                # Count reads with hits
                reads_with_hits = len(hits)
                logger.info(f"Reads with hits to {primer_name}: {reads_with_hits}")
                
                # Save list of reads that will be split
                if reads_with_hits > 0:
                    reads_to_split_file = round_outputs_dir / f"round_{i+1}_reads_to_split_{primer_name}.txt"
                    with open(reads_to_split_file, 'w') as f:
                        f.write(f"Reads to be split by primer {primer_name}:\n")
                        f.write("=" * 50 + "\n")
                        for read_id in hits:
                            f.write(f"{read_id}\n")
                            for hit in hits[read_id]:
                                f.write(f"  Hit: {hit['start']}-{hit['end']}, identity={hit['identity']}%, overlap={hit['overlap_fraction']:.2f}\n")
                    logger.info(f"Saved list of reads to split to: {reads_to_split_file}")
                
                if reads_with_hits == 0:
                    logger.info(f"No hits found for primer {primer_name}, continuing...")
                    # Save unchanged file for this round
                    after_primer_file = round_outputs_dir / f"round_{i+1}_after_{primer_name}_no_hits.fasta"
                    shutil.copy2(current_input, after_primer_file)
                    logger.info(f"Saved unchanged reads to: {after_primer_file}")
                    continue
                
                # Split reads and create new file
                new_records = []
                splits_this_round = 0
                split_details = []
                
                logger.info(f"Processing {len(list(SeqIO.parse(current_input, 'fasta')))} reads for splitting")
                if self.verbose:
                    logger.info(f"Hits dictionary keys: {list(hits.keys())}")
                
                for record in SeqIO.parse(current_input, 'fasta'):
                    if self.verbose:
                        logger.debug(f"Processing read: {record.id} (length: {len(record.seq)})")
                    
                    if record.id in hits:
                        # Split this read
                        if self.verbose:
                            logger.info(f"Splitting read {record.id} with {len(hits[record.id])} hits")
                            logger.info(f"Hit details for {record.id}: {hits[record.id]}")
                        
                        fragments = self._split_read_at_primer(record, hits[record.id], primer_name)
                        
                        if self.verbose:
                            logger.info(f"Created {len(fragments)} fragments for read {record.id}")
                        
                        if len(fragments) > 0:
                            new_records.extend(fragments)
                            splits_this_round += 1
                            split_details.append({
                                'original_id': record.id,
                                'original_length': len(record.seq),
                                'fragments_created': len(fragments),
                                'fragment_ids': [f.id for f in fragments]
                            })
                            if self.verbose:
                                logger.debug(f"Split read {record.id} into {len(fragments)} fragments")
                        else:
                            if self.verbose:
                                logger.warning(f"No fragments created for read {record.id}, keeping original")
                            new_records.append(record)
                    else:
                        # Keep read unchanged
                        if self.verbose:
                            logger.debug(f"Read {record.id} has no hits, keeping unchanged")
                        new_records.append(record)
                
                # Save splitting details
                split_details_file = round_outputs_dir / f"round_{i+1}_split_details_{primer_name}.txt"
                with open(split_details_file, 'w') as f:
                    f.write(f"Splitting details for primer {primer_name}:\n")
                    f.write("=" * 50 + "\n")
                    f.write(f"Total reads processed: {len(new_records)}\n")
                    f.write(f"Reads split this round: {splits_this_round}\n\n")
                    for detail in split_details:
                        f.write(f"Original read: {detail['original_id']} (length: {detail['original_length']})\n")
                        f.write(f"  Fragments created: {detail['fragments_created']}\n")
                        for frag_id in detail['fragment_ids']:
                            f.write(f"    {frag_id}\n")
                        f.write("\n")
                logger.info(f"Saved splitting details to: {split_details_file}")
                
                logger.info(f"Total reads processed: {len(new_records)}")
                logger.info(f"Reads split this round: {splits_this_round}")
                
                # Write new file (always FASTA)
                new_input = temp_dir / f"round_{i+1}_output.fasta"
                SeqIO.write(new_records, new_input, 'fasta')
                
                # Save after-primer file for debugging
                after_primer_file = round_outputs_dir / f"round_{i+1}_after_{primer_name}.fasta"
                shutil.copy2(new_input, after_primer_file)
                logger.info(f"Saved reads after {primer_name} processing to: {after_primer_file}")
                
                # Update current input for next round
                current_input = new_input
                
                # Count reads after this round
                new_count = len(new_records)
                logger.info(f"Reads split with {primer_name}: {splits_this_round}")
                logger.info(f"Total reads after {primer_name}: {new_count}")
                
                total_splits += splits_this_round
            
            # Copy final result to output (always FASTA)
            shutil.copy2(current_input, self.output_file)
            
            # Save final round file
            final_round_file = round_outputs_dir / f"round_final_output.fasta"
            shutil.copy2(current_input, final_round_file)
            logger.info(f"Saved final output to: {final_round_file}")
            
            # Final statistics
            final_count = self._count_reads(self.output_file)
            logger.info(f"Final read count: {final_count}")
            logger.info(f"Total reads split: {total_splits}")
            logger.info(f"Net increase in reads: {final_count - initial_count}")
            logger.info(f"Output written to: {self.output_file}")
            logger.info(f"Temp files available in: {temp_dir}")
            logger.info(f"Round output files available in: {round_outputs_dir}")
            
        except Exception as e:
            logger.error(f"Error during processing: {e}")
            logger.info(f"Temp files are preserved in: {temp_dir} for debugging")
            logger.info(f"Round output files are preserved in: {round_outputs_dir} for debugging")
            raise


def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Split nanopore sequencing reads at primer positions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python adapter_splitter.py input.fastq primers.fasta output.fastq
  python adapter_splitter.py input.fasta primers.fasta output.fasta --blast-db-dir /tmp/blast_db
        """
    )
    
    parser.add_argument('input_file', help='Input FASTA/FASTQ file')
    parser.add_argument('primer_file', help='Primer sequences in FASTA format')
    parser.add_argument('output_file', help='Output file (always FASTA format)')
    parser.add_argument('--blast-db-dir', help='Directory for BLAST databases (default: current directory)')
    parser.add_argument('--identity', type=int, default=80, help='Minimum identity percentage (default: 80)')
    parser.add_argument('--coverage', type=int, default=90, help='Minimum query coverage percentage (default: 90)')
    parser.add_argument('--evalue', type=float, default=1e-10, help='Maximum E-value (default: 1e-10)')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    parser.add_argument('--verbose', action='store_true', help='Print detailed output for each split (default: summary only)')
    parser.add_argument('--save-blast-outputs', action='store_true', help='Save BLAST output files for debugging')
    parser.add_argument('--use-reversed-blast', action='store_true', help='Use reversed BLAST (reads as DB, primer as query) (default: False)')
    
    args = parser.parse_args()
    
    # Set up logging level
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
    
    try:
        # Create splitter instance
        splitter = NanoporeConcatemerSplitter(
            input_file=args.input_file,
            primer_file=args.primer_file,
            output_file=args.output_file,
            blast_db_dir=args.blast_db_dir,
            save_blast_outputs=args.save_blast_outputs,
            use_reversed_blast=args.use_reversed_blast,
            verbose=args.verbose
        )
        
        # Update BLAST parameters if provided via command line
        if args.identity != 80:
            splitter.blast_params['identity'] = args.identity
        if args.coverage != 90:
            splitter.blast_params['query_coverage'] = args.coverage
        if args.evalue != 1e-10:
            splitter.blast_params['evalue'] = args.evalue
        
        # Process the data
        splitter.process()
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
