# Nanopore Concatemer Splitter

A Python script for splitting nanopore sequencing reads containing concatemers by identifying primer sequences and splitting reads before the primer matches. This tool is specifically optimized for the lower accuracy of nanopore reads.

## Features

- **Strand-aware splitting**: Automatically detects forward vs reverse primer hits and splits accordingly
- **Multiple primer support**: Process multiple primers in sequence
- **Optimized for nanopore reads**: Uses relaxed BLAST parameters suitable for nanopore sequencing accuracy
- **Flexible filtering**: Uses primer overlap and identity thresholds instead of problematic coverage thresholds
- **Comprehensive debugging**: Optional verbose output and intermediate file saving
- **FASTQ/FASTA support**: Handles both input formats, always outputs FASTA

## Installation

### Prerequisites

- Python 3.6+
- BLAST+ (makeblastdb, blastn)
- Biopython

### Install BLAST+

**macOS:**
```bash
brew install blast
```

**Ubuntu/Debian:**
```bash
sudo apt-get install ncbi-blast+
```

**CentOS/RHEL:**
```bash
sudo yum install blast
```

### Install Python dependencies

```bash
pip install biopython
```

## Usage

### Basic Usage

```bash
python adapter_splitter.py input.fastq primers.fasta output.fasta
```

### Examples

**Example 1: Basic splitting with debug output**
```bash
python adapter_splitter.py first_20_reads.fastq adapter_primers.fasta output_test.fasta --debug --save-blast-outputs
```

**Example 2: Verbose output for detailed analysis**
```bash
python adapter_splitter.py first_20_reads.fastq adapter_primers.fasta output_verbose.fasta --verbose
```

**Example 3: Production run on large dataset**
```bash
python adapter_splitter.py input_sequences/trimmed-barcode01_ctxm.fastq adapter_primers.fasta output_split_fr.fasta
```

**Example 4: Custom BLAST parameters**
```bash
python adapter_splitter.py input.fastq primers.fasta output.fasta --identity 85 --evalue 1e-6
```

## Command Line Options

### Required Arguments
- `input_file`: Input FASTA/FASTQ file
- `primer_file`: Primer sequences in FASTA format
- `output_file`: Output file (always FASTA format)

### Optional Arguments
- `--blast-db-dir`: Directory for BLAST databases (default: current directory)
- `--identity`: Minimum identity percentage (default: 80)
- `--coverage`: Minimum query coverage percentage (default: 90)
- `--evalue`: Maximum E-value (default: 1e-10)
- `--debug`: Enable debug logging
- `--verbose`: Print detailed output for each split (default: summary only)
- `--save-blast-outputs`: Save BLAST output files for debugging
- `--use-reversed-blast`: Use reversed BLAST (reads as DB, primer as query) (default: False)

## Primer File Format

Create a FASTA file containing your primer sequences:

```fasta
>y_adapter
CCGTCACGCTGTTGTTAGG
>ctx-m(+)
CCGTCACGCTGTTGTTAGGCCGTCACGCTGTTGTTAGGCCGTCACGCTGTTGTTAGG
>KPC(-)
CTCATTCAAGGGCTTTCTTGCTGCCGCTGTGCTGGCTCGCAGCCAGCAGCAGGCCGGCTTGCTGGACACACCCATCCGTTACGGCAAAAATGCGCTGGT
```

## Output Format

The script outputs FASTA files with headers indicating the split positions and primers used:

```
>original_read_id[start:end][primer_name]
>original_read_id[primer_name][start:end]
```

**Examples:**
```
>50428d5d-cae7-4a52-a0d6-bafeae860a7c_part1[1:100][y_adapter]
>50428d5d-cae7-4a52-a0d6-bafeae860a7c_part1[y_adapter][101:200]
```

## Tunable Parameters

All parameters can be modified at the top of the script:

### BLAST Parameters
```python
BLAST_PARAMS = {
    'identity': 80,  # 80% identity for nanopore reads
    'query_coverage': 90,  # 90% query coverage
    'evalue': 1e-5,  # More permissive for short sequences
    'word_size': 4,  # Smaller word size for short primers
    'task': 'blastn-short'  # Optimized for short sequences
}
```

### Parsing Parameters
```python
PARSING_PARAMS = {
    'identity_threshold': 80,  # Minimum identity percentage for hits
    'min_primer_overlap': 0.9,  # Minimum fraction of primer that must overlap
}
```

### Splitting Parameters
```python
SPLITTING_PARAMS = {
    'min_fragment_length': 100,  # Minimum fragment length to keep
    'split_before_primer': True,  # Split before primer (True) or after primer (False)
}
```

### Debug Parameters
```python
DEBUG_PARAMS = {
    'save_blast_outputs': False,  # Save BLAST output files for debugging
    'save_round_outputs': True,  # Save intermediate files for each round
    'use_reversed_blast': False,  # Use reversed BLAST approach
    'verbose': False,  # Print detailed output for each split
}
```

## How It Works

1. **BLAST Search**: Uses BLAST to find primer sequences in reads
2. **Hit Filtering**: Filters hits based on identity and primer overlap thresholds
3. **Strand Detection**: Determines if primer hits are forward or reverse
4. **Splitting Logic**:
   - **Forward hits**: Split BEFORE the primer
   - **Reverse hits**: Split AFTER the primer
5. **Fragment Creation**: Creates new fragments with appropriate headers
6. **Iterative Processing**: Processes each primer sequentially

## Debugging

### Verbose Output
Use the `--verbose` flag to see detailed information about each split:

```bash
python adapter_splitter.py input.fastq primers.fasta output.fasta --verbose
```

This will show:
- Each read being processed
- Hit details for each primer
- Split positions and fragment creation
- Strand information for each hit

### Debug Files
Use `--save-blast-outputs` to save BLAST output files:

```bash
python adapter_splitter.py input.fastq primers.fasta output.fasta --save-blast-outputs
```

This creates:
- `blast_outputs/`: Directory containing BLAST output files
- `round_outputs/`: Directory containing intermediate files for each round
- `tmp_nanopore_splitter/`: Temporary files for debugging

### Debug Logging
Use `--debug` to enable detailed logging:

```bash
python adapter_splitter.py input.fastq primers.fasta output.fasta --debug
```

## Troubleshooting

### Common Issues

1. **No hits found**: Check primer sequences and BLAST parameters
2. **Empty output**: Verify input file format and primer sequences
3. **BLAST errors**: Ensure BLAST+ is properly installed
4. **Memory issues**: For large files, consider processing in chunks

### Performance Tips

- For large files (100,000+ reads), use `--verbose=False` to reduce output
- Use `--save-blast-outputs=False` to save disk space
- Adjust `min_fragment_length` based on your needs

## Output Files

- **Main output**: FASTA file with split fragments
- **Round outputs**: Intermediate files showing progress through each primer
- **BLAST outputs**: Raw BLAST results for debugging
- **Split details**: Detailed information about each split operation

## License

This project is open source and available under the MIT License.

## Citation

If you use this tool in your research, please cite:

```
Nanopore Concatemer Splitter
Author: Assistant
Date: 2024
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. 