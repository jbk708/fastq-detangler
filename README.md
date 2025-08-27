# FASTQ Detangler

A Python tool to detangle interweaved FASTQ files containing mixed R1 and R2 reads into separate, organized output files.

## Overview

The FASTQ Detangler processes interweaved FASTQ files where R1 and R2 reads are mixed together in random order and separates them into four organized output files:

1. **R1 reads with missing R2 pairs** - R1 reads that don't have matching R2 reads
2. **R2 reads with missing R1 pairs** - R2 reads that don't have matching R1 reads
3. **Paired R1 reads** - R1 reads that have matching R2 reads
4. **Paired R2 reads** - R2 reads that have matching R1 reads

## Features

- **Automatic read identification** - Identifies R1/R2 reads based on header patterns (`/1` and `/2` suffixes)
- **Missing pair detection** - Finds reads that don't have matching pairs
- **Ordered output** - All output files are sorted by read name for consistency
- **Efficient processing** - Single-pass parsing of input files
- **No external dependencies** - Uses only Python standard library
- **Comprehensive testing** - Full test suite with example data

## Installation

### From source

```bash
git clone https://github.com/yourusername/fastq-detangler.git
cd fastq-detangler
pip install -e .
```

### Development installation

```bash
pip install -e ".[dev]"
```

## Usage

### Command Line Interface

```bash
# Basic usage
python -m fastq_detangler input.fastq output_prefix

# Example
python -m fastq_detangler interweaved_reads.fastq sample_001
```

This will create four output files:
- `sample_001_R1_ordered_with_missing_R2.fastq`
- `sample_001_R2_ordered_with_missing_R1.fastq`
- `sample_001_R1_paired.fastq`
- `sample_001_R2_paired.fastq`

### Python API

```python
from fastq_detangler import FastqDetangler
from pathlib import Path

# Create detangler instance
detangler = FastqDetangler()

# Process input file
input_file = Path("interweaved_reads.fastq")
detangler.detangle(input_file, "output_prefix")
```

## Input Format

The tool expects FASTQ files with read headers that end in `/1` for R1 reads and `/2` for R2 reads:

```
@read1/1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read1/2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

## Output Files

### R1 reads with missing R2 pairs
Contains R1 reads that don't have corresponding R2 reads in the input file.

### R2 reads with missing R1 pairs
Contains R2 reads that don't have corresponding R1 reads in the input file.

### Paired R1 reads
Contains R1 reads that have matching R2 reads, sorted by read name.

### Paired R2 reads
Contains R2 reads that have matching R1 reads, sorted by read name.

## Testing

Run the test suite:

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=fastq_detangler

# Run specific test file
pytest tests/test_fastq_detangler.py
```

## Development

### Code Style

The project follows Google Python Style Guide with:
- Google docstring format
- Minimal inline comments
- KISS (Keep It Simple, Stupid) programming principles

### Running Tests

```bash
# Install development dependencies
pip install -r requirements-dev.txt

# Run tests
pytest

# Run linting
flake8 fastq_detangler tests

# Format code
black fastq_detangler tests
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass
6. Submit a pull request

## Example

The repository includes a test dataset with 25 reads (15 R1, 10 R2) where some reads are missing pairs. You can test the tool with:

```bash
python -m fastq_detangler tests/test_data/interweaved_input.fastq test_output
```

This will create the four output files demonstrating the tool's functionality.
