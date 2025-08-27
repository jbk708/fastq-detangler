"""FASTQ detangler module for separating interweaved paired-end reads."""

import re
from pathlib import Path
from typing import Dict, Tuple, List


class FastqDetangler:
    """Detangles interweaved FASTQ files into separate R1/R2 files.
    
    This class processes interweaved FASTQ files containing mixed R1 and R2 reads
    and separates them into four output files:
    - R1 reads with missing R2 pairs
    - R2 reads with missing R1 pairs  
    - Paired R1 reads
    - Paired R2 reads
    
    Attributes:
        r1_pattern: Regex pattern to identify R1 read headers
        r2_pattern: Regex pattern to identify R2 read headers
    """
    
    def __init__(self):
        """Initialize the FastqDetangler."""
        self.r1_pattern = re.compile(r'^@(.+)/1$')
        self.r2_pattern = re.compile(r'^@(.+)/2$')
    
    def detangle(self, input_file: Path, output_prefix: str) -> None:
        """Detangle interweaved FASTQ file into separate output files.
        
        Args:
            input_file: Path to input interweaved FASTQ file
            output_prefix: Prefix for output file names
            
        Raises:
            FileNotFoundError: If input file doesn't exist
            ValueError: If input file is empty or invalid
        """
        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        r1_reads, r2_reads = self._parse_fastq_file(input_file)
        
        if not r1_reads and not r2_reads:
            raise ValueError("Input file contains no valid reads")
        
        missing_r1, missing_r2 = self._identify_missing_pairs(r1_reads, r2_reads)
        paired_r1, paired_r2 = self._identify_paired_reads(r1_reads, r2_reads)
        
        self._write_fastq_file(missing_r1, f"{output_prefix}_R1_ordered_with_missing_R2.fastq")
        self._write_fastq_file(missing_r2, f"{output_prefix}_R2_ordered_with_missing_R1.fastq")
        self._write_fastq_file(paired_r1, f"{output_prefix}_R1_paired.fastq")
        self._write_fastq_file(paired_r2, f"{output_prefix}_R2_paired.fastq")
    
    def _identify_read_type(self, header: str) -> str:
        """Identify if a read header is R1 or R2.
        
        Args:
            header: Read header line
            
        Returns:
            "R1" for R1 reads, "R2" for R2 reads, None if invalid
        """
        if self.r1_pattern.match(header):
            return "R1"
        elif self.r2_pattern.match(header):
            return "R2"
        return None
    
    def _extract_base_name(self, header: str) -> str:
        """Extract base read name from header.
        
        Args:
            header: Read header line
            
        Returns:
            Base read name without /1 or /2 suffix
        """
        match = self.r1_pattern.match(header) or self.r2_pattern.match(header)
        return match.group(1) if match else header[1:]
    
    def _parse_fastq_file(self, input_file: Path) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
        """Parse FASTQ file into R1 and R2 read collections.
        
        Args:
            input_file: Path to input FASTQ file
            
        Returns:
            Tuple of (r1_reads, r2_reads) where each is a dict mapping
            base names to lists of 4-line FASTQ records
        """
        r1_reads = {}
        r2_reads = {}
        
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        i = 0
        while i < len(lines):
            if lines[i].startswith('@'):
                header = lines[i].strip()
                read_type = self._identify_read_type(header)
                
                if read_type and i + 3 < len(lines):
                    base_name = self._extract_base_name(header)
                    record = lines[i:i+4]
                    
                    if read_type == "R1":
                        r1_reads[base_name] = record
                    else:
                        r2_reads[base_name] = record
                    
                    i += 4
                else:
                    i += 1
            else:
                i += 1
        
        return r1_reads, r2_reads
    
    def _identify_missing_pairs(self, r1_reads: Dict[str, List[str]], 
                               r2_reads: Dict[str, List[str]]) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
        """Identify reads with missing pairs.
        
        Args:
            r1_reads: Dictionary of R1 reads
            r2_reads: Dictionary of R2 reads
            
        Returns:
            Tuple of (missing_r1, missing_r2) where each contains reads
            that don't have matching pairs
        """
        missing_r1 = {name: reads for name, reads in r1_reads.items() 
                     if name not in r2_reads}
        missing_r2 = {name: reads for name, reads in r2_reads.items() 
                     if name not in r1_reads}
        
        return missing_r1, missing_r2
    
    def _identify_paired_reads(self, r1_reads: Dict[str, List[str]], 
                              r2_reads: Dict[str, List[str]]) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
        """Identify reads with matching pairs.
        
        Args:
            r1_reads: Dictionary of R1 reads
            r2_reads: Dictionary of R2 reads
            
        Returns:
            Tuple of (paired_r1, paired_r2) where each contains reads
            that have matching pairs
        """
        paired_names = set(r1_reads.keys()) & set(r2_reads.keys())
        
        paired_r1 = {name: r1_reads[name] for name in paired_names}
        paired_r2 = {name: r2_reads[name] for name in paired_names}
        
        return paired_r1, paired_r2
    
    def _write_fastq_file(self, reads: Dict[str, List[str]], output_file: str) -> None:
        """Write reads to FASTQ file in sorted order.
        
        Args:
            reads: Dictionary of reads to write
            output_file: Path to output file
        """
        with open(output_file, 'w') as f:
            for base_name in sorted(reads.keys()):
                f.writelines(reads[base_name]) 