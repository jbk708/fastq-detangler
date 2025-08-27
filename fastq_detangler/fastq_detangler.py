"""FASTQ detangler module for separating interweaved paired-end reads."""

import logging
import re
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple


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
        logger: Logger instance for tracking operations
    """

    def __init__(self, log_level: str = "INFO"):
        """Initialize the FastqDetangler.

        Args:
            log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        self.r1_pattern = re.compile(r"^@(.+)/1$")
        self.r2_pattern = re.compile(r"^@(.+)/2$")

        # Set up logging
        self.logger = logging.getLogger(__name__)
        try:
            self.logger.setLevel(getattr(logging, log_level.upper()))
        except AttributeError:
            self.logger.warning(f"Invalid log level '{log_level}', defaulting to INFO")
            self.logger.setLevel(logging.INFO)

        # Create console handler if none exists
        if not self.logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            )
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)

    def detangle(self, input_file: Path, output_prefix: str) -> None:
        """Detangle interweaved FASTQ file into separate output files.

        Args:
            input_file: Path to input interweaved FASTQ file
            output_prefix: Prefix for output file names

        Raises:
            FileNotFoundError: If input file doesn't exist
            ValueError: If input file is empty or invalid
            OSError: If there are issues writing output files
        """
        start_time = time.time()
        self.logger.info("Starting FASTQ detangling process")
        self.logger.info(f"Input file: {input_file}")
        self.logger.info(f"Output prefix: {output_prefix}")

        # Input validation
        if not input_file.exists():
            error_msg = f"Input file not found: {input_file}"
            self.logger.error(error_msg)
            raise FileNotFoundError(error_msg)

        if not input_file.is_file():
            error_msg = f"Input path is not a file: {input_file}"
            self.logger.error(error_msg)
            raise ValueError(error_msg)

        if input_file.stat().st_size == 0:
            error_msg = f"Input file is empty: {input_file}"
            self.logger.error(error_msg)
            raise ValueError(error_msg)

        self.logger.info(f"Input file size: {input_file.stat().st_size:,} bytes")

        # Parse FASTQ file
        parse_start = time.time()
        self.logger.info("Parsing FASTQ file...")
        r1_reads, r2_reads = self._parse_fastq_file(input_file)
        parse_time = time.time() - parse_start

        if not r1_reads and not r2_reads:
            error_msg = "Input file contains no valid reads"
            self.logger.error(error_msg)
            raise ValueError(error_msg)

        self.logger.info(f"Parsing completed in {parse_time:.2f} seconds")
        self.logger.info(f"Total R1 reads found: {len(r1_reads):,}")
        self.logger.info(f"Total R2 reads found: {len(r2_reads):,}")

        # Identify missing and paired reads
        analysis_start = time.time()
        self.logger.info("Analyzing read pairs...")
        missing_r1, missing_r2 = self._identify_missing_pairs(r1_reads, r2_reads)
        paired_r1, paired_r2 = self._identify_paired_reads(r1_reads, r2_reads)
        analysis_time = time.time() - analysis_start

        self.logger.info(f"Analysis completed in {analysis_time:.2f} seconds")
        self.logger.info(f"R1 reads with missing pairs: {len(missing_r1):,}")
        self.logger.info(f"R2 reads with missing pairs: {len(missing_r2):,}")
        self.logger.info(f"Paired R1 reads: {len(paired_r1):,}")
        self.logger.info(f"Paired R2 reads: {len(paired_r2):,}")

        # Write output files
        write_start = time.time()
        self.logger.info("Writing output files...")

        output_files = []
        try:
            r1_missing_file = f"{output_prefix}_R1_ordered_with_missing_R2.fastq"
            r2_missing_file = f"{output_prefix}_R2_ordered_with_missing_R1.fastq"
            r1_paired_file = f"{output_prefix}_R1_paired.fastq"
            r2_paired_file = f"{output_prefix}_R2_paired.fastq"

            self._write_fastq_file(missing_r1, r1_missing_file)
            output_files.append(r1_missing_file)
            self.logger.info(f"Written: {r1_missing_file} ({len(missing_r1):,} reads)")

            self._write_fastq_file(missing_r2, r2_missing_file)
            output_files.append(r2_missing_file)
            self.logger.info(f"Written: {r2_missing_file} ({len(missing_r2):,} reads)")

            self._write_fastq_file(paired_r1, r1_paired_file)
            output_files.append(r1_paired_file)
            self.logger.info(f"Written: {r1_paired_file} ({len(paired_r1):,} reads)")

            self._write_fastq_file(paired_r2, r2_paired_file)
            output_files.append(r2_paired_file)
            self.logger.info(f"Written: {r2_paired_file} ({len(paired_r2):,} reads)")

        except Exception as e:
            error_msg = f"Error writing output files: {e}"
            self.logger.error(error_msg)
            # Clean up any files that were created
            for file_path in output_files:
                try:
                    Path(file_path).unlink(missing_ok=True)
                except Exception:
                    pass
            raise OSError(error_msg) from e

        write_time = time.time() - write_start
        total_time = time.time() - start_time

        self.logger.info(f"File writing completed in {write_time:.2f} seconds")
        self.logger.info(f"Total processing time: {total_time:.2f} seconds")
        self.logger.info(
            f"Processing speed: {len(r1_reads) + len(r2_reads) / total_time:.0f} reads/second"
        )
        self.logger.info("FASTQ detangling completed successfully")

    def _identify_read_type(self, header: str) -> Optional[str]:
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

    def _parse_fastq_file(
        self, input_file: Path
    ) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
        """Parse FASTQ file into R1 and R2 read collections.

        Args:
            input_file: Path to input FASTQ file

        Returns:
            Tuple of (r1_reads, r2_reads) where each is a dict mapping
            base names to lists of 4-line FASTQ records
        """
        r1_reads = {}
        r2_reads = {}
        total_lines = 0
        valid_reads = 0
        invalid_reads = 0
        malformed_records = 0

        try:
            with open(input_file, "r") as f:
                lines = f.readlines()

            total_lines = len(lines)
            self.logger.debug(f"Processing {total_lines:,} lines from input file")

            i = 0
            while i < len(lines):
                if lines[i].startswith("@"):
                    header = lines[i].strip()
                    read_type = self._identify_read_type(header)

                    if read_type and i + 3 < len(lines):
                        # Check if we have a complete 4-line record
                        record = lines[i : i + 4]
                        if len(record) == 4 and record[2].startswith("+"):
                            base_name = self._extract_base_name(header)

                            if read_type == "R1":
                                r1_reads[base_name] = record
                            else:
                                r2_reads[base_name] = record

                            valid_reads += 1
                            i += 4
                        else:
                            self.logger.warning(
                                f"Malformed record at line {i+1}: incomplete or invalid FASTQ format"
                            )
                            malformed_records += 1
                            i += 1
                    else:
                        if not read_type:
                            self.logger.warning(
                                f"Invalid read header at line {i+1}: {header}"
                            )
                        else:
                            self.logger.warning(
                                f"Incomplete record at line {i+1}: expected 4 lines, got {min(len(lines) - i, 4)}"
                            )
                        invalid_reads += 1
                        i += 1
                else:
                    i += 1

            self.logger.debug(
                f"Parsing summary: {valid_reads:,} valid reads, {invalid_reads:,} invalid headers, {malformed_records:,} malformed records"
            )

        except Exception as e:
            error_msg = f"Error reading input file: {e}"
            self.logger.error(error_msg)
            raise OSError(error_msg) from e

        return r1_reads, r2_reads

    def _identify_missing_pairs(
        self, r1_reads: Dict[str, List[str]], r2_reads: Dict[str, List[str]]
    ) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
        """Identify reads with missing pairs.

        Args:
            r1_reads: Dictionary of R1 reads
            r2_reads: Dictionary of R2 reads

        Returns:
            Tuple of (missing_r1, missing_r2) where each contains reads
            that don't have matching pairs
        """
        missing_r1 = {
            name: reads for name, reads in r1_reads.items() if name not in r2_reads
        }
        missing_r2 = {
            name: reads for name, reads in r2_reads.items() if name not in r1_reads
        }

        self.logger.debug(
            f"Missing pairs identified: {len(missing_r1):,} R1, {len(missing_r2):,} R2"
        )

        return missing_r1, missing_r2

    def _identify_paired_reads(
        self, r1_reads: Dict[str, List[str]], r2_reads: Dict[str, List[str]]
    ) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
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

        self.logger.debug(f"Paired reads identified: {len(paired_r1):,} pairs")

        return paired_r1, paired_r2

    def _write_fastq_file(self, reads: Dict[str, List[str]], output_file: str) -> None:
        """Write reads to FASTQ file in sorted order.

        Args:
            reads: Dictionary of reads to write
            output_file: Path to output file

        Raises:
            OSError: If there are issues writing the file
        """
        try:
            with open(output_file, "w") as f:
                for base_name in sorted(reads.keys()):
                    f.writelines(reads[base_name])

            # Verify file was written and has content
            output_path = Path(output_file)
            if not output_path.exists():
                raise OSError(f"Output file was not created: {output_file}")

            if output_path.stat().st_size == 0:
                self.logger.warning(f"Output file is empty: {output_file}")

        except Exception as e:
            error_msg = f"Error writing to {output_file}: {e}"
            self.logger.error(error_msg)
            raise OSError(error_msg) from e
