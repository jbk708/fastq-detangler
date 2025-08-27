"""Command-line interface for FASTQ detangler."""

import argparse
import sys
import logging
from pathlib import Path

from .fastq_detangler import FastqDetangler


def main():
    """Main entry point for command-line interface."""
    parser = argparse.ArgumentParser(
        description="Detangle interweaved FASTQ files into separate R1/R2 files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m fastq_detangler input.fastq output_prefix
  python -m fastq_detangler /path/to/input.fastq /path/to/output
  python -m fastq_detangler --log-level DEBUG input.fastq output_prefix
        """
    )
    
    parser.add_argument(
        "input_file",
        type=Path,
        help="Input interweaved FASTQ file"
    )
    
    parser.add_argument(
        "output_prefix",
        type=str,
        help="Output file prefix (will create 4 files with this prefix)"
    )
    
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set logging level (default: INFO)"
    )
    
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress all output except errors"
    )
    
    args = parser.parse_args()
    
    # Set up logging level
    if args.quiet:
        log_level = "ERROR"
    else:
        log_level = args.log_level
    
    try:
        # Create detangler with specified log level
        detangler = FastqDetangler(log_level=log_level)
        
        # Run detangling
        detangler.detangle(args.input_file, args.output_prefix)
        
        # Print summary if not quiet
        if not args.quiet:
            print(f"\n✅ Successfully detangled {args.input_file}")
            print(f"📁 Output files created with prefix: {args.output_prefix}")
            print("\n📊 Files created:")
            print(f"  • {args.output_prefix}_R1_ordered_with_missing_R2.fastq")
            print(f"  • {args.output_prefix}_R2_ordered_with_missing_R1.fastq")
            print(f"  • {args.output_prefix}_R1_paired.fastq")
            print(f"  • {args.output_prefix}_R2_paired.fastq")
            print("\n🎉 Processing completed successfully!")
        
    except FileNotFoundError as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)
    except OSError as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\n⚠️  Operation cancelled by user", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Unexpected error: {e}", file=sys.stderr)
        if log_level == "DEBUG":
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main() 