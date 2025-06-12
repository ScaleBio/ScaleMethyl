#!/usr/bin/env python
"""
Generate chrom sizes file from fasta
"""
import argparse
from pathlib import Path


def parse_fasta(fasta_file: Path):
    """
    Parse a fasta file and yield (seq_id, sequence) tuples

    Args:
        fasta_file: path to the fasta file
    """
    with open(fasta_file, "r") as f:
        seq_id = None
        sequence = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    yield (seq_id, "".join(sequence))
                seq_id = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
        if seq_id:
            yield (seq_id, "".join(sequence))


def generate_chrom_sizes(fasta_file: Path, output_file: Path):
    """
    Generate a chrom sizes file from a fasta file

    Args:
        fasta_file: path to the fasta file
        output_file: path to the output chrom sizes file
    """
    with open(output_file, "w") as out_f:
        for seq_id, sequence in parse_fasta(fasta_file):
            out_f.write(f"{seq_id}\t{len(sequence)}\n")


def main():
    parser = argparse.ArgumentParser(description="Generate chrom sizes file from fasta")
    parser.add_argument("--fastaFile", type=Path, help="Input fasta file")
    parser.add_argument("--output", type=Path, help="Output chrom sizes file")
    args = parser.parse_args()

    generate_chrom_sizes(args.fastaFile, args.output)


if __name__ == "__main__":
    main()
