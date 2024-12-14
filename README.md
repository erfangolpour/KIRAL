# KIRAL
KIRAL is a fast and efficient tool for aligning sequencing reads to a database of over 1600 Killer Immunoglobulin-like Receptor (KIR) allele sequences. Designed for genomic research and immunological studies, KIRAL provides powerful options for aligning reads, selecting alignment methods, and generating detailed reports.

## Features
- **Multiple Alignment Methods**: Choose between `naive`, `regional`, and `categorical` alignment strategies to suit your analysis.
- **High Customizability**: Adjust parameters such as the number of representative alleles, error tolerance, and paired-read alignment for maximum flexibility.
- **Threading Support**: Take full advantage of multi-core systems with configurable threading for faster processing.
- **Comprehensive Reporting**: Generate detailed reports filtered by read, KIR, or allele identifiers.

---

## Installation
To use KIRAL, clone this repository, then clone and build [minimap2](https://github.com/lh3/minimap2) in the cloned directory, and build the tool:
```bash
git clone https://github.com/erfangolpour/KIRAL
cd KIRAL
git clone https://github.com/lh3/minimap2
cd minimap2
make
cd ..
make
```

---

## Usage

### 1. Align Reads
Use the `align` command to align sequencing reads to the KIR allele database.

#### Command:
```bash
./main align <database> <reads> [--method <method_name>] [-r <num_representatives>] [--pair] [-t <threads>] [-o <output_file>]
```

#### Options:
- **`<database>`**: Path to the KIR allele database.
- **`<reads>`**: Path to the sequencing reads file.
- **`--method <method_name>`**: Alignment method to use. Options are:
  - `naive` (simple alignment to all alleles),
  - `regional` (align to representatives for each gene),
  - `categorical` (align by categorical grouping).  
  Default: `regional`.
- **`-r <num_representatives>`**: Number of representative alleles per gene for `regional` or `categorical` alignment. Default: 1.
- **`--pair`**: Include paired reads in the second pass for `regional` or `categorical` alignment.
- **`-t <threads>`**: Number of threads to use. Default: Number of hardware threads.
- **`-o <output_file>`**: Path to save the alignment results.

### 2. Analyze Reports
Use the `report` command to analyze and filter results from a previously generated alignment file.

#### Command:
```bash
./main report <alignments_file> [--head <num_results>] [-r <read_id>] [-k <KIR read_id>] [-a <allele read_id>]
```

#### Options:
- **`<alignments_file>`**: Path to the file with alignment results.
- **`--head <num_results>`**: Display only the first `<num_results>` alignments.
- **`-r <read_id>`**: Show results for a specific read ID.
- **`-k <KIR read_id>`**: Show results for a specific KIR ID.
- **`-a <allele read_id>`**: Show results for a specific allele ID.

---

## Examples

### Example 1: Align Reads with Default Settings
```bash
./main align KIR_database.fasta reads.fastq -o alignments.txt
```

### Example 2: Use 4 Threads and Pair Alignment
```bash
./main align KIR_database.fasta reads.fastq --method regional --pair -t 4 -o paired_alignments.txt
```

### Example 3: Generate a Report for Specific Read ID
```bash
./main report alignments.txt -r 12345
```

### Example 4: Show Top 10 Results
```bash
./main report alignments.txt --head 10
```

---

## Contributing
Contributions are welcome! Feel free to submit issues or pull requests to improve KIRAL.

---

## License
KIRAL is released under the [GNU General Public License v3.0](LICENSE).

---

## Acknowledgments

Special thanks to:
- The developers of [minimap2](https://github.com/lh3/minimap2) for their amazing project. minimap2’s exceptional performance and versatility in sequence alignment have been instrumental in the development of KIRAL.
