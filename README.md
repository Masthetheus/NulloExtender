# NulloExtender

NulloExtender: Constrained Seed Extension from Nullomers and Minimal Absent Words

**Version 1.0.0**  | **Experimental Release**

** The current version is extremely experimental, and its results  shouldn't be taken as face value for final usage. **

## 1. Description

NulloExtender is a command-line utility written in C for selecting and extending nullomers or minimal absent words from binary genomic indices. The program supports inputting of sequence constraints, up to this moments specifically maximum GC content and maximum homopolymer run length, to generate longer absent oligomers. The tool is intended for applications in barcoding design and other practical usage of nullomers where long sequences are needed but time and computational resources are limited.

This implementation prioritises simplicity and correct adjustment to the filters, with further work intended to tailor for memory efficiency and deterministic execution. It uses an in-place Fisher-Yates swap to guarantee non-redundant seed selection without auxiliary data structures.

## 2. Compilation and Dependencies

### 2.1 Requirements

- A C99-compliant compiler (GCC, Clang, or equivalent).
- Standard C library (libc). No external dependencies are required.

### 2.2 Build Instructions

Currently, given the whole code being present on a single file, clone the repository and compile the source file directly:
```
git clone https://github.com/Masthetheus/NulloExtender.git
cd NulloExtender
gcc -O3 -std=c99 -o nullo_extender src/main.c -lm
```
An optional Makefile is currently being worked on, given the future objective of organizing the scrip in proper files.

## 3. Usage

### 3.1 Command-line syntax
```
./nullo_extender <nullomer_file> <gc_max> <homopolymer_max> <number_of_seeds> <target_length>
```
### 3.2 Arguments

Argument: nullomer_file
Type: string
Description: Path to the binary nullomer index produced by the NulloRetriever pipeline.
Obs: NulloRetrieve has a tool for translating nullomers into the correct format. In the future accepting oher type of nullomeric files may be added.

Argument: gc_max
Type: float
Description: Maximum allowed GC percentage, expressed as a value between 0 and 100 (e.g., 60 denotes 60 percent).

Argument: homopolymer_max
Type: int
Description: Maximum allowable length of consecutive identical bases (e.g., 3 forbids runs of four or more A, C, T, or G).

Argument: number_of_seeds
Type: int
Description: The number of distinct seed sequences to select and extend. Must not exceed the total number of nullomers present in the index.

Agument: target_length
Type: int
Description: Size of the final extended sequence.

## 3.3 Example
```
./nullo_extender data/subtilis_result 55 3 100 20
```
This invocation reads the subtilis_result index, selects 100 seeds with a GC content not exceeding 55 percent, and rejects any seed containing homopolymers longer than three bases. The extended sequences, independent of the original k value, will have length equal to k = 20 and are written to standard output. Writing to output file is currently at work.


## 4. Algorithmic Overview

The program currently proceeds through four consecutive phases. Due its alpha state, a lot of code is up to enhancement.

### 4.1 Index parsing

The binary index is read according to a fixed header format. The header stores the k-mer length (k), the half-k parameter, the byte sizes for the prefix and counter fields, and the total number of records. The nullomer list is loaded into a dynamically resized array. Reallocation follows a doubling strategy starting from an initial capacity of 100,000 elements.

### 4.2 Seed selection

A seed is selected uniformly at random from the remaining pool. To guarantee uniqueness across all n selections, the selected element is swapped with the last unselected element, and the pool size is decremented. This procedure is equivalent to a partial Fisher-Yates shuffle.

### 4.3 Constraint validation

For each candidate seed, the program computes:
- The total number of G and C bases.
- The length of the longest contiguous run of identical bases.

If the GC count exceeds (gc_max / 100) * k, or if the longest homopolymer run exceeds homopolymer_max, the seed is rejected. A replacement seed is drawn from the remaining pool, and the validation is repeated.

### 4.4 Extension

Accepted seeds are extended from length k to length target_length by sequential random base addition. At each extension step, the algorithm attempts to choose a base that does not violate the GC and homopolymer constraints for the extended sequence. If the initially drawn base violates a constraint, alternative bases are tested.

## 5. Input File Format

The program expects the binary format produced by the NulloRetriever pipeline. A tool aimed to transcribing raw txt to such format is currently available on NulloRetriever and shall be ported to NulloExtender in due time.

RAM usage will directly relate to the nullomeric set size, complexity and seeds targeted.

## 6. Output

The program writes extended sequences to standard output in plain text format, one sequence per line. No FASTA headers are added by default; however, redirection may be used to capture the output:
```
./maw_extension data/hg38.nullomer 55 3 100 20 > selected_sequences.txt
```
Diagnostic messages, including the selection of replacement seeds, are printed to stderr to avoid interfering with the primary output stream.
Some DEBUG prints are currently available and shall be patched soon.

## 7. Limitations and Future Work

- Orthogonality: The current version does not enforce pairwise Hamming distance between extended sequences. Two seeds may converge to similar sequences after extension. A future release will incorporate post-extension distance checks and, if necessary, automatic reseeding.
- Parallelisation: The seed selection and extension loops are currently sequential. Multi-threading support is planned for large-scale designs (n greater than 10,000).
- Output formats: Additional output modes (FASTA, CSV, and BED) are under consideration to facilitate integration with downstream alignment tools.

## 8. Citation

If you use this software in a publication, please cite it as:

Cassol, Matheus Pedron Cassol, "NulloExtender" GitHub repository, https://github.com/Masthetheus/NulloExtender.git, 2026.

## 9. License

This project is distributed under the GPL-3.0 License. See the LICENSE file in the repository root for full terms.


## 10. Contact

Correspondence and bug reports should be directed to:

Matheus Pedron Cassol
matheuspedroncassol@gmail.com
LBCM
Universidade Federal do Rio Grande do Sul

Project repository: https://github.com/Masthetheus/NulloExtender.git.


## 11. Acknowledgements

This software was developed at LBCM-UFRGS and supported in spirit and mind by my colleagues.
