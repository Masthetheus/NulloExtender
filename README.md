# NulloExtender

NulloExtender: Constrained Seed Extension from Nullomers and Minimal Absent Words

**Version 1.1.0**  | **Experimental Release**

**The current version is extremely experimental, and its results shouldn't be taken as face value for final usage.**

## 1. Description

NulloExtender is a command-line utility written in C for selecting and extending nullomers or minimal absent words from binary genomic indices. The program supports sequence constraints on GC content (minimum and maximum), homopolymer run length, melting temperature (Tm), and pairwise orthogonality (Hamming distance), to generate longer absent oligomers. The tool is intended for applications in barcoding design and other practical usage of nullomers where long sequences are needed but time and computational resources are limited.

This implementation prioritises simplicity and correct adjustment to the filters, with further work intended to tailor for memory efficiency, deterministic execution, and support for sequences longer than 32 bases (see Section 7). It uses an in-place Fisher-Yates swap to guarantee non-redundant seed selection without auxiliary data structures.

## 2. Compilation and Dependencies

### 2.1 Requirements

- A C99-compliant compiler (GCC, Clang, or equivalent).
- Standard C library (libc), including `<math.h>` (link with `-lm`).
- No external dependencies are required.

### 2.2 Build Instructions

The source is organized into separate translation units under `src/`, with corresponding headers under `include/`. Build using the provided Makefile:
```
git clone https://github.com/Masthetheus/NulloExtender.git
cd NulloExtender
make
```
The resulting binary is placed at `bin/nullo_extender`. To rebuild from scratch:
```
make clean && make
```

## 3. Usage

### 3.1 Command-line syntax
```
./bin/nullo_extender <nullomer_file> <gc_max> <gc_min> <homopolymer_max> <number_of_seeds> <target_length> <tm_min> <tm_max>
```

### 3.2 Arguments

| Argument | Type | Description |
|---|---|---|
| `nullomer_file` | string | Path to the binary nullomer index produced by the NulloRetriever pipeline. |
| `gc_max` | float | Maximum allowed GC percentage (0–100). |
| `gc_min` | float | Minimum allowed GC percentage (0–100). Must be strictly lower than `gc_max`. |
| `homopolymer_max` | int | Maximum allowable length of consecutive identical bases (e.g., 3 forbids runs of four or more identical bases). |
| `number_of_seeds` | int | Number of distinct seed sequences to select and extend. Must not exceed the total number of nullomers present in the index. |
| `target_length` | int | Length of the final extended sequence. **Currently limited to 32 bases** due to the internal `uint64_t` bit-packed representation (see Section 7). |
| `tm_min` | float | Minimum acceptable melting temperature, in °C. |
| `tm_max` | float | Maximum acceptable melting temperature, in °C. Must be strictly greater than `tm_min`. |

Sodium concentration used in the salt-correction term of the Tm calculation is currently fixed internally (`NA_CONC = 0.05`, i.e. 50mM) and is not user-configurable in this release.

### 3.3 Example
```
./bin/nullo_extender data/subtilis_result 55 30 3 20 32 55 65
```
This invocation reads the `subtilis_result` index, selects 20 seeds with GC content between 30% and 55%, rejects any candidate containing homopolymer runs longer than 3 bases, extends each accepted seed to a final length of 32 bases, and requires a melting temperature between 55°C and 65°C. Extended sequences are written to standard output.

## 4. Algorithmic Overview

The program proceeds through five consecutive phases. Due to its alpha state, a lot of code is up for enhancement.

### 4.1 Index parsing

The binary index is read according to a fixed header format. The header stores the k-mer length (k), the half-k parameter, and the byte sizes for the prefix and counter fields. The nullomer list is loaded into a dynamically resized array, growing by doubling from an initial capacity of 100,000 elements.

### 4.2 Seed selection

A seed is selected uniformly at random from the remaining pool. To guarantee uniqueness across all `n` selections, the selected element is swapped with the last unselected element of the pool, and the pool size is decremented — a partial Fisher-Yates shuffle.

### 4.3 Seed validation

Each candidate seed is checked, base by base, against the GC maximum and homopolymer constraints, while simultaneously accumulating nearest-neighbor thermodynamic parameters (ΔH, ΔS) for the eventual Tm calculation. Seeds failing GC or homopolymer constraints are discarded and replaced.

### 4.4 Extension

Accepted seeds are extended from length k to `target_length` by sequential random base selection. At each step, a base is proposed and checked against the running GC and homopolymer constraints; if rejected, an alternative base is drawn until one satisfies both. Nearest-neighbor thermodynamic contributions are accumulated throughout the extension.

Once a sequence reaches its target length, two additional checks are performed on the completed sequence as a whole (not incrementally, since both properties are only meaningful over the full sequence):
- **Minimum GC content**, computed against `target_length`.
- **Melting temperature (Tm)**, computed via the nearest-neighbor method (unified SantaLucia 1998 parameters) with a fixed sodium concentration salt correction.

If either check fails, the sequence is discarded, a new seed is drawn, and the seed-validation and extension steps are repeated for that slot.

### 4.5 Orthogonality

Once all `n` sequences have been generated and individually validated, each pair is compared using Hamming distance (computed via XOR and popcount over the bit-packed representation). Pairs with a Hamming distance smaller than 3 are considered insufficiently orthogonal; in that case, one of the two sequences is discarded, a new seed is drawn and extended in its place, and the comparison is repeated.

## 5. Input File Format

The program expects the binary format produced by the NulloRetriever pipeline. A tool for converting raw text nullomer lists into this format is available in NulloRetriever and is planned to be ported to NulloExtender in a future release.

RAM usage scales with the size and complexity of the nullomer set and the number of seeds requested.

## 6. Output

The program writes extended sequences to standard output in plain text format, one sequence per line. No FASTA headers are added by default; redirection may be used to capture the output:
```
./bin/nullo_extender data/hg38.nullomer 55 30 3 100 32 55 65 > selected_sequences.txt
```
**Known issue:** diagnostic and status messages (e.g. "All seeds checked", seed-replacement notices) are currently printed to standard output rather than standard error, and are therefore mixed with the primary sequence output. Redirecting stdout to a file will currently capture these messages as well. This will be corrected in a future release.

## 7. Limitations and Future Work

- **Sequence length ceiling:** the internal representation packs each sequence into a single `uint64_t` (2 bits/base), limiting `target_length` to 32 bases. Support for longer sequences (up to ~250 bases, as required by some barcoding applications) requires migrating to a multi-block representation and is planned for a future release.
- **Fixed salt concentration:** sodium concentration for the Tm salt correction is hardcoded and not currently exposed as a parameter.
- **Output formats:** additional output modes (FASTA, CSV, BED) are under consideration to facilitate integration with downstream tools.
- **Diagnostic/output stream separation:** see Section 6.

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

Project repository: https://github.com/Masthetheus/NulloExtender.git

## 11. Acknowledgements

This software was developed at LBCM-UFRGS and supported in spirit and mind by my colleagues.
