# Project Report: Minimizer-Based Mutation Finder

Status as of 25 Oct 2026. Nothing is implemented in C yet. This file is written for a teammate catching up.

## 1. What we are building

A small tool, written in C with no machine learning, that finds where a virus's genetic text differs from a known correct copy.

- **Input:** a known reference genome (a long text over the letters A, C, G, T) and many short, shuffled, slightly erroneous pieces of the patient's virus (called reads).
- **Output:** the positions where the patient's virus differs from the reference (mutations), with sequencing errors (typos made by the machine) filtered out.
- **Why it suits the course:** every stage uses a classic data structure from Data Structures and Applications.
- **On top of the tool:** a benchmarking study (how typo rate and settings affect accuracy, speed and memory), to be written up as a conference paper later. We do not claim any algorithm is new.

## 2. The two tracks

| Track | Sample data | Needed by |
|---|---|---|
| **Project submission (option B)** | Our program copies the reference, plants known substitutions, and makes reads with injected typos. The answer key is exact because we planted it. | Partial output 4 Nov; full deadline not known yet |
| **Paper (option C), later** | A real Indian SARS-CoV-2 genome is the sample. Its real differences from the reference are the answer key. | After the project |

Also due next week: a literature survey and a dataset list.

## 3. The nine stages (data structure in brackets)

1. **Load:** read the reference into a compact 2-bit array [packed array].
2. **Pick words:** cut the sequence into short words and keep a small, repeatable subset called minimizers [queue built from two stacks].
3. **Index:** table from each kept word to the reference positions where it appears [hash table with linked lists].
4. **Pre-check:** quickly skip words that cannot be in the reference [Bloom filter]. Speed only.
5. **Vote:** each read word suggests a start position on the reference; count the suggestions [binary search tree, required].
6. **Short-list:** keep the best 3 candidate starts [min-heap]. Speed only.
7. **Verify:** compare all letters of the read at a candidate start and accept if close enough [dynamic-programming grid, edit distance].
8. **Pile up:** after all reads are placed, count the letters at every position and call a mutation where most reads agree [counting arrays].
9. **Score:** compare calls with the answer key [sorted array and binary search]. Test data only.

A JavaScript reference version of the toy example exists (`minimizer_reference_core.js`). The C version must reproduce its numbers: 53 words, 24 kept pairs, 21 index entries, 12 of 12 reads placed, 2 true positives, 0 false positives, 0 false negatives.

## 4. Decisions made since the original plan

- Project uses option B; paper uses option C.
- For option C, stage 8 will be extended to use the alignment path from stage 7 (the traceback through the grid), so that deletions and insertions do not shift positions.
- Any word that contains an N (even one letter) is skipped. We believe existing mappers do the same, but this is **not verified** and must be checked before the paper says so.
- Stage 1 reads the FASTA file one character at a time and writes straight into the 2-bit array. No intermediate string is built.
- N handling in stage 1: the 2-bit array has one cell per position, with a placeholder code (for example `00`) under each N, so a position's index equals its position. A separate sorted list stores each N run as (start, length). The list is always checked before an array letter is used.

## 5. Terms in plain words

- **Reference / sample:** the known correct genome, and the patient's genome.
- **Read:** a short piece (100 to 300 letters in real data, 20 in our toy) of the sample, in random order, with a few wrong letters.
- **Mutation:** a position where the sample's real letter differs from the reference.
- **Sequencing error (typo):** a wrong letter the machine wrote into a read. Must not be reported as a mutation.
- **Coverage / depth:** how many reads cover a position. High depth lets us outvote typos.
- **Word (k-mer):** a short fixed-length piece of text, cut by sliding one letter at a time. Toy length 4, real about 15.
- **Window:** a fixed number of consecutive words. We keep the one with the smallest scrambled number in each window.
- **Minimizer:** the kept word of a window.
- **Scrambled number (hash):** the word's code run through a fixed arithmetic recipe so that "smallest" is spread evenly. It does not reorder words.
- **Vote vs pile up:** the vote (stage 5) decides where one read belongs. The pile up (stage 8) decides which letter belongs at one position using all reads. These are different.
- **2-bit packing:** A, C, G, T are stored as 00, 01, 10, 11, so 16 letters fit in a 32-bit integer.
- **Substitution:** one letter replaced by another. This is what option B plants.
- **Indel:** an insertion or deletion. Letters are added or missing. After one, every later letter sits shifted against the reference, so stage 8's rule "read letter i sits at start + i" breaks unless the alignment path is used.
- **N region:** a stretch where the sequencer could not decide the letter and wrote `N` for "unknown". Common in real genomes, and `N` has no 2-bit code.
- **FASTA:** a plain text format for storing sequences. A header line starts with `>` and holds the name, then the sequence follows, split across lines. The line breaks are not part of the sequence. The reference and the Indian genomes come in this format.
- **Minimap2:** an existing, widely used program that finds where one sequence aligns to another. We use it only in the paper phase: to build the answer key from a real genome, and as a baseline to compare against. It is not part of our tool.
- **SAM / CIGAR:** SAM is the text format aligners write. CIGAR is a string inside it such as `50M2D30M` (50 aligned, 2 deleted, 30 aligned). SAM counts positions from 1; our code counts from 0, so conversions need care.

## 6. Datasets and tools

| What | Link | Status |
|---|---|---|
| Reference genome NC_045512.2 (SARS-CoV-2 Wuhan-Hu-1, about 29.9 kb) | https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2 | Proposed. Confirm accession and length on NCBI. |
| Indian genomes (paper phase) | https://www.ncbi.nlm.nih.gov/labs/virus/vssi/ then filter Geographic Region to India | No specific genome chosen yet. Needs accession numbers. |
| GISAID (most Indian sequences) | https://gisaid.org | Needs registration and data terms. |
| Minimap2 | https://github.com/lh3/minimap2 | From memory; verify. |

## 7. Existing work (unverified, from background knowledge)

Minimap2 uses minimizer seeds. BWA-MEM and Bowtie2 use an FM-index. GATK, bcftools, FreeBayes, LoFreq and iVar call variants from piled-up reads. DeepVariant and Clair3 use neural networks, so they are outside our no-ML rule. Industry tools handle indels by gapped alignment, reassembly and left-alignment, and handle N by treating it as missing information. We found no IEEE paper from 2024 or later on this topic through our search tool, so the IEEE literature search must be done directly on IEEE Xplore.

## 8. Open items

- Rule for when a window overlaps an N run (Anand is working this out).
- Verify that mappers skip N-containing words.
- Pick the reference accession and any Indian genomes.
- Full project deadline.
- Who writes which part of stage 1.
- Whether a faculty mentor reviews the paper draft.
