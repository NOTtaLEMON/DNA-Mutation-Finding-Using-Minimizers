# Project Map

## Purpose
Find mutations (differences from a reference genome) from short noisy reads, in C, no ML. Source: PROJECT_REPORT.md (status 25 Oct 2026; nothing implemented in C yet).

## Requirements
- Nine stages, each using a course data structure (packed array, two-stack queue, hash table, Bloom filter, BST, min-heap, DP grid, counting arrays, sorted array).
- Project track uses planted substitutions (option B); paper track uses a real Indian SARS-CoV-2 genome (option C).
- Paper explanations in this session relate papers to the project (e.g. paper 4 on minimizers, stage 2).

## Components
- PROJECT_REPORT.md: project description and decisions.
- paper 4.pdf: Verbeek et al., "Minimizing the Minimizers via Alphabet Reordering".
- Other files in folder: explanation PDFs, pitch pptx, sequence.fasta, Datasets_and_Literature.pdf (not yet read in this session).

## Main Flow
Load -> Pick words (minimizers) -> Index -> Pre-check -> Vote -> Short-list -> Verify -> Pile up -> Score (from PROJECT_REPORT.md; not yet implemented).

## Data and Trust Boundaries
Reference NC_045512.2 proposed; Indian genomes not yet chosen. Unverified.

## Build and Deployment
None yet; no C code exists.

## Unknowns
- Full project deadline; N-run window rule; which paper sections matter most to the project.
