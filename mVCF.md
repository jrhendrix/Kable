# The Modification Variant Call Format (mVCF) Version 1.0 Specification

May 25, 2023

## The mVCF specification
mVCF is a text file format that follows the format specifications for a general Variant Call Format (VCF) version 4.2. To find more information on the base file format, see the full specification document [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf). 

## Understanding the mVCF haplotype representation
For each locus in the file, haplotype is represented in the REF and ALT fields. The REF field indicates the allele found in the refernce while the ALT field contains a list of alternative alleles found in the other sequences being compared. A general VCF is limited to reporting simple haplotypes of the four canonical bases (A, C, G, and T). 

The mVCF expands the allele alphabet to include base modifications where each type of base modification is given a single letter representation.

## Base modification alphabet

| Modification | Short name | Symbol |
| ------------ | ---------- | ------ |
| 5-methylcytosine | 5mC    | m      |
| 5-hydroxymethylcytosine | 5hmC | h |
| 5-formylcytosine | 5fC    | f      |
| 4-methylcytosine | 4mC    | v      |
| 6-methyladenine | 6mA     | 6      |
