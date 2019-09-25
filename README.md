This pipeline maps 10X Genomics reads to selected genome (by using Tophat2, Hisat2 or STAR), measures digital gene expression (by using ESAT) and finally creates UMI distributions table for expression analysis. 

Steps:
  1. Umitools (https://github.com/CGATOxford/UMI-tools) is used to extract valid reads by checking cell barcode list. By default, interested sequence is found on only the R2 side of the reads,  therefore only these files are send to split process.
  2. Splitted fastq files are aligned by Tophat2, STAR or HISAT2 and converted to bam files.
  3. Mapped reads are merged by samtools.
  4. Merged bams are sorted and indexed by samtools. 
  5. Python script (countUniqueAlignedBarcodes_fromFile.py) used to count reads aligned to a single cell for filtering purposes. It uses the cutoff_for_reads_per_cell to print valid cells.(eg. 3000 reads per cell). 
  6. Python script (filter_lowCountBC_bam_print.py) removes cells from bam files that are below cutoff value (cutoff_for_reads_per_cell) (eg. 3000 reads per cell). 
  7. ESAT(http://garberlab.umassmed.edu/software/esat/) create UMI distributions table.
  8. Python script (cleanLowEndUmis.py) merges low count UMIs with high count UMIs and created output file (*_umiClean.txt). 

Outputs:

*UMI table: The output file (*_umiClean.txt) is tab separated gene/transcript vs cell_Barcode matrix filled with count data as shown at the example below.

    | gene  | ATCAATCGCGAACCGA | ACCCTCAACTCAAACA | ACTCATACCCGGAAAT |
    |-------|------------------|------------------|------------------|
    | RNF14 | 0                | 0                | 0                |
    | MZT2B | 0                | 12               | 0                |
    | SPN   | 0                | 2                | 8                |

