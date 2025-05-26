# Antibiotic-Resistance-Detection-Tool
This tool detects antibiotic resistance genes (ARGs) and pathogenicity islands (PAIs) in bacterial genome sequences. It uses sequence alignment against curated databases to identify ARGs and analyzes genomic features such as mobility genes and tRNA sites to predict PAIs.

Detection Tool: Antibiotic Resistance and PAIs


You can find in here a simple brief on how to interact with the mentioned Detection Tool.


Features
1. Detection of Antibiotic Resistance Genes
2. Gene Presence/Absence across genomes
3. Jaccard distance matrix
4. Most common antibiotic resistance gene across genomes
5. Phylogenetic Tree 
6. Detection of Pathogeinicty Islands
7. Applying sliding window to get highest GC-regions
8. Visualizing whole circular genome.

Requirements
!pip install biopython
!sudo apt-get install ncbi-blast+
!pip install -q condacolab
import condacolab
condacolab.install()
!conda install -c bioconda blast -y
!makeblastdb -in reference.fasta -dbtype nucl -out reference_db
!pip install biopython
!sudo apt-get install ncbi-blast+
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import subprocess
from Bio.Blast.Applications import NcbiblastnCommandline
from matplotlib.patches import Wedge
from Bio import SeqIO

####Note: While importing the above requirements, it takes a few minutes, so please be patient while it runs.
##For running the driver codes, the names of the files are commeneted next to each input cell.
#The output of the tree is the branch length between nodes/genomes, and they are redundant (for an unknown reason). But it doesn't affect the results in any way.


1. For Antibiotic Resistance Genes

The antibiotic resistance part functions by getting two files from the user first file is the CARD genes that is parsed 
and used to search for genes in the second file, which is a file containing the reference genomes (27 genomes) that is used to 
make a local database to run blast alignment. 


2. For PAIs

Tools and Libraries Used
BLAST+: Used for sequence alignment and construction of BLAST databases for genome analysis.
Biopython: It is used for parsing genomes and for interaction with BLAST through the NcbiblastnCommandline class.
Matplotlib: Used to generate circular genome plots to visualize alignments and GC-rich regions.
Subprocess: This module allows you to run entirely independent processes. For example, making the blast db.
Python Standard Libraries: are used mainly for file handling and data management and implementing core functionalities. 

Input Files
1- Genome File: A FASTA file containing the bacterial genome sequence. For this demo, the file MRSA-LUX10.fna is used. If you want to try others just change the genome name in the two needed inputs for genome and go on.
2- Query Files: FASTA files containing sequences of interest. For the demo:
VFDB_setA_nt.fasta: Database of virulence factors.
trna-motifs.fasta: Database of tRNA motifs.
ANTI_resis.fasta: Sequences related to antibiotic resistance genes.

Outputs
1- Alignment Results:
virulence_blast.txt: Results for virulence factors.
trna_blast.txt: Results for tRNA motifs.
antibiotic_blast.txt: Results for antibiotic resistance genes.
2-Top GC Content Regions
3-Filtered Alignments
4-Circular Genome Plot


Contact Info
s-ahmed.eid@zewailcity.edu.eg
s-haitham.mohamed@zewailcity.edu.eg
s-shorouq.elzentahy@zewailcity.edu.eg
