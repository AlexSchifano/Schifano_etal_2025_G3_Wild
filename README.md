# Schifano_etal_2025_G3_Wild
Scripts used for the data production (from DNA sequences to mutation table and in-between steps) and for data analysis (statistical tests, figures and necessary data formatting) of the article:
Schifano Alexandre, Bergthorsson Ulfar, Katju Vaishali, 2025. Mutational Biases and Selection in Mitochondrial Genomes: Insights from a Comparative Analysis of Natural and Laboratory Populations of Caenorhabditis elegans. https://doi.org/10.1101/2025.09.03.674070

### DATA PRODUCTION SCRIPTS ###
Scripts used for data production encapsulate all steps used in the study from acquiring data from the CaeNDR database to calling variants, annotate them and adjusting positional information.
They are numbered in order of use.
The file architecture is not created by the scripts.

### DATA ANALYSIS SCRIPTS ###
Those scripts are essentially in R and include all analysis done on the dataset produced by the set of "Data production scripts".

# Details for each script #
## DATA PRODUCTION SCRIPTS ##

*S0_Pipeline.sh*  --> Contains a pipeline executing scripts in order.

*S1_get_BAMS_MtDNA_CeNDR.sh* --> Downloads the sequences of a previously made isotype.list (at this stage already we remove the sequences we think to have signatures of contamination, the list is in the script. It then extracts only the reads mapping to the mtDNA and creates a consensus sequence.

*S2_merge_fasta.sh* --> Merges the consensus sequences into one fasta file. Done twice, once including only the wild isolates and once including C. briggsae and C. tropicalis in order to root the tree.

*S3_align_mafft.sh* --> Aligns both multiple sequence fasta files using MAFFT (https://doi.org/10.1093/molbev/mst010).

*S4_blast_REF_trim.sh* + *R4_cut_MSA.R* --> These two scripts in conjunction are used to cut the MSA file in the right places. It works by first blasting the reference sequences of atp-6, tRNA-K and the last tRNA. The R script then cuts between atp-6 and tRNA-K and everything after the last tRNA.

*S5_CD-HIT.sh* --> Using CD-HIT (https://doi.org/10.1093/bioinformatics/17.3.282) to cluster identical mtDNA sequences into one sequence (using the name of one of the lines clustered that way).

*S6_IQ-TREE.sh* --> Using IQ-TREE (https://doi.org/10.1093/molbev/msaa015), computes the phylogenies based on both MSA. The tree without the other species is rooted with the most basal cluster from the tree rooted with the sister species. Only one line from this cluster is selected as an outgroup.

*R7_root_tree.R* --> Adds a propper root to the tree file for use later

*S7_ARPIP.sh* --> Using ARPIP (https://doi.org/10.1093/sysbio/syac050) to reconstruc ancestral sequences with indels along the tree

*R7.5_Mod_Node_rel.R* --> ARPIP outputs a *Node_rel.txt* table with the relations between the sequences and the nodes (a sort of tree in a table format). This R script is used to take this as an input and re-format it in order to create 3 tables. First, *Table_Node.txt* listing each node and its next node ancestor. Second, *Table_VC.txt* listing all mitotypes and there closest node ancestor. And finally, *Table_VC_ALL.txt* listing all isotypes (mitotypes + lines clustered in the mitotypes) and there closest node ancestor.

*S8_count_mut_tree.sh* + *R8_count_mut_tree.R* --> Uses the *Table_Node.txt* file in order to count the mutations between the nodes along the tree. Also adds mutation details such as codon position, synonymous or nonsynonymous, site degeneracy, radical or conservative, transition or transversion, gene id etc.

*S8.5_blast_GTF.sh* + *R8.5_blast2gtf.sh* --> Blasts the reference mtDNA against the sequence of CB4932 (mitotype which includes N2) from *Anc.fasta*, this file is a MSA outputed by ARPIP which includes mtDNA sequences of all mitotypes and of all reconstructed ancestors. Due to the reconstruction of indels along with SNPs (on top of the trimming steps), the gene positions might not match from the reference positions. These scripts create a new gtf file with gene posittions after blastn.

*S8.9_ref_ASR.sh* --> Blasts the region surrounding each previously identified mutation against the reference in order to provide positions matching the reference positions (might seem convoluted with the previous step but the above step is necessary to account for frameshift mutations and other indels while this process here is necessary in order to compare mutations between lines and present the results in the necessary framework of the reference sequence).

*S9_Variant_Calling.sh* --> Script using SAMtools (https://doi.org/10.1093/bioinformatics/btp352), GATK (Van der Auwera GA, O’Connor BD. 2020. Genomics in the cloud: using Docker, GATK, and WDL in Terra. O’Reilly Media, Inc.), BCFtools (https://doi.org/10.1093/bioinformatics/btr509), BWA (https://doi.org/10.1093/bioinformatics/btp324) and Seqkit (https://doi.org/10.1371/journal.pone.0163962). Calls variants using the bam files downloaded in *S1_get_BAMS_MtDNA_CeNDR.sh* by reading the *Table_VC_ALL.txt* file (which lists all lines and there corresponding direct ancestor). For each pair of Line + Ancestor, the script does: (1) Extract the Ancestor's reconstructed sequence to serve as a reference sequence to call mutations from; (2) Trimms the BAM file to remove the D-loop; (3) Converts the BAM into FASTQ; (4) Aligns the FASTQ against the ancestor's sequence; (5) Mark duplicates; (6) Call variants; (7) Add VAF and TYPE info to the VCF.

*S9.6_gtf_ref.sh* + *R9.6_blast2gtf.R* --> Makes a new GTF file for each Ancestor (named REF in the script). This is due to the INDEL reconstruction, leading to different gtf after removing gaps from the alignement.

*S9.5_vcf_results.sh* + *R9.5_vcf_results.R* --> Uses the VCF file of a strain and GTF file of its ancestor to add mutation details. It also formats the VCF results into a table usuable in R.

*S9.9_ref_pos.sh* --> Similar to *S8.5_blast_GTF.sh* to add positions in the reference (N2) perspective.

*R9.9_Annotation.R* --> Adds mutation details such as codon position, synonymous or nonsynonymous, site degeneracy, radical or conservative, transition or transversion, gene id etc.


## DATA ANALYSIS SCRIPTS ##















