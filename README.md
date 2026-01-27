# MYC screening paper  

Scripts used in MYC Cas9-SpRY screening paper  
Preprint: [Decoding the MYC locus reveals a druggable ultraconserved RNA element]()  
Feiyue Yang, Jan 27th, 2026  


<pre>
.  
├── data  
│   └── MYC-lib-for-mageck.type.myc2_indel.organized.tsv	# Screening sgRNA library  
├── README.md  
└── script  
    ├── ASO_rnaseq	# ASO RNA-Seq analysis (Figure 5)  
    │   ├── bam_2_bw.sh
    │   ├── deseq2_final.R
    │   ├── enrichment_analysis.R
    │   ├── fastp.sh
    │   ├── featureCounts.sh
    │   ├── gene_TPM.ipynb
    │   ├── GSEA_barplot.ipynb
    │   └── star.sh
    ├── coding_vs_noncoding	# Coding / Non-coding comparison (Figure 2)  
    │   ├── coding_vs_noncoding_whole_genome_conservation.myc.ipynb
    │   ├── coding_vs_noncoding_whole_genome_conservation.py
    │   └── grna_lfc_nucleotide_top_select.ipynb
    ├── conservation	# Screening phenotype & conservation (Figure 3)  
    │   ├── conservation_HL.ipynb
    │   ├── conservation_prefilter_RPS19.sh
    │   ├── conservation_prefilter.sh
    │   ├── grna_conservation_nucleotide.ipynb
    │   ├── grna_conservation_nucleotide_sig.ipynb
    │   ├── grna_conservation_PMID_38889719.ipynb
    │   └── jarvis_sig.ipynb
    ├── ENCODE	# Screening phenotype & ENCODE publica data (Figure 2)  
    │   └── ENCODE_cCRE.ipynb
    ├── neg_control_select	# Screening data preprocess: perfect negative controls  
    │   ├── grna_nc_bowtie.sh
    │   └── grna_nc_select_bowtie.ipynb
    ├── off_targets	# Screening pehnotype & off-targets (Figure 1)  
    │   ├── grna_bowtie_process.ipynb
    │   ├── grna_bowtie.sh
    │   └── grna_offtargets_viz.ipynb
    ├── pam	# Screening phenotype & PAM motif (Figure 1)  
    │   └── grna_lfc_pam_enrich.ipynb
    ├── PMID_35561311_STING	# ASO RNA-Seq & Public STING reexpression RNA-Seq (Figure 5)  
    │   ├── ASO_corr.ipynb
    │   ├── deseq2_gsea.R
    │   ├── fastp.sh
    │   ├── featureCounts.sh
    │   └── star.sh
    ├── screening_data_preprocess	# Screening data preprocess with MAGeCK  
    │   ├── mageck_Cas9-HMOI_D0_unified_baseline.sh
    │   ├── mageck_duplicate_sgRNA.ipynb
    │   └── mageck_duplicate_sgRNA.sh
    ├── screening_data_qc	# Screening data QC (Figure 1, Figure 2)  
    │   ├── grna_d8_d20_corr_plot.ipynb
    │   ├── grna_plot_cons_ridge_nucleotide.ipynb
    │   └── grna_plot_lfc_seq_wt_nucleotide.ipynb
    └── sgRNA_on_target	# Screening phenotype & predicted efficiency (Supplementary Figure)  
        └── on_target_LFC.ipynb
</pre>