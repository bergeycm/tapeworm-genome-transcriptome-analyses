
# ========================================================================================
# --- Run analyses for phylogenomics and functional overrep portion of tapeworm project
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Variables
# ----------------------------------------------------------------------------------------

# Get IDs of WGS individuals (also some throw away variables with parts of read filenames)
(IDS_WGS, SSTR, LSTR) = glob_wildcards("/gpfs/group/ghp3/default/Bergey/tapeworm_WGS/tapeworms_hiseq_2017-02/GP01272017/{ind_id}_S{S_str}_L{L_str}_R1_001.fastq.gz")

# For now, ignore non-Taenia samples
import re
regex = re.compile('^T')
IDS_WGS = [x for x in IDS_WGS if regex.match(x)]

# Make lists that include reference genome from Tsai et al (coded as T00)
IDS_WGS_WITH_REF = IDS_WGS[:]
IDS_WGS_WITH_REF.append("T00")

PROTEOMES=("echinococcus_multilocularis/echinococcus_multilocularis.PRJEB122.WBPS12")

PROT_HUNKS_ECHINOCOCCUS, = glob_wildcards("genomes/echinococcus_multilocularis/echinococcus_multilocularis.PRJEB122.WBPS12.protein.oneline.fa_{prot_part}")

# Remove run marker files
regex = re.compile('.*atram')
PROT_HUNKS_ECHINOCOCCUS = [x for x in PROT_HUNKS_ECHINOCOCCUS if not regex.match(x)]

GENE_TREE_PARTS=list(map(lambda x: str(x).zfill(2), range(0,24)))

print ("Individual IDs for WGS:")
print (IDS_WGS)

# ----------------------------------------------------------------------------------------
# --- Make all
# ----------------------------------------------------------------------------------------

rule all:
    input:
        # symlink_to_reads:
        expand("data/WGS/{ind}.{read}.fastq.gz", ind=IDS_WGS, read=['R1','R2']),
        # download_ref_genomes:
        "genomes/Tsolium_Mexico_v1/Tsolium_Mexico_v1.fa",
        "genomes/TGD_Tsolium/TGD_Tsolium.fa",
        "genomes/TGD_Tsaginata/TGD_Tsaginata.fa",
        "genomes/TGD_Tasiatica/TGD_Tasiatica.fa",
        # trim_reads:
        expand("data/WGS/{ind}_TRIM.{read}.fastq.gz", ind=IDS_WGS, read=['R1','R2']),
        # remove_dups:
        expand("data/WGS/{ind}_TRIM.{read}.rmdup.fastq.gz", ind=IDS_WGS, read=['R1','R2']),
        # estimate_kmer:
        expand("results/kmer_estimation/{ind}.bestkmer.txt", ind=IDS_WGS),
        # assemble_mt_dna:
        expand("results/MITObim_{ind}/MITObim_log_{ind}.txt", ind=IDS_WGS),
        expand("results/MITObim_{ind}/MITObim_mt_{ind}-mt-final_noIUPAC.fasta", ind=IDS_WGS),
        # concat_mt_dna:
        "results/MITObim_all.fasta",
        # download_mt_genomes:
        "data/Guo2017_tapeworm_mt_genomes.fa",
        # align_mt_genomes:
        "results/mt_genomes.clade_AH.aln.fa",
        "results/mt_genomes.non_clade_AH.aln.fa",
        "results/mt_genes/gene_aln_list.txt",
        # infer_mt_tree:
        expand("results/phylogeny/mt_phylogeny/mt_genes_concat.{ending}",
            ending=["best_scheme.nex","best_scheme","model.gz","mldist","bionj",
                    "best_model.nex","splits.nex","contree","treefile","iqtree", "log"]),
        expand("results/phylogeny/mt_phylogeny/mt_genes_loci.{ending}",
            ending=["parstree", "best_scheme.nex", "best_scheme", "model.gz",
                    "best_model.nex", "treefile", "iqtree", "log"]),
        expand("results/phylogeny/mt_phylogeny/mt_genes_concord.{ending}",
            ending=["cf.tree", "cf.tree.nex", "cf.branch", "cf.stat", "log"]),
        # download_annotations:
        "genomes/annotations/taenia_solium.PRJNA170813.WBPS9.annotations.gff3",
        "genomes/annotations/taenia_solium.PRJNA170813.WBPS9.annotations.gff3.gz",
        "genomes/annotations/Tas.v1.gff3",
        "genomes/annotations/Tsa.v1.gff3",
        # download_ref_proteome_echinococcus:
        expand("genomes/{proteome}.protein{ending}",
            proteome=PROTEOMES, ending=[".fa", ".oneline.fa"]),
        expand("genomes/{proteome}.CDS_transcripts.fa",
            proteome=PROTEOMES),
        # split_proteome:
        expand("genomes/{proteome}.protein.oneline.fa-split_marker",
            proteome=PROTEOMES),
        # build_atram_libs:
        expand("results/atram/atram_libraries/{ind}.sqlite.db", ind=IDS_WGS),
        # run_atram_echinococcus:
        expand("genomes/echinococcus_multilocularis/echinococcus_multilocularis.PRJEB122.WBPS12.protein.oneline.fa_{prot_part}.{ind}.atram_run_marker.txt",
            prot_part=PROT_HUNKS_ECHINOCOCCUS, ind=IDS_WGS),
        # copy_combined_atram_output_echinococcus:
        "results/atram/atram_output_echinococcus/combined/copy_done_marker.txt",
        # do_exon_stitching_echinococcus:
        "exon_stitching_echinococcus/final/summary_stats.per.gene.csv",
        "exon_stitching_echinococcus/final/summary_stats.per.taxon.csv",
        # make_gene_trees_cmd_lists:
        "results/phylogeny/gene_tree_commands/gene_tree_cmds.sh",
        expand("results/phylogeny/gene_tree_commands/gene_tree_cmds_part{gene_tree_part}",
            gene_tree_part=GENE_TREE_PARTS),
        # infer_all_gene_trees:
        expand("results/phylogeny/gene_tree_commands/gene_tree_cmds_part{gene_tree_part}.log",
            gene_tree_part=GENE_TREE_PARTS),
        # concat_loci:
        "results/phylogeny/concat/combined.phy",
        "results/phylogeny/concat/combined.nex",
        # infer_species_tree:
        "results/phylogeny/concat/trees/RAxML_bipartitions.concat_atram_boot.sp.tre",
        # infer_species_tree_iqtree:
        expand("results/phylogeny/concat/wgs_iq_tree/wgs_concat.{ending}",
            ending=["iqtree", "treefile", "mldist", "splits.nex", "contree", "log"]),
        expand("results/phylogeny/concat/wgs_iq_tree/wgs_loci.{ending}",
            ending=["iqtree", "treefile", "log"]),
        expand("results/phylogeny/concat/wgs_iq_tree/wgs_concord.{ending}",
            ending=["cf.tree", "cf.tree.nex", "cf.branch", "cf.stat"]),
        # infer_species_tree_mb:
        "results/phylogeny/concat/combined.mb.nex.t",
        # make_ML_gene_trees_cmd_lists:
        expand("results/phylogeny/gene_tree_commands/gene_tree_ML_cmds_part{gene_tree_part}.cmd_list",
            gene_tree_part=GENE_TREE_PARTS),
        # infer_all_ML_gene_trees:
        expand("results/phylogeny/gene_tree_commands/gene_tree_ML_cmds_part{gene_tree_part}.log",
            gene_tree_part=GENE_TREE_PARTS),
        # compute_brlen_ratios:
        "results/gene_tree_distance_ratios_human_nonhuman.txt",
        # do_go_overrep_brlen:
        "results/gene_tree_distance_ratios_human_nonhuman.tbl.txt",
        # overrep_test_drug_targets:
        "results/human_accel_drug_targets.txt",
        "results/human_accel_drug_targets_known.txt",
        # overrep_test_heat:
        "results/human_accel_heat_DE.allheated.brlen_perc0.99.chisq.txt",
        "results/human_accel_heat_DE.allheated.brlen_perc0.99.hypergeometric.txt",
        "results/human_accel_heat_DE.allheated.brlen_perc0.99.wilcox.txt",
        "results/human_slowed_heat_DE.allheated.brlen_perc0.01.chisq.txt",
        "results/human_slowed_heat_DE.allheated.brlen_perc0.01.hypergeometric.txt",
        "results/human_slowed_heat_DE.allheated.brlen_perc0.01.wilcox.txt",
        "results/human_accel_heat_DE.highestheated.brlen_perc0.99.chisq.txt",
        "results/human_accel_heat_DE.highestheated.brlen_perc0.99.hypergeometric.txt",
        "results/human_accel_heat_DE.highestheated.brlen_perc0.99.wilcox.txt",
        "results/human_slowed_heat_DE.highestheated.brlen_perc0.01.chisq.txt",
        "results/human_slowed_heat_DE.highestheated.brlen_perc0.01.hypergeometric.txt",
        "results/human_slowed_heat_DE.highestheated.brlen_perc0.01.wilcox.txt",
        # count_reads:
        "reports/ind_read_count.txt",

localrules: all

# ========================================================================================
# --- Assemble genomes
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Symlink to reads for WGS
# ----------------------------------------------------------------------------------------

rule symlink_to_reads:
    output:
        expand("data/WGS/{ind}.{read}.fastq.gz", ind=IDS_WGS, read=['R1','R2']),
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/symlink_to_reads.sh;"

# ----------------------------------------------------------------------------------------
# --- Download reference genomes
# ----------------------------------------------------------------------------------------

rule download_ref_genomes:
    output:
        "genomes/Tsolium_Mexico_v1/Tsolium_Mexico_v1.fa",
        "genomes/TGD_Tsolium/TGD_Tsolium.fa",
        "genomes/TGD_Tsaginata/TGD_Tsaginata.fa",
        "genomes/TGD_Tasiatica/TGD_Tasiatica.fa",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/download_ref_genomes.sh;"

# ----------------------------------------------------------------------------------------
# --- Trim reads
# ----------------------------------------------------------------------------------------

rule trim_reads:
    input:
        expand("data/WGS/{{ind}}.{read}.fastq.gz", read=['R1','R2'])
    output:
        expand("data/WGS/{{ind}}_TRIM.{read}.fastq.gz", read=['R1','R2'])
    threads: 8
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/trim_reads.sh {wildcards.ind}"

# ----------------------------------------------------------------------------------------
# --- Remove PCR duplicates
# ----------------------------------------------------------------------------------------

rule remove_dups:
    input:
        expand("data/WGS/{{ind}}_TRIM.{read}.fastq.gz", read=['R1','R2'])
    output:
        expand("data/WGS/{{ind}}_TRIM.{read}.rmdup.fastq.gz", read=['R1','R2'])
    threads: 1
    params: runtime="96",
            mem=",mem=24gb"
    shell:
        "sh scripts/remove_pcr_dups.sh {wildcards.ind}"

# ----------------------------------------------------------------------------------------
# --- Estimate ideal kmer size
# ----------------------------------------------------------------------------------------

rule estimate_kmer:
    input:
        expand("data/WGS/{{ind}}_TRIM.{read}.fastq.gz", read=['R1','R2'])
    output:
        "results/kmer_estimation/{ind}.bestkmer.txt"
    threads: 8
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/estimate_kmer.sh {wildcards.ind}"

# ----------------------------------------------------------------------------------------
# --- Assemble mtDNA
# ----------------------------------------------------------------------------------------

rule assemble_mt_dna:
    input:
        expand("data/WGS/{{ind}}_TRIM.{read}.fastq.gz", read=['R1','R2'])
    output:
        "results/MITObim_{ind}/MITObim_log_{ind}.txt",
        "results/MITObim_{ind}/MITObim_mt_{ind}-mt-final_noIUPAC.fasta"
    threads: 8
    params: runtime="24",
            mem=",mem=96gb"
    shell:
        "sh scripts/assemble_mt_genomes.sh {wildcards.ind}"

# ----------------------------------------------------------------------------------------
# --- Concatenate mt genomes
# ----------------------------------------------------------------------------------------

rule concat_mt_dna:
    input:
        expand("results/MITObim_{ind}/MITObim_mt_{ind}-mt-final_noIUPAC.fasta", ind=IDS_WGS)
    output:
        "results/MITObim_all.fasta",
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "cat {input} > {output}"

# ----------------------------------------------------------------------------------------
# --- Download reference mt genomes
# ----------------------------------------------------------------------------------------

rule download_mt_genomes:
    input:
        "data/Guo2017_S2_tapeworm_mt_genome_accessions.txt",
    output:
        "data/Guo2017_tapeworm_mt_genomes.fa",
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "Rscript scripts/download_mt_genomes.R"

# ----------------------------------------------------------------------------------------
# --- Align Guo 2017 and our mt genomes
# ----------------------------------------------------------------------------------------

rule align_mt_genomes:
    input:
        "data/Guo2017_tapeworm_mt_genomes.fa",
        "results/MITObim_all.fasta",
    output:
        "results/mt_genomes.clade_AH.aln.fa",
        "results/mt_genomes.non_clade_AH.aln.fa",
        "results/mt_genes/gene_aln_list.txt"
    threads: 1
    params: runtime="48",
            mem=",mem=25gb"
    shell:
        "sh scripts/align_mt.sh"

# ----------------------------------------------------------------------------------------
# --- Infer mitochondrial tree
# ----------------------------------------------------------------------------------------

rule infer_mt_tree:
    input:
        "results/mt_genes/gene_aln_list.txt"
    output:
        expand("results/phylogeny/mt_phylogeny/mt_genes_concat.{ending}",
            ending=["best_scheme.nex","best_scheme","model.gz","mldist","bionj",
                    "best_model.nex","splits.nex","contree","treefile","iqtree", "log"]),
        expand("results/phylogeny/mt_phylogeny/mt_genes_loci.{ending}",
            ending=["parstree", "best_scheme.nex", "best_scheme", "model.gz",
                    "best_model.nex", "treefile", "iqtree", "log"]),
        expand("results/phylogeny/mt_phylogeny/mt_genes_concord.{ending}",
            ending=["cf.tree", "cf.tree.nex", "cf.branch", "cf.stat", "log"]),
    threads: 8
    params: runtime="24",
            mem=",mem=5gb"
    shell:
        "sh scripts/infer_mt_tree.sh"

# ----------------------------------------------------------------------------------------
# --- Make mt phylogeny figure
# ----------------------------------------------------------------------------------------

       

# ----------------------------------------------------------------------------------------
# --- Download annotations
# ----------------------------------------------------------------------------------------

rule download_annotations:
    output:
        "genomes/annotations/taenia_solium.PRJNA170813.WBPS9.annotations.gff3",
        "genomes/annotations/taenia_solium.PRJNA170813.WBPS9.annotations.gff3.gz",
        "genomes/annotations/Tas.v1.gff3",
        "genomes/annotations/Tsa.v1.gff3",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/download_annotations.sh"

# ----------------------------------------------------------------------------------------
# --- Download reference proteomes
# ----------------------------------------------------------------------------------------

rule download_ref_proteome_echinococcus:
    output:
        "genomes/echinococcus_multilocularis/echinococcus_multilocularis.PRJEB122.WBPS12.protein.fa",
        "genomes/echinococcus_multilocularis/echinococcus_multilocularis.PRJEB122.WBPS12.protein.oneline.fa",
        "genomes/echinococcus_multilocularis/echinococcus_multilocularis.PRJEB122.WBPS12.CDS_transcripts.fa",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/download_proteomes.sh echinococcus_multilocularis"

# ----------------------------------------------------------------------------------------
# --- Split proteome into hunks to be run through aTRAM in parallel
# ----------------------------------------------------------------------------------------

rule split_proteome:
    input:
        "genomes/{proteome}.protein.oneline.fa"
    output:
        "genomes/{proteome}.protein.oneline.fa-split_marker"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/split_proteomes.sh {input}"

# ----------------------------------------------------------------------------------------
# --- Build aTRAM libraries
# ----------------------------------------------------------------------------------------

rule build_atram_libs:
    input:
        expand("data/WGS/{{ind}}_TRIM.{read}.rmdup.fastq.gz", read=['R1','R2'])
    output:
        "results/atram/atram_libraries/{ind}.sqlite.db"
    threads: 8
    params: runtime="24",
            mem=",mem=5gb"
    shell:
        "sh scripts/build_atram_libraries.sh {wildcards.ind}"

# ----------------------------------------------------------------------------------------
# --- Run aTRAM
# ----------------------------------------------------------------------------------------

rule run_atram_echinococcus:
    input:
        split_marker="genomes/echinococcus_multilocularis/echinococcus_multilocularis.PRJEB122.WBPS12.protein.oneline.fa-split_marker",
        prot_fasta="genomes/echinococcus_multilocularis/echinococcus_multilocularis.PRJEB122.WBPS12.protein.oneline.fa_{prot_part}"
    output:
        "genomes/echinococcus_multilocularis/echinococcus_multilocularis.PRJEB122.WBPS12.protein.oneline.fa_{prot_part}.{ind}.atram_run_marker.txt"
    threads: 8
    params: runtime="24",
            mem=",mem=5gb"
    shell:
        #"rm -f results/atram/atram_output_echinococcus/combined/copy_done_marker.txt; "
        "sh scripts/run_atram.sh {wildcards.ind} {input.prot_fasta} echinococcus"

# ----------------------------------------------------------------------------------------
# --- Copy combined aTRAM results
# ----------------------------------------------------------------------------------------

rule copy_combined_atram_output_echinococcus:
    input:
        expand("genomes/echinococcus_multilocularis/echinococcus_multilocularis.PRJEB122.WBPS12.protein.oneline.fa_{prot_part}.{ind}.atram_run_marker.txt",
            ind=IDS_WGS, prot_part=PROT_HUNKS_ECHINOCOCCUS)
    output:
        "results/atram/atram_output_echinococcus/combined/copy_done_marker.txt",
    threads: 1
    params: runtime="8",
            mem=",mem=5gb"
    shell:
        "sh scripts/copy_combined_atram_output.sh echinococcus"

# ----------------------------------------------------------------------------------------
# --- Do exon stitching
# ----------------------------------------------------------------------------------------

rule do_exon_stitching_echinococcus:
    input:
        cp_marker="results/atram/atram_output_echinococcus/combined/copy_done_marker.txt",
        prot_fasta="genomes/echinococcus_multilocularis/echinococcus_multilocularis.PRJEB122.WBPS12.protein.oneline.fa"
    output:
        "exon_stitching_echinococcus/final/summary_stats.per.gene.csv",
        "exon_stitching_echinococcus/final/summary_stats.per.taxon.csv",
    threads: 20
    params: runtime="48",
            mem=",mem=5gb"
    shell:
        "sh scripts/do_exon_stitching.sh echinococcus {input.prot_fasta}"

# ----------------------------------------------------------------------------------------
# --- Make lists of commands for alignment and Bayesian gene tree inference in parallel
# ----------------------------------------------------------------------------------------

rule make_gene_trees_cmd_lists:
    input:
        "exon_stitching_echinococcus/final/summary_stats.per.gene.csv",
        "exon_stitching_echinococcus/final/summary_stats.per.taxon.csv",
    output:
        "results/phylogeny/gene_tree_commands/gene_tree_cmds.sh",
        expand("results/phylogeny/gene_tree_commands/gene_tree_cmds_part{gene_tree_part}",
            gene_tree_part=GENE_TREE_PARTS)
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/make_gene_tree_cmd_lists.sh"

# ----------------------------------------------------------------------------------------
# --- Infer all gene trees with Bayesian methods (and do alignment along the way)
# ----------------------------------------------------------------------------------------

rule infer_all_gene_trees:
    input:
        cmd_list="results/phylogeny/gene_tree_commands/gene_tree_cmds_part{gene_tree_part}"
    output:
        "results/phylogeny/gene_tree_commands/gene_tree_cmds_part{gene_tree_part}.log"
    threads: 20
    params: runtime="30",
            mem=",mem=5gb"
    shell:
        "module load parallel/20170422;"
        "echo 'Beginning parallel processing of commands in {input.cmd_list}.' > {output};"
        # Run parallel, tolerating failures
        "cat {input.cmd_list} | parallel -j 32 || true;"
        "echo 'Finished processing commands in parallel.' >> {output};"
        "echo `date` >> {output}"

# ----------------------------------------------------------------------------------------
# --- Concatenate all loci
# ----------------------------------------------------------------------------------------

rule concat_loci:
    input:
        expand("results/phylogeny/gene_tree_commands/gene_tree_cmds_part{gene_tree_part}.log",
            gene_tree_part=GENE_TREE_PARTS)
    output:
        "results/phylogeny/concat/combined.phy",
        "results/phylogeny/concat/combined.nex"
    threads: 1
    params: runtime="48",
            mem=",mem=5gb"
    shell:
        "sh scripts/concatenate_loci.sh"

# ----------------------------------------------------------------------------------------
# --- Infer species tree with RAxML
# ----------------------------------------------------------------------------------------

rule infer_species_tree:
    input:
        "results/phylogeny/concat/combined.phy"
    output:
        "results/phylogeny/concat/trees/RAxML_bipartitions.concat_atram_boot.sp.tre"
    threads: 1
    params: runtime="48",
            mem=",mem=48gb"
    shell:
        "sh scripts/infer_species_tree.sh"

# ----------------------------------------------------------------------------------------
# --- Infer species tree with IQ-TREE
# ----------------------------------------------------------------------------------------

rule infer_species_tree_iqtree:
    input:
        expand("results/phylogeny/gene_tree_commands/gene_tree_cmds_part{gene_tree_part}.log",
            gene_tree_part=GENE_TREE_PARTS)
    output:
        expand("results/phylogeny/concat/wgs_iq_tree/wgs_concat.{ending}",
            ending=["iqtree", "treefile", "mldist", "splits.nex", "contree", "log"]),
        expand("results/phylogeny/concat/wgs_iq_tree/wgs_loci.{ending}",
            ending=["iqtree", "treefile", "log"]),
        expand("results/phylogeny/concat/wgs_iq_tree/wgs_concord.{ending}",
            ending=["cf.tree", "cf.tree.nex", "cf.branch", "cf.stat"])
    threads: 1
    params: runtime="96",
            mem=",mem=48gb"
    shell:
        "sh scripts/infer_species_tree_iqtree.sh"

# ----------------------------------------------------------------------------------------
# --- Infer species tree with MrBayes
# ----------------------------------------------------------------------------------------

rule infer_species_tree_mb:
    input:
        "results/phylogeny/concat/combined.nex"
    output:
        "results/phylogeny/concat/combined.mb.nex.t"
    threads: 1
    params: runtime="96",
            mem=",mem=48gb"
    shell:
        "sh scripts/infer_species_tree_mrbayes.sh"

# ----------------------------------------------------------------------------------------
# --- Make WGS phylogeny figure
# ----------------------------------------------------------------------------------------

       

# ----------------------------------------------------------------------------------------
# --- Make lists of commands for ML gene tree inference in parallel
# ----------------------------------------------------------------------------------------

rule make_ML_gene_trees_cmd_lists:
    input:
        "results/phylogeny/gene_tree_commands/gene_tree_cmds_part{gene_tree_part}"
    output:
        "results/phylogeny/gene_tree_commands/gene_tree_ML_cmds_part{gene_tree_part}.cmd_list"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sed -e \"s/align_and_infer/infer_ML/\" {input} > {output}"

# ----------------------------------------------------------------------------------------
# --- Infer all gene trees with Bayesian methods (and do alignment along the way)
# ----------------------------------------------------------------------------------------

rule infer_all_ML_gene_trees:
    input:
        cmd_list="results/phylogeny/gene_tree_commands/gene_tree_ML_cmds_part{gene_tree_part}.cmd_list"
    output:
        "results/phylogeny/gene_tree_commands/gene_tree_ML_cmds_part{gene_tree_part}.log"
    threads: 20
    params: runtime="30",
            mem=",mem=5gb"
    shell:
        "module load parallel/20170422;"
        "echo 'Beginning parallel processing of commands in {input.cmd_list}.' > {output};"
        # Run parallel, tolerating failures
        "cat {input.cmd_list} | parallel -j 32 || true;"
        "echo 'Finished processing commands in parallel.' >> {output};"
        "echo `date` >> {output}"

# ----------------------------------------------------------------------------------------
# --- Compute human/non-human branch length ratio for all gene trees
# ----------------------------------------------------------------------------------------

rule compute_brlen_ratios:
    input:
        expand("results/phylogeny/gene_tree_commands/gene_tree_cmds_part{gene_tree_part}.log",
            gene_tree_part=GENE_TREE_PARTS)
    output:
        "results/gene_tree_distance_ratios_human_nonhuman.txt",
    threads: 1
    params: runtime="48",
            mem=",mem=5gb"
    shell:
        "sh scripts/compute_all_brlen_ratios.sh"

# ----------------------------------------------------------------------------------------
# --- Do GO overrep tests for accelerated (or decelerated) ev genes
# ----------------------------------------------------------------------------------------

rule do_go_overrep_brlen:
    input:
        "results/gene_tree_distance_ratios_human_nonhuman.txt",
    output:
        "results/gene_tree_distance_ratios_human_nonhuman.tbl.txt",
        "reports/outlier_high_RHR.gProfiler_results.csv",
        "reports/outlier_low_RHR.gProfiler_results.csv"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "Rscript scripts/do_GO_enrichment_for_RHR_genes.R"

# ----------------------------------------------------------------------------------------
# --- Test if accelerated ev genes are overrep'd for predicted or known drug targets
# ----------------------------------------------------------------------------------------

rule overrep_test_drug_targets:
    input:
        "results/gene_tree_distance_ratios_human_nonhuman.txt",
    output:
        "results/human_accel_drug_targets.txt",
        "results/human_accel_drug_targets_known.txt",
        "reports/outlier_high_RHR_drug_targets.txt",
        "reports/outlier_low_RHR_drug_targets.txt",
        "reports/outlier_high_RHR_drug_targets.gProfiler_results.csv",
        "reports/outlier_low_RHR_drug_targets.gProfiler_results.csv"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "Rscript scripts/overrep_test_accelerated_pred_drug_targets.R"

# ----------------------------------------------------------------------------------------
# --- Test if accelerated ev genes are overrep'd for heat DE genes
# ----------------------------------------------------------------------------------------

rule overrep_test_heat:
    input:
        "results/gene_tree_distance_ratios_human_nonhuman.txt",
    output:
        "results/human_accel_heat_DE.allheated.brlen_perc0.99.chisq.txt",
        "results/human_accel_heat_DE.allheated.brlen_perc0.99.hypergeometric.txt",
        "results/human_accel_heat_DE.allheated.brlen_perc0.99.wilcox.txt",
        "results/human_slowed_heat_DE.allheated.brlen_perc0.01.chisq.txt",
        "results/human_slowed_heat_DE.allheated.brlen_perc0.01.hypergeometric.txt",
        "results/human_slowed_heat_DE.allheated.brlen_perc0.01.wilcox.txt",
        "results/human_accel_heat_DE.highestheated.brlen_perc0.99.chisq.txt",
        "results/human_accel_heat_DE.highestheated.brlen_perc0.99.hypergeometric.txt",
        "results/human_accel_heat_DE.highestheated.brlen_perc0.99.wilcox.txt",
        "results/human_slowed_heat_DE.highestheated.brlen_perc0.01.chisq.txt",
        "results/human_slowed_heat_DE.highestheated.brlen_perc0.01.hypergeometric.txt",
        "results/human_slowed_heat_DE.highestheated.brlen_perc0.01.wilcox.txt",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "Rscript scripts/overrep_test_accelerated_heat.R"

# ----------------------------------------------------------------------------------------
# --- Count up reads
# ----------------------------------------------------------------------------------------

rule count_reads:
    input:
        expand("data/WGS/{ind}.{read}.fastq.gz", ind=IDS_WGS, read=['R1','R2']),
    output:
        "reports/ind_read_count.txt"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/compute_read_count.sh > {output}"

# ========================================================================================
# --- Set default for optional mem parameter
# ========================================================================================

default = ""
name = 'mem'
for r in workflow.rules:
    try:
        getattr(r.params, name)
    except AttributeError:
        r.params.append(default)
        r.params.add_name(name)
