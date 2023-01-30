from pathlib import Path

#export AUGUR_RECURSION_LIMIT=10000

# input
input_fasta = Path("data/tc2_bvic_na.fasta")
input_metadata = Path("data/tc2_bvic.tsv")
lineage = 'vic' #vic h1n1pdm h3n2
segment = 'na' #ha na


#### Do not touch anything below this line ####

reference_name = {
    "h3n2"    : "3C.2a1b.1_A/Wisconsin/327/2017",
    "h1n1pdm" : "6B.1A.2_A/Townsville/21/2018", 
    "vic"     : "V1A_B/Brisbane/60/2008"}

color_taxa = 'month_abbr'
color_tips_shapes = 'type'

subtype_name = {
    'h3n2na' : 'n2',
    'h3n2ha' : 'h3',
    'h1n1pdmha' : 'h1',
    'h1n1pdmna' : 'n1',
    'vicha' : 'bvic_ha',
    'vicna' : 'bvic_na'
}

rootby = reference_name[lineage]

# output file names
nexus_out = "results" + "/" + input_fasta.stem + ".tree"
auspice_json = "results" + "/" + input_fasta.stem + ".json"
ggtree_out = "results" + "/" + input_fasta.stem + ".pdf"
alignment_out = "results" + "/" + input_fasta.stem + ".fasta"
clade_data = "data" + "/" + segment + "_clades.json"

# references
# must be nextstrain references!
reference = "config/reference/reference_%s_%s.gb" % (lineage, segment)
clades = "config/clades/clades_%s_%s.mod.tsv" % (lineage, segment)
auspice_config = "config/auspice/auspice_config_%s.json" % (lineage)

#config - not used
dropped_strains = "config/bvic_meta_19Aug21.tsv"
colors = "config/auspice/colors.tsv"
lat_longs = "config/auspice/lat_longs.tsv"


def clock_rate(lineage, segment):
    # these rates are from 12y runs on 2019-10-18   
    rate = {
     ('h1n1pdm', 'ha'): 0.00329,
     ('h1n1pdm', 'na'): 0.00326,
     ('h1n1pdm', 'np'): 0.00221,
     ('h1n1pdm', 'pa'): 0.00217,
     ('h1n1pdm', 'pb1'): 0.00205,
     ('h1n1pdm', 'pb2'): 0.00277,
     ('h3n2', 'ha'): 0.00382,
     ('h3n2', 'na'): 0.00267,
     ('h3n2', 'np'): 0.00157,
     ('h3n2', 'pa'): 0.00178,
     ('h3n2', 'pb1'): 0.00139,
     ('h3n2', 'pb2'): 0.00218,
     ('vic', 'ha'): 0.00145,
     ('vic', 'na'): 0.00133,
     ('vic', 'np'): 0.00132,
     ('vic', 'pa'): 0.00178,
     ('vic', 'pb1'): 0.00114,
     ('vic', 'pb2'): 0.00106,
     ('yam', 'ha'): 0.00176,
     ('yam', 'na'): 0.00177,
     ('yam', 'np'): 0.00133,
     ('yam', 'pa'): 0.00112,
     ('yam', 'pb1'): 0.00092,
     ('yam', 'pb2'): 0.00113}
    return rate.get((lineage, segment), 0.001)

def clock_std_dev(lineage, segement):
    return 0.2*clock_rate(lineage, segement)

rule all:
    input:
        #auspice_json = auspice_json,
        nexus_out = nexus_out,
        ggtree_out = ggtree_out,
        clade_data = clade_data
        
rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering.
        """
    input:
        sequences = input_fasta
    output:
        sequence_index = "results/sequenceindex.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule align:
    message:
        """
        Aligning sequences to tmp
          - filling gaps with N
        """
    input:
        sequences = input_fasta
    output:
        alignment = alignment_out
    params:
        refseq = reference
    threads: 15
    shell:
        """
        augur align \
            --sequences {input.sequences} \
        	--reference-sequence {params.refseq} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference \
            --nthreads {threads} \
            2>&1> results/align.log
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    # params:
    #     tree_builder_args = "-B 1000 -alrt 1000 -bnni"
    #     tree_builder_args = "-b 100"
            #     --tree-builder-args="{params.tree_builder_args}" \
            # --override-default-args \

    threads: 15
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {threads} \
            2>&1> results/tree.log
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:  
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = input_metadata
    output:
        tree = "results/tree_refined.nwk",
        node_data = "results/branch_lengths.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_rate = clock_rate(lineage, segment),
        clock_std_dev = clock_std_dev(lineage, segment),
        root = rootby,
        clock_filter_iqd = 4,
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --root {params.root} \
             2>&1> results/refine.log
        """
            # --timetree \
            # --clock-std-dev {params.clock_std_dev} \
            # --coalescent {params.coalescent} \
            # --date-confidence \
            # --date-inference {params.date_inference} \
            # --clock-filter-iqd {params.clock_filter_iqd} 2>&1> results/refine.log


rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = reference
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """

if segment == 'ha':
    rule clades:
        message: "Assigning and Annotating clades - ha"
        input:
            tree = rules.refine.output.tree,
            aa_muts = rules.translate.output.node_data,
            nuc_muts = rules.ancestral.output.node_data,
            clades = clades
        output:
            clade_data = "data/ha_clades.json"
        shell:
            """
            augur clades \
            --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
            """

if segment == 'na':
    rule clades:
        message: "Assigning and Annotating clades - na"
        input:
            tree = rules.refine.output.tree,
            clades = "data/ha_clades.json"
        output:
            clade_data = "data/na_clades.json"
        shell:
            """
            python3 scripts/import_tip_clades.py \
            --tree {input.tree} \
            --clades {input.clades} \
            --output {output.clade_data}
            """

rule recency:
    input:
        metadata = input_metadata
    params: 
        submission_date_field = "date", # the date which recency gets calculated off
        date_bins = [7, 14, 21, 29],
        date_bin_labels = ["last week", "last 2wks", "last 3wks", "last 4wks"],
        upper_bin_label = "Older"
    output:
        node_data = "results/recency.json",
    shell:
        """
        python3 scripts/construct_receny.py \
            --metadata {input.metadata} \
            --submission-date-field {params.submission_date_field} \
            --date-bins {params.date_bins} \
            --date-bin-labels {params.date_bin_labels:q} \
            --upper-bin-label {params.upper_bin_label} \
            --output {output.node_data}
        """

# rule export_auspice:
#     message: "Exporting data files for for auspice"
#     input:
#         tree = rules.refine.output.tree,
#         metadata = input_metadata,
#         branch_lengths = rules.refine.output.node_data,
#         #clades = rules.clades.output.clade_data,  {input.clades}
#         nt_muts = rules.ancestral.output.node_data,
#         aa_muts = rules.translate.output.node_data,
#         recency = rules.recency.output.node_data,
#         colors = colors,
#         lat_longs = lat_longs,
#         auspice_config = auspice_config
#     output:
#         tmp = "results/status.txt",
#         auspice_json = auspice_json
#     shell:
#         """
#         augur export v2 \
#             --tree {input.tree} \
#             --metadata {input.metadata} \
#             --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts}  {input.recency} \
#             --colors {input.colors} \
#             --lat-longs {input.lat_longs} \
#             --auspice-config {input.auspice_config} \
#             --output {output.auspice_json}
#         touch {output.tmp}
#         """

# export is not working for mroe than 1 trait
        # --traits-name {params.recency_name} \
        # --traits-json {input.recency} \
rule export:    
    message: "Exporting data files for for ggtree"
    input:
        tree = rules.refine.output.tree,
        aa = rules.translate.output.node_data,
        nt = rules.ancestral.output.node_data,
        recency = rules.recency.output.node_data,
        clades = rules.clades.output.clade_data
    output:
        annotated_tree = nexus_out
    params:
        add_gene_names = 'no',
        clade_name = 'clade_membership',
        recency_name = "recency"
    shell: 
        """
        python3 scripts/export.py \
        --tree {input.tree} \
        --aa {input.aa} \
        --nt {input.nt} \
        --add-gene-name {params.add_gene_names} \
        --traits-name {params.clade_name} \
        --traits-json {input.clades} \
        --output {output.annotated_tree}
        """

rule plot_ggtree:
    message: "Plotting tree to pdf"
    input:
        tree = rules.export.output.annotated_tree,
        meta = input_metadata
    output:
        ggtree = ggtree_out
    params:
        color_taxa = color_taxa,
        color_tip = color_tips_shapes,
        title = "%s-%s Phylogenetic Analysis" % (lineage, segment),
        size = "A3p"
    shell:
        """
        scripts/treeme.R  \
        -t {input.tree} \
        -m {input.meta} \
        -o {output.ggtree} \
        -c config/clades/clades_muts.txt \
        -l {params.color_taxa} \
        -p {params.color_tip} \
        -g {params.title} \
        -s {params.size}
        """