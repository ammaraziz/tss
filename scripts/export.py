#! /usr/bin/env python

import Bio
import Bio.Phylo
import argparse
import sys

sys.setrecursionlimit(10000)

def read_newick(filepath):
    with open(filepath, "r") as filehandle:
        tree = Bio.Phylo.read(filehandle, "newick")
    return tree


def read_json(filepath):
    import json

    with open(filepath, "r") as filehandle:
        contents = json.load(filehandle)

    return contents


def write_nexus(tree, filepath):
    with open(filepath, "w+") as filehandle:
        Bio.Phylo.write(tree, filehandle, "nexus")

def add_aa_muts(tree_node, muts, add_gene_name=True):
    """
    Add mutations (aa or nt) to branches on a tree
    tree_node           :  output of tree.find_clades()
    muts                :  json file read in via read_json() containing mutations
    add_gene_name       :  bool if gene name should be added 
    """

    muts_type_long = "aa_muts"
    muts_name_out = "aa_muts"


    node_muts = muts["nodes"][tree_node.name][muts_type_long]
    all_genes_muts = []

    for gene_name, mut_names in node_muts.items():
        if len(mut_names) > 0:
            # add gene names
            if add_gene_name:
                one_gene_muts_string = "{}|{}".format(
                    gene_name, " ".join(mut_names))
                all_genes_muts.append(one_gene_muts_string)
            # no gene
            else:
                all_genes_muts.append(" ".join(mut_names))

    if add_gene_name:
        all_genes_muts_string = "--".join(all_genes_muts)
    else:
        all_genes_muts_string = " ".join(all_genes_muts)
        
    node_comment = '{}={{"{}"}}'.format(muts_name_out, all_genes_muts_string)
    return(node_comment)

def add_nt_muts(tree_node, muts, add_gene_name=True):
    """
    Add nt mutations to branches on a tree
    tree_node           :  output of tree.find_clades()
    muts                :  json file read in via read_json() containing mutations
    add_gene_name       :  bool if gene name should be added 
    """

    muts_type_long = "muts"
    muts_name_out = "nt_muts"

    node_nt_muts = muts["nodes"][tree_node.name][muts_type_long]
    all_genes_muts = " ".join(node_nt_muts)

    if add_gene_name:
        all_genes_muts_string = "".join(all_genes_muts)
    else:
        all_genes_muts_string = "".join(all_genes_muts)
        
    node_comment = '{}={{"{}"}}'.format(muts_name_out, all_genes_muts_string)
    return(node_comment)


def add_trait(tree_node, trait_json, trait_name):
    
    if tree_node.name in trait_json["nodes"]:
        node_trait = trait_json["nodes"][tree_node.name][trait_name]
        node_comment = '{}={{"{}"}}'.format(trait_name, node_trait)

        return(node_comment)

def run(newick_tree_filepath, json_aa_filepath, json_nt_filepath, json_trait,
    json_traits_name, nexus_output_filepath, add_gene_name):

    tree = read_newick(newick_tree_filepath)
    aa_muts = read_json(json_aa_filepath)
    nt_muts = read_json(json_nt_filepath)
    if json_trait is not None:
        traits = read_json(json_trait)

    all_tree_nodes = tree.find_clades()
    for tree_node in all_tree_nodes:
        aa_node_comment = add_aa_muts(tree_node = tree_node, muts = aa_muts,
            add_gene_name = add_gene_name)
        nt_node_comment = add_nt_muts(tree_node = tree_node, muts = nt_muts,
            add_gene_name = add_gene_name)
        
        if json_trait is not None:
            trait_node_comment = add_trait(tree_node, traits, json_traits_name)
            complete_comment = "&{},{},{}".format(aa_node_comment, nt_node_comment, trait_node_comment)
        if json_trait is None:
            complete_comment = "&{},{}".format(aa_node_comment, nt_node_comment)
        tree_node.comment = complete_comment

    write_nexus(tree, nexus_output_filepath)

    return


if __name__ == "__main__":

    def collect_args():
        import argparse

        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description="Combine json files and a newick file into a nexus file",
            epilog="""
            export.py --tree nexus.tree 
                      --output nexus_new.tree 
                      --aa aa.json 
                      --nt nt.json 
                      --traits_name clades 
                      --traits_json clades.json
            """)

        if len(sys.argv) == 1:
            parser.print_help()

        parser.add_argument("-t", "--tree", required=True,
            type=str, help="Newick tree file")
        parser.add_argument("-o", "--output", required=True,
            type=str, help="Output nexus filename")
        parser.add_argument("-a", "--aa", required=True, 
            type=str, help="JSON file with amino acid mutations")
        parser.add_argument("-n", "--nt", required=True, 
            type=str, help="JSON with nucleotide mutations")
        parser.add_argument("--traits-name", required=False, 
            type=str, help="Name of the trait to add. eg: clade_membership. must match json id")
        parser.add_argument("--traits-json", required=False, 
            type=str, help="JSON with traits to add")
        parser.add_argument("--add-gene-name", required=False,
            default="yes", type=str, help="Should gene names be added to aa/nt labels? [yes]/no")

        params = parser.parse_args()

        if params.add_gene_name in ('yes', 'Yes'):
            params.add_gene_name = True
        elif params.add_gene_name in ('no', 'No'):
            params.add_gene_name = False

        # required as they are parsed together 
        required_together = ("traits_name", "traits_json")
        if not all([getattr(params, x) for x in required_together]):
            raise RuntimeError("Please supply both --traits-name and --traits-json together. ")

        return params

    global_params = collect_args()

    run(
        global_params.tree,
        global_params.aa,
        global_params.nt,
        global_params.traits_json,
        global_params.traits_name,
        global_params.output,
        global_params.add_gene_name
    )
