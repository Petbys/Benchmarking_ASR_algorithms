import phylopandas as pd
import dendropy as d
import pyasr
import toytree
import argparse




if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="testing pyasr")
    parser.add_argument(
        "--fasta",
        required=True,
    )
    parser.add_argument(
        "--tree",
        required=True,
    )
    parser.add_argument(
        "--dir",
        required=True,
    )
    args = parser.parse_args()
    fasta_path= args.fasta
    tree_path = args.tree
    work_dir = args.dir
    df = pd.read_fasta(fasta_path)
    df = df.phylo.read_newick(tree_path, combine_on="id")
    df = pyasr.reconstruct(df, working_dir=work_dir, alpha=1.235)

    # Slice out ancestors
    ancestors = df[df.type == 'node']

    # Write out Ancestors CSV
    ancestors.to_csv('ancestors.csv')

    # Preview some ancestors
    # Get a newick string to feed into ToyTree
    newick = df.phylo.to_newick(taxon_col='id', node_col='reconstruct_label')

    # Draw tree.
    tree_to_draw = toytree.tree(newick)
    tree_to_draw.draw(width=400, height=400,
        tip_labels_align=True,
        use_edge_lengths=True,
        node_labels=True)

    with open("output_aa.fasta", 'w') as output_file:
            output_file.write(tree_to_draw)
            output_file.close()

