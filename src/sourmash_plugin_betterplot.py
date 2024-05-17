"""betterplot plugin implementation"""

usage="""
   sourmash scripts plot2 <@CTB>
"""

epilog="""

Need help? Have questions? Ask at http://github.com/sourmash-bio/sourmash/issues!
"""

import argparse
import sourmash
from sourmash import sourmash_args
import os

import numpy
import pylab
import scipy.cluster.hierarchy as sch
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, csr_matrix

from collections import defaultdict

from sourmash.logging import debug_literal, error, notify, print_results
from sourmash.plugins import CommandLinePlugin


###

#
# CLI plugin - supports 'sourmash scripts plot2'
#

class Command_Plot2(CommandLinePlugin):
    command = 'plot2'             # 'scripts <command>'
    description = "@CTB"       # output with -h
    usage = "@CTB"               # output with no args/bad args as well as -h
    epilog = epilog             # output with -h
    formatter_class = argparse.RawTextHelpFormatter # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)
        # add argparse arguments here.
        subparser.add_argument("distances",
                               help='output from "sourmash compare"')
        subparser.add_argument("labels_from",
                               help='output from "sourmash compare --labels-to"')
        subparser.add_argument(
            "--vmin",
            default=0.0,
            type=float,
            help="lower limit of heatmap scale; default=%(default)f",
        )
        subparser.add_argument(
            "--vmax",
            default=1.0,
            type=float,
            help="upper limit of heatmap scale; default=%(default)f",
        )
        subparser.add_argument("--figsize-x", type=int, default=11)
        subparser.add_argument("--figsize-y", type=int, default=8)
        subparser.add_argument(
            "--subsample",
            type=int,
            metavar="N",
            help="randomly downsample to this many samples, max",
        )
        subparser.add_argument(
            "--subsample-seed",
            type=int,
            default=1,
            metavar="S",
            help="random seed for --subsample; default=1",
        )
        subparser.add_argument(
            "-f", "--force", action="store_true", help="forcibly plot non-distance matrices"
        )
        subparser.add_argument(
            "-o", "--output-figure", help="output figure to this file",
            required=True
        )
        subparser.add_argument(
            "--cut-point", type=float,
            help="cut point for dendrogram, to produce clusters"
        )
        subparser.add_argument(
            "--cluster-out", action='store_true',
            help="output clusters"
        )

    def main(self, args):
        # code that we actually run.
        super().main(args)
        plot2(args)


def plot_composite_matrix(
    D, labelinfo, show_labels=True, vmax=1.0, vmin=0.0, force=False,
        cut_point=None, figsize_x=11, figsize_y=8,
):
    """Build a composite plot showing dendrogram + distance matrix/heatmap.

    Returns a matplotlib figure.

    If show_labels is True, display labels. Otherwise, no labels are
    shown on the plot.
    """
    if D.max() > 1.0 or D.min() < 0.0:
        error(
            "This matrix doesn't look like a distance matrix - min value {}, max value {}",
            D.min(),
            D.max(),
        )
        if not force:
            raise ValueError("not a distance matrix")
        else:
            notify("force is set; scaling to [0, 1]")
            D -= D.min()
            D /= D.max()

    fig = pylab.figure(figsize=(figsize_x, figsize_y))
    ax1 = fig.add_axes([0.09, 0.1, 0.2, 0.6])

    # plot dendrogram
    Y = sch.linkage(D, method="single")  # centroid

    dend_kwargs = {}
    if cut_point is not None:
        cut_point = float(cut_point)
        dend_kwargs = dict(color_threshold=float(cut_point))
        
    labeltext = [ row['label'] for row in labelinfo ]
    Z1 = sch.dendrogram(
        Y,
        orientation="left",
        labels=labeltext,
        no_labels=not show_labels,
        get_leaves=True,
    )
    #ax1.set_xticks([])

    if cut_point is not None:
        ax1.axvline(x=cut_point, c='red', linestyle='dashed')

    xstart = 0.45
    width = 0.45
    if not show_labels:
        xstart = 0.315
    scale_xstart = xstart + width + 0.01

    idx1 = Z1["leaves"]

    # reorder D by the clustering in the dendrogram
    D = D[idx1, :]
    D = D[:, idx1]

    # show matrix
    axmatrix = fig.add_axes([xstart, 0.1, width, 0.6])

    im = axmatrix.matshow(
        D, aspect="auto", origin="lower", cmap=pylab.cm.YlGnBu, vmin=vmin, vmax=vmax
    )
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar.
    axcolor = fig.add_axes([scale_xstart, 0.1, 0.02, 0.6])
    pylab.colorbar(im, cax=axcolor)

    return fig, Y, D


def plot2(args):
    "Produce a clustering matrix and plot."
    import matplotlib as mpl

    #mpl.use("Agg")

    # load files
    D_filename = args.distances

    notify(f"loading comparison matrix from {D_filename}...")
    with open(D_filename, "rb") as f:
        D = numpy.load(f)
    notify(f"...got {D.shape[0]} x {D.shape[1]} matrix.", *D.shape)

    display_labels = True
    labelfilename = args.labels_from
    notify(f"loading labels from CSV file '{labelfilename}'")

    labelinfo = []
    with sourmash_args.FileInputCSV(labelfilename) as r:
        labelinfo = list(r)
        labelinfo.sort(key=lambda row: int(row['sort_order']))

    if len(labelinfo) != D.shape[0]:
        error("{} labels != matrix size, exiting", len(labeltext))
        sys.exit(-1)

    ### make the dendrogram:
    fig = pylab.figure(figsize=(8, 5))
    ax1 = fig.add_axes([0.1, 0.1, 0.7, 0.8])
    ax1.set_xticks([])
    ax1.set_yticks([])

    # subsample?
    if args.subsample:
        numpy.random.seed(args.subsample_seed)

        sample_idx = list(range(len(labelinfo)))
        numpy.random.shuffle(sample_idx)
        sample_idx = sample_idx[: args.subsample]

        np_idx = numpy.array(sample_idx)
        D = D[numpy.ix_(np_idx, np_idx)]
        labelinfo = [labelinfo[idx] for idx in sample_idx]

    ### make the dendrogram+matrix:
    (fig, linkage_Z, rmat) = plot_composite_matrix(
        D,
        labelinfo,
        show_labels=display_labels,
        vmin=args.vmin,
        vmax=args.vmax,
        force=args.force,
        cut_point=args.cut_point,
        figsize_x=args.figsize_x,
        figsize_y=args.figsize_y,
    )
    fig.savefig(args.output_figure, bbox_inches='tight')
    notify(f"wrote numpy distance matrix to: {args.output_figure}")

    # re-order labels along rows, top to bottom
    #reordered_labels = [labelinfo[i] for i in idx1]

    # output reordered labels with their clusters?
    if args.cut_point is not None and args.cluster_out:
        cut_point = float(args.cut_point)
        # @CTB 'distance'...
        assignments = sch.fcluster(linkage_Z, cut_point, 'distance')

        print(assignments)

        cluster_d = defaultdict(list)
        for cluster_n, label_row in zip(assignments, labelinfo):
            cluster_d[cluster_n].append(label_row)

        for k, v in cluster_d.items():
            print(f"cluster {k}")
            for row in v:
                print(f"\t{row['label']}")

            with open(f"cluster.{k}.list", "w") as fp:
                for row in v:
                    fp.write(f"{row['signature_file']}\n")


class Command_MDS(CommandLinePlugin):
    command = 'mds'             # 'scripts <command>'
    description = "@CTB"       # output with -h
    usage = "@CTB"               # output with no args/bad args as well as -h
    epilog = epilog             # output with -h
    formatter_class = argparse.RawTextHelpFormatter # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)

        subparser.add_argument('comparison_matrix')
        subparser.add_argument('-o', '--output-figure', required=True)

    def main(self, args):
        # code that we actually run.
        super().main(args)

        with open(args.comparison_matrix, 'rb') as f:
            mat = numpy.load(f)

        # Example usage
        # Assume object indices instead of names for simplicity
        #similarity_tuples = [(0, 1, 0.7), (0, 2, 0.4), (1, 2, 0.5)]
        #num_objects = 3  # You should know the total number of objects
        #sparse_matrix = create_sparse_similarity_matrix(similarity_tuples, num_objects)
        plot_mds_sparse(mat)

        plt.savefig(args.output_figure)


def create_sparse_similarity_matrix(tuples, num_objects):
    # Initialize matrix in LIL format for efficient setup
    similarity_matrix = lil_matrix((num_objects, num_objects))

    for obj1, obj2, similarity in tuples:
        similarity_matrix[obj1, obj2] = similarity
        if obj1 != obj2:
            similarity_matrix[obj2, obj1] = similarity

    # Ensure diagonal elements are 1
    similarity_matrix.setdiag(1)

    # Convert to CSR format for efficient operations later
    return similarity_matrix.tocsr()


def plot_mds_sparse(matrix):
    # Convert sparse similarity to dense dissimilarity matrix
    #dissimilarities = 1 - matrix.toarray()
    dissimilarities = 1 - matrix
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    mds_coords = mds.fit_transform(dissimilarities)
    plt.scatter(mds_coords[:, 0], mds_coords[:, 1])
    plt.title('MDS Plot')
    plt.xlabel('Dimension 1')
    plt.ylabel('Dimension 2')
