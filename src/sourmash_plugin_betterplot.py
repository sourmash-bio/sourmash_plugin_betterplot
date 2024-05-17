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
        plot(args)


def plot_composite_matrix(
    D, labelinfo, show_labels=True, vmax=1.0, vmin=0.0, force=False,
    cut_point=None,    
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

    fig = pylab.figure(figsize=(11, 8))
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


def plot(args):
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
