
"""betterplot plugin implementation"""

epilog = """

Need help? Have questions? Ask at http://github.com/sourmash-bio/sourmash/issues!
"""

import sys
import argparse
import os
import csv
from collections import defaultdict, Counter
from itertools import chain, combinations
import pickle

import numpy
import pylab
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib_venn import venn2, venn3
import scipy.cluster.hierarchy as sch
from sklearn.manifold import MDS, TSNE
from scipy.sparse import lil_matrix, csr_matrix
from matplotlib.lines import Line2D
import seaborn as sns
import upsetplot
import pandas as pd
import plotly.graph_objects as go
import squarify

# this turns off a warning in presence_filter, but results in an error in
# upsetplot :sweat_smile:
#pd.options.mode.copy_on_write = True

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix

import sourmash
from sourmash import sourmash_args
from sourmash.logging import debug_literal, error, notify, print_results
from sourmash.plugins import CommandLinePlugin
import sourmash_utils
from sourmash.cli.utils import (add_ksize_arg, add_moltype_args, add_scaled_arg)


### utility functions

def load_labelinfo_csv(filename):
    "Load file output by 'sourmash compare --labels-to'"
    with sourmash_args.FileInputCSV(filename) as r:
        labelinfo = list(r)

    if len(labelinfo) == 0:
        raise Exception("ERROR: no labels found!?")

    if not 'sort_order' in labelinfo[0]:
        raise Exception("ERROR: this doesn't look like a 'labels' file produced by 'sourmash compare --labels-to'")

    labelinfo.sort(key=lambda row: int(row["sort_order"]))
    return labelinfo


def load_categories_csv(filename, labelinfo):
    "Load categories file, integrate with labelinfo => colors"
    with sourmash_args.FileInputCSV(filename) as r:
        categories = list(r)

    category_to_color = {}
    colors = None
    if categories:
        # first, figure out which column is matching between labelinfo
        # and categories file.
        assert labelinfo
        keys = set(categories[0].keys())
        if "category" not in keys:
            raise ValueError(f"no 'category' column found in keys: only {keys}")
        keys -= {"category"}

        key = None
        for k in keys:
            if k in labelinfo[0].keys():
                notify(f"found category key: {k}")
                key = k
                break

        # found one? awesome. load in all the categories & assign colors.

        if key:
            # get distinct categories
            category_values = set([row["category"] for row in categories])
            category_values = list(sorted(category_values))

            # map categories to colormap colors
            cat_colors = list(map(plt.cm.tab10, range(len(category_values))))

            # build map of category => color
            category_to_color = {}
            for v, color in zip(category_values, cat_colors):
                category_to_color[v] = color

            # build map of key => color
            label_to_color = {}
            for row in categories:
                value = row[key]
                category = row["category"]
                color = category_to_color[category]
                label_to_color[value] = color

            # build list of colors in sample order
            colors = []
            missing_values = []
            for row in labelinfo:
                value = row[key]
                if value not in label_to_color:
                    missing_values.append(value)
                    continue
                color = label_to_color[value]
                colors.append(color)

            if missing_values:
                raise ValueError(f"values {missing_values} are missing in categories file")

        else:
            notify(f"no valid key column found in categories file '{filename}'.")
    else:
        notify(f"nothing in categories file '{filename}'?!")

    return category_to_color, colors


def load_categories_csv_for_labels(filename, samples_d):
    "Load a categories CSV that uses the 'label' column."
    with sourmash_args.FileInputCSV(filename) as r:
        categories = list(r)

    category_to_color = {}
    colors = None
    if categories:
        key = "label"

        # load distinct categories
        category_values = list(set([row["category"] for row in categories]))
        category_values.sort()

        # map categories to color
        cat_colors = list(map(plt.cm.tab10, range(len(category_values))))
        category_to_color = {}
        for v, color in zip(category_values, cat_colors):
            category_to_color[v] = color

        # map label to color
        label_to_color = {}
        for row in categories:
            label = row[key]
            cat = row["category"]
            label_to_color[label] = category_to_color[cat]

        # build list of colors
        colors = []
        for label, idx in samples_d.items():
            color = label_to_color[label]
            colors.append(color)
    else:
        notify(f"nothing in categories file '{filename}'?!")

    return category_to_color, colors


def manysearch_rows_to_index(rows, *, column_name='query_name'):
    """Extract # of samples and build name -> sample_index map from manysearch.

    Note, column names are "query_name", "match_name", or "both".
    """
    if column_name in ('query_name', 'match_name'):
        samples = set(( row[column_name] for row in rows ))
    elif column_name == 'both':
        samples = set()
        for col in ('query_name', 'match_name'):
            samples.update(( row[col] for row in rows ))
    else:
        raise ValueError(f"unknown column_name '{column_name}'")

    samples = list(sorted(samples))
    sample_d = {}
    for n, sample_name in enumerate(samples):
        sample_d[sample_name] = n

    return sample_d


def labelinfo_to_idents(labelinfo):
    "Make x/yticklabels from a list of labelinfo rows"
    xx = []
    for row in labelinfo:
        ident = row["label"].split(' ')
        ident = ident[0]
        xx.append(ident)

    return xx


def sample_d_to_idents(sample_d):
    "Make x/yticklabels from a list of (k, v) from sample_d.items()."
    xx = []
    for k, v in sample_d:
        ident = k.split(' ')
        ident = ident[0]
        xx.append(ident)

    return xx

#
# CLI plugin code
#

class Command_Plot2(CommandLinePlugin):
    command = "plot2"  # 'scripts <command>'
    description = (
        "plot a distance matrix produced by 'sourmash compare'"  # output with -h
    )
    usage = "sourmash scripts plot <matrix> <labels_csv> -o <output.png>"  # output with no args/bad args as well as -h
    epilog = epilog  # output with -h
    formatter_class = argparse.RawTextHelpFormatter  # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)
        subparser.add_argument("distances", help='output from "sourmash compare"')
        subparser.add_argument(
            "labels_from", help='output from "sourmash compare --labels-to"'
        )
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
            "-f",
            "--force",
            action="store_true",
            help="forcibly plot non-distance matrices",
        )
        subparser.add_argument(
            "-o", "--output-figure", help="output figure to this file", required=True
        )
        subparser.add_argument(
            "--cut-point",
            type=float,
            help="cut point for dendrogram, to produce clusters",
        )
        subparser.add_argument(
            "--cluster-out", action="store_true", help="output clusters"
        )
        subparser.add_argument(
            "--cluster-prefix",
            default=None,
            help="prefix to prepend to cluster names; default is cmp file",
        )
        subparser.add_argument(
            "--dendrogram-only",
            "--no-matrix",
            action="store_true",
            help="plot only the dendrogram",
        )

    def main(self, args):
        super().main(args)
        plot2(args)


def plot_composite_matrix(
    D,
    labelinfo,
    show_labels=True,
    vmax=1.0,
    vmin=0.0,
    force=False,
    cut_point=None,
    figsize_x=11,
    figsize_y=8,
    dendrogram_only=False,
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

    labeltext = [row["label"] for row in labelinfo]

    truncate_name = lambda x: x[:30-3] + '...' if len(x) >= 30 else x
    labeltext = [truncate_name(label) for label in labeltext]

    Z1 = sch.dendrogram(
        Y,
        orientation="left",
        labels=labeltext,
        no_labels=not show_labels,
        get_leaves=True,
    )

    # draw cut point
    if cut_point is not None:
        ax1.axvline(x=cut_point, c="red", linestyle="dashed")

    # draw matrix
    xstart = 0.45
    width = 0.45
    if not show_labels:
        xstart = 0.315
    scale_xstart = xstart + width + 0.01

    idx1 = Z1["leaves"]

    if not dendrogram_only:
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
    # load files
    D_filename = args.distances

    notify(f"loading comparison matrix from {D_filename}...")
    with open(D_filename, "rb") as f:
        D = numpy.load(f)
    notify(f"...got {D.shape[0]} x {D.shape[1]} matrix.", *D.shape)

    display_labels = True
    labelfilename = args.labels_from
    notify(f"loading labels from CSV file '{labelfilename}'")

    labelinfo = load_labelinfo_csv(labelfilename)

    if len(labelinfo) != D.shape[0]:
        error("{} labels != matrix size, exiting", len(labelinfo))
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
        dendrogram_only=args.dendrogram_only,
    )
    fig.savefig(args.output_figure, bbox_inches="tight")
    notify(f"wrote plot to: {args.output_figure}")

    # re-order labels along rows, top to bottom
    # reordered_labels = [labelinfo[i] for i in idx1]

    # output reordered labels with their clusters?
    if args.cut_point is not None and args.cluster_out:
        cut_point = float(args.cut_point)
        prefix = args.cluster_prefix or os.path.basename(D_filename)
        notify(f"outputting clusters with prefix '{prefix}'")

        # generate clusters using 'fcluster'
        # @CTB 'distance'...? does this conflict with 'linkage'?
        assignments = sch.fcluster(linkage_Z, cut_point, "distance")

        # reorganize labelinfo by cluster
        cluster_d = defaultdict(list)
        for cluster_n, label_row in zip(assignments, labelinfo):
            cluster_d[cluster_n].append(label_row)

        # output labelinfo rows.
        notify(f"writing {len(cluster_d)} clusters.")
        for k, v in cluster_d.items():
            filename = f"{prefix}.{k}.csv"
            with sourmash_args.FileOutputCSV(filename) as fp:
                w = csv.DictWriter(fp, fieldnames=labelinfo[0].keys())
                w.writeheader()
                for row in v:
                    w.writerow(row)


def plot_mds(matrix, *, colors=None, category_map=None, metric=True):
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42,
              metric=metric)
    mds_coords = mds.fit_transform(matrix)
    plt.scatter(mds_coords[:, 0], mds_coords[:, 1], color=colors)
    plt.xlabel("Dimension 1")
    plt.ylabel("Dimension 2")

    if colors and category_map:
        # create a custom legend of just the categories
        legend_elements = []
        for k, v in category_map.items():
            legend_elements.append(Line2D([0], [0], color=v, label=k, marker="o", lw=0))
        plt.legend(handles=legend_elements)


class Command_MDS(CommandLinePlugin):
    command = "mds"  # 'scripts <command>'
    description = "plot a 2-D multidimensional scaling plot from 'sourmash compare' output"  # output with -h
    usage = "sourmash scripts mds <matrix> <labels_csv> -o <figure.png>"  # output with no args/bad args as well as -h
    epilog = epilog  # output with -h
    formatter_class = argparse.RawTextHelpFormatter  # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)

        subparser.add_argument(
            "comparison_matrix", help="output from 'sourmash compare'"
        )
        subparser.add_argument(
            "labels_from", help="output from 'sourmash compare --labels-to'"
        )
        subparser.add_argument(
            "-C", "--categories-csv", help="CSV mapping label columns to categories"
        )
        subparser.add_argument("-o", "--output-figure", required=True)
        subparser.add_argument("--metric", dest="metric", default=True,
                               action="store_true",
                               help="compute MDS (metric) - the default")
        subparser.add_argument("--nmds", dest="metric",
                               action="store_false",
                               help="compute NMDS (non-metric)")

    def main(self, args):
        super().main(args)

        if args.metric:
            notify(f"building metric plot (MDS)")
        else:
            notify(f"building non-metric plot (NMDS)")

        with open(args.comparison_matrix, "rb") as f:
            mat = numpy.load(f)

        labelinfo = load_labelinfo_csv(args.labels_from)

        # load categories?
        category_map = None
        colors = None
        if args.categories_csv:
            category_map, colors = load_categories_csv(args.categories_csv, labelinfo)

        dissim = 1 - mat
        plot_mds(dissim, colors=colors, category_map=category_map,
                 metric=args.metric)

        notify(f"writing figure to '{args.output_figure}'")
        plt.savefig(args.output_figure)


class Command_PairwiseToMatrix(CommandLinePlugin):
    command = "pairwise_to_matrix"  # 'scripts <command>'
    description = "convert pairwise CSV output to a 'compare' matrix"  # output with -h
    usage = "sourmash scripts pairwise_to_matrix <pairwise_csv> -o <matrix_cmp>"  # output with no args/bad args as well as -h
    epilog = epilog  # output with -h
    formatter_class = argparse.RawTextHelpFormatter  # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)

        subparser.add_argument(
            "pairwise_csv", help="output from 'sourmash scripts pairwise'"
        )
        subparser.add_argument("-o", "--output-matrix", required=True)
        subparser.add_argument("--labels-to")
        subparser.add_argument(
            "-u",
            "--use-column",
            default="jaccard",
            help="column name to use in matrix (default: jaccard)",
        )

    def main(self, args):
        super().main(args)

        notify(f"loading '{args.pairwise_csv}'")
        with sourmash_args.FileInputCSV(args.pairwise_csv) as r:
            rows = list(r)

        sample_d = manysearch_rows_to_index(rows)
        notify(f"loaded {len(rows)} rows containing {len(sample_d)} distinct samples")

        mat = numpy.zeros((len(sample_d), len(sample_d)))
        colname = args.use_column

        for row in rows:
            # get unique indices for each query/match pair.
            q = row["query_name"]
            qi = sample_d[q]
            m = row["match_name"]
            mi = sample_d[m]
            jaccard = float(row[colname])

            mat[qi, mi] = jaccard
            mat[mi, qi] = jaccard

        numpy.fill_diagonal(mat, 1)

        notify(f"writing output matrix to '{args.output_matrix}'")
        with open(args.output_matrix, "wb") as fp:
            numpy.save(fp, mat)

        notify(f"writing output labels.txt to '{args.output_matrix}.labels.txt'")
        with open(args.output_matrix + ".labels.txt", "wt") as fp:
            for label, n in sample_d.items():
                fp.write(label + "\n")

        if args.labels_to:
            notify(f"writing output labels csv to '{args.labels_to}'")
            with open(args.labels_to, "w", newline="") as fp:
                w = csv.writer(fp)
                w.writerow(["sort_order", "label"])
                for label, n in sample_d.items():
                    w.writerow([n, label])


class Command_MDS2(CommandLinePlugin):
    command = "mds2"  # 'scripts <command>'
    description = "plot a 2-D multidimensional scaling plot from branchwater plugin's 'pairwise' output"  # output with -h
    usage = "sourmash scripts mds2 <matrix> -o <figure.png>"  # output with no args/bad args as well as -h
    epilog = epilog  # output with -h
    formatter_class = argparse.RawTextHelpFormatter  # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)

        subparser.add_argument(
            "pairwise_csv", help="output from 'sourmash scripts pairwise'"
        )
        subparser.add_argument(
            "-C", "--categories-csv", help="CSV mapping label columns to categories"
        )
        subparser.add_argument("-o", "--output-figure", required=True)
        subparser.add_argument("--metric", dest="metric", default=True,
                               action="store_true",
                               help="compute MDS (metric) - the default")
        subparser.add_argument("--nmds", dest="metric",
                               action="store_false",
                               help="compute NMDS (non-metric)")

    def main(self, args):
        super().main(args)

        if args.metric:
            notify(f"building metric plot (MDS)")
        else:
            notify(f"building non-metric plot (NMDS)")

        with sourmash_args.FileInputCSV(args.pairwise_csv) as r:
            rows = list(r)

        # pick out all the distinct queries/matches.
        notify(f"loaded {len(rows)} rows from '{args.pairwise_csv}'")
        sample_d = manysearch_rows_to_index(rows, column_name='both')
        notify(f"loaded {len(sample_d)} total elements")

        mat = numpy.zeros((len(sample_d), len(sample_d)))

        for row in rows:
            # get unique indices for each query/match pair.
            q = row["query_name"]
            qi = sample_d[q]
            m = row["match_name"]
            mi = sample_d[m]
            jaccard = float(row["jaccard"])

            mat[qi, mi] = jaccard
            mat[mi, qi] = jaccard

        numpy.fill_diagonal(mat, 1)

        # load categories?
        category_map = None
        colors = None
        if args.categories_csv:
            category_map, colors = load_categories_csv_for_labels(
                args.categories_csv, sample_d
            )

        dissim = 1 - mat
        plot_mds(dissim, colors=colors, category_map=category_map,
                 metric=args.metric)

        notify(f"writing figure to '{args.output_figure}'")
        plt.savefig(args.output_figure)


# @CTB unused code for sparse matrix foo. Revisit!
def create_sparse_dissimilarity_matrix(tuples, num_objects):
    # Initialize matrix in LIL format for efficient setup
    similarity_matrix = lil_matrix((num_objects, num_objects))

    for obj1, obj2, similarity in tuples:
        similarity_matrix[obj1, obj2] = 1 - similarity
        if obj1 != obj2:
            similarity_matrix[obj2, obj1] = 1 - similarity

    # Ensure diagonal elements are 1
    similarity_matrix.setdiag(1)

    # Convert to array format
    # @CTB use tocsr or tocoo instead??
    return similarity_matrix.toarray()


class Command_Plot3(CommandLinePlugin):
    command = "plot3"  # 'scripts <command>'
    description = (
        "plot a distance matrix produced by 'sourmash compare'"  # output with -h
    )
    usage = "sourmash scripts plot <matrix> <labels_csv> -o <output.png>"  # output with no args/bad args as well as -h
    epilog = epilog  # output with -h
    formatter_class = argparse.RawTextHelpFormatter  # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)
        subparser.add_argument("distances", help='output from "sourmash compare"')
        subparser.add_argument(
            "labels_from", help='output from "sourmash compare --labels-to"'
        )
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
            "-o", "--output-figure", help="output figure to this file", required=True
        )
        subparser.add_argument(
            "-C", "--categories-csv", help="CSV mapping label columns to categories"
        )
        subparser.add_argument(
            "--no-labels", action="store_true",
            help="disable X & Y axis labels"
        )

    def main(self, args):
        super().main(args)
        D_filename = args.distances

        notify(f"loading comparison matrix from {D_filename}...")
        with open(D_filename, "rb") as f:
            D = numpy.load(f)
        notify(f"...got {D.shape[0]} x {D.shape[1]} matrix.", *D.shape)

        labelfilename = args.labels_from
        notify(f"loading labels from CSV file '{labelfilename}'")
        labelinfo = load_labelinfo_csv(labelfilename)

        if len(labelinfo) != D.shape[0]:
            error("{} labels != matrix size, exiting", len(labelinfo))
            sys.exit(-1)

        # load categories?
        category_map = None
        colors = None
        if args.categories_csv:
            category_map, colors = load_categories_csv(args.categories_csv, labelinfo)
        # subsample?
        if args.subsample:
            numpy.random.seed(args.subsample_seed)

            sample_idx = list(range(len(labelinfo)))
            numpy.random.shuffle(sample_idx)
            sample_idx = sample_idx[: args.subsample]

            np_idx = numpy.array(sample_idx)
            D = D[numpy.ix_(np_idx, np_idx)]
            labelinfo = [labelinfo[idx] for idx in sample_idx]

        # turn into dissimilarity matrix
        # dissim = 1 - D
        # numpy.fill_diagonal(dissim, 1)
        dissim = D

        if args.no_labels:
            yticklabels=[]
        else:
            yticklabels=labelinfo_to_idents(labelinfo)

        # plot!
        fig = sns.clustermap(
            dissim,
            figsize=(args.figsize_x, args.figsize_y),
            vmin=args.vmin,
            vmax=args.vmax,
            col_colors=colors,
            yticklabels=yticklabels,
            xticklabels=[],
            cmap="flare",
        )

        if colors and category_map:
            # create a custom legend of just the categories
            legend_elements = []
            for k, v in category_map.items():
                legend_elements.append(
                    Line2D([0], [0], color=v, label=k, marker="o", lw=0)
                )
            fig.ax_col_dendrogram.legend(handles=legend_elements)

        # turn off column dendrogram
        fig.ax_row_dendrogram.set_visible(False)

        fig.savefig(args.output_figure, bbox_inches="tight")
        notify(f"wrote plot to: {args.output_figure}")


class Command_Clustermap1(CommandLinePlugin):
    command = "clustermap1"  # 'scripts <command>'
    description = "plot the results of 'manysearch'"  # output with -h
    usage = "sourmash scripts clustermap1 <manysearch_csv> -o <output.png>"  # output with no args/bad args as well as -h
    epilog = epilog  # output with -h
    formatter_class = argparse.RawTextHelpFormatter  # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)
        subparser.add_argument("manysearch_csv", help='output from "sourmash compare"')

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
            "-o", "--output-figure", help="output figure to this file", required=True
        )
        subparser.add_argument(
            "-R",
            "--row-categories-csv",
            help="CSV mapping labels @CTB query or against? to categories",
        )
        subparser.add_argument(
            "-C",
            "--col-categories-csv",
            help="CSV mapping labels @CTB query or against? to categories",
        )
        subparser.add_argument(
            "-u",
            "--use-column",
            default="jaccard",
            help="column name to use in matrix (default: jaccard)",
        )
        subparser.add_argument(
            "--boolean", action="store_true", help="convert values into 0/1"
        )
        subparser.add_argument(
            "--no-labels", action="store_true",
            help="disable X & Y axis labels"
        )
        subparser.add_argument(
            "--no-x-labels", action="store_true",
            help="disable X axis labels"
        )
        subparser.add_argument(
            "--no-y-labels", action="store_true",
            help="disable Y axis labels"
        )

    def main(self, args):
        super().main(args)
        with sourmash_args.FileInputCSV(args.manysearch_csv) as r:
            rows = list(r)

        # pick out all the distinct queries/matches.
        notify(f"loaded {len(rows)} rows from '{args.manysearch_csv}'")

        query_d = manysearch_rows_to_index(rows, column_name='query_name')
        against_d = manysearch_rows_to_index(rows, column_name='match_name')

        notify(f"loaded {len(query_d)} x {len(against_d)} total elements")

        query_d_items = list(sorted(query_d.items(), key=lambda x: x[1]))
        against_d_items = list(sorted(against_d.items(), key=lambda x: x[1]))

        mat = numpy.zeros((len(query_d), len(against_d)))

        colname = args.use_column
        notify(f"using column '{colname}'")
        make_bool = args.boolean
        if make_bool:
            notify(f"forcing values to 0 / 1 and disabling color bar because of --boolean")

        for row in rows:
            q = row["query_name"]
            qi = query_d[q]
            m = row["match_name"]
            mi = against_d[m]
            value = float(row[colname])
            if make_bool:
                value = 1 if value else 0

            mat[qi, mi] = value

        # load categories?
        row_category_map = None
        row_colors = None
        if args.row_categories_csv:
            row_category_map, row_colors = load_categories_csv_for_labels(
                args.row_categories_csv, query_d
            )

        col_category_map = None
        col_colors = None
        if args.col_categories_csv:
            col_category_map, col_colors = load_categories_csv_for_labels(
                args.col_categories_csv, against_d
            )

        kw_args = {}
        if args.boolean:        # turn off colorbar if boolean.
            kw_args['cbar_pos'] = None

        yticklabels=sample_d_to_idents(query_d_items)
        xticklabels=sample_d_to_idents(against_d_items)
        if args.no_labels:
            xticklabels = []
            yticklabels = []
        elif args.no_x_labels:
            xticklabels = []
        elif args.no_y_labels:
            yticklabels = []

        # turn into dissimilarity matrix
        # plot!
        fig = sns.clustermap(
            mat,
            figsize=(args.figsize_x, args.figsize_y),
            vmin=args.vmin,
            vmax=args.vmax,
            col_colors=col_colors,
            row_colors=row_colors,
            xticklabels=xticklabels,
            yticklabels=yticklabels,
            cmap="flare",
            **kw_args
        )

        if col_colors and col_category_map:
            # create a custom legend of just the categories
            legend_elements = []
            for k, v in col_category_map.items():
                legend_elements.append(
                    Line2D([0], [0], color=v, label=k, marker="o", lw=0)
                )
            fig.ax_col_dendrogram.legend(handles=legend_elements)

        if row_colors and row_category_map:
            # create a custom legend of just the categories
            legend_elements = []
            for k, v in row_category_map.items():
                legend_elements.append(
                    Line2D([0], [0], color=v, label=k, marker="o", lw=0)
                )
            fig.ax_row_dendrogram.legend(handles=legend_elements)

        fig.savefig(args.output_figure, bbox_inches="tight")
        notify(f"wrote plot to: {args.output_figure}")


class Command_Upset(CommandLinePlugin):
    command = "upset"  # 'scripts <command>'
    description = "visualize intersections of sketches using upsetplot"  # output with -h
    usage = "sourmash scripts upset <sketches> [<sketches> ...] -o <output>.png"  # output with no args/bad args as well as -h
    epilog = epilog  # output with -h
    formatter_class = argparse.RawTextHelpFormatter  # do not reformat multiline

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('sketches', nargs='+')
        p.add_argument('--show-singletons', action='store_true',
                       help='show membership of single sketches as well')
        p.add_argument('-o', '--output-figure', required=True)
        p.add_argument('--truncate-labels-at', default=30, type=int,
                       help="limit labels to this length (default: 30)")
        sourmash_utils.add_standard_minhash_args(p)
        p.add_argument('--sort-by', default='cardinality',
                       choices=['cardinality', 'degree', '-cardinality', '-degree'],
                       help='sort display by size of intersection, or number of categories intersected')
        p.add_argument('--min-subset-size', default="0%",
                       type=str,
                       help="omit sets below this size or percentage (default: '0%%')")
        p.add_argument('--show-percentages', action="store_true",
                       help='show percentages on plot')
        p.add_argument('--save-intersections-to-file', default=None,
                       help='save intersections to a file, to avoid expensive recalculations')
        p.add_argument('--load-intersections-from-file', default=None,
                       help='load precalculated intersections from a file')

#        p.add_argument('--save-names-to-file', default=None,
#                       help='save set names to a file, for editing & customization')
#        p.add_argument('--load-names-from-file', default=None,
#                       help='load set names from a file to customize plots')

    def main(self, args):
        super().main(args)

        # https://docs.python.org/3/library/itertools.html
        def powerset(iterable, *, start=2):
            "powerset([1,2,3]) → () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
            s = list(iterable)
            return chain.from_iterable(combinations(s, r) for r in range(start, len(s)+1))

        select_mh = sourmash_utils.create_minhash_from_args(args)
        print(f"selecting sketches: {select_mh}")
        scaled = select_mh.scaled

        siglist = []
        for filename in args.sketches:
            print(f"loading sketches from file {filename}")
            db = sourmash_utils.load_index_and_select(filename, select_mh)

            # downsample?
            for ss in db.signatures():
                if ss.minhash.scaled != scaled:
                    with ss.update() as ss:
                        ss.minhash = ss.minhash.downsample(scaled=scaled)
                siglist.append(ss)

        notify(f"Loaded {len(siglist)} signatures & downsampled to scaled={scaled}")

        names_check = [ ss.name for ss in siglist ]
        if len(set(names_check)) != len(names_check):
            notify("ERROR: duplicate names or sketches; please fix!!")
            cnt = Counter(names_check)
            for k, v in cnt.most_common():
                if v > 1:
                    print(f"\t* {k} shows up {v} times")
            sys.exit(-1)

        # @CTB: check scaled, ksize, etc.

        if not siglist:
            notify(f"ERROR: found no sketches. Exiting!")
            sys.exit(-1)

        if len(siglist) > 10:
            notify(f"WARNING: this is probably too many sketches.")

        start = 2
        if args.show_singletons:
            notify(f"Showing individual sketch membership b/c of --show-singletons")
            start = 1
        else:
            notify(f"Omitting individual sketch membership; use --show-singletons to see.")

        pset = list(powerset(siglist, start=start))
        pset.sort(key=lambda x: -len(x))
        #get_name = lambda x: [ ss.name.split(' ')[0] for ss in x ]
        truncate_at = args.truncate_labels_at
        truncate_name = lambda x: x[:truncate_at-3] + '...' if len(x) >= truncate_at else x
        get_name = lambda x: [ truncate_name(ss.name) for ss in x ]
        names = [ get_name(combo) for combo in pset ]

        notify(f"powerset of distinct combinations: {len(pset)}")

        # CTB: maybe turn the intersection code below into a class?

        if args.load_intersections_from_file:
            notify(f"loading intersections from '{args.load_intersections_from_file}'")
            with open(args.load_intersections_from_file, 'rb') as fp:
                check_names, nonzero_names, counts = pickle.load(fp)

            # confirm!
            if check_names != names:
                error("ERROR: saved intersections do not match provided sketches!?")
                sys.exit(-1)
        else:
            notify(f"generating intersections...")
            counts = []
            nonzero_names = []
            subtract_me = set()
            for n, combo in enumerate(pset):
                if n and n % 10 == 0:
                    notify(f"...{n} of {len(pset)}", end="\r")

                combo = list(combo)
                ss = combo.pop()
                hashes = set(ss.minhash.hashes) - subtract_me

                while combo and hashes:
                    ss = combo.pop()
                    hashes.intersection_update(ss.minhash.hashes)

                if hashes:
                    counts.append(len(hashes) * scaled)
                    nonzero_names.append(names[n])
                    subtract_me.update(hashes)
            notify(f"\n...done! {len(nonzero_names)} non-empty intersections of {len(names)} total.")

            # maybe decrease memory, but also prevent re/mis-use of these :)
            del subtract_me
            del hashes
#            del names

            if args.save_intersections_to_file:
                notify(f"saving intersections to '{args.save_intersections_to_file}'")
                with open(args.save_intersections_to_file, 'wb') as fp:
                    pickle.dump((names, nonzero_names, counts), fp)

#        if args.save_names_to_file:
#            with open(args.save_names_to_file, 'w', newline='') as fp:
#                w = csv.writer(fp)
#                w.writerow(['sort_order', 'name'])
#                for n, name in enumerate(names):
#                    w.writerow([n, name])
#            notify(f"saved {len(names)} names to '{args.save_names_to_file}'")

#        if args.load_names_from_file:
#            with open(args.load_names_from_file, 'r', newline='') as fp:
#                r = csv.DictReader(fp)
#                rows = list(r)
#                if 'sort_order' not in rows[0].keys():
#                    error("'sort_order' must be a column in names file '{args.load_names_from_file}'")
#                if 'name' not in rows[0].keys():
#                    error("'name' must be a column in names file '{args.load_names_from_file}'")
#
#                rows.sort(key=lambda x: int(x["sort_order"]))
#                names = [ row["name"] for row in rows ]
#            notify("loaded {len(names)} names from '{args.load_names_from_file}'")

        ## now! calculate actual data for upsetplot...

        data = upsetplot.from_memberships(nonzero_names, counts)

        try:
            min_subset_size = float(args.min_subset_size)
            notify(f"setting min_subset_size={min_subset_size:g} (number)")
        except ValueError:
            min_subset_size = args.min_subset_size
            notify(f"setting min_subset_size='{min_subset_size}' (percentage)")

        print(data)
        upsetplot.plot(data, sort_by=args.sort_by,
                       min_subset_size=min_subset_size,
                       show_percentages=args.show_percentages)

        notify(f"saving upsetr figure to '{args.output_figure}'")
        plt.savefig(args.output_figure, bbox_inches="tight")
        # @CTB use 'notify'
        

def plot_tsne(matrix, *, colors=None, category_map=None):
    perplexity = min(len(matrix) - 1, 50)
    tsne = TSNE(n_components=2, random_state=42, perplexity=perplexity)
    tsne_coords = tsne.fit_transform(matrix)
    plt.scatter(tsne_coords[:, 0], tsne_coords[:, 1], color=colors)
    plt.xlabel("Dimension 1")
    plt.ylabel("Dimension 2")

    if colors and category_map:
        # create a custom legend of just the categories
        legend_elements = []
        for k, v in category_map.items():
            legend_elements.append(Line2D([0], [0], color=v, label=k, marker="o", lw=0))
        plt.legend(handles=legend_elements)


class Command_TSNE(CommandLinePlugin):
    command = "tsne"  # 'scripts <command>'
    description = "plot a tSNE from 'sourmash compare' output"  # output with -h
    usage = "sourmash scripts tsne <matrix> <labels_csv> -o <figure.png>"  # output with no args/bad args as well as -h
    epilog = epilog  # output with -h
    formatter_class = argparse.RawTextHelpFormatter  # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)

        subparser.add_argument(
            "comparison_matrix", help="output from 'sourmash compare'"
        )
        subparser.add_argument(
            "labels_from", help="output from 'sourmash compare --labels-to'"
        )
        subparser.add_argument(
            "-C", "--categories-csv", help="CSV mapping label columns to categories"
        )
        subparser.add_argument("-o", "--output-figure", required=True)

    def main(self, args):
        super().main(args)

        with open(args.comparison_matrix, "rb") as f:
            mat = numpy.load(f)

        labelinfo = load_labelinfo_csv(args.labels_from)

        # load categories?
        category_map = None
        colors = None
        if args.categories_csv:
            category_map, colors = load_categories_csv(args.categories_csv, labelinfo)

        dissim = 1 - mat
        plot_tsne(dissim, colors=colors, category_map=category_map)

        notify(f"writing figure to '{args.output_figure}'")
        plt.savefig(args.output_figure)


class Command_TSNE2(CommandLinePlugin):
    command = "tsne2"  # 'scripts <command>'
    description = "plot a 2-D multidimensional scaling plot from branchwater plugin's 'pairwise' output"  # output with -h
    usage = "sourmash scripts tsne2 <pairwise_csv> -o <figure.png>"  # output with no args/bad args as well as -h
    epilog = epilog  # output with -h
    formatter_class = argparse.RawTextHelpFormatter  # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)

        subparser.add_argument(
            "pairwise_csv", help="output from 'sourmash scripts pairwise'"
        )
        subparser.add_argument(
            "-C", "--categories-csv", help="CSV mapping label columns to categories"
        )
        subparser.add_argument("-o", "--output-figure", required=True)
        subparser.add_argument("--save-matrix", help="save a numpy matrix")
        subparser.add_argument("--save-labels-to", help="save a labels_to csv")

    def main(self, args):
        super().main(args)

        with sourmash_args.FileInputCSV(args.pairwise_csv) as r:
            rows = list(r)

        # pick out all the distinct queries/matches.
        notify(f"loaded {len(rows)} rows from '{args.pairwise_csv}'")
        sample_d = manysearch_rows_to_index(rows, column_name='both')
        notify(f"loaded {len(sample_d)} total elements")

        mat = numpy.zeros((len(sample_d), len(sample_d)))

        for row in rows:
            # get unique indices for each query/match pair.
            q = row["query_name"]
            qi = sample_d[q]
            m = row["match_name"]
            mi = sample_d[m]
            jaccard = float(row["jaccard"])

            mat[qi, mi] = jaccard
            mat[mi, qi] = jaccard

        numpy.fill_diagonal(mat, 1)

        if args.save_matrix:
            notify(f"writing numpy matrix to '{args.save_matrix}'")
            with open(args.save_matrix, "wb") as fp:
                numpy.save(fp, mat)

        if args.save_labels_to:
            notify(f"writing output labels csv to '{args.save_labels_to}'")
            with open(args.save_labels_to, "w", newline="") as fp:
                w = csv.writer(fp)
                w.writerow(["sort_order", "label"])
                for label, n in sample_d.items():
                    w.writerow([n, label])

        # load categories?
        category_map = None
        colors = None
        if args.categories_csv:
            category_map, colors = load_categories_csv_for_labels(
                args.categories_csv, sample_d
            )

        dissim = 1 - mat
        plot_tsne(dissim, colors=colors, category_map=category_map)

        notify(f"writing figure to '{args.output_figure}'")
        plt.savefig(args.output_figure)


class Command_ClusterToCategories(CommandLinePlugin):
    command = "cluster_to_categories"  # 'scripts <command>'
    description = "convert branchwater plugin 'cluster' output to a betterplot categories CSV"  # output with -h
    usage = "sourmash scripts cluster_to_categories <manysearch_csv> <cluster_csv> -o <categories_csv>"  # output with no args/bad args as well as -h
    epilog = epilog  # output with -h
    formatter_class = argparse.RawTextHelpFormatter  # do not reformat multiline

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('manysearch_csv')
        p.add_argument('cluster_csv')
        p.add_argument('-o', '--output-categories-csv', required=True)

    def main(self, args):
        super().main(args)

        # load samples
        with open(args.manysearch_csv, newline='') as fp:
            r = csv.DictReader(fp)
            rows = list(r)

        samples_d = manysearch_rows_to_index(rows, column_name='both')
        notify(f"loaded {len(samples_d)} samples from '{args.manysearch_csv}'")

        ident_d = {}
        for name, sample_idx in samples_d.items():
            ident = name.split(' ')[0]
            ident_d[ident] = name

        with open(args.cluster_csv, newline='') as fp:
            r = csv.DictReader(fp)
            rows = list(r)

        cluster_to_idents = defaultdict(set)
        n_samples_clustered = 0
        for row in rows:
            cluster = row['cluster']
            nodes = row['nodes'].split(';')
            if len(nodes) == 1:
                cluster = 'unclustered'
            cluster_to_idents[cluster].update(nodes)
            n_samples_clustered += len(nodes)

        notify(f"loaded {len(cluster_to_idents)} clusters containing {n_samples_clustered} members total")
        notify(f"{len(cluster_to_idents['unclustered'])} singletons => 'unclustered'")

        notfound = set(ident_d)

        with open(args.output_categories_csv, 'w', newline='') as fp:
            w = csv.writer(fp)
            w.writerow(['label', 'category'])
            for cluster_name, idents in cluster_to_idents.items():
                for ident in idents:
                    name = ident_d[ident]
                    w.writerow([name, cluster_name])
                notfound -= idents

            if notfound:
                notify(f"{len(notfound)} unmentioned samples => 'unclustered'")
                for ident in notfound:
                    name = ident_d[ident]
                    w.writerow([name, 'unclustered'])


def set_size(x):
    return len(x)


def _venn2_sizes(a, b):
    return (set_size(a - b), set_size(b - a), set_size(a & b))

def _venn3_sizes(a, b, c):
    return (
        set_size(
            a - (b | c)
        ),  # TODO: This is certainly not the most efficient way to compute.
        set_size(b - (a | c)),
        set_size((a & b) - c),
        set_size(c - (a | b)),
        set_size((a & c) - b),
        set_size((b & c) - a),
        set_size(a & b & c),
    )

def set_venn_label(v, loc, label):
    x = v.get_label_by_id(loc)
    if x is not None:
        x.set_text(label)


def format_bp(bp):
    "Pretty-print bp information."
    bp = float(bp)
    if bp < 500:
        return f"{bp:.0f} bp"
    elif bp <= 500e3:
        return f"{round(bp / 1e3, 1):.1f} kbp"
    elif bp < 500e6:
        return f"{round(bp / 1e6, 1):.1f} Mbp"
    elif bp < 500e9:
        return f"{round(bp / 1e9, 1):.1f} Gbp"
    return "???"


class Command_Venn(CommandLinePlugin):
    command = 'venn'
    description = """\
create and write out a pairwise or three-way Venn set overlap diagram.

Calculate and display overlaps between two or three sourmash sketches.
Abundances are ignored.
"""

    usage = """
   sourmash scripts venn <sketches> -o fig.png
"""
    epilog = epilog
    formatter_class = argparse.RawTextHelpFormatter

    def __init__(self, subparser):
        super().__init__(subparser)
        # add argparse arguments here.
        debug_literal('RUNNING cmd venn __init__')
        subparser.add_argument('sketches', nargs='+',
                               help="file(s) containing two or three sketches")
        subparser.add_argument('-o', '--output', default=None,
                               help="save Venn diagram image to this file",
                               required=True)
        subparser.add_argument('--name1', default=None,
                               help="override name for first sketch")
        subparser.add_argument('--name2', default=None,
                               help="override name for second sketch")
        subparser.add_argument('--name3', default=None,
                               help="override name for (optional) third sketch")
        subparser.add_argument('--ident', action='store_true', dest='ident',
                               help="use first space-separated identifier for sequence name")
        add_ksize_arg(subparser, default=31)
        add_moltype_args(subparser)
        add_scaled_arg(subparser)

    def main(self, args):
        # code that we actually run.
        super().main(args)
        moltype = sourmash_args.calculate_moltype(args)

        debug_literal(f'RUNNING cmd {self} {args}')

        sketch_files = list(args.sketches)

        sketches = []
        for filename in sketch_files:
            print_moltype = moltype
            if print_moltype is None:
                print_moltype = '*'
            notify(f"Loading sketches from {filename} with k={args.ksize} moltype={print_moltype}")
            x = list(sourmash.load_file_as_signatures(filename,
                                                      ksize=args.ksize,
                                                      select_moltype=moltype))
            notify(f"...loaded {len(x)} sketches from {filename}.")
            sketches.extend(x)

        if not len(sketches):
            error("ERROR: no sketches found. Must supply 2 or 3.")
            sys.exit(-1)
        elif len(sketches) == 1:
            error("ERROR: only found one sketch. Must supply 2 or 3.")
            sys.exit(-1)
        elif len(sketches) > 3:
            error("ERROR: found more than three sketches. Must supply 2 or 3.")
            sys.exit(-1)

        mh1 = sketches[0].minhash

        scaled = args.scaled
        if scaled is None:
            scaled = mh1.scaled

        mh1 = mh1.downsample(scaled=scaled)
        mh2 = sketches[1].minhash.downsample(scaled=scaled)
        mh1.jaccard(mh2)        # test for general compatibility :)

        hashes1 = set(mh1.hashes)
        hashes2 = set(mh2.hashes)

        label1 = args.name1
        if not label1:
            label1 = sketches[0].name
            if args.ident:
                label1 = sketches[0].name.split(' ')[0]

        label2 = args.name2
        if not label2:
            label2 = sketches[1].name
            if args.ident:
                label2 = sketches[1].name.split(' ')[0]

        if len(sketches) == 2:
            notify("found two sketches - outputting a 2-part Venn diagram.")
            if mh1.track_abundance or mh2.track_abundance:
                notify("NOTE: abundances detected, but not used; try weighted_venn")
            sizes = _venn2_sizes(hashes1, hashes2)
            sizes = [ size * mh1.scaled for size in sizes ]
            v = venn2(sizes, set_labels=(label1, label2))
            set_venn_label(v, '10', format_bp(sizes[0]))
            set_venn_label(v, '01', format_bp(sizes[1]))
            set_venn_label(v, '11', format_bp(sizes[2]))

        elif len(sketches) == 3:
            notify("found three sketches - outputting a 3-part Venn diagram.")
            mh3 = sketches[2].minhash.downsample(scaled=scaled)
            mh1.jaccard(mh3)    # again, test for compatibility

            if mh1.track_abundance or mh2.track_abundance or mh3.track_abundance:
                notify("NOTE: abundances detected, but not used; try weighted_venn")

            hashes3 = set(mh3.hashes)
            label3 = args.name3
            if not label3:
                label3 = sketches[2].name
                if args.ident:
                    label3 = sketches[2].name.split(' ')[0]

            sizes = _venn3_sizes(hashes1, hashes2, hashes3)
            sizes = [ size * mh1.scaled for size in sizes ]
            v = venn3(sizes, set_labels=(label1, label2, label3))
            set_venn_label(v, '100', format_bp(sizes[0]))
            set_venn_label(v, '010', format_bp(sizes[1]))
            set_venn_label(v, '110', format_bp(sizes[2]))
            set_venn_label(v, '001', format_bp(sizes[3]))
            set_venn_label(v, '101', format_bp(sizes[4]))
            set_venn_label(v, '011', format_bp(sizes[5]))
            set_venn_label(v, '111', format_bp(sizes[6]))

        notify(f"saving to '{args.output}'")
        pylab.savefig(args.output)


class Command_WeightedVenn(CommandLinePlugin):
    command = 'weighted_venn'
    description = """\
create and write out a pairwise Venn diagram, weighted by the abundance of
the first sketch.
"""

    usage = """
   sourmash scripts weighted_venn <sketches> -o fig.png
"""
    epilog = epilog
    formatter_class = argparse.RawTextHelpFormatter

    def __init__(self, subparser):
        super().__init__(subparser)
        # add argparse arguments here.
        debug_literal('RUNNING cmd weighted_venn __init__')
        subparser.add_argument('sketches', nargs='+',
                               help="file(s) containing two sketches")
        subparser.add_argument('-o', '--output', default=None,
                               help="save Venn diagram image to this file",
                               required=True)
        subparser.add_argument('--name1', default=None,
                               help="override name for first sketch")
        subparser.add_argument('--name2', default=None,
                               help="override name for second sketch")

        subparser.add_argument('--ident', action='store_true', dest='ident',
                               help="use first space-separated identifier for sequence name")
        add_ksize_arg(subparser, default=31)
        add_moltype_args(subparser)
        add_scaled_arg(subparser)

    def main(self, args):
        # code that we actually run.
        super().main(args)
        moltype = sourmash_args.calculate_moltype(args)

        debug_literal(f'RUNNING cmd {self} {args}')

        sketch_files = list(args.sketches)

        sketches = []
        for filename in sketch_files:
            print_moltype = moltype
            if print_moltype is None:
                print_moltype = '*'
            notify(f"Loading sketches from {filename} with k={args.ksize} moltype={print_moltype}")
            x = list(sourmash.load_file_as_signatures(filename,
                                                      ksize=args.ksize,
                                                      select_moltype=moltype))
            notify(f"...loaded {len(x)} sketches from {filename}.")
            sketches.extend(x)

        if len(sketches) != 2:
            error(f"ERROR: {len(sketches)} sketches found. Must supply exactly 2.")
            sys.exit(-1)

        mh1 = sketches[0].minhash

        scaled = args.scaled
        if scaled is None:
            scaled = mh1.scaled

        mh1 = mh1.downsample(scaled=scaled)
        mh2 = sketches[1].minhash.downsample(scaled=scaled)
        mh1.jaccard(mh2)        # test for general compatibility :)

        label1 = args.name1
        if not label1:
            label1 = sketches[0].name
            if args.ident:
                label1 = sketches[0].name.split(' ')[0]

        label2 = args.name2
        if not label2:
            label2 = sketches[1].name
            if args.ident:
                label2 = sketches[1].name.split(' ')[0]

        notify("found two sketches - outputting a 2-part Venn diagram.")
        notify("Venn diagram will be weighted by abundances in first sketch.")
        if not mh1.track_abundance:
            notify("ERROR: first sketch MUST have abundances.")
            sys.exit(-1)
        if mh2.track_abundance:
            notify("WARNING: abundances on second sketch will be ignored.")

        abunds = mh1.hashes     # dictionary w/abunds
        mh2_hashes = set(mh2.hashes)
        a_sub_b = 0
        b_sub_a = 0
        isect_count = 0

        # count overlapping, weighted by abund of first
        for h in mh2_hashes:
            isect_count += abunds.get(h, 0)

        # count disjoint, second only.
        b_sub_a = len(mh2_hashes - set(abunds))

        # count disjoint, first only, weighted by abund
        for h in set(abunds) - mh2_hashes:
            a_sub_b += abunds[h]

        sizes = [ a_sub_b, b_sub_a, isect_count ]
        sizes = [ size * mh1.scaled for size in sizes ]

        v = venn2(sizes, set_labels=(label1, label2))
        set_venn_label(v, '10', format_bp(sizes[0]))
        set_venn_label(v, '01', format_bp(sizes[1]))
        set_venn_label(v, '11', format_bp(sizes[2]))

        notify(f"saving to '{args.output}'")
        pylab.savefig(args.output)


class Command_PresenceFilter(CommandLinePlugin):
    command = 'presence_filter'
    description = """\
Provide a filtered view of 'gather' output, plotting detection or ANI
against average abund for significant matches.
"""

    usage = """
   sourmash scripts presence_filter gather.csv -o presence.png
"""
    epilog = epilog
    formatter_class = argparse.RawTextHelpFormatter

    def __init__(self, subparser):
        super().__init__(subparser)
        # add argparse arguments here.
        subparser.add_argument('gather_csv')
        subparser.add_argument('-o', '--output', default=None,
                               help="save image to this file",
                               required=True)
        subparser.add_argument('-N', '--min-num-hashes', type=int,
                               default=3, help='threshold (default: 3)')
        subparser.add_argument('--min-ani', type=float, default=0.0,
                               help='ANI threshold (default: None)')
        subparser.add_argument('--min-fraction', type=float, default=0.0,
                               help='detection threshold (default: None)')
        subparser.add_argument('--detection', action="store_true",
                               default=True)
        subparser.add_argument('--detection-column-name',
                               default='f_match_orig')
        subparser.add_argument('--ani', dest='detection',
                               action="store_false")
        subparser.add_argument('--green-color',
                               help="color genomes with matching names green")
        subparser.add_argument('--red-color',
                               help="color genomes with matching names red")
        subparser.add_argument('--blue-color',
                               help="color genomes with matching names blue")

    def main(self, args):
        df = pd.read_csv(args.gather_csv)
        notify(f"loaded {len(df)} rows from '{args.gather_csv}'")

        scaled = set(df['scaled'])
        assert len(scaled) == 1
        scaled = list(scaled)[0]

        if 'name' in df.columns:
            df['match_name'] = df['name'] # correct for gather/fastgather column names

        notify(f"loaded {len(df)} rows.")
        if args.min_num_hashes:
            threshold = args.min_num_hashes * scaled
            df = df[df['unique_intersect_bp'] >= threshold]
            notify(f"filtered down to {len(df)} rows with unique_intersect_bp >= {threshold}")

        if args.min_ani:
            df = df[df['match_containment_ani'] >= args.min_ani]
            notify(f"filtered down to {len(df)} rows with match_containment_ani >= {args.min_ani} (--min-ani)")

        if args.min_fraction:
            df = df[df[args.detection_column_name] >= args.min_fraction]
            notify(f"filtered down to {len(df)} rows with {args.detection_column_name} >= {args.min_fraction} (--min-fraction)")

        if args.detection:
            plt.plot(df[args.detection_column_name], df.average_abund, 'k.')
        else:
            plt.plot(df.match_containment_ani, df.average_abund, 'k.')

        dfs = []
        colors = []
        if args.green_color:
            df2 = df[df['match_name'].str.contains(args.green_color)]
            notify(f"{len(df2)} matches to {args.green_color} => green circles")
            dfs.append(df2)
            colors.append('go')
        if args.red_color:
            df2 = df[df['match_name'].str.contains(args.red_color)]
            notify(f"{len(df2)} matches to {args.red_color} => red crosses")

            dfs.append(df2)
            colors.append('r+')
        if args.blue_color:
            df2 = df[df['match_name'].str.contains(args.blue_color)]
            notify(f"{len(df2)} matches to {args.blue_color} => blue triangles")
            dfs.append(df2)
            colors.append('bv')

        for (df2, color) in zip(dfs, colors):
            if args.detection:
                plt.plot(df2.f_match_orig, df2.average_abund, color)
            else:               # ANI!
                plt.plot(df2.match_containment_ani, df2.average_abund, color)

        ax = plt.gca()
        ax.set_yscale('log')

        if args.detection:
            ax.set_xlabel('fraction of genome detected')
            ax.set_xlim((0, 1))
        else:
            ax.set_xlabel('cANI of match')

        ax.set_ylabel('log abundance (copy number)')

        notify(f"saving figure to '{args.output}'")
        plt.tight_layout()
        plt.savefig(args.output)


def save_sankey_diagram(fig, output_file):
        if output_file:
            if output_file.endswith(".html"):
                fig.write_html(output_file)
                notify(f"Saved interactive HTML: {output_file}")
            elif output_file.endswith((".png", ".jpg", ".jpeg", ".pdf", ".svg")):
                fig.write_image(output_file)
                notify(f"Saved image file: {output_file}")
            else:
                notify("Unsupported file format. Use .html, .png, .jpg, .jpeg, .pdf, or .svg.")
        else:
            fig.show()  # Show the plot if no output file is specified
    
def process_csv_for_sankey(input_csv, csv_type):
    nodes = []  # List of unique taxonomy nodes
    node_map = {}  # Map taxonomic label to index
    links = []  # List of link connections with flow values
    hover_texts = []  # Custom hover text for percentages
    processed_lineages = set()  # Tracks added lineage links

    # Read CSV file
    with open(input_csv, 'r') as inF:
        reader = csv.DictReader(inF)
        data = list(reader)

    # Determine the appropriate headers based on csv_type
    if csv_type == "csv_summary":
        fraction_key = "f_weighted_at_rank"
    elif csv_type == "with-lineages":
        fraction_key = "f_unique_weighted"
    else:
        raise ValueError("Invalid csv_type. Use 'csv_summary' or 'with-lineages'.")

    # Process each row in the dataset
    for n, row in enumerate(data):
        fraction = float(row[fraction_key]) * 100  # Convert to percentage
        lineage_parts = row["lineage"].split(";")  # Taxonomic hierarchy

        # Iterate through lineage levels and create source-target links
        for i in range(len(lineage_parts) - 1):
            source_label = lineage_parts[i].strip()
            target_label = lineage_parts[i + 1].strip()

            # Since 'tax metagenome' is already summarized, skip duplicates to prevent overcounting
            if csv_type == "csv_summary" and (source_label, target_label) in processed_lineages:
                continue

            # Assign indices to nodes
            if source_label not in node_map:
                node_map[source_label] = len(nodes)
                nodes.append(source_label)

            if target_label not in node_map:
                node_map[target_label] = len(nodes)
                nodes.append(target_label)

            # Create a link between source and target
            links.append({
                "source": node_map[source_label],
                "target": node_map[target_label],
                "value": fraction
            })
            processed_lineages.add((source_label, target_label))  # Track added links
            hover_texts.append(f"{source_label} → {target_label}<br>{fraction:.2f}%")
    notify(f"loaded {n+1} rows from '{input_csv}'")

    return nodes, links, hover_texts

class Command_Sankey(CommandLinePlugin):
    command = 'sankey'
    description = """\
Build a sankey plot to visualize taxonomic profiling. Uses sourmash 'gather' --> 'tax' output ('tax metagenome' csv_summary format or 'tax annotate' output).
"""

    usage = """
   sourmash scripts sankey --summary-csv gather.summarized.csv -o gather.html
"""
    epilog = epilog
    formatter_class = argparse.RawTextHelpFormatter

    def __init__(self, subparser):
        super().__init__(subparser)
        # add argparse arguments here.
        group = subparser.add_mutually_exclusive_group(required=True)
        group.add_argument("--summary-csv", type=str, help="Path to csv_summary generated by running 'sourmash tax metagenome' on a sourmash gather csv")
        group.add_argument("--annotate-csv", type=str, help="Path to 'with-lineages' file generated by running 'sourmash tax annotate' on a sourmash gather csv")
        
        subparser.add_argument("-o", "--output", type=str, help="output file for alluvial flow diagram")
        subparser.add_argument("--title", type=str, help="Plot title (default: use input filename)")

        subparser.epilog = "You must provide either --summary-csv or --annotate-csv, but not both."

    def main(self, args):
        # Build info appropriately based on input file type
        if args.summary_csv:
            input_csv = args.summary_csv
            csv_type = "csv_summary"
            required_headers = ["f_weighted_at_rank", "lineage"]
        else:
            input_csv = args.annotate_csv
            csv_type = "with-lineages"
            required_headers = ["f_unique_weighted", "lineage"]

        # Check if the required headers are present
        with open(input_csv, 'r') as file:
            reader = csv.DictReader(file)
            if not all(header in reader.fieldnames for header in required_headers):
                raise ValueError(f"Expected headers {required_headers} not found. Is this a correct file for '{csv_type}' type?")

        # process csv
        nodes, links, hover_texts = process_csv_for_sankey(input_csv, csv_type)
        base_title = os.path.basename(input_csv.rsplit(".csv")[0])

        # Create Sankey diagram
        fig = go.Figure(go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                label=nodes
            ),
            link=dict(
                source=[link["source"] for link in links],
                target=[link["target"] for link in links],
                value=[link["value"] for link in links],
                customdata=hover_texts,
                hovertemplate="%{customdata}<extra></extra>"  # Use custom hover text
            )
        ))

        if args.title:
            title = args.title
        else:
            title = base_title 
        fig.update_layout(title_text=f"{title}",
                        font_size=10,
                        autosize=False,
                        width=1500,  # Increase width
                        height=900   # Increase height
                        )

        # Save output based on file extension
        save_sankey_diagram(fig, args.output)


def read_compare_matrix(input_matrix, labels_from, matrix_type="similarity"):
    """Read and process a similarity/distance matrix and labels file."""
    with open(input_matrix, "rb") as f:
        D = numpy.load(f)
    notify(f"...got {D.shape[0]} x {D.shape[1]} matrix.", *D.shape)

    labelfilename = labels_from
    notify(f"loading labels from CSV file '{labelfilename}'")

    labelinfo = load_labelinfo_csv(labelfilename)

    if len(labelinfo) != D.shape[0]:
        error("{} labels != matrix size, exiting", len(labelinfo))
        sys.exit(-1)

    if matrix_type == "similarity":
        matrix = 1 - D  # Convert similarity (e.g. jaccard) to distance

    labels = [row["label"] for row in labelinfo]

    # replace commas with '_' because commas are treated differently for plotting
    labels = [x.replace(',', '_') for x in labels]

    return matrix, labels

def pairwise_to_matrix(input_csv, use_column='jaccard'):
    """Convert pairwise CSV to a distance matrix."""
    with open(input_csv, 'r') as file:
        reader = csv.DictReader(file, quotechar='"', delimiter=',')
        rows = list(reader)

    sample_names = sorted(set(row['query_name'].strip() for row in rows) |
                          set(row['match_name'].strip() for row in rows))
    sample_d = {name: idx for idx, name in enumerate(sample_names)}
    notify(f"...found {len(sample_d.keys())} (total {len(rows)} pairwise comparisons).")

    # Initialize distance matrix
    matrix = numpy.full((len(sample_d), len(sample_d)), numpy.float64(1.0))
    numpy.fill_diagonal(matrix, 0.0)

    for row in rows:
        q, m = row['query_name'], row['match_name']
        qi, mi = sample_d[q], sample_d[m]
        value = numpy.float64(row[use_column])
        distance = numpy.float64(1 - value)  # Convert similarity to distance

        matrix[qi, mi] = distance
        matrix[mi, qi] = distance

    labels = [x.replace(',', '_') for x in sample_d.keys()]  # Replace commas for safe plotting
    return matrix, labels


def build_nj_tree(distance_matrix, labels):
    """Build a Neighbor-Joining tree from a matrix using BioPython."""
    # Convert full matrix to lower triangle format
    lower_triangle = []
    for i in range(len(labels)):
        lower_triangle.append([distance_matrix[i, j] for j in range(i + 1)])

    # Create BioPython DistanceMatrix
    dm = DistanceMatrix(names=labels, matrix=lower_triangle)

    # Construct NJ tree using BioPython
    constructor = DistanceTreeConstructor()
    nj_tree = constructor.nj(dm)

    return nj_tree  # Return BioPython tree directly


def save_tree(tree, output_file):
    """Save tree to file in Newick format using BioPython."""
    if output_file.endswith(".nwk"):
        Phylo.write(tree, output_file, "newick")
        print(f"Saved Newick tree: {output_file}")
    else:
        print("Unsupported file format. Use .nwk for Newick format.")


def plot_tree_ete(tree, layout, output_image=None, show=False):
    """Render and save tree image using ete3."""
    try:
        from ete3 import Tree, TreeStyle
    except ImportError:
        print("** WARNING: could not import TreeStyle; maybe PyQT5 is not installed?",
              file=sys.stderr)
        print("** Will not be able to output trees. About to fail in 1... 2... 3...",
               file=sys.stderr)

    ete_tree = Tree(tree.format('newick'), format=1)
    ts = TreeStyle()
    ts.show_leaf_name = True
    if layout == "circular":
        ts.mode = "c"
    else:
        ts.mode = "r"  # Default to rectangular

    if output_image:
        if output_image.endswith((".png", ".jpg", ".jpeg", ".svg", ".pdf")):
            ete_tree.render(output_image, tree_style=ts)
            print(f"Saved tree image: {output_image}")
        else:
            print("Unsupported file format. Use .png, .jpg, .jpeg, .svg, or .pdf.")
    if show:
        ete_tree.show(tree_style=ts)

class Command_DistTree(CommandLinePlugin):
    command = 'tree'
    description = """\
Build a neighbor-joining tree from 'sourmash compare' or 'sourmash scripts pairwise' output.
"""

    usage = """
   sourmash scripts tree --compare-matrix matrix.np -o output.nwk --image tree.png
"""
    epilog = epilog
    formatter_class = argparse.RawTextHelpFormatter

    def __init__(self, subparser):
        super().__init__(subparser)
        # add argparse arguments here.
        group = subparser.add_mutually_exclusive_group(required=True)
        group.add_argument("--compare-matrix", type=str, help="Path to the distance matrix (numpy .npy file)")
        group.add_argument("--pairwise-csv", type=str, help="Path to a pairwise CSV file")
        subparser.add_argument(
            "--labels-from", help="output from 'sourmash compare --labels-to'"
        )
        subparser.add_argument("--matrix-type", type=str, choices=["similarity", "distance"], default="similarity", help="Are matrix values 'similarity' (e.g. jaccard, containment) or 'distance' (1-jaccard, 1-containment, etc)? default: 'similarity' (also default output for sourmash compare).")
        subparser.add_argument("--use-column", type=str, default="jaccard", choices=["jaccard", "max_containment"], help="column name to use in pairwise CSV (default: jaccard)",)
        subparser.add_argument("--newick", type=str, help="Output tree in  Newick format. File must end in '.nwk'")
        subparser.add_argument("-o", "--output", type=str, help="Output file for tree image (.png, .jpg, .jpeg, .svg, .pdf)")
        subparser.add_argument("--show", action="store_true", default=False, help="Open tree image in ETE3 browser window (default=False)")
        subparser.add_argument("--tree-layout", type=str, choices=["rectangular", "circular"], default="rectangular", help="Tree layout (rectangular or circular)")

        subparser.epilog = "You must provide either --compare-matrix and --labels-from or --pairwise-csv, but not both."

    def main(self, args):
        if args.compare_matrix:
            if not args.labels_from:
                notify("Must provide --labels-from when using --compare-matrix")
                sys.exit(-1)
            matrix, labels = read_compare_matrix(args.compare_matrix, args.labels_from, args.matrix_type)
        elif args.pairwise_csv:
            matrix, labels = pairwise_to_matrix(args.pairwise_csv, args.use_column)

        tree = build_nj_tree(matrix, labels)
        if args.newick:
            save_tree(tree, args.newick)
        if not args.show:
            os.environ['QT_QPA_PLATFORM']='offscreen'
        if args.show or args.output:
            plot_tree_ete(tree, args.tree_layout, output_image=args.output, show=args.show)


class Command_TreeMap(CommandLinePlugin):
    command = 'treemap'
    description = """\
Build a treemap with proportional representation of a metagenome taxonomy.
"""

    usage = """
   sourmash scripts treemap csv_summary -o treemap.png [ -r <rank> ]
"""
    epilog = epilog
    formatter_class = argparse.RawTextHelpFormatter

    def __init__(self, subparser):
        super().__init__(subparser)
        subparser.add_argument('csvfile', help='csv_summary output from tax metagenome')
        subparser.add_argument('-o', '--output', required=True,
                               help='output figure to this file')
        subparser.add_argument('-r', '--rank', default='phylum',
                               help='display at this rank')
        subparser.add_argument('-n', '--num-to-display', type=int,
                               default=25,
                               help="display at most these many taxa; aggregate the remainder (default: 25; 0 to display all)")

        
    def main(self, args):
        super().main(args)
        plot_treemap(args)


def plot_treemap(args):
    from string import ascii_lowercase
    import itertools
    cmap = colormaps['viridis']

    df = pd.read_csv(args.csvfile)

    print(f"reading input file '{args.csvfile}'")
    for colname in ('query_name', 'rank', 'f_weighted_at_rank', 'lineage'):
        if colname not in df.columns:
            print(f"input is missing column '{colname}'; is this a csv_summary file?")
            sys.exit(-1)

    df = df.sort_values(by='f_weighted_at_rank')

    # select rank
    df2 = df[df['rank'] == args.rank]
    df2['name'] = df2['lineage'].apply(lambda x: x.split(';')[-1])

    fractions = list(df2['f_weighted_at_rank'].tolist())
    names = list(df2['name'].tolist())
    fractions.reverse()
    names.reverse()

    num = max(args.num_to_display, 0) # non-negative
    num = min(args.num_to_display, len(names)) # list of names
    if num:
        display_fractions = fractions[:num]
        display_names = names[:num]

        print(f'treemap: displaying {num} taxa of {len(fractions)} total')
        if fractions[num:]:
            remaining = fractions[num:]
            remainder = sum(remaining)
            display_fractions.append(remainder)
            display_names.append(f'{len(remaining)} remaining taxa')
            print(f'aggregating {len(remaining)} remaining taxa into one box')

        fractions, names = display_fractions, display_names

    # use to build labels: a, b, c..., aa, ab, ac...
    def iter_all_strings():
        for size in itertools.count(1):
            for s in itertools.product(ascii_lowercase, repeat=size):
                yield "".join(s)

    labels = []
    for name, label in zip(names, iter_all_strings()):
        labels.append(label)

    # Create treemap
    fig, ax = plt.subplots(figsize=(8, 6))

    squarify.plot(sizes=fractions, label=labels,
                  alpha=0.7, ax=ax,
                  color=cmap(numpy.linspace(1, 0, len(labels)))
    )
    plt.axis('off')

    # Position the text on the right side
    # The x-coordinate (1.05) places the text slightly outside the right edge of the axes.
    # The y-coordinate is calculated to distribute the text vertically.
    # 'ha' (horizontal alignment) is set to 'left' to align the text from its left edge.
    # 'va' (vertical alignment) is set to 'center' to center the text vertically.
    # 'transform=ax.transAxes' ensures coordinates are relative to the axes.

    for i, (name, label, f) in enumerate(zip(names, labels, fractions)):
        y_position = 0.95 - (i * 0.035)  # Adjust for desired spacing and starting point
        ax.text(1.03, y_position, f'{label}: {name} ({f*100:.1f}%)',
                ha='left', va='center', transform=ax.transAxes, fontsize=10)

    plt.tight_layout()

    print(f"saving output to '{args.output}'")
    plt.savefig(args.output)
