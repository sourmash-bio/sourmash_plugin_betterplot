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
from matplotlib_venn import venn2, venn3
import scipy.cluster.hierarchy as sch
from sklearn.manifold import MDS, TSNE
from scipy.sparse import lil_matrix, csr_matrix
from matplotlib.lines import Line2D
import seaborn as sns
import upsetplot
import pandas as pd

import sourmash
from sourmash import sourmash_args
from sourmash.logging import debug_literal, error, notify, print_results
from sourmash.plugins import CommandLinePlugin
import sourmash_utils
from sourmash.cli.utils import (add_ksize_arg, add_moltype_args)


### utility functions

def load_labelinfo_csv(filename):
    "Load file output by 'sourmash compare --labels-to'"
    with sourmash_args.FileInputCSV(filename) as r:
        labelinfo = list(r)

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
            error("{} labels != matrix size, exiting", len(labeltext))
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

'venn' provides figure output.
"""

    usage = """
   sourmash scripts venn -k 31 <sketches>
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
        subparser.add_argument('--ident', action='store_true',
                               help="use first space-separated identifier for sequence name")
        add_ksize_arg(subparser)
        add_moltype_args(subparser)

    def main(self, args):
        # code that we actually run.
        super().main(args)
        moltype = sourmash_args.calculate_moltype(args)

        debug_literal(f'RUNNING cmd {self} {args}')

        sketch_files = list(args.sketches)

        sketches = []
        for filename in sketch_files:
            notify(f"Loading sketches from {filename}")
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
        mh2 = sketches[1].minhash
        assert mh1.scaled == mh2.scaled

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
            sizes = _venn2_sizes(hashes1, hashes2)
            sizes = [ size * mh1.scaled for size in sizes ]
            v = venn2(sizes, set_labels=(label1, label2))
            set_venn_label(v, '10', format_bp(sizes[0]))
            set_venn_label(v, '01', format_bp(sizes[1]))
            set_venn_label(v, '11', format_bp(sizes[2]))

        elif len(sketches) == 3:
            notify("found three sketches - outputting a 3-part Venn diagram.")
            mh3 = sketches[2].minhash
            assert mh1.scaled == mh3.scaled

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

        if args.output:
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
        subparser.add_argument('-N', '--min-num-hashes',
                               default=3, help='threshold (default: 3)')
        subparser.add_argument('--detection', action="store_true",
                               default=True)
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

        threshold = args.min_num_hashes * scaled
        df = df[df['unique_intersect_bp'] >= threshold]
        notify(f"filtered down to {len(df)} rows with unique_intersect_bp >= {threshold}")

        if args.detection:
            plt.plot(df.f_match_orig, df.average_abund, 'k.')
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
            else:
                plt.plot(df2.match_containment_ani, df2.average_abund, color)

        ax = plt.gca()
        ax.set_ylabel('number of copies')
        ax.set_yscale('log')

        if args.detection:
            ax.set_xlabel('fraction of genome detected')
        else:
            ax.set_xlabel('cANI of match')

        notify(f"saving figure to '{args.output}'")
        plt.savefig(args.output)
