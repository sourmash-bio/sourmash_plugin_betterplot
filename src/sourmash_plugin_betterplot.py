"""betterplot plugin implementation"""

epilog = """

Need help? Have questions? Ask at http://github.com/sourmash-bio/sourmash/issues!
"""

import argparse
import sourmash
from sourmash import sourmash_args
import os
import csv
from collections import defaultdict


import numpy
import pylab
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from sklearn.manifold import MDS
from scipy.sparse import lil_matrix, csr_matrix
from matplotlib.lines import Line2D

from sourmash.logging import debug_literal, error, notify, print_results
from sourmash.plugins import CommandLinePlugin


###


def load_labelinfo_csv(filename):
    with sourmash_args.FileInputCSV(filename) as r:
        labelinfo = list(r)

    labelinfo.sort(key=lambda row: int(row["sort_order"]))
    return labelinfo


def load_categories_csv(filename, labelinfo):
    with sourmash_args.FileInputCSV(filename) as r:
        categories = list(r)

    category_map = {}
    colors = None
    if categories:
        assert labelinfo
        keys = set(categories[0].keys())
        keys -= {"category"}

        key = None
        for k in keys:
            if k in labelinfo[0].keys():
                notify(f"found category key: {k}")
                key = k
                break

        if key:
            category_values = list(set([row["category"] for row in categories]))
            category_values.sort()

            cat_colors = list(map(plt.cm.tab10, range(len(category_values))))
            category_map = {}
            for v, color in zip(category_values, cat_colors):
                category_map[v] = color

            category_map2 = {}
            for row in categories:
                category_map2[row[key]] = category_map[row["category"]]

            colors = []
            for row in labelinfo:
                value = row[key]
                color = category_map2[value]
                colors.append(color)

        else:
            notify(f"no valid key column found in categories file '{filename}'.")
    else:
        notify(f"nothing in categories file '{filename}'?!")

    return category_map, colors


def load_categories_csv_for_labels(filename, queries):
    "Load a categories CSV that must use label name."
    with sourmash_args.FileInputCSV(filename) as r:
        categories = list(r)

    category_map = {}
    colors = None
    if categories:
        key = "label"

        category_values = list(set([row["category"] for row in categories]))
        category_values.sort()

        cat_colors = list(map(plt.cm.tab10, range(len(category_values))))
        category_map = {}
        for v, color in zip(category_values, cat_colors):
            category_map[v] = color

        category_map2 = {}
        for row in categories:
            label = row[key]
            cat = row["category"]
            category_map2[label] = category_map[cat]

        colors = []
        for label, idx in queries:
            color = category_map2[label]
            colors.append(color)
    else:
        notify(f"nothing in categories file '{filename}'?!")

    return category_map, colors


#
# CLI plugin - supports 'sourmash scripts plot2'
#


class Command_Plot2(CommandLinePlugin):
    command = "plot2"  # 'scripts <command>'
    description = "plot a distance matrix produced by 'sourmash compare'"  # output with -h
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
    Z1 = sch.dendrogram(
        Y,
        orientation="left",
        labels=labeltext,
        no_labels=not show_labels,
        get_leaves=True,
    )
    # ax1.set_xticks([])

    if cut_point is not None:
        ax1.axvline(x=cut_point, c="red", linestyle="dashed")

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
    notify(f"wrote numpy distance matrix to: {args.output_figure}")

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
        plot_mds(dissim, colors=colors, category_map=category_map)

        plt.savefig(args.output_figure)


class Command_PairwiseToCompare(CommandLinePlugin):
    command = "pairwise_to_compare"  # 'scripts <command>'
    description = "convert pairwise CSV output to a 'compare' matrix"  # output with -h
    usage = "sourmash scripts pairwise_to_compare <pairwise_csv> -o <matrix_cmp>"  # output with no args/bad args as well as -h
    epilog = epilog  # output with -h
    formatter_class = argparse.RawTextHelpFormatter  # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)

        subparser.add_argument(
            "pairwise_csv", help="output from 'sourmash scripts pairwise'"
        )
        subparser.add_argument('-o', '--output-matrix', required=True)
        subparser.add_argument('--labels-to')

    def main(self, args):
        super().main(args)

        with sourmash_args.FileInputCSV(args.pairwise_csv) as r:
            rows = list(r)

        # pick out all the distinct queries/matches.
        print(f"loaded {len(rows)} rows from '{args.pairwise_csv}'")
        queries = set( [ row['query_name'] for row in rows ] )
        queries.update(set( [ row['match_name'] for row in rows ] ))
        print(f"loaded {len(queries)} total elements")

        queries = list(sorted(queries))

        sample_d = {}
        for n, sample_name in enumerate(queries):
            sample_d[sample_name] = n

        assert n == len(queries) - 1

        mat = numpy.zeros((len(queries), len(queries)))

        pairs = []
        for row in rows:
            # get unique indices for each query/match pair.
            q = row['query_name']
            qi = sample_d[q]
            m = row['match_name']
            mi = sample_d[m]
            jaccard = float(row['jaccard'])

            mat[qi, mi] = jaccard
            mat[mi, qi] = jaccard

        numpy.fill_diagonal(mat, 1)

        with open(args.output_matrix, 'wb') as fp:
            numpy.save(fp, mat)

        with open(args.output_matrix + '.labels.txt', 'wt') as fp:
            for label, n in sample_d.items():
                fp.write(label + "\n")

        if args.labels_to:
            with open(args.labels_to, 'w', newline="") as fp:
                w = csv.writer(fp)
                w.writerow(['sort_order', 'label'])
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

    def main(self, args):
        super().main(args)

        with sourmash_args.FileInputCSV(args.pairwise_csv) as r:
            rows = list(r)

        # pick out all the distinct queries/matches.
        print(f"loaded {len(rows)} rows from '{args.pairwise_csv}'")
        queries = set( [ row['query_name'] for row in rows ] )
        queries.update(set( [ row['match_name'] for row in rows ] ))
        print(f"loaded {len(queries)} total elements")

        queries = list(sorted(queries))

        sample_d = {}
        for n, sample_name in enumerate(queries):
            sample_d[sample_name] = n

        assert n == len(queries) - 1

        mat = numpy.zeros((len(queries), len(queries)))

        pairs = []
        for row in rows:
            # get unique indices for each query/match pair.
            q = row['query_name']
            qi = sample_d[q]
            m = row['match_name']
            mi = sample_d[m]
            jaccard = float(row['jaccard'])

            mat[qi, mi] = jaccard
            mat[mi, qi] = jaccard

        numpy.fill_diagonal(mat, 1)

        # load categories?
        category_map = None
        colors = None
        if args.categories_csv:
            category_map, colors = load_categories_csv_for_labels(args.categories_csv, sample_d.items())

        dissim = 1 - mat
        plot_mds(dissim, colors=colors, category_map=category_map)

        plt.savefig(args.output_figure)


#@CTB unused again...
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


def plot_mds(matrix, *, colors=None, category_map=None):
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
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
