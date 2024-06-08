"""
Tests for sourmash_plugin_betterplot.
"""
import os
import pytest

import sourmash
import sourmash_tst_utils as utils
from sourmash_tst_utils import SourmashCommandFailed
from sourmash import sourmash_args

from sourmash_plugin_betterplot import (load_labelinfo_csv,
                                        load_categories_csv_for_labels,
                                        load_categories_csv,
                                        manysearch_rows_to_index)


def test_run_sourmash(runtmp):
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('', fail_ok=True)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert runtmp.last_result.status != 0                    # no args provided, ok ;)


def test_load_categories():
    labels_csv = utils.get_test_data('10sketches.cmp.labels_to.csv')
    cat_csv = utils.get_test_data('10sketches-categories.csv')

    labelinfo = load_labelinfo_csv(labels_csv)
    category_map, colors = load_categories_csv(cat_csv, labelinfo)

    assert len(labelinfo) == 10    # list: number of labels
    assert len(category_map) == 5  # dict: categories -> colors mapping
    assert len(colors) == 10       # list: labels => colors, in sort order


def test_load_categories_for_labels():
    pairwise_csv = utils.get_test_data('10sketches.pairwise.csv')
    cat_csv = utils.get_test_data('10sketches-categories.csv')

    with sourmash_args.FileInputCSV(pairwise_csv) as r:
        rows = list(r)

    sample_d = manysearch_rows_to_index(rows, column_name='query_name')

    category_map, colors = load_categories_csv_for_labels(cat_csv, sample_d)

    assert len(sample_d) == 10     # dict: sample name to sample index
    assert len(category_map) == 5  # dict: categories -> colors mapping
    assert len(colors) == 10       # list: labels => colors, in sort order
