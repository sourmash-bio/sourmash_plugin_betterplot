"""
Tests for sourmash_plugin_betterplot.
"""
import os
import pytest

import sourmash
import sourmash_tst_utils as utils
from sourmash_tst_utils import SourmashCommandFailed
from sourmash import sourmash_args
from collections import defaultdict

from sourmash_plugin_betterplot import (load_labelinfo_csv,
                                        load_categories_csv_for_labels,
                                        load_categories_csv,
                                        manysearch_rows_to_index,
                                        expand_with_ancestors_sum,
                                        is_lins_lineage,
                                        detect_lins,
                                        rows_to_edges,
                                        edges_to_links,
                                        build_links_taxonomy,
                                        strip_prefix,
                                        path_to_display,
                                        make_hover,
                                        process_csv_for_sankey)


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


def test_expand_with_ancestors_sum_simple_taxonomy():
    rows = [{"lineage": "a;b;c", "frac": 1.0}]
    expanded = expand_with_ancestors_sum(rows, "frac")

    # Should include all ancestors: a, a;b, a;b;c
    labels = {r["lineage"] for r in expanded}
    assert {"a", "a;b", "a;b;c"} <= labels

    # Fractions should propagate upward
    lineage_map = {r["lineage"]: r["frac"] for r in expanded}
    assert lineage_map["a;b;c"] == 1.0
    assert lineage_map["a;b"] == 1.0
    assert lineage_map["a"] == 1.0

def test_expand_with_ancestors_sum_multiple_children():
    rows = [
        {"lineage": "a;b;c1", "frac": 2.0},
        {"lineage": "a;b;c2", "frac": 3.0},
    ]
    expanded = expand_with_ancestors_sum(rows, "frac")
    lineage_map = {r["lineage"]: r["frac"] for r in expanded}

    # parent a;b should sum c1 + c2
    assert lineage_map["a;b"] == pytest.approx(5.0)
    # root a should also be 5.0
    assert lineage_map["a"] == pytest.approx(5.0)

def test_is_and_detect_lins():
    assert is_lins_lineage("0;1;2")
    assert not is_lins_lineage("a;1;2")

    rows = [{"lineage": "0;1;1"}, {"lineage": "0;2;3"}]
    assert detect_lins(rows)

    rows = [{"lineage": "a;b;c"}]
    assert not detect_lins(rows)

def test_rows_to_edges_basic():
    rows = [
        {"lineage": "a;b", "frac": 0.2},
        {"lineage": "a;b;c", "frac": 0.15},
        {"lineage": "a;b;d", "frac": 0.05},
    ]
    edges = rows_to_edges(rows, "frac", lins=False)
    # Expect edges (a->a;b, a;b->a;b;c)
    edge_paths = {(src, tgt) for src, tgt, _ in edges}
    assert ("a;b", "a;b") in edge_paths
    assert ("a;b", "a;b;c") in edge_paths
    assert ("a;b", "a;b;d") in edge_paths
    assert len(edges) == 3
    # Check percents for each edge
    edge_map = {(src, tgt): frac for src, tgt, frac in edges}
    assert edge_map[("a;b", "a;b")] == pytest.approx(20.0)
    assert edge_map[("a;b", "a;b;c")] == pytest.approx(15.0)
    assert edge_map[("a;b", "a;b;d")] == pytest.approx(5.0)


def test_rows_to_edges_with_lins_prefixing():
    rows = [{"lineage": "0;1", "frac": 1.0}]
    edges = rows_to_edges(rows, "frac", lins=True)
    print(edges)
    # should prefix with p0:/p1:
    src, tgt, frac = edges[0]
    assert src.startswith("p0:")
    assert tgt.startswith("p0:")
    assert frac == 100.0  # fraction converted to percent


def test_edges_to_links_roundtrip():
    edges = [("a", "a;b", 50.0), ("a;b", "a;b;c", 25.0)]
    labels, links, hovers = edges_to_links(edges)

    assert "a" in "".join(labels)
    assert any("a → a;b" in h for h in hovers)
    # Links use numeric indices
    assert all(isinstance(link["source"], int) for link in links)


def test_build_links_taxonomy_simple():
    rows = [
        {"lineage": "a;b;c", "f_unique_weighted": 0.5},
        {"lineage": "a;b;d", "f_unique_weighted": 0.5},
    ]
    nodes, links, hovers = build_links_taxonomy(rows, "f_unique_weighted", csv_type="with-lineages")
    assert "a" in nodes and "b" in nodes
    assert any("a → b" in h for h in hovers)


def test_strip_prefix_and_display():
    path = "p0:1;p1:2"
    stripped = strip_prefix(path)
    assert stripped == "1;2"

    # test display mapping
    lin2name = {"1;2": "GroupX"}
    disp = path_to_display("1;2", lin2name)
    assert "GroupX" in disp

def test_make_hover_with_lin2name():
    lin2name = {"a;b": "GroupAB"}
    hover = make_hover("a", "a;b", 12.3456, lin2name)
    assert "GroupAB" in hover
    assert "12.35%" in hover


def test_rows_to_edges_collapses_passthrough():
    rows = [
        {"lineage": "a;b", "frac": 0.5},
        {"lineage": "a;b;c", "frac": 0.5},
    ]
    edges = rows_to_edges(rows, "frac", lins=False, lin2name={})
    edge_paths = {(src, tgt) for src, tgt, _ in edges}
    # 'a;b' collapsed, so only "a;b;c → a;b;c"
    assert ("a;b;c", "a;b;c") in edge_paths


def test_rows_to_edges_collapses_passthrough_2():
    rows = [
        {"lineage": "a;b", "frac": 0.5},
        {"lineage": "a;b;c", "frac": 0.5},
        {"lineage": "a;b;c;d", "frac": 0.2},
        {"lineage": "a;b;c;e", "frac": 0.3},
        {"lineage": "a;b;c;e;f", "frac": 0.3},
    ]
    edges = rows_to_edges(rows, "frac", lins=False, lin2name={})
    edge_paths = {(src, tgt) for src, tgt, _ in edges}
    # 'a;b' collapsed, so root is "a;b;c → a;b;c"
    assert ("a;b;c", "a;b;c") in edge_paths
    # 'a;b;c;e' collapsed, so "a;b;c → a;b;c;e;f"
    assert ("a;b;c", "a;b;c;e;f") in edge_paths
    # 'a;b;c;d' NOT collapsed, so "a;b;c → a;b;c;d"
    assert ("a;b;c", "a;b;c;d") in edge_paths
    assert len(edges) == 4  # 4 edges total
    # note: self edges are ignored for display
    # edges we have: [('a;b;c', 'a;b;c', 50.0), ('a;b;c', 'a;b;c;d', 20.0), ('a;b;c', 'a;b;c;e;f', 30.0), ('a;b;c;e;f', 'a;b;c;e;f', 30.0)]


def test_process_csv_for_sankey_multiple_query_names_annotate(runtmp):
    input_csv = utils.get_example_data('tax/test.gather.with-lineages.csv')
    mult_query = runtmp.output('mult-query.csv')

    # modify input to have multiple query_name values
    with open(input_csv, newline="") as inF, open(mult_query, 'w', newline='') as outF:
        # write directly to outF
        for i, line in enumerate(inF):
            if i == 0:
                outF.write(line)
            else:
                if i % 2 == 0:
                    line = line.replace("test1", "test2")
                outF.write(line)

    # now try to run process_csv_for_sankey, should fail
    with pytest.raises(ValueError) as exc:
        process_csv_for_sankey(mult_query, csv_type="with-lineages")
    print(str(exc.value))
    assert "Multiple query_name values detected:" in str(exc.value)
    assert "test2" in str(exc.value)


def test_process_csv_for_sankey_multiple_query_names_summarized(runtmp):
    input_csv = utils.get_example_data('tax/test.tax-mg.summarized.csv')
    mult_query = runtmp.output('mult-query.csv')

    # modify input to have multiple query_name values
    with open(input_csv, newline="") as inF, open(mult_query, 'w', newline='') as outF:
        # write directly to outF
        for i, line in enumerate(inF):
            if i == 0:
                outF.write(line)
            else:
                if i % 2 == 0:
                    line = line.replace("test1", "test2")
                outF.write(line)

    # now try to run process_csv_for_sankey, should fail
    with pytest.raises(ValueError) as exc:
        process_csv_for_sankey(mult_query, csv_type="csv_summary")
    print(str(exc.value))
    assert "Multiple query_name values detected:" in str(exc.value)
    assert "test2" in str(exc.value)

def test_process_csv_for_sankey_annotate():
    input_csv = utils.get_example_data('tax/test.gather.with-lineages.csv')

    nodes, links, hover_texts = process_csv_for_sankey(input_csv, csv_type="with-lineages")
    print(nodes)
    assert nodes == ['d__Bacteria', 'p__Proteobacteria', 'c__Gammaproteobacteria', 'o__Enterobacterales', 'f__Enterobacteriaceae', 'g__Escherichia', 's__Escherichia coli', 'p__Bacteroidota', 'c__Bacteroidia', 'o__Bacteroidales', 'f__Bacteroidaceae', 'g__Prevotella', 's__Prevotella copri', 'g__Phocaeicola', 's__Phocaeicola vulgatus']
    print(links)
    assert links == [{'source': 0, 'target': 1, 'value': 5.815279361459521}, {'source': 1, 'target': 2, 'value': 5.815279361459521}, {'source': 2, 'target': 3, 'value': 5.815279361459521}, {'source': 3, 'target': 4, 'value': 5.815279361459521}, {'source': 4, 'target': 5, 'value': 5.815279361459521}, {'source': 5, 'target': 6, 'value': 5.815279361459521}, {'source': 0, 'target': 7, 'value': 5.04968235869034}, {'source': 7, 'target': 8, 'value': 5.04968235869034}, {'source': 8, 'target': 9, 'value': 5.04968235869034}, {'source': 9, 'target': 10, 'value': 5.04968235869034}, {'source': 10, 'target': 11, 'value': 5.04968235869034}, {'source': 11, 'target': 12, 'value': 5.04968235869034}, {'source': 0, 'target': 7, 'value': 1.5637726014008795}, {'source': 7, 'target': 8, 'value': 1.5637726014008795}, {'source': 8, 'target': 9, 'value': 1.5637726014008795}, {'source': 9, 'target': 10, 'value': 1.5637726014008795}, {'source': 10, 'target': 13, 'value': 1.5637726014008795}, {'source': 13, 'target': 14, 'value': 1.5637726014008795}, {'source': 0, 'target': 7, 'value': 0.6515719172503665}, {'source': 7, 'target': 8, 'value': 0.6515719172503665}, {'source': 8, 'target': 9, 'value': 0.6515719172503665}, {'source': 9, 'target': 10, 'value': 0.6515719172503665}, {'source': 10, 'target': 11, 'value': 0.6515719172503665}, {'source': 11, 'target': 12, 'value': 0.6515719172503665}]
    print(hover_texts)
    assert hover_texts == ['d__Bacteria → p__Proteobacteria<br>5.82%', 'p__Proteobacteria → c__Gammaproteobacteria<br>5.82%', 'c__Gammaproteobacteria → o__Enterobacterales<br>5.82%', 'o__Enterobacterales → f__Enterobacteriaceae<br>5.82%', 'f__Enterobacteriaceae → g__Escherichia<br>5.82%', 'g__Escherichia → s__Escherichia coli<br>5.82%', 'd__Bacteria → p__Bacteroidota<br>5.05%', 'p__Bacteroidota → c__Bacteroidia<br>5.05%', 'c__Bacteroidia → o__Bacteroidales<br>5.05%', 'o__Bacteroidales → f__Bacteroidaceae<br>5.05%', 'f__Bacteroidaceae → g__Prevotella<br>5.05%', 'g__Prevotella → s__Prevotella copri<br>5.05%', 'd__Bacteria → p__Bacteroidota<br>1.56%', 'p__Bacteroidota → c__Bacteroidia<br>1.56%', 'c__Bacteroidia → o__Bacteroidales<br>1.56%', 'o__Bacteroidales → f__Bacteroidaceae<br>1.56%', 'f__Bacteroidaceae → g__Phocaeicola<br>1.56%', 'g__Phocaeicola → s__Phocaeicola vulgatus<br>1.56%', 'd__Bacteria → p__Bacteroidota<br>0.65%', 'p__Bacteroidota → c__Bacteroidia<br>0.65%', 'c__Bacteroidia → o__Bacteroidales<br>0.65%', 'o__Bacteroidales → f__Bacteroidaceae<br>0.65%', 'f__Bacteroidaceae → g__Prevotella<br>0.65%', 'g__Prevotella → s__Prevotella copri<br>0.65%']

