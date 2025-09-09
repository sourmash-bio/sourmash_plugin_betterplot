"""
Tests for sourmash_plugin_betterplot.
"""
import os
import pytest
import math

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


def _links_by_label(nodes, links):
    # turn index-based links into label-keyed dict
    pairs = {}
    for L in links:
        s = nodes[L["source"]]
        t = nodes[L["target"]]
        pairs[(s, t)] = pairs.get((s, t), 0.0) + float(L["value"])
    return pairs

def _approx_equal_maps(a, b, tol=1e-9):
    if set(a.keys()) != set(b.keys()):
        return False
    return all(math.isclose(a[k], b[k], rel_tol=1e-9, abs_tol=tol) for k in a.keys())

def test_process_csv_annotate_vs_summary():
    ann_csv = utils.get_example_data('tax/test.gather.with-lineages.csv')
    sum_csv = utils.get_example_data('tax/test.tax-mg.summarized.csv')

    nodes_a, links_a, _ = process_csv_for_sankey(ann_csv, csv_type="with-lineages")
    nodes_b, links_b, _ = process_csv_for_sankey(sum_csv, csv_type="csv_summary")

    expect = {
        ("g__Escherichia", "s__Escherichia coli"):      5.815279361459521,
        ("g__Prevotella",  "s__Prevotella copri"):      5.701254275940707,  # 5.04968235869034 + 0.6515719172503665
        ("g__Phocaeicola", "s__Phocaeicola vulgatus"):  1.5637726014008795,
    }

    pairs_a = _links_by_label(nodes_a, links_a)
    print("annotate:", pairs_a)
    pairs_b = _links_by_label(nodes_b, links_b)
    print("summary_csv:", pairs_b)

    # check for expected values
    for k in expect.keys():
        assert k in pairs_a
        assert k in pairs_b
        assert math.isclose(pairs_a[k], expect[k], rel_tol=1e-9, abs_tol=1e-9)
        assert math.isclose(pairs_b[k], expect[k], rel_tol=1e-9, abs_tol=1e-9)

    # check that annotate and summary give same results
    assert _approx_equal_maps(pairs_a, pairs_b)


def test_lins_annotate_vs_summary_equivalence_no_lingroup():
    # Adjust paths to how you vend these fixtures in your test utils
    ann_csv = utils.get_test_data('SRR29654720_k31_gather.with-lineages.csv')
    sum_csv = utils.get_test_data('SRR29654720_k31_gather.summarized.no-lingroups.csv')

    # Both files should hit the LINS codepath (detect_lins=True)
    nodes_a, links_a, _ = process_csv_for_sankey(ann_csv, csv_type="with-lineages")
    nodes_b, links_b, _ = process_csv_for_sankey(sum_csv, csv_type="csv_summary")

    # Node labels should be raw integer paths (prefixes stripped), not pN:*
    assert not any(("p0:" in n or "p1:" in n) for n in nodes_a + nodes_b)
    assert all(all(tok.isdigit() for tok in n.split(";")) for n in nodes_a + nodes_b)

    # Pick a known leaf from the sample (first row of the annotate file)
    leaf = "864;0;0;1;0;1;0;0;0;0;2;1;0;2;0;0;0;0;0;0"
    # It should appear as a node in both results
    assert leaf in nodes_a
    assert leaf in nodes_b

    # check leaf value
    leaf_idx_a = nodes_a.index(leaf)
    leaf_idx_b = nodes_b.index(leaf)
    leaf_value_a = sum(L["value"] for L in links_a if L["target"] == leaf_idx_a)
    leaf_value_b = sum(L["value"] for L in links_b if L["target"] == leaf_idx_b)
    print(f"leaf {leaf} values: annotate {leaf_value_a}, summary {leaf_value_b}")
    assert math.isclose(leaf_value_a, leaf_value_b, rel_tol=1e-9, abs_tol=1e-9)

    # Compare the flows keyed by (src_label, tgt_label)
    pairs_a = _links_by_label(nodes_a, links_a)
    pairs_b = _links_by_label(nodes_b, links_b)

    # They should be the same (within tiny tolerance; values are percentages)
    assert _approx_equal_maps(pairs_a, pairs_b, tol=1e-9)

def _is_unlabeled(label: str) -> bool:
    # LINGROUP labels look like "name (lin)"; raw LINS have no '('
    return "(" not in label

def _assert_summary_subset_allow_unlabeled_bridge(summary_pairs, annotate_pairs, tol=1e-9):
    # Build adjacency from annotate
    children = defaultdict(set)
    for (s, t) in annotate_pairs.keys():
        children[s].add(t)

    missing = []
    mismatched = []

    for (src, tgt), v in summary_pairs.items():
        # direct match?
        if (src, tgt) in annotate_pairs:
            if not math.isclose(annotate_pairs[(src, tgt)], v, rel_tol=1e-9, abs_tol=tol):
                mismatched.append(((src, tgt), annotate_pairs[(src, tgt)], v))
            continue

        # try a one-hop unlabeled bridge: src -> U (unlabeled) -> tgt
        bridged_ok = False
        for mid in children.get(src, ()):
            if _is_unlabeled(mid) and (mid, tgt) in annotate_pairs:
                if math.isclose(annotate_pairs[(mid, tgt)], v, rel_tol=1e-9, abs_tol=tol):
                    bridged_ok = True
                    break

        if not bridged_ok:
            missing.append((src, tgt))

    msg = []
    if missing:
        msg.append(f"Missing {len(missing)} links (even with unlabeled bridge): {missing[:5]}...")
    if mismatched:
        msg.append("Value mismatches (annotate vs summary) e.g. " +
                   "; ".join([f"{k}: {a} vs {b}" for k, a, b in mismatched[:5]]) + " ...")
    assert not missing and not mismatched, "\n".join(msg)


def test_lins_summary_with_lingroup_is_subset_of_annotate_with_lingroup():
    ann_csv = utils.get_test_data('SRR29654720_k31_gather.with-lineages.csv')
    lingroups_file = utils.get_test_data('ralstonia_lingroups.csv')
    sum_csv = utils.get_test_data('SRR29654720_k31_gather.summarized-lingroups.csv')

    # Both files should hit the LINS codepath (detect_lins=True)
    nodes_a, links_a, _ = process_csv_for_sankey(ann_csv, csv_type="with-lineages", lingroup_map=lingroups_file)
    nodes_b, links_b, _ = process_csv_for_sankey(sum_csv, csv_type="csv_summary")
    print("nodes_a:", nodes_a)
    print("nodes_b:", nodes_b)
    # Node labels should be raw integer paths (prefixes stripped), not pN:*
    assert not any(("p0:" in n or "p1:" in n) for n in nodes_a + nodes_b)
    # Pick a known named lingroup from the sample
    leaf_lin = "864;0;0;1;0;2;0;0;0;1;0"
    leaf_display = f"phylotype iv ({leaf_lin})"
    # It should appear as a node in both results
    assert leaf_display in nodes_a
    assert leaf_display in nodes_b
    # check leaf value
    leaf_idx_a = nodes_a.index(leaf_display)
    leaf_idx_b = nodes_b.index(leaf_display)
    leaf_value_a = sum(L["value"] for L in links_a if L["target"] == leaf_idx_a)
    leaf_value_b = sum(L["value"] for L in links_b if L["target"] == leaf_idx_b)
    print(f"leaf {leaf_display} values: annotate {leaf_value_a}, summary {leaf_value_b}")
    assert math.isclose(leaf_value_a, leaf_value_b, rel_tol=1e-9, abs_tol=1e-9)
    # check higher-level node too
    higher_lin = "864;0;0;1;0;0"
    higher_display = f"Ralstonia solanacearum ({higher_lin})"
    assert higher_display in nodes_a
    assert higher_display in nodes_b
    higher_idx_a = nodes_a.index(higher_display)
    higher_idx_b = nodes_b.index(higher_display)
    higher_value_a = sum(L["value"] for L in links_a if L["target"] == higher_idx_a)
    higher_value_b = sum(L["value"] for L in links_b if L["target"] == higher_idx_b)
    print(f"higher {higher_display} values: annotate {higher_value_a}, summary {higher_value_b}")
    assert math.isclose(higher_value_a, higher_value_b, rel_tol=1e-9, abs_tol=1e-9)

    # # Compare the flows keyed by (src_label, tgt_label)
    pairs_a = _links_by_label(nodes_a, links_a)
    pairs_b = _links_by_label(nodes_b, links_b)
    print("annotate pairs:", pairs_a)
    print("summary pairs:", pairs_b)
    # summary should be a subset of annotate (some paths may be missing. Allow unlabeled bridges)
    _assert_summary_subset_allow_unlabeled_bridge(pairs_b, pairs_a, tol=1e-9)
