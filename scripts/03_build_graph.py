#!/usr/bin/env python3
"""
Step 3: Assemble the Cytoscape.js-ready graph JSON from all data sources.

Node types:  gene | chd_phenotype | ndd_phenotype | variant
Edge types:  gene_in_panel | gene_causes_phenotype | gene_has_variant
             gene_in_multiple_panels (overlap edges)

Outputs: data/processed/graph.json
"""

import json, re
from pathlib import Path
from collections import defaultdict

PANELAPP = Path(__file__).parent.parent / "data" / "raw" / "panelapp_genes.json"
CLINVAR  = Path(__file__).parent.parent / "data" / "raw" / "clinvar_variants.json"
OUT      = Path(__file__).parent.parent / "data" / "processed" / "graph.json"

# Panel display colours
PANEL_COLOURS = {
    "Congenital Heart Defect": "#e74c3c",   # red
    "Intellectual Disability":  "#3498db",  # blue
    "Autism":                   "#2ecc71",  # green
}

CONFIDENCE_COLOURS = {
    "Green": "#27ae60",
    "Amber": "#f39c12",
}


def slugify(s):
    return re.sub(r'[^a-z0-9_]', '_', s.lower().strip())[:80]


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    pa = json.loads(PANELAPP.read_text())
    cv = json.loads(CLINVAR.read_text()) if CLINVAR.exists() else {"variants": []}

    nodes, edges = [], []
    node_ids = set()

    def add_node(nid, data):
        if nid not in node_ids:
            nodes.append({"data": {"id": nid, **data}})
            node_ids.add(nid)

    def add_edge(eid, source, target, data):
        edges.append({"data": {"id": eid, "source": source, "target": target, **data}})

    # ── Gene nodes ───────────────────────────────────────────────────────────
    # Track which panels each gene appears in
    gene_panels  = defaultdict(set)
    gene_entries = defaultdict(list)
    for g in pa["genes"]:
        sym = g["gene_symbol"]
        gene_panels[sym].add(g["panel_name"])
        gene_entries[sym].append(g)

    for sym, entries in gene_entries.items():
        panels = sorted(gene_panels[sym])
        is_overlap = len(panels) > 1
        has_chd    = "Congenital Heart Defect" in panels
        has_ndd    = any(p in panels for p in ["Intellectual Disability","Autism"])

        # Colour: overlap genes get special colour
        if has_chd and has_ndd:
            colour = "#9b59b6"       # purple — CHD+NDD overlap
        elif has_chd:
            colour = "#e74c3c"       # red — CHD only
        else:
            colour = "#3498db"       # blue — NDD only

        # Best confidence across entries
        best_conf = max(e["confidence"] for e in entries)
        conf_label = "Green" if best_conf == 3 else "Amber"

        # Collect all phenotypes across entries
        all_phenotypes = []
        for e in entries:
            all_phenotypes.extend(e.get("phenotypes", []))

        # MOI from CHD entry if available, else first
        moi = next((e["mode_of_inheritance"] for e in entries
                    if e["panel_name"] == "Congenital Heart Defect"), entries[0].get("mode_of_inheritance",""))

        add_node(sym, {
            "label":          sym,
            "type":           "gene",
            "full_name":      entries[0].get("gene_name",""),
            "hgnc_id":        entries[0].get("hgnc_id",""),
            "omim_gene":      entries[0].get("omim_gene",[]),
            "panels":         panels,
            "is_overlap":     is_overlap,
            "has_chd":        has_chd,
            "has_ndd":        has_ndd,
            "confidence":     conf_label,
            "moi":            moi,
            "colour":         colour,
            "phenotypes":     all_phenotypes[:10],  # cap for perf
        })

    # ── Phenotype nodes + gene→phenotype edges ────────────────────────────────
    for g in pa["genes"]:
        sym   = g["gene_symbol"]
        ptype = "chd_phenotype" if g["panel_name"] == "Congenital Heart Defect" else "ndd_phenotype"
        pcol  = "#e8a49c" if ptype == "chd_phenotype" else "#a8c8e8"

        for ph in g.get("phenotypes", []):
            label = ph.get("label","").strip()
            if not label:
                continue
            pid = "ph_" + slugify(label)
            add_node(pid, {
                "label":  label[:60] + ("…" if len(label) > 60 else ""),
                "label_full": label,
                "type":   ptype,
                "mim":    ph.get("mim"),
                "colour": pcol,
            })
            eid = f"e_{sym}_{pid}"
            add_edge(eid, sym, pid, {
                "type":   "gene_causes_phenotype",
                "panel":  g["panel_name"],
                "colour": "#999",
            })

    # ── Variant nodes + gene→variant edges ───────────────────────────────────
    variant_counts = defaultdict(int)
    for v in cv.get("variants", []):
        sym = v.get("gene_symbol","")
        if not sym or sym not in node_ids:
            continue
        variant_counts[sym] += 1
        # Only add individual variant nodes for genes with ≤10 variants
        # (avoids overwhelming the graph for common genes)
        if variant_counts[sym] > 10:
            continue
        vid = f"var_{v['clinvar_id']}"
        add_node(vid, {
            "label":    v.get("title","")[:50],
            "type":     "variant",
            "clinsig":  v.get("clinical_significance",""),
            "accession":v.get("accession",""),
            "conditions": v.get("conditions",[]),
            "colour":   "#f1c40f",
        })
        add_edge(f"e_{sym}_{vid}", sym, vid, {
            "type":   "gene_has_variant",
            "colour": "#ccc",
        })

    # Annotate gene nodes with total variant count
    for node in nodes:
        if node["data"]["type"] == "gene":
            sym = node["data"]["id"]
            node["data"]["variant_count"] = variant_counts.get(sym, 0)

    # ── Stats ────────────────────────────────────────────────────────────────
    n_genes    = sum(1 for n in nodes if n["data"]["type"] == "gene")
    n_pheno    = sum(1 for n in nodes if "phenotype" in n["data"]["type"])
    n_variants = sum(1 for n in nodes if n["data"]["type"] == "variant")
    n_overlap  = sum(1 for n in nodes if n["data"].get("is_overlap"))

    print(f"Nodes:  {len(nodes):,}  ({n_genes} genes, {n_pheno} phenotypes, {n_variants} variants)")
    print(f"Edges:  {len(edges):,}")
    print(f"Overlap genes (multi-panel): {n_overlap}")

    graph = {
        "stats": {
            "total_nodes": len(nodes), "total_edges": len(edges),
            "genes": n_genes, "phenotypes": n_pheno,
            "variants": n_variants, "overlap_genes": n_overlap,
        },
        "elements": {"nodes": nodes, "edges": edges}
    }
    OUT.write_text(json.dumps(graph, indent=2))
    print(f"Saved → {OUT}")


if __name__ == "__main__":
    main()
