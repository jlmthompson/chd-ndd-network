#!/usr/bin/env python3
"""
Step 1: Pull green and amber genes from PanelApp Australia panels.
Panels:
  250 - Intellectual disability syndromic and non-syndromic
   51 - Autism
   76 - Congenital Heart Defect

Outputs: data/raw/panelapp_genes.json
"""

import requests, json, time, re
from pathlib import Path
from collections import defaultdict

PANELS = {
    "250": "Intellectual Disability",
    "51":  "Autism",
    "76":  "Congenital Heart Defect"
}

INCLUDE_CONFIDENCE = {2, 3}  # 3=Green, 2=Amber

BASE = "https://panelapp-aus.org/api/v1/panels"
OUT  = Path(__file__).parent.parent / "data" / "raw" / "panelapp_genes.json"


def fetch_all_genes(panel_id):
    genes, url = [], f"{BASE}/{panel_id}/genes/?page=1&page_size=200"
    while url:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        data = r.json()
        genes.extend(data.get("results", []))
        url = data.get("next")
        if url:
            time.sleep(0.3)
    return genes


def parse_gene(g, panel_id, panel_name):
    gd   = g.get("gene_data", {})
    conf = int(g.get("confidence_level", 0))
    if conf not in INCLUDE_CONFIDENCE:
        return None

    raw_phenotypes = g.get("phenotypes", [])
    phenotypes = []
    for p in raw_phenotypes:
        if not p.strip():
            continue
        mim_match = re.search(r'MIM#?\s*(\d+)', p)
        phenotypes.append({
            "label": p.strip(),
            "mim":   mim_match.group(1) if mim_match else None
        })

    return {
        "gene_symbol":          gd.get("gene_symbol"),
        "gene_name":            gd.get("gene_name"),
        "hgnc_id":              gd.get("hgnc_id"),
        "omim_gene":            gd.get("omim_gene", []),
        "confidence":           conf,
        "confidence_label":     "Green" if conf == 3 else "Amber",
        "mode_of_inheritance":  g.get("mode_of_inheritance"),
        "mode_of_pathogenicity":g.get("mode_of_pathogenicity"),
        "penetrance":           g.get("penetrance"),
        "phenotypes":           phenotypes,
        "publications":         g.get("publications", []),
        "panel_id":             panel_id,
        "panel_name":           panel_name,
    }


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    all_genes = []

    for panel_id, panel_name in PANELS.items():
        print(f"Fetching panel {panel_id}: {panel_name}...")
        raw    = fetch_all_genes(panel_id)
        parsed = [parse_gene(g, panel_id, panel_name) for g in raw]
        parsed = [g for g in parsed if g and g["gene_symbol"]]
        print(f"  {len(parsed)} green/amber genes")
        all_genes.extend(parsed)

    unique_symbols = {g["gene_symbol"] for g in all_genes}
    print(f"\nTotal entries:  {len(all_genes)}")
    print(f"Unique genes:   {len(unique_symbols)}")

    gene_panels = defaultdict(set)
    for g in all_genes:
        gene_panels[g["gene_symbol"]].add(g["panel_name"])
    overlap = {sym: sorted(panels) for sym, panels in gene_panels.items() if len(panels) > 1}
    print(f"Cross-panel overlaps: {len(overlap)}")

    print("\nTop overlap genes (CHD + NDD/Autism):")
    for sym, panels in list(overlap.items())[:20]:
        print(f"  {sym}: {', '.join(panels)}")

    OUT.write_text(json.dumps({
        "panels":        PANELS,
        "total_entries": len(all_genes),
        "unique_genes":  len(unique_symbols),
        "overlap_count": len(overlap),
        "overlap_genes": overlap,
        "genes":         all_genes
    }, indent=2))
    print(f"\nSaved → {OUT}")


if __name__ == "__main__":
    main()
