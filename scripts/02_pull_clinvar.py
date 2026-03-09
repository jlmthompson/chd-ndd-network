#!/usr/bin/env python3
"""
Step 2: Pull ClinVar pathogenic/likely-pathogenic variants for all panel genes.
Uses NCBI eutils — no API key required (rate-limited to 3 req/s).

Outputs: data/raw/clinvar_variants.json
"""

import requests, json, time
from pathlib import Path

PANELAPP = Path(__file__).parent.parent / "data" / "raw" / "panelapp_genes.json"
OUT       = Path(__file__).parent.parent / "data" / "raw" / "clinvar_variants.json"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CLINSIG = '("pathogenic"[Clinical significance] OR "likely pathogenic"[Clinical significance])'


def search_clinvar(gene_symbol):
    query = f'{gene_symbol}[Gene Name] AND {CLINSIG} AND "single nucleotide variant"[Type]'
    r = requests.get(f"{EUTILS}/esearch.fcgi", params={
        "db": "clinvar", "term": query,
        "retmax": 50, "retmode": "json", "usehistory": "y"
    }, timeout=20)
    r.raise_for_status()
    data = r.json().get("esearchresult", {})
    ids  = data.get("idlist", [])
    wh   = data.get("webenv")
    qk   = data.get("querykey")
    return ids, wh, qk


def fetch_summaries(ids, webenv, querykey):
    if not ids:
        return []
    r = requests.get(f"{EUTILS}/esummary.fcgi", params={
        "db": "clinvar", "webenv": webenv, "query_key": querykey,
        "retmax": 50, "retmode": "json"
    }, timeout=30)
    r.raise_for_status()
    result = r.json().get("result", {})
    uids   = result.get("uids", [])
    variants = []
    for uid in uids:
        v = result.get(uid, {})
        # Extract germline classification
        germ = v.get("germline_classification", {})
        clnsig = germ.get("description", v.get("clinical_significance", {}).get("description", ""))
        if not any(s in clnsig.lower() for s in ["pathogenic", "likely pathogenic"]):
            continue
        # Preferred name / HGVS
        title = v.get("title", "")
        # Gene info
        genes_info = v.get("genes", [{}])
        gene_sym   = genes_info[0].get("symbol", "") if genes_info else ""
        variants.append({
            "clinvar_id":    uid,
            "title":         title,
            "gene_symbol":   gene_sym,
            "clinical_significance": clnsig,
            "variation_type": v.get("obj_type", ""),
            "accession":     v.get("accession", ""),
            "last_evaluated": germ.get("last_evaluated", ""),
            "conditions":    [c.get("name","") for c in v.get("trait_set", [])],
            "review_status": germ.get("review_status", ""),
        })
    return variants


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    pa = json.loads(PANELAPP.read_text())

    # Get unique gene symbols
    gene_symbols = sorted({g["gene_symbol"] for g in pa["genes"] if g.get("gene_symbol")})
    print(f"Pulling ClinVar variants for {len(gene_symbols)} genes...")

    all_variants = []
    errors = []
    for i, sym in enumerate(gene_symbols):
        if i % 50 == 0:
            print(f"  [{i}/{len(gene_symbols)}] {sym}")
        try:
            ids, wh, qk = search_clinvar(sym)
            if ids:
                variants = fetch_summaries(ids, wh, qk)
                all_variants.extend(variants)
            time.sleep(0.35)  # Stay under 3 req/s
        except Exception as e:
            errors.append({"gene": sym, "error": str(e)})
            time.sleep(1)

    print(f"\nTotal P/LP variants retrieved: {len(all_variants)}")
    print(f"Genes with variants: {len({v['gene_symbol'] for v in all_variants})}")
    if errors:
        print(f"Errors: {len(errors)}")

    OUT.write_text(json.dumps({
        "total_variants": len(all_variants),
        "variants": all_variants,
        "errors": errors
    }, indent=2))
    print(f"Saved → {OUT}")


if __name__ == "__main__":
    main()
