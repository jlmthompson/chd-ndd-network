# CHD-NDD Network Explorer

An interactive gene-phenotype network connecting congenital heart disease (CHD) and neurodevelopmental disorders (NDD).

## Data Sources
- [PanelApp Australia](https://panelapp-aus.org) — gene panels 250 (Intellectual Disability), 51 (Autism), 76 (Congenital Heart Defect)
- ClinGen — gene-disease dosage sensitivity
- ClinVar — pathogenic/likely pathogenic variants
- HPO — phenotype ontology
- gnomAD — constraint scores

## Structure
- `data/raw/` — source data from APIs
- `data/processed/` — graph-ready JSON
- `scripts/` — data assembly pipeline
- `docs/` — GitHub Pages site

## Development
Built with [Cytoscape.js](https://js.cytoscape.org/).

