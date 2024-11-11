## Data Download

### download.sh
**Purpose**: Downloads all required data files from the Alliance Genome Database.

**Usage**:
```bash
bash data/download.sh
```

**Actions**:
1. Creates necessary directory structure under `data/raw/`
2. Downloads TSV files from Alliance Genome Database
3. Renames HUMAN gene description file to HGNC format

**Dependencies**:
- wget
- bash shell

## Data Organization

### Directory Structure
```
data/
├── raw/                                  # Original files from Alliance Genome Database
│   ├── Disease/
│   │   └── DISEASE-ALLIANCE_COMBINED.tsv
│   ├── GeneDescriptions/
│   │   ├── GENE-DESCRIPTION-TSV_FB.tsv    # Fly Base
│   │   ├── GENE-DESCRIPTION-TSV_HGNC.tsv  # Human
│   │   ├── GENE-DESCRIPTION-TSV_MGI.tsv   # Mouse
│   │   ├── GENE-DESCRIPTION-TSV_SGD.tsv   # Yeast
│   │   ├── GENE-DESCRIPTION-TSV_WB.tsv    # Worm Base
│   │   ├── GENE-DESCRIPTION-TSV_XBXL.tsv  # Xenopus
│   │   ├── GENE-DESCRIPTION-TSV_XBXT.tsv  # Xenopus
│   │   └── GENE-DESCRIPTION-TSV_ZFIN.tsv  # Zebrafish
│   ├── GeneticInteractions/
│   │   └── INTERACTION-GEN_COMBINED.tsv
│   ├── MolecularInteractions/
│   │   └── INTERACTION-MOL_COMBINED.tsv
│   └── Orthology/
│       └── ORTHOLOGY-ALLIANCE_COMBINED.tsv
│
└── processed/                            # Generated files from processing pipeline
    ├── GeneDescriptions/
    │   ├── gene_nodes.csv               # Combined gene descriptions
    │   └── species_metadata.json        # Database and species metadata
    ├── GeneticInteractions/
    │   ├── extracted_genetic_interactions.csv
    │   ├── valid_interactions.csv
    │   ├── invalid_interactions.csv
    │   ├── interactions_stats.csv
    │   └── species_metadata.json
    └── gene_synonyms.json               # Dictionary of gene synonyms
```

## Processing Pipeline

### 1. CombineAllGeneDescription.py
**Purpose**: Combines gene descriptions from multiple species into a single standardized format.

**Input**:
- `data/raw/GeneDescriptions/*.tsv` files containing gene descriptions for each species
  - Format: Tab-separated files with columns: geneId, Symbol, Description

**Output**:
- `data/processed/GeneDescriptions/gene_nodes.csv`
  - Columns: database, geneId, Symbol, Description, Species
  - All text fields are quoted
- `data/processed/GeneDescriptions/species_metadata.json`
  - Contains list of unique databases and species

### 2. getSynym.py
**Purpose**: Processes molecular interaction data to build a dictionary of gene synonyms across species.

**Input**:
- `data/raw/MolecularInteractions/INTERACTION-MOL_COMBINED.tsv`
  - Tab-separated file containing molecular interactions
- Gene descriptions from step 1

**Output**:
- `data/processed/gene_synonyms.json`
  - Hierarchical dictionary organized by:
    - taxon_id → gene_id → synonym → database_name

### 3. GeneInteractionProcessor.py
**Purpose**: Processes genetic interaction data, mapping gene IDs using synonyms dictionary.

**Input**:
- `data/raw/GeneticInteractions/INTERACTION-GEN_COMBINED.tsv`
  - Tab-separated file containing genetic interactions
- `data/processed/gene_synonyms.json` from step 2

**Output**:
- `data/processed/GeneticInteractions/extracted_genetic_interactions.csv`
  - Columns: database, taxonId, fromGeneId, toGeneId
- `data/processed/GeneticInteractions/species_metadata.json`
  - Contains validation summary and examples per species

### 4. validate_gene_interactions.py
**Purpose**: Validates processed genetic interactions against gene descriptions.

**Input**:
- `data/processed/GeneticInteractions/extracted_genetic_interactions.csv` from step 3
- `data/processed/GeneDescriptions/gene_nodes.csv` from step 1

**Output**:
- `data/processed/GeneticInteractions/valid_interactions.csv`
  - Contains only interactions where both genes exist in descriptions
- `data/processed/GeneticInteractions/invalid_interactions.csv`
  - Contains interactions where at least one gene is missing
- `data/processed/GeneticInteractions/interactions_stats.csv`
  - Statistical summary of valid/invalid interactions by species and database

## Execution Order
1. Run `CombineAllGeneDescription.py` first to create unified gene descriptions
2. Run `getSynym.py` to build gene synonyms dictionary
3. Run `GeneInteractionProcessor.py` to process genetic interactions
4. Run `validate_gene_interactions.py` to validate and filter interactions

## Dependencies
- pandas
- pathlib
- json
- utils.species_utils (custom utility module)




