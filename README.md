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

### Purpose and Validation Strategy
The processing pipeline serves multiple critical functions:

1. **Data Standardization**: The scripts consolidate gene information from multiple species databases into a single file.

2. **ID Validation Challenge**: A key challenge is that molecular interactions often reference genes using various ID formats rather than the standardized IDs from gene descriptions. For example:
   - Gene descriptions use database-specific IDs (e.g., Xenbase)
   - Molecular interactions may reference Entrez Gene/LocusLink IDs
   - Matching IDs from gene descriptions may appear in the synonym section of the molecular interaction file

3. **Synonym Resolution**: To address this:
   - `getSynonym.py` builds a comprehensive synonym dictionary from molecular interaction data
   - This maps alternative gene identifiers to their canonical IDs from gene descriptions
   - Structure: `taxon_id → gene_id → [list of known synonyms] → source database`

4. **Validation Process**:
   - All interactions are validated against  gene IDs from gene descriptions
   - Interactions are classified as:
     - Valid: Both interacting genes have matching  IDs from gene descriptions
     - Invalid: One or both genes cannot be mapped to canonical IDs
   - Statistics are generated to track validation success rates per species

### Script Workflow
1. `CombineAllGeneDescription.py`: Establishes canonical gene IDs and descriptions
2. `getSynonym.py`: Creates ID mapping infrastructure for validation
3. **GeneInteractionProcessor.py**: Processes and validates genetic interactions
   - Extracts genetic interactions from INTERACTION-GEN_COMBINED.tsv
   - Maps gene IDs to their canonical forms using the synonym dictionary
   - Validates interaction pairs against gene descriptions
   - Outputs:
     - extracted_genetic_interactions.csv: Raw extracted interactions
     - valid_interactions.csv: Interactions where both genes are validated
     - invalid_interactions.csv: Interactions with unmappable gene IDs
     - interactions_stats.csv: Validation statistics per species
4. `validate_gene_interactions.py`: Performs final validation and generates reports

This validation ensures data integrity and provides transparency about data quality across different species databases.

## Execution Order
1. Run `CombineAllGeneDescription.py` first to create unified gene descriptions
2. Run `getSynonym.py` to build gene synonyms dictionary
3. Run `GeneInteractionProcessor.py` to process genetic interactions
4. Run `validate_gene_interactions.py` to validate and filter interactions

## Dependencies
- pandas
- pathlib
- json
- utils.species_utils (custom utility module)




