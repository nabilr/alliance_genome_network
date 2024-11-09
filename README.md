## Data Organization

### Data Folder Structure

### Description
- **raw/**: Original TSV format files from [Alliance Genome Database](https://www.alliancegenome.org/downloads), organized into categories:
  - **Disease/**: Disease annotation data
    - DISEASE-ALLIANCE_COMBINED.tsv
  - **GeneDescriptions/**: Gene descriptions across multiple species
    - GENE-DESCRIPTION-TSV_FB.tsv (Fly Base)
    - GENE-DESCRIPTION-TSV_HUMAN.tsv (Human)
    - GENE-DESCRIPTION-TSV_MGI.tsv (Mouse)
    - GENE-DESCRIPTION-TSV_SGD.tsv (Yeast)
    - GENE-DESCRIPTION-TSV_WB.tsv (Worm Base)
    - GENE-DESCRIPTION-TSV_XBXL.tsv (Xenopus)
    - GENE-DESCRIPTION-TSV_XBXT.tsv (Xenopus)
    - GENE-DESCRIPTION-TSV_ZFIN.tsv (Zebrafish)
  - **GeneticInteractions/**: Gene-to-gene interaction data
    - INTERACTION-GEN_COMBINED.tsv
  - **MolecularInteractions/**: Molecular interaction data
    - INTERACTION-MOL_COMBINED.tsv
  - **Orthology/**: Cross-species gene orthology data
    - ORTHOLOGY-ALLIANCE_COMBINED.tsv
- **processed/**: Directory for cleaned and analyzed datasets derived from the raw files

