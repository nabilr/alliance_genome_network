# Create directories first
mkdir -p data/raw/{Disease,GeneDescriptions,GeneticInteractions,MolecularInteractions,Orthology}

# Create processed directory structure
mkdir -p data/processed/{GeneDescriptions,GeneticInteractions}

# Download files (using latest release URL pattern)
wget -P data/raw/Disease/ https://fms.alliancegenome.org/download/DISEASE-ALLIANCE_COMBINED.tsv

wget -P data/raw/GeneDescriptions/ https://fms.alliancegenome.org/download/GENE-DESCRIPTION-TSV_FB.tsv
wget -P data/raw/GeneDescriptions/ https://fms.alliancegenome.org/download/GENE-DESCRIPTION-TSV_HUMAN.tsv
mv data/raw/GeneDescriptions/GENE-DESCRIPTION-TSV_HUMAN.tsv data/raw/GeneDescriptions/GENE-DESCRIPTION-TSV_HGNC.tsv
wget -P data/raw/GeneDescriptions/ https://fms.alliancegenome.org/download/GENE-DESCRIPTION-TSV_MGI.tsv
wget -P data/raw/GeneDescriptions/ https://fms.alliancegenome.org/download/GENE-DESCRIPTION-TSV_SGD.tsv
wget -P data/raw/GeneDescriptions/ https://fms.alliancegenome.org/download/GENE-DESCRIPTION-TSV_WB.tsv
wget -P data/raw/GeneDescriptions/ https://fms.alliancegenome.org/download/GENE-DESCRIPTION-TSV_XBXL.tsv
wget -P data/raw/GeneDescriptions/ https://fms.alliancegenome.org/download/GENE-DESCRIPTION-TSV_XBXT.tsv
wget -P data/raw/GeneDescriptions/ https://fms.alliancegenome.org/download/GENE-DESCRIPTION-TSV_ZFIN.tsv

wget -P data/raw/GeneticInteractions/ https://fms.alliancegenome.org/download/INTERACTION-GEN_COMBINED.tsv
wget -P data/raw/MolecularInteractions/ https://fms.alliancegenome.org/download/INTERACTION-MOL_COMBINED.tsv
wget -P data/raw/Orthology/ https://fms.alliancegenome.org/download/ORTHOLOGY-ALLIANCE_COMBINED.tsva