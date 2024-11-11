import pandas as pd
from pathlib import Path
from utils.species_utils import load_species_map
import csv

def validate_gene_interactions():
    """
    Validate that all genes referenced in interactions exist in gene descriptions.
    Returns a DataFrame with only valid interactions where both genes exist.
    """
    # Load data files and species map
    species_map = load_species_map()
    # Fix: species_map now returns Dict[str, Dict[str, str]], so we need to get names differently
    taxon_names = {taxon_id: info['name'] for taxon_id, info in species_map.items()}
    
    # Read gene_nodes with double quotes (since they're quoted in the file)
    gene_nodes = pd.read_csv('data/processed/GeneDescriptions/gene_nodes.csv', 
                            low_memory=False,
                            quoting=csv.QUOTE_ALL)
    
    # Read interactions with no special quoting (since they're plain CSV)
    interactions = pd.read_csv('data/processed/GeneticInteractions/extracted_genetic_interactions.csv', 
                             low_memory=False,
                             names=['database', 'taxonId', 'fromGeneId', 'toGeneId'])  # Specify column names
    
    
    
    
   
    
    # Add species names to interactions and drop any rows with NA
    interactions['species_name'] = interactions['taxonId'].map(taxon_names)
    interactions = interactions.dropna()
    print("\nFirst few rows of interactions:")
    print(interactions.head())
    
    # Create composite key for genes (database + geneId + taxonId)
    # First, print the columns to debug
    print("\nColumns in gene_nodes:")
    print(gene_nodes.columns)
    
    # Debug: Print column names
    print("\nColumns in gene_nodes DataFrame:")
    print(gene_nodes.columns.tolist())
    
    # Map species names to taxon IDs using db_name instead of name
    species_to_taxon = {}
    for taxon_id, info in species_map.items():
        db_name = info['db_name']
        species_name = info['short_name'].lower()
        species_to_taxon[species_name] = taxon_id
        species_to_taxon[db_name] = taxon_id  # Also map database names to taxon IDs

    print(species_to_taxon)
    input("Press Enter to continue...")
    # Map taxon IDs for gene_nodes, trying both Species and database fields
    gene_nodes['taxonId'] = gene_nodes['Species'].str.lower().map(species_to_taxon)
    # If taxonId is NA, try using the database field
    mask = gene_nodes['taxonId'].isna()
    gene_nodes.loc[mask, 'taxonId'] = gene_nodes.loc[mask, 'database'].map(species_to_taxon)

    # Map database names using species map
    def get_db_name(row):
       
        species_info = species_map.get(str(row['taxonId']), {})
        # print(species_info,row)
        return species_info.get('db_name', row['database'].lower()).lower()

    # Create composite key using mapped database names
    gene_nodes['gene_key'] = (gene_nodes.apply(get_db_name, axis=1) + ':' + 
                             gene_nodes['geneId'] + ':' + 
                             gene_nodes['taxonId'].astype(str))
    valid_genes = set(gene_nodes['gene_key'])

    print("\nFirst few rows of gene_nodes:")
    print(gene_nodes[gene_nodes['Species']=='HGNC'].head())

    input("Press Enter to continue...")
    
    # Create composite keys for interaction genes
    interactions['from_key'] = (interactions['database'].str.lower() + ':' + 
                              interactions['fromGeneId'] + ':' + 
                              interactions['taxonId'].astype(str))
    interactions['to_key'] = (interactions['database'].str.lower() + ':' + 
                            interactions['toGeneId'] + ':' + 
                            interactions['taxonId'].astype(str))
    
        # Print head of interactions DataFrame
    print("\nFirst few rows of interactions:")
    print(interactions.head())
    
    # Print head of gene_nodes DataFrame 
    print("\nFirst few rows of gene_nodes:")
    print(gene_nodes.head())

    # print(valid_genes)
    
    # Check which genes are missing
    missing_from = ~interactions['from_key'].isin(valid_genes)
    missing_to = ~interactions['to_key'].isin(valid_genes)
    
    # Filter to valid/invalid interactions
    valid_interactions = interactions[~(missing_from | missing_to)].copy()
    invalid_interactions = interactions[missing_from | missing_to].copy()

    print(f"Valid interactions: {len(valid_interactions)}")
    print(f"Invalid interactions: {len(invalid_interactions)}")
    
    # Create statistics DataFrame
    stats_rows = []
    
    # Calculate overall statistics
    total = len(interactions)
    
    # Add overall totals
    stats_rows.append({
        'level': 'overall',
        'taxon_id': 'all',
        'species_name': 'all',
        'database': 'all',
        'interaction_type': 'all',
        'count': total
    })
    stats_rows.append({
        'level': 'overall',
        'taxon_id': 'all',
        'species_name': 'all',
        'database': 'all',
        'interaction_type': 'valid',
        'count': len(valid_interactions)
    })
    stats_rows.append({
        'level': 'overall',
        'taxon_id': 'all',
        'species_name': 'all',
        'database': 'all',
        'interaction_type': 'invalid',
        'count': len(invalid_interactions)
    })
    
    # Calculate statistics by taxon_id and database
    for taxon_id in sorted(interactions['taxonId'].unique()):
        species_name = taxon_names.get(taxon_id, 'Unknown')
        taxon_mask = interactions['taxonId'] == taxon_id
        taxon_data = interactions[taxon_mask]
        
        # Fix: Use boolean indexing with loc for proper alignment
        taxon_valid = taxon_data.loc[~(missing_from[taxon_mask] | missing_to[taxon_mask])]
        taxon_invalid = taxon_data.loc[missing_from[taxon_mask] | missing_to[taxon_mask]]
        
        stats_rows.append({
            'level': 'taxon',
            'taxon_id': taxon_id,
            'species_name': species_name,
            'database': 'all',
            'interaction_type': 'all',
            'count': len(taxon_data)
        })
        stats_rows.append({
            'level': 'taxon',
            'taxon_id': taxon_id,
            'species_name': species_name,
            'database': 'all',
            'interaction_type': 'valid',
            'count': len(taxon_valid)  # Fixed: Count valid interactions
        })
        stats_rows.append({
            'level': 'taxon',
            'taxon_id': taxon_id,
            'species_name': species_name,
            'database': 'all',
            'interaction_type': 'invalid',
            'count': len(taxon_invalid)  # Fixed: Count invalid interactions
        })
        
        # Calculate database-level statistics for this taxon
        for database in sorted(taxon_data['database'].unique()):
            db_mask = taxon_data['database'] == database
            db_data = taxon_data[db_mask]
            
            # Fix: Use boolean indexing with loc for proper alignment
            db_mask_indices = db_data.index
            valid_db = db_data.loc[~(missing_from[db_mask_indices] | missing_to[db_mask_indices])]
            invalid_db = db_data.loc[missing_from[db_mask_indices] | missing_to[db_mask_indices]]

            stats_rows.extend([
                {
                    'level': 'database',
                    'taxon_id': taxon_id,
                    'species_name': species_name,
                    'database': database,
                    'interaction_type': 'all',
                    'count': len(db_data)
                },
                {
                    'level': 'database',
                    'taxon_id': taxon_id,
                    'species_name': species_name,
                    'database': database,
                    'interaction_type': 'valid',
                    'count': len(valid_db)
                },
                {
                    'level': 'database',
                    'taxon_id': taxon_id,
                    'species_name': species_name,
                    'database': database,
                    'interaction_type': 'invalid',
                    'count': len(invalid_db)
                }
            ])
    
    # Create and save statistics DataFrame
    stats_df = pd.DataFrame(stats_rows)
    stats_df = stats_df[['level', 'taxon_id', 'species_name', 'database', 'interaction_type', 'count']]
    stats_df = stats_df.sort_values(['level', 'taxon_id', 'database', 'interaction_type'])
    
    # Save files with same format as input
    base_path = 'data/processed/GeneticInteractions'
    valid_output = f'{base_path}/valid_interactions.csv'
    invalid_output = f'{base_path}/invalid_interactions.csv'
    stats_output = f'{base_path}/interactions_stats.csv'
    
    # Save without any special quoting to match input format
    valid_interactions.to_csv(valid_output, index=False)
    invalid_interactions.to_csv(invalid_output, index=False)
    stats_df.to_csv(stats_output, index=False)
    
    print(f"\nFiles saved:")
    print(f"Valid interactions: {valid_output}")
    print(f"Invalid interactions: {invalid_output}")
    print(f"Statistics: {stats_output}")
    
    # Sample of invalid interactions for debugging
    if len(invalid_interactions) > 0:
        print("\nSample of invalid interactions:")
        sample = invalid_interactions.head()
        for _, row in sample.iterrows():
            from_exists = row['from_key'] in valid_genes
            to_exists = row['to_key'] in valid_genes
            print(f"\nInteraction in {row['species_name']} (taxon {row['taxonId']}, {row['database']}):")
            print(f"  {row['from_key']} -> {row['to_key']}")
            print(f"  From gene exists: {from_exists}")
            print(f"  To gene exists: {to_exists}")
    
    return valid_interactions

if __name__ == "__main__":
    validate_gene_interactions()