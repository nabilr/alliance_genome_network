import pandas as pd
import json
from pathlib import Path
from utils.species_utils import (
    load_species_map, 
    get_species_name, 
    get_species_db_name,
    get_species_shortname,
    get_valid_databases,
    is_valid_database
)
from typing import Dict

# Define constants
RAW_DIR = Path('data/raw/GeneticInteractions')
PROCESSED_DIR = Path('data/processed/GeneticInteractions')
INPUT_FILE = RAW_DIR / 'INTERACTION-GEN_COMBINED.tsv'
OUTPUT_FILE = PROCESSED_DIR / 'extracted_genetic_interactions.csv'
METADATA_FILE = PROCESSED_DIR / 'species_metadata.json'
SYNONYMS_FILE = Path('data/processed/gene_synonyms.json')

INTERACTOR_COLS = [
    'ID(s) interactor A', 'ID(s) interactor B',
    'Taxid interactor A', 'Taxid interactor B'
]

# Load species map from config file
SPECIES_MAP = load_species_map()

class DataValidationError(Exception):
    """Custom exception for data validation errors."""
    pass

def validate_taxon_id(taxon_id: str, species_map: Dict[str, Dict[str, str]]) -> bool:
    """
    Validate if taxon ID exists in species map.
    
    Args:
        taxon_id: NCBI taxonomy ID
        species_map: Dictionary of species information
    
    Returns:
        bool: True if taxon ID is valid
    """
    if taxon_id not in species_map:
        raise DataValidationError(f"Taxon ID {taxon_id} not found in species map")
    return True

def validate_species_data(species_map: Dict[str, Dict[str, str]], taxon_id: str) -> None:
    """
    Validate species data exists and is complete.
    
    Args:
        species_map: Dictionary of species information
        taxon_id: NCBI taxonomy ID
        
    Raises:
        DataValidationError: If species data is invalid or incomplete
    """
    species_data = species_map.get(taxon_id, {})
    if not species_data:
        raise DataValidationError(f"No species data found for taxon ID: {taxon_id}")
    
    required_fields = ['name', 'db_name', 'short_name']
    for field in required_fields:
        if field not in species_data:
            raise DataValidationError(f"No {field} found for taxon ID: {taxon_id}")

# Replace hardcoded VALID_DATABASES with function call
VALID_DATABASES = set(get_valid_databases())

# 1. Move database aliases to constants at module level
DATABASE_ALIASES = {
    'entrez': {'entrezgene', 'entrez gene/locuslink', 'geneid', 'gene/locuslink'},
    'sgd': {'sgd', 'sgdid'},
    'wb': {'wormbase', 'wormbaseid', 'wb'},
    'fb': {'flybase', 'fbgn'}
}

def get_remapped_database(database: str, taxon_id: str) -> str:
    """Optimized database remapping function"""
    if not database:
        return ''
    
    database = database.lower()
    species_db = get_species_db_name(SPECIES_MAP, taxon_id)
    
    if species_db != 'unknown' and is_valid_database(species_db):
        # Check if database matches any aliases for the species_db
        if database in DATABASE_ALIASES.get(species_db, set()):
            return species_db
    
    return database

def split_gene_id(gene_id: str, taxon_id: str = '') -> tuple:
    """
    Split gene ID into database and ID components.
    
    Args:
        gene_id: String in format 'database:id'
        taxon_id: NCBI taxonomy ID for database remapping
        
    Returns:
        tuple: (database, gene_id)
    """
    if pd.isna(gene_id):
        return ('', '')
    
    parts = gene_id.split(':')
    if len(parts) != 2:
        return ('', gene_id)
    
    database, gene_id = parts
    database = get_remapped_database(database, taxon_id)
    
    return (database, gene_id)

def load_gene_synonyms():
    """Load gene synonyms from JSON file."""
    try:
        with open(SYNONYMS_FILE) as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"Warning: Gene synonyms file not found at {SYNONYMS_FILE}")
        return {}

def map_gene_id(gene_id: str, taxon_id: str, gene_synonyms: dict) -> str:
    """Map gene ID to canonical form using synonyms dictionary."""
    if not taxon_id or taxon_id not in gene_synonyms:
        return gene_id
        
    # Look through all genes in this taxon
    for canonical_id, synonyms in gene_synonyms[taxon_id].items():
        if gene_id in synonyms:
            return canonical_id
    return gene_id

def main():
    # Load gene synonyms at start of main
    gene_synonyms = load_gene_synonyms()
    
    # Change set to list for JSON serialization
    metadata = {
        'species': {},
        'validation_summary': {},
        'invalid_taxons': [],
        'unmatched_databases': []  # Changed from set() to list
    }

    # Read only required columns from the input file
    genetic_interactions = pd.read_csv(
        INPUT_FILE, 
        sep='\t', 
        comment='#', 
        usecols=INTERACTOR_COLS
    )

    # Process data in chunks for better memory usage
    chunk_size = 10000
    processed_chunks = []
    
    for chunk_start in range(0, len(genetic_interactions), chunk_size):
        chunk = genetic_interactions[chunk_start:chunk_start + chunk_size].copy()
        
        # Extract taxon ID and process gene IDs
        chunk['taxonId'] = chunk['Taxid interactor A'].str.extract(r'taxid:(\d+)')
        
        # Process gene IDs and look up synonyms
        chunk[['database', 'fromGeneId']] = pd.DataFrame(
            [split_gene_id(x, t) for x, t in zip(chunk['ID(s) interactor A'], chunk['taxonId'])],
            index=chunk.index
        )
        
        chunk[['database', 'toGeneId']] = pd.DataFrame(
            [split_gene_id(x, t) for x, t in zip(chunk['ID(s) interactor B'], chunk['taxonId'])],
            index=chunk.index
        )
        
        # Map gene IDs using synonyms with taxon ID
        chunk['fromGeneId'] = chunk.apply(
            lambda row: map_gene_id(row['fromGeneId'], row['taxonId'], gene_synonyms),
            axis=1
        )
        
        chunk['toGeneId'] = chunk.apply(
            lambda row: map_gene_id(row['toGeneId'], row['taxonId'], gene_synonyms),
            axis=1
        )
        
        processed_chunks.append(chunk[['database', 'taxonId', 'fromGeneId', 'toGeneId']])

    # Combine processed chunks
    interactions_subset = pd.concat(processed_chunks, ignore_index=True)


    # Display processed data
    print("\nProcessed data:")
    print(interactions_subset)

    print(SPECIES_MAP)

    # Validate all taxon IDs exist in species map
    unique_taxon_ids = interactions_subset['taxonId'].unique()
    for taxon_id in unique_taxon_ids:
        if pd.notna(taxon_id):  # Skip NA values
            try:
                validate_taxon_id(taxon_id, SPECIES_MAP)
                validate_species_data(SPECIES_MAP, taxon_id)
            except DataValidationError as e:
                print(f"Warning: {e}")
                metadata['invalid_taxons'].append(taxon_id)  # Track invalid taxons
                interactions_subset = interactions_subset[interactions_subset['taxonId'] != taxon_id]

    # Ensure we still have valid data after filtering
    assert not interactions_subset.empty, "No valid interactions remaining after taxon ID validation"

    # Optimize species processing
    species_data_dict = {}
    for taxon_id in interactions_subset['taxonId'].unique():
        print(taxon_id)
        if pd.isna(taxon_id):
            continue
            
        species_data = interactions_subset[interactions_subset['taxonId'] == taxon_id]
        all_databases = set(species_data['database'].unique())
        valid_dbs = all_databases & VALID_DATABASES
        
        # Process examples more efficiently
        examples = {
            db: species_data[species_data['database'] == db]
            .head(5)[['fromGeneId', 'toGeneId']]
            .to_dict('records')
            for db in valid_dbs
            if not species_data[species_data['database'] == db].empty
        }
        
        species_data_dict[taxon_id] = {
            'name': get_species_name(SPECIES_MAP, taxon_id),
            'db_name': get_species_db_name(SPECIES_MAP, taxon_id),
            'shortname': get_species_shortname(SPECIES_MAP, taxon_id),
            'examples': examples
        }

    metadata['species'] = species_data_dict

    # Validate we have species data
    assert metadata['species'], "No valid species data found after validation"

    # Add validation summary to metadata
    metadata['validation_summary'] = {
        'total_taxon_ids_found': len(unique_taxon_ids),
        'valid_taxon_ids': len(metadata['species']),
        'invalid_taxon_ids': len(metadata['invalid_taxons']),
        'invalid_taxon_list': metadata['invalid_taxons'],  # Add list of invalid taxons
        'processed_interactions': len(interactions_subset),
        'unmatched_databases': sorted(list(set(metadata['unmatched_databases'])))  # Deduplicate and sort
    }

    # Print summary after processing
    print("\nProcessing Summary:")
    print(f"Total Taxon IDs found: {len(unique_taxon_ids)}")
    print(f"Valid Taxon IDs: {len(metadata['species'])}")
    print(f"Invalid Taxon IDs: {len(metadata['invalid_taxons'])}")
    if metadata['invalid_taxons']:
        print(f"Invalid Taxon IDs list: {metadata['invalid_taxons']}")
    
    if metadata['unmatched_databases']:
        print("\nUnmatched databases found:")
        print(sorted(list(metadata['unmatched_databases'])))
    
    # print("\nDatabase usage by species:")
    # for taxon_id, species_info in metadata['species'].items():
    #     print(f"\n{species_info['name']} (Taxon ID: {taxon_id}):")
    #     print(f"Databases used: {list(species_info['examples'].keys())}")
    #     print("Example gene IDs:")
    #     for db, examples in species_info['examples'].items():
    #         print(f"\n{db}:")
    #         for ex in examples[:3]:  # Show first 3 examples
    #             print(f"  {ex['fromGeneId']} → {ex['toGeneId']}")

    # Save metadata to JSON file
    with open(METADATA_FILE, 'w') as f:
        json.dump(metadata, f, indent=2)


    # Save processed data
    interactions_subset.to_csv(OUTPUT_FILE, index=False)

if __name__ == "__main__":
    try:
        main()
    except (AssertionError, DataValidationError) as e:
        print(f"Error: {e}")
        exit(1)
