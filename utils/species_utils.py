import json
from pathlib import Path
from typing import Dict, Union, List

def load_species_map() -> Dict[str, Dict[str, str]]:
    """
    Load species mapping from JSON file.
    
    Returns:
        Dict[str, Dict[str, str]]: Dictionary mapping taxon IDs to species information
    """
    config_file = Path('data/config/species_map.json')
    try:
        with open(config_file, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"Species map file not found at {config_file}")

def get_species_name(species_map: Dict[str, Dict[str, str]], taxon_id: str) -> str:
    """
    Get species scientific name from taxon ID.
    
    Args:
        species_map: The loaded species map dictionary
        taxon_id: NCBI taxonomy ID
    
    Returns:
        str: Scientific name of the species or 'Unknown species'
    """
    return species_map.get(taxon_id, {}).get('name', 'Unknown species')

def get_species_db_name(species_map: Dict[str, Dict[str, str]], taxon_id: str) -> str:
    """
    Get species database short name from taxon ID.
    
    Args:
        species_map: The loaded species map dictionary
        taxon_id: NCBI taxonomy ID
    
    Returns:
        str: Database short name for the species or 'unknown'
    """
    return species_map.get(taxon_id, {}).get('db_name', 'unknown')

def get_species_shortname(species_map: Dict[str, Dict[str, str]], taxon_id: str) -> str:
    """
    Get species shortname from taxon ID.
    
    Args:
        species_map: The loaded species map dictionary
        taxon_id: NCBI taxonomy ID
    
    Returns:
        str: Shortname of the species or 'unknown'
    """
    return species_map.get(taxon_id, {}).get('short_name', 'unknown')

def is_valid_database(database: str) -> bool:
    """
    Check if a database name is valid according to the species map.
    
    Args:
        database: Database name to check
    
    Returns:
        bool: True if database is valid, False otherwise
    """
    species_map = load_species_map()
    unique_dbs = {species_info['db_name'] for species_info in species_map.values()}
    return database in unique_dbs

def is_valid_species_code(taxon_id: str) -> bool:
    """
    Check if a taxon ID exists in the species map.
    
    Args:
        taxon_id: NCBI taxonomy ID to check
    
    Returns:
        bool: True if taxon ID is valid, False otherwise
    """
    species_map = load_species_map()
    return taxon_id in species_map

def get_valid_databases() -> List[str]:
    """
    Get list of valid database names from species map.
    
    Returns:
        List[str]: List of unique database names
    """
    species_map = load_species_map()
    return sorted(list({species_info['db_name'] for species_info in species_map.values()}))

def get_valid_species_codes() -> List[str]:
    """
    Get list of valid taxon IDs.
    
    Returns:
        List[str]: List of valid taxon IDs
    """
    species_map = load_species_map()
    return sorted(list(species_map.keys()))