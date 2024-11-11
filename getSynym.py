from utils.species_utils import is_valid_species_code, load_species_map
def parse_synonyms(line, field_index):
    # Split the line by tabs
    fields = line.strip().split('\t')
    
    # Modified to store tuples of (synonym, db_name)
    synonyms = []
    # Check fields[4] for interactor A aliases and fields[5] for interactor B aliases
    for field in [fields[field_index]]:  # Changed to specifically look at alias fields
        if "(" in field and ")" in field:  # Check for any bracketed text
            # Split by pipe character for multiple synonyms
            syns = field.split('|')
            for syn in syns:
                if "(" in syn and ")" in syn:  # Check each synonym for bracketed text
                    # Extract the database name from within brackets
                    db_in_brackets = syn[syn.find("(")+1:syn.find(")")]
                    
                    # Skip if it contains "public_name"
                    if "public_name" in db_in_brackets:
                        continue
                        
                    # Extract the part after the colon if it exists, otherwise use the whole name
                    name = syn.split('(')[0]
                    # Extract database name if present
                    db_name = None
                    if ':' in name:
                        parts = name.split(':')
                        db_name = parts[0]
                        name = parts[1]
                    
                    synonyms.append((name.strip(), db_name))
    
    return synonyms

def get_taxon_id(field):
    # Extract taxon ID from field containing format "taxid:6239(caeel)|taxid:6239(...)"
    if field and 'taxid:' in field:
        return field.split('|')[0].split('(')[0].split(':')[1]
    return None

def is_valid_species_code(taxon):
    # Import the utility function to check valid species codes
    from utils.species_utils import is_valid_species_code
    return is_valid_species_code(taxon)

def get_database_name(gene_id):
    # Extract database name from gene ID like "wormbase:WBGene00002996"
    if gene_id and ':' in gene_id:
        return gene_id.split(':')[0]
    return None

def find_gene_in_descriptions(gene_id, synonyms, gene_descriptions, taxon):
    """
    Check if gene_id or any of its synonyms match with gene descriptions for a specific taxon.
    Returns (found_gene_id, original_term) if found, (None, None) if not found
    """
    # First check if taxon exists in descriptions
    if taxon not in gene_descriptions:
        return None, None
        
    # Check the gene_id itself
    if gene_id in gene_descriptions[taxon]:
        return gene_id, None
    
    # Then check all synonyms
    for synonym, db_name in synonyms:
        if synonym in gene_descriptions[taxon]:
            return synonym, gene_id
    
    return None, None

def load_gene_descriptions():
    gene_desc_files = [
        "data/raw/GeneDescriptions/GENE-DESCRIPTION-TSV_FB.tsv",
        "data/raw/GeneDescriptions/GENE-DESCRIPTION-TSV_HUMAN.tsv",
        "data/raw/GeneDescriptions/GENE-DESCRIPTION-TSV_MGI.tsv",
        "data/raw/GeneDescriptions/GENE-DESCRIPTION-TSV_SGD.tsv",
        "data/raw/GeneDescriptions/GENE-DESCRIPTION-TSV_WB.tsv",
        "data/raw/GeneDescriptions/GENE-DESCRIPTION-TSV_XBXL.tsv",
        "data/raw/GeneDescriptions/GENE-DESCRIPTION-TSV_XBXT.tsv",
        "data/raw/GeneDescriptions/GENE-DESCRIPTION-TSV_ZFIN.tsv"
    ]
    
    # Load species map
    species_map = load_species_map()
    # Create reverse lookup from short_name to taxon_id
    short_name_to_taxon = {info['short_name']: taxon_id 
                          for taxon_id, info in species_map.items()}
    
    gene_descriptions = {}
    for file_path in gene_desc_files:
        try:
            # Extract short name from filename (e.g., GENE-DESCRIPTION-TSV_FB.tsv -> fb)
            short_name = file_path.split('_')[-1].split('.')[0].lower()
            # Get actual taxon ID from species map
            taxon_id = short_name_to_taxon.get(short_name)
            
            if not taxon_id:
                print(f"Warning: Could not find taxon ID for {short_name} from {file_path}")
                continue

            print(taxon_id)    
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    fields = line.strip().split('\t')
                    gene_id = fields[0].split(':')[1] if ':' in fields[0] else fields[0]
                    
                    if taxon_id not in gene_descriptions:
                        gene_descriptions[taxon_id] = set()  # Using a set to avoid duplicates
                    gene_descriptions[taxon_id].add(gene_id)
        except FileNotFoundError:
            print(f"Warning: Could not find file {file_path}")
    return gene_descriptions
    
def process_interaction_file(filename, gene_descriptions):
    # Modified to organize by taxon ID
    gene_synonyms = {}
    formatted_synonyms_dict = {}
    taxon_db_pairs = set()
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            gene_a = fields[0].split(':')[1] if ':' in fields[0] else fields[0]
            gene_b = fields[1].split(':')[1] if ':' in fields[1] else fields[1]
            taxon_a = get_taxon_id(fields[9])
            taxon_b = get_taxon_id(fields[10])
            
            # Get database names
            db_a = get_database_name(fields[0])
            db_b = get_database_name(fields[1])
            
            # Store unique combinations
            if taxon_a and db_a:
                taxon_db_pairs.add((taxon_a, db_a))
            if taxon_b and db_b:
                taxon_db_pairs.add((taxon_b, db_b))
            
            # Initialize taxon dictionaries if they don't exist
            if taxon_a and taxon_a not in formatted_synonyms_dict:
                formatted_synonyms_dict[taxon_a] = {}
            if taxon_b and taxon_b not in formatted_synonyms_dict:
                formatted_synonyms_dict[taxon_b] = {}
            
            # Process interactor A
            if taxon_a:
                synonyms_a = parse_synonyms(line, 4)
                found_id_a, original_a = find_gene_in_descriptions(gene_a, synonyms_a, gene_descriptions, taxon_a)
                if found_id_a:
                    if found_id_a not in formatted_synonyms_dict[taxon_a]:
                        formatted_synonyms_dict[taxon_a][found_id_a] = {}
                    
                    # Add synonyms to the dictionary
                    for syn, db in synonyms_a:
                        formatted_synonyms_dict[taxon_a][found_id_a][syn] = db
                    
                    # Add original ID if it exists
                    if original_a:
                        formatted_synonyms_dict[taxon_a][found_id_a][original_a] = db_a
            
            # Process interactor B
            if taxon_b:
                synonyms_b = parse_synonyms(line, 5)
                found_id_b, original_b = find_gene_in_descriptions(gene_b, synonyms_b, gene_descriptions, taxon_b)
                if found_id_b:
                    if found_id_b not in formatted_synonyms_dict[taxon_b]:
                        formatted_synonyms_dict[taxon_b][found_id_b] = {}
                    
                    # Add synonyms to the dictionary
                    for syn, db in synonyms_b:
                        formatted_synonyms_dict[taxon_b][found_id_b][syn] = db
                    
                    # Add original ID if it exists
                    if original_b:
                        formatted_synonyms_dict[taxon_b][found_id_b][original_b] = db_b
    
    return taxon_db_pairs, formatted_synonyms_dict

# Example usage
if __name__ == "__main__":
    filename = "data/raw/MolecularInteractions/INTERACTION-MOL_COMBINED.tsv"
    try:
        # First load all gene descriptions
        gene_descriptions = load_gene_descriptions()
        # Load species map for db_name lookup
        species_map = load_species_map()
        print(f"Loaded {len(gene_descriptions)} gene descriptions")
        
        # Then process interaction file
        taxon_db_pairs, synonyms_dict = process_interaction_file(filename, gene_descriptions)
        
        # Save synonyms dictionary to file
        import json
        output_file = "data/processed/gene_synonyms.json"
        with open(output_file, 'w') as f:
            json.dump(synonyms_dict, f, indent=2)
        print(f"\nSaved synonyms dictionary to {output_file}")
        
        # Print sample of the saved format
        print("\nSample of saved synonym format:")
        sample_items = list(synonyms_dict.items())[:3]
        for gene, syn_dict in sample_items:
            print(f"{gene}: {syn_dict}")
            
    except Exception as e:
        print(f"Error processing file: {str(e)}")