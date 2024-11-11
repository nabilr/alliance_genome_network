import pandas as pd
import os
import json

def process_description_file(filepath):
    """Process a single gene description TSV file."""
    # Skip the header comments that start with #
    with open(filepath, 'r') as f:
        lines = f.readlines()
        start_line = 0
        for i, line in enumerate(lines):
            if not line.startswith('#'):
                start_line = i
                break
    
    # Read TSV file starting after comments with correct column names
    df = pd.read_csv(filepath, sep='\t', skiprows=start_line, header=None,
                     names=['fullGeneId', 'Symbol', 'Description'])
    
    # Split database and geneId
    df[['database', 'geneId']] = df['fullGeneId'].str.split(':', n=1, expand=True)
    
    # Standardize database names
    df['database'] = df['database'].str.lower()
    
    # Replace specific database names
    df['database'] = df['database'].replace({
        'fb': 'fb',
        'wb': 'wb',
        'entrezgene': 'entrez',
        'gene/locuslink': 'entrez',
        'geneid': 'entrez',
        'sgdid': 'sgd'
    })
    
    # Extract species from filename and convert to standard format
    species = filepath.split('_')[-1].split('.')[0]
    df['Species'] = species.upper()  # Matches "HUMAN" format seen in the image
    
    # Select and reorder columns
    df = df[['database', 'geneId', 'Symbol', 'Description', 'Species']]
    
    return df

def combine_descriptions(input_dir, output_dir):
    """Combine all gene description files in the directory."""
    all_data = []
    metadata = {}
    species_set = set()  # New set to collect species
    
    # Process each TSV file in the directory
    for filename in os.listdir(input_dir):
        if filename.startswith('GENE-DESCRIPTION-TSV_') and filename.endswith('.tsv'):
            # Extract species from filename
            species = filename.split('_')[-1].split('.')[0].lower()
            species_set.add(species)
            
            filepath = os.path.join(input_dir, filename)
            df = process_description_file(filepath)
            all_data.append(df)
    
    # Combine all dataframes
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Print and collect unique databases and species
    unique_databases = sorted(combined_df['database'].unique())
    unique_species = sorted(species_set)
    print("\nFound databases:", unique_databases)
    print("\nFound species:", unique_species)
    
    metadata['databases'] = unique_databases
    metadata['species'] = unique_species  # Add species to metadata
    
    # Save metadata to JSON file
    metadata_file = os.path.join(output_dir, 'species_metadata.json')
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"\nMetadata saved to: {metadata_file}")
    
    # Clean descriptions
    combined_df['Description'] = combined_df['Description'].str.replace('"', "'")
    
    # Create output CSV file
    output_file = os.path.join(output_dir, 'gene_nodes.csv')
    combined_df.to_csv(output_file, index=False, quoting=1)
    
    return output_file

# Usage
input_directory = "data/raw/GeneDescriptions"
output_directory = "data/processed/GeneDescriptions"
output_file = combine_descriptions(input_directory, output_directory)
print(f"\nCombined descriptions saved to: {output_file}")