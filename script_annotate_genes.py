#!/usr/bin/env python3
"""
==============================================================================
COMPREHENSIVE GENE ANNOTATION SCRIPT
==============================================================================

This script annotates genes with pathways from multiple databases:
- KEGG (Kyoto Encyclopedia of Genes and Genomes)
- REACTOME (Reactome Pathway Database)
- Hallmark (MSigDB Hallmark Gene Sets)
- WikiPathways (Community-curated pathway database)
- GO (Gene Ontology): Biological Process, Molecular Function, Cellular Component

INSTALLATION INSTRUCTIONS:
==============================================================================
1. Install Python 3.8 or higher
2. Install required packages:
   pip install pandas mygene gseapy requests

USAGE:
==============================================================================
python gene_annotator_standalone.py <input_file.csv> [output_file.csv]

Example:
python gene_annotator_standalone.py my_genes.csv annotated_genes.csv

INPUT FILE FORMAT:
==============================================================================
The script expects a CSV file with gene symbols in the first column.
Header row is expected. Example:

Gene,Other_Column1,Other_Column2
BRCA1,value1,value2
TP53,value3,value4
...

OUTPUT FILE FORMAT:
==============================================================================
Gene,KEGG_Count,KEGG_Pathways,REACTOME_Count,REACTOME_Pathways,Hallmark_Count,Hallmark_Pathways,WikiPathways_Count,WikiPathways,GO_BP_Count,GO_Biological_Process,GO_MF_Count,GO_Molecular_Function,GO_CC_Count,GO_Cellular_Component

==============================================================================
"""

import pandas as pd
import sys
import os

# Try to import required libraries
try:
    import mygene
    MYGENE_AVAILABLE = True
except ImportError:
    print("WARNING: mygene not installed. Install with: pip install mygene")
    MYGENE_AVAILABLE = False

try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    print("WARNING: gseapy not installed. Install with: pip install gseapy")
    GSEAPY_AVAILABLE = False

def get_gene_list(input_file):
    """Extract unique gene names from input file"""
    print(f"\n{'='*70}")
    print(f"Reading gene list from: {input_file}")
    print(f"{'='*70}")
    
    df = pd.read_csv(input_file)
    
    # Assuming first column contains gene names
    gene_column = df.columns[0]
    genes = df[gene_column].dropna().unique().tolist()
    
    print(f"✓ Found {len(genes)} unique genes in column '{gene_column}'")
    return genes

def annotate_with_mygene(genes, batch_size=500):
    """Annotate genes using MyGene.info API"""
    print(f"\n{'='*70}")
    print("METHOD 1: MyGene.info Annotation")
    print(f"{'='*70}")
    
    if not MYGENE_AVAILABLE:
        print("✗ Skipping: mygene package not installed")
        return {}
    
    print(f"Querying MyGene.info for {len(genes)} genes...")
    mg = mygene.MyGeneInfo()
    
    # Query in batches
    all_results = []
    for i in range(0, len(genes), batch_size):
        batch = genes[i:i+batch_size]
        batch_num = i//batch_size + 1
        total_batches = (len(genes)-1)//batch_size + 1
        
        print(f"  Batch {batch_num}/{total_batches}: Processing {len(batch)} genes...")
        
        try:
            results = mg.querymany(
                batch,
                scopes='symbol',
                fields='go,pathway.kegg,pathway.reactome,pathway.wikipathways',
                species='human',
                returnall=True,
                as_dataframe=False
            )
            all_results.extend(results['out'])
        except Exception as e:
            print(f"  ✗ Error in batch {batch_num}: {str(e)[:60]}")
            continue
    
    # Create lookup dictionary
    gene_lookup = {}
    for result in all_results:
        if 'query' in result and not result.get('notfound', False):
            gene_lookup[result['query']] = result
    
    print(f"✓ Successfully annotated {len(gene_lookup)} genes via MyGene.info")
    return gene_lookup

def get_msigdb_libraries():
    """Download gene sets from MSigDB via gseapy"""
    print(f"\n{'='*70}")
    print("METHOD 2: MSigDB/Enrichr Annotation")
    print(f"{'='*70}")
    
    if not GSEAPY_AVAILABLE:
        print("✗ Skipping: gseapy package not installed")
        return {}
    
    libraries = {}
    db_names = {
        'KEGG_2021': 'KEGG_2021_Human',
        'Reactome_2022': 'Reactome_2022',
        'Hallmark': 'MSigDB_Hallmark_2020',
        'WikiPathways': 'WikiPathway_2023_Human',
        'GO_BP': 'GO_Biological_Process_2023',
        'GO_MF': 'GO_Molecular_Function_2023',
        'GO_CC': 'GO_Cellular_Component_2023'
    }
    
    print("Downloading gene set libraries from Enrichr/MSigDB...")
    for db_key, db_name in db_names.items():
        try:
            print(f"  {db_key}...", end=' ')
            library = gp.get_library(name=db_name, organism='Human')
            libraries[db_key] = library
            print(f"✓ ({len(library)} gene sets)")
        except Exception as e:
            print(f"✗ Error: {str(e)[:40]}")
            libraries[db_key] = {}
    
    return libraries

def extract_pathways_from_mygene(gene_data, gene):
    """Extract pathway information from MyGene result"""
    result = {
        'kegg': [],
        'reactome': [],
        'wikipathways': []
    }
    
    if not gene_data:
        return result
    
    # KEGG pathways
    if 'pathway' in gene_data and gene_data['pathway']:
        pathway_data = gene_data['pathway']
        
        # KEGG
        if 'kegg' in pathway_data:
            kegg = pathway_data['kegg']
            if isinstance(kegg, list):
                result['kegg'] = [p.get('name', '') for p in kegg if isinstance(p, dict)]
            elif isinstance(kegg, dict):
                result['kegg'] = [kegg.get('name', '')]
        
        # Reactome
        if 'reactome' in pathway_data:
            reactome = pathway_data['reactome']
            if isinstance(reactome, list):
                result['reactome'] = [p.get('name', '') for p in reactome if isinstance(p, dict)]
            elif isinstance(reactome, dict):
                result['reactome'] = [reactome.get('name', '')]
        
        # WikiPathways
        if 'wikipathways' in pathway_data:
            wiki = pathway_data['wikipathways']
            if isinstance(wiki, list):
                result['wikipathways'] = [p.get('name', '') for p in wiki if isinstance(p, dict)]
            elif isinstance(wiki, dict):
                result['wikipathways'] = [wiki.get('name', '')]
    
    return result

def extract_go_from_mygene(gene_data):
    """Extract GO terms from MyGene result"""
    go_result = {
        'BP': [],
        'MF': [],
        'CC': []
    }
    
    if not gene_data or 'go' not in gene_data:
        return go_result
    
    go_data = gene_data['go']
    
    for category in ['BP', 'MF', 'CC']:
        if category in go_data:
            terms = go_data[category]
            if isinstance(terms, list):
                go_result[category] = [t.get('term', '') for t in terms if isinstance(t, dict)]
            elif isinstance(terms, dict):
                go_result[category] = [terms.get('term', '')]
    
    return go_result

def find_gene_in_libraries(gene, libraries):
    """Find pathways containing the gene in MSigDB libraries"""
    result = {
        'kegg': [],
        'reactome': [],
        'hallmark': [],
        'wikipathways': [],
        'go_bp': [],
        'go_mf': [],
        'go_cc': []
    }
    
    gene_upper = gene.upper()
    
    # Map library keys to result keys
    lib_mapping = {
        'KEGG_2021': 'kegg',
        'Reactome_2022': 'reactome',
        'Hallmark': 'hallmark',
        'WikiPathways': 'wikipathways',
        'GO_BP': 'go_bp',
        'GO_MF': 'go_mf',
        'GO_CC': 'go_cc'
    }
    
    for lib_key, result_key in lib_mapping.items():
        if lib_key not in libraries:
            continue
        
        library = libraries[lib_key]
        for pathway_name, gene_list in library.items():
            gene_list_upper = [g.upper() for g in gene_list]
            if gene_upper in gene_list_upper:
                # Clean pathway name
                clean_name = pathway_name
                for prefix in ['KEGG_', 'REACTOME_', 'HALLMARK_', 'WP_', 'GOBP_', 'GOMF_', 'GOCC_']:
                    if clean_name.startswith(prefix):
                        clean_name = clean_name[len(prefix):]
                        break
                
                clean_name = clean_name.replace('_', ' ')
                result[result_key].append(clean_name)
    
    return result

def merge_annotations(mygene_annotation, msigdb_annotation):
    """Merge annotations from both sources"""
    merged = {}
    
    # For each database, combine unique entries
    for db in ['kegg', 'reactome', 'wikipathways']:
        combined = set(mygene_annotation.get(db, []) + msigdb_annotation.get(db, []))
        merged[db] = sorted(list(combined))
    
    # Hallmark only from MSigDB
    merged['hallmark'] = msigdb_annotation.get('hallmark', [])
    
    # GO terms
    merged['go_bp'] = list(set(mygene_annotation.get('go_bp', []) + msigdb_annotation.get('go_bp', [])))
    merged['go_mf'] = list(set(mygene_annotation.get('go_mf', []) + msigdb_annotation.get('go_mf', [])))
    merged['go_cc'] = list(set(mygene_annotation.get('go_cc', []) + msigdb_annotation.get('go_cc', [])))
    
    return merged

def create_annotations(genes, mygene_lookup, msigdb_libraries):
    """Create comprehensive annotations for all genes"""
    print(f"\n{'='*70}")
    print("Creating Comprehensive Annotations")
    print(f"{'='*70}")
    
    annotations = []
    total = len(genes)
    
    for idx, gene in enumerate(genes, 1):
        if idx % 100 == 0 or idx == 1:
            print(f"  Progress: {idx}/{total} genes processed ({100*idx//total}%)")
        
        # Get MyGene annotation
        mygene_data = mygene_lookup.get(gene, {})
        pathway_mygene = extract_pathways_from_mygene(mygene_data, gene)
        go_mygene = extract_go_from_mygene(mygene_data)
        
        mygene_annotation = {
            'kegg': pathway_mygene['kegg'],
            'reactome': pathway_mygene['reactome'],
            'wikipathways': pathway_mygene['wikipathways'],
            'go_bp': go_mygene['BP'],
            'go_mf': go_mygene['MF'],
            'go_cc': go_mygene['CC']
        }
        
        # Get MSigDB annotation
        msigdb_annotation = find_gene_in_libraries(gene, msigdb_libraries)
        
        # Merge annotations
        final_annotation = merge_annotations(mygene_annotation, msigdb_annotation)
        
        # Create annotation row
        annotation = {
            'Gene': gene,
            'KEGG_Count': len(final_annotation['kegg']),
            'KEGG_Pathways': '; '.join(final_annotation['kegg']) if final_annotation['kegg'] else 'None',
            'REACTOME_Count': len(final_annotation['reactome']),
            'REACTOME_Pathways': '; '.join(final_annotation['reactome']) if final_annotation['reactome'] else 'None',
            'Hallmark_Count': len(final_annotation['hallmark']),
            'Hallmark_Pathways': '; '.join(final_annotation['hallmark']) if final_annotation['hallmark'] else 'None',
            'WikiPathways_Count': len(final_annotation['wikipathways']),
            'WikiPathways': '; '.join(final_annotation['wikipathways']) if final_annotation['wikipathways'] else 'None',
            'GO_BP_Count': len(final_annotation['go_bp']),
            'GO_Biological_Process': '; '.join(final_annotation['go_bp'][:100]) if final_annotation['go_bp'] else 'None',
            'GO_MF_Count': len(final_annotation['go_mf']),
            'GO_Molecular_Function': '; '.join(final_annotation['go_mf'][:100]) if final_annotation['go_mf'] else 'None',
            'GO_CC_Count': len(final_annotation['go_cc']),
            'GO_Cellular_Component': '; '.join(final_annotation['go_cc'][:100]) if final_annotation['go_cc'] else 'None',
        }
        
        annotations.append(annotation)
    
    print(f"✓ Completed annotations for all {total} genes")
    
    df = pd.DataFrame(annotations)
    return df

def print_summary(df):
    """Print annotation summary statistics"""
    print(f"\n{'='*70}")
    print("ANNOTATION SUMMARY")
    print(f"{'='*70}")
    print(f"Total genes annotated: {len(df)}")
    print(f"\nPathway Coverage:")
    print(f"  KEGG pathways:      {(df['KEGG_Count'] > 0).sum():4d} genes ({100*(df['KEGG_Count'] > 0).sum()/len(df):.1f}%)")
    print(f"  REACTOME pathways:  {(df['REACTOME_Count'] > 0).sum():4d} genes ({100*(df['REACTOME_Count'] > 0).sum()/len(df):.1f}%)")
    print(f"  Hallmark pathways:  {(df['Hallmark_Count'] > 0).sum():4d} genes ({100*(df['Hallmark_Count'] > 0).sum()/len(df):.1f}%)")
    print(f"  WikiPathways:       {(df['WikiPathways_Count'] > 0).sum():4d} genes ({100*(df['WikiPathways_Count'] > 0).sum()/len(df):.1f}%)")
    print(f"\nGO Term Coverage:")
    print(f"  Biological Process: {(df['GO_BP_Count'] > 0).sum():4d} genes ({100*(df['GO_BP_Count'] > 0).sum()/len(df):.1f}%)")
    print(f"  Molecular Function: {(df['GO_MF_Count'] > 0).sum():4d} genes ({100*(df['GO_MF_Count'] > 0).sum()/len(df):.1f}%)")
    print(f"  Cellular Component: {(df['GO_CC_Count'] > 0).sum():4d} genes ({100*(df['GO_CC_Count'] > 0).sum()/len(df):.1f}%)")
    print(f"{'='*70}")

def main():
    """Main execution function"""
    if len(sys.argv) < 2:
        print(__doc__)
        print("\nERROR: No input file specified")
        print("\nUsage: python gene_annotator_standalone.py <input_file.csv> [output_file.csv]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'annotated_genes.csv'
    
    if not os.path.exists(input_file):
        print(f"ERROR: Input file '{input_file}' not found")
        sys.exit(1)
    
    print("\n" + "="*70)
    print("COMPREHENSIVE GENE ANNOTATION PIPELINE")
    print("="*70)
    
    # Step 1: Read genes
    genes = get_gene_list(input_file)
    
    # Step 2: Get MyGene annotations
    mygene_lookup = annotate_with_mygene(genes)
    
    # Step 3: Get MSigDB libraries
    msigdb_libraries = get_msigdb_libraries()
    
    # Step 4: Create comprehensive annotations
    annotated_df = create_annotations(genes, mygene_lookup, msigdb_libraries)
    
    # Step 5: Save results
    print(f"\n{'='*70}")
    print(f"Saving results to: {output_file}")
    print(f"{'='*70}")
    annotated_df.to_csv(output_file, index=False)
    print(f"✓ Output file saved successfully")
    
    # Step 6: Print summary
    print_summary(annotated_df)
    
    # Show examples
    print(f"\nFirst 5 genes (preview):")
    print(annotated_df[['Gene', 'KEGG_Count', 'REACTOME_Count', 'Hallmark_Count', 'GO_BP_Count']].head().to_string(index=False))
    
    print(f"\n{'='*70}")
    print("ANNOTATION COMPLETE!")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
