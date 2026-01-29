import argparse
import pandas as pd
from collections import defaultdict
from Bio import SearchIO
import json

#######################
# Auxiliary functions #
#######################

# Function to select the row with the lowest evalue or randomly if there is a tie
def select_lowest_evalue(group):
    # Sort the group by 'evalue' and take the first row of the sorted group
    return group.sort_values(by='evalue').head(1)

# Function to select the row with the lowest evalue or randomly if there is a tie
def select_highest_confidence(group):
    # Sort the group by 'evalue' and take the first row of the sorted group
    return group.sort_values(by='confidence', ascending=False).head(1)

# Function to extract the part after the underscore in ID and append it to seqid
def append_suffix_to_seqid(row):
    # Extract the part after 'ID=' and before the first ';'
    id_part = row['attributes'].split(';')[0].replace('ID=', '')
    # Extract the part after the underscore
    suffix = id_part.split('_')[-1]
    return f"{row['seqid']}_{suffix}"

#################
# Main function #
#################

def merge_annotations(gff_file, kegg_file, keggdb_file, pfam_file, ec_file, cazy_file, vf_file, vfdb_file, amr_file, amrdb_file, signalp_file, output_file, defense_file=None):

    ##############
    # Load genes #
    ##############

    annotations = pd.read_csv(gff_file, sep='\t', comment='#', header=None,
                     names=['seqid', 'source', 'type', 'start', 'end',
                            'score', 'strand', 'phase', 'attributes'])
    annotations['seqid'] = annotations.apply(append_suffix_to_seqid, axis=1)
    annotations = annotations.drop(columns=['attributes', 'source', 'score', 'type', 'phase'])
    annotations = annotations.rename(columns={'seqid': 'gene'})

    ######################
    # Load mapping files #
    ######################

    #PFAM to EC
    pfam_to_ec = pd.read_csv(ec_file, sep='\t', comment='#', header=0)
    pfam_to_ec = pfam_to_ec[pfam_to_ec['Type'] == 'GOLD']
    pfam_to_ec = pfam_to_ec.rename(columns={'Confidence-Score': 'confidence'})
    pfam_to_ec = pfam_to_ec.rename(columns={'Pfam-Domain': 'pfam'})
    pfam_to_ec = pfam_to_ec.rename(columns={'EC-Number': 'ec'})
    pfam_to_ec['confidence'] = pd.to_numeric(pfam_to_ec['confidence'], errors='coerce')
    pfam_to_ec = pfam_to_ec.groupby('pfam', group_keys=False)[['pfam','ec','confidence']].apply(select_highest_confidence, include_groups=False)

    #Entry to VF
    entry_to_vf = pd.read_csv(vfdb_file, sep='\t', comment='#', header=0)

    #AMR to class
    amr_to_class = pd.read_csv(amrdb_file, sep='\t', header=0)
    amr_to_class = amr_to_class.rename(columns={'#hmm_accession': 'accession'})

    #KEGG hierarchy
    kegg_json=json.load(open(keggdb_file))
    kegg_df = []
    for main in kegg_json['children']:
        for broad in main['children']:
            for sub in broad['children']:
                for gn in sub.get('children', [None]):
                    if gn:
                        # Extract the relevant parts of the 'name' field
                        name = gn['name']
                        kegg = name.split(' ')[0] if len(name.split(' ')) > 0 else ''
                        description_and_ec = ' '.join(name.split(' ')[1:])

                        # Split the description and EC number
                        description_parts = description_and_ec.split(' [')
                        description = description_parts[0]
                        ec = description_parts[1][:-1] if len(description_parts) > 1 else ''  # Remove the trailing ']'

                        # Append data to the list
                        kegg_df.append([kegg, description, ec, sub['name'], broad['name'], main['name']])
                    else:
                        # Handle terms with no children if necessary
                        pass
    kegg_hierachy = pd.DataFrame(kegg_df, columns=['kegg', 'Description', 'ec', 'Subcategory', 'Broad Category', 'Main Category'])

    #####################
    # Parse annotations #
    #####################

    hmm_attribs = ['accession', 'bitscore', 'evalue', 'id', 'overlap_num', 'region_num']
    evalue_threshold=0.00001

    # Parse keggdb
    keggdb_hits = defaultdict(list)
    query_ids = []

    with open(kegg_file) as handle:
        for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
            query_id = queryresult.id  # Capture the query result id
            query_ids.extend([query_id] * len(queryresult.hits))  # Extend query_ids list to match the number of hits

            for hit in queryresult.hits:
                for attrib in hmm_attribs:
                    # Use `getattr` to fetch the attribute value from the hit object
                    value = getattr(hit, attrib, None)  # Use `None` as a default if the attribute does not exist
                    keggdb_hits[attrib].append(value)

    keggdb_hits['query_id'] = query_ids
    keggdb_df = pd.DataFrame.from_dict(keggdb_hits)
    keggdb_df = keggdb_df.rename(columns={'query_id': 'gene'})
    keggdb_df['evalue'] = pd.to_numeric(keggdb_df['evalue'], errors='coerce')
    keggdb_df = keggdb_df[keggdb_df['evalue'] < evalue_threshold]
    keggdb_df = keggdb_df.rename(columns={'id': 'kegg'})
    keggdb_df = pd.merge(keggdb_df, kegg_hierachy[['kegg', 'ec']], on='kegg', how='left')
    keggdb_df = keggdb_df.groupby('gene', group_keys=False)[['gene','kegg','ec','evalue']].apply(select_lowest_evalue, include_groups=False).reset_index(drop=True)
    keggdb_df['ec'] = keggdb_df['ec'].str.replace('EC:', '', regex=False)

    # Parse PFAM
    pfam_hits = defaultdict(list)
    query_ids = []

    with open(pfam_file) as handle:
        for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
            query_id = queryresult.id  # Capture the query result id
            query_ids.extend([query_id] * len(queryresult.hits))  # Extend query_ids list to match the number of hits

            for hit in queryresult.hits:
                for attrib in hmm_attribs:
                    # Use `getattr` to fetch the attribute value from the hit object
                    value = getattr(hit, attrib, None)  # Use `None` as a default if the attribute does not exist
                    pfam_hits[attrib].append(value)

    pfam_hits['query_id'] = query_ids
    pfam_df = pd.DataFrame.from_dict(pfam_hits)
    pfam_df = pfam_df.rename(columns={'query_id': 'gene'})
    pfam_df = pfam_df.rename(columns={'accession': 'pfam'})
    pfam_df['pfam'] = pfam_df['pfam'].str.split('.').str[0]
    pfam_df['evalue'] = pd.to_numeric(pfam_df['evalue'], errors='coerce')
    pfam_df = pfam_df[pfam_df['evalue'] < evalue_threshold]
    pfam_df = pfam_df.groupby('gene', group_keys=False)[['gene','pfam','evalue']].apply(select_lowest_evalue, include_groups=False).reset_index(drop=True)
    pfam_df = pd.merge(pfam_df, pfam_to_ec[['pfam', 'ec']], on='pfam', how='left')

    # Parse CAZY
    cazy_hits = defaultdict(list)
    query_ids = []

    with open(cazy_file) as handle:
        for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
            query_id = queryresult.id  # Capture the query result id
            query_ids.extend([query_id] * len(queryresult.hits))  # Extend query_ids list to match the number of hits

            for hit in queryresult.hits:
                for attrib in hmm_attribs:
                    # Use `getattr` to fetch the attribute value from the hit object
                    value = getattr(hit, attrib, None)  # Use `None` as a default if the attribute does not exist
                    cazy_hits[attrib].append(value)

    cazy_hits['query_id'] = query_ids
    cazy_df = pd.DataFrame.from_dict(cazy_hits)
    cazy_df = cazy_df.rename(columns={'query_id': 'gene'})
    cazy_df['id'] = cazy_df['id'].str.replace('.hmm', '', regex=False)
    cazy_df['evalue'] = pd.to_numeric(cazy_df['evalue'], errors='coerce')
    cazy_df = cazy_df[cazy_df['evalue'] < evalue_threshold]
    cazy_df = cazy_df.rename(columns={'id': 'cazy'})
    cazy_df = cazy_df.groupby('gene', group_keys=False)[['gene','cazy','evalue']].apply(select_lowest_evalue, include_groups=False).reset_index(drop=True)

    # Parse AMR
    amr_hits = defaultdict(list)
    query_ids = []

    with open(amr_file) as handle:
        for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
            query_id = queryresult.id  # Capture the query result id
            query_ids.extend([query_id] * len(queryresult.hits))  # Extend query_ids list to match the number of hits

            for hit in queryresult.hits:
                for attrib in hmm_attribs:
                    # Use `getattr` to fetch the attribute value from the hit object
                    value = getattr(hit, attrib, None)  # Use `None` as a default if the attribute does not exist
                    amr_hits[attrib].append(value)

    amr_hits['query_id'] = query_ids
    amr_df = pd.DataFrame.from_dict(amr_hits)
    amr_df = amr_df.rename(columns={'query_id': 'gene'})
    amr_df['evalue'] = pd.to_numeric(amr_df['evalue'], errors='coerce')
    amr_df = amr_df[amr_df['evalue'] < evalue_threshold]
    amr_df = amr_df.rename(columns={'id': 'amr'})
    amr_df = amr_df.groupby('gene', group_keys=False)[['gene','amr','accession','evalue']].apply(select_lowest_evalue, include_groups=False).reset_index(drop=True)
    amr_df = pd.merge(amr_df, amr_to_class[['accession','subtype','subclass']], on='accession', how='left')
    amr_df = amr_df.rename(columns={'subtype': 'resistance_type'})
    amr_df = amr_df.rename(columns={'subclass': 'resistance_target'})

    # Parse VFDB
    vfdb_df = pd.read_csv(vf_file, sep='\t', comment='#', header=None,
                     names=['gene', 'entry', 'identity', 'length', 'mismatches',
                            'gaps', 'query_start', 'query_end', 'target_start', 'target_end', 'evalue', 'bitscore'])
    vfdb_df['evalue'] = pd.to_numeric(vfdb_df['evalue'], errors='coerce')
    vfdb_df = vfdb_df[vfdb_df['evalue'] < evalue_threshold]
    vfdb_df = vfdb_df.groupby('gene', group_keys=False)[['gene','entry','evalue']].apply(select_lowest_evalue, include_groups=False).reset_index(drop=True)
    vfdb_df = pd.merge(vfdb_df, entry_to_vf[['entry','vf','vfc','vf_type']], on='entry', how='left')

    # Parse SIGNALP
    signalp_df = pd.read_csv(signalp_file, sep='\t', comment='#', header=None, names=['gene', 'signalp', 'confidence'])
    signalp_df['confidence'] = pd.to_numeric(signalp_df['confidence'], errors='coerce')
    signalp_df = signalp_df.groupby('gene', group_keys=False)[['gene','signalp','confidence']].apply(select_highest_confidence).reset_index(drop=True)

    #####################
    # Merge annotations #
    #####################

    annotations = pd.merge(annotations, keggdb_df[['gene', 'kegg', 'ec']], on='gene', how='left')
    # Perform the merge, adding 'ec' from pfam_df as 'pfam_ec' to avoid conflict
    annotations = pd.merge(annotations, pfam_df[['gene', 'pfam', 'ec']], on='gene', how='left', suffixes=('', '_pfam'))
    # Update 'ec' column in annotations only where it is empty
    annotations['ec'] = annotations.apply(
        lambda row: row['ec_pfam'] if pd.isna(row['ec']) or row['ec'] == '' else row['ec'],
        axis=1)
    # Drop the temporary 'ec_pfam' column
    annotations.drop(columns=['ec_pfam'], inplace=True)
    annotations = pd.merge(annotations, cazy_df[['gene', 'cazy']], on='gene', how='left')
    annotations = pd.merge(annotations, amr_df[['gene','resistance_type','resistance_target']], on='gene', how='left')
    annotations = pd.merge(annotations, vfdb_df[['gene', 'vf', 'vf_type']], on='gene', how='left')
    annotations = pd.merge(annotations, signalp_df[['gene', 'signalp']], on='gene', how='left')

    # Parse DefenseFinder (optional)
    if defense_file:
        defense_df = pd.read_csv(defense_file, sep='\t')
        defense_df = defense_df.rename(columns={'hit_id': 'gene'})
        defense_df['activity'] = defense_df['activity'].fillna('')
        defense_df['gene_name'] = defense_df['gene_name'].fillna('')
        defense_df['type'] = defense_df['type'].fillna('')

        defense_hits = defense_df[defense_df['activity'] == 'Defense'][
            ['gene', 'gene_name', 'type']
        ].rename(columns={'gene_name': 'defense', 'type': 'defense_type'})
        antidefense_hits = defense_df[defense_df['activity'] == 'Antidefense'][
            ['gene', 'gene_name', 'type']
        ].rename(columns={'gene_name': 'antidefense', 'type': 'antidefense_type'})

        defense_hits = defense_hits.groupby('gene', group_keys=False).head(1).reset_index(drop=True)
        antidefense_hits = antidefense_hits.groupby('gene', group_keys=False).head(1).reset_index(drop=True)

        annotations = pd.merge(annotations, defense_hits, on='gene', how='left')
        annotations = pd.merge(annotations, antidefense_hits, on='gene', how='left')
    else:
        annotations['defense'] = pd.NA
        annotations['defense_type'] = pd.NA
        annotations['antidefense'] = pd.NA
        annotations['antidefense_type'] = pd.NA

    # Output the final DataFrame to the output file
    annotations.to_csv(output_file, sep='\t', index=False)

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Compare GFF and PFAM files and output the results.')
    parser.add_argument('-gff', required=True, type=str, help='Path to the GFF file')
    parser.add_argument('-kegg', required=True, type=str, help='Path to the kegg file')
    parser.add_argument('-keggdb', required=True, type=str, help='Path to the keggdb file')
    parser.add_argument('-pfam', required=True, type=str, help='Path to the PFAM file')
    parser.add_argument('-ec', required=True, type=str, help='Path to the EC file')
    parser.add_argument('-cazy', required=True, type=str, help='Path to the CAZY file')
    parser.add_argument('-vf', required=True, type=str, help='Path to the VF file')
    parser.add_argument('-vfdb', required=True, type=str, help='Path to the VFDB file')
    parser.add_argument('-amr', required=True, type=str, help='Path to the AMR file')
    parser.add_argument('-amrdb', required=True, type=str, help='Path to the AMRDB file')
    parser.add_argument('-signalp', required=True, type=str, help='Path to the SIGNALP file')
    parser.add_argument('-o', required=True, type=str, help='Path to the OUTPUT file')
    parser.add_argument('-defense', required=False, type=str, help='Path to DefenseFinder genes TSV (optional)')

    # Parse the arguments
    args = parser.parse_args()

    # Process the files
    merge_annotations(args.gff, args.kegg, args.keggdb, args.pfam, args.ec, args.cazy, args.vf, args.vfdb, args.amr, args.amrdb, args.signalp, args.o, args.defense)

if __name__ == '__main__':
    main()
