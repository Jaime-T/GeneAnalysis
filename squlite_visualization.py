#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wednesday 18 Sep

@author: Pablo Acera Mateos. Jaime Taitz added GUI functionality

Description: Added GUI with Plotly Dash

Notes: Change the SQLite connection depending on where you store the dataset 
"""

import pandas as pd
import requests
import sys
import re
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
from dash import Dash, html, dcc, Input, Output


# get gene info, use either ensembl id or gene symbol at input
def get_ensembl_gene_info(goi='KRAS'):
    if goi[:4].lower() == 'ensg':
        url = f'http://rest.ensembl.org/lookup/id/{goi}'
        r = requests.get(url, headers={ "Content-Type" : "application/json"})
    else:
        url = f'http://rest.ensembl.org/lookup/symbol/homo_sapiens/{goi}'
        r = requests.get(url, headers={ "Content-Type" : "application/json"})
        
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    output = r.json()
    
    return output


# get transcript's intron positions
def get_enst_intron_df(enst='ENST00000311936'):
    
    url = f'http://rest.ensembl.org/lookup/id/{enst}?expand=1'
    r = requests.get(url, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    enst_annot = r.json()
    
    intron_start_list = [k['end']+1 for k in enst_annot['Exon']]
    intron_end_list = [k['start']-1 for k in enst_annot['Exon']]
    start_tmp = sorted(intron_start_list)[:-1][::enst_annot['strand']]
    end_tmp = sorted(intron_end_list)[1:][::enst_annot['strand']]

    output = pd.DataFrame({'seqnames':enst_annot['seq_region_name'], 
                           'start':start_tmp, 
                           'end':end_tmp, 
                           'strand':{1:'+',-1:'-'}[enst_annot['strand']],
                           'transcript_name':enst_annot['display_name'],
                           'ensg':enst_annot['Parent'],
                           'enst':enst_annot['id'],
                           'object_type':'intron',
                           'annotation':list(range(1,len(start_tmp)+1))})
    
    return output


## get gene stats and wildtype sequence (from ensembl)
def ensembl_trxseq(enst = 'ENST00000311936'):
    if enst is not None:
        # CDNA = requests.get('http://rest.ensembl.org/sequence/id/'+ensg+'?', headers={ "Content-Type" : "application/json"})
        CDNA = requests.get('http://rest.ensembl.org/sequence/id/'+enst+'?mask_feature=1', headers={ "Content-Type" : "application/json"})
        data = CDNA.json()

        trx_stat = dict()
        _, _, chromosome, start, end, strand = data['desc'].split(':')
        trx_stat['chr'] = f'chr{chromosome}'
        trx_stat['start'] = int(start)
        trx_stat['end'] = int(end)
        trx_stat['strand'] = int(strand)
        trx_stat['length'] = len(data['seq'])

        WT_sequence = data['seq']

        return trx_stat, WT_sequence
    else:
        return None, None

def merge_exons(exon_list):
    merged = []
    sorted_exons = sorted(exon_list, key=lambda x: x.start)
    current_exon = sorted_exons[0]
    for exon in sorted_exons[1:]:
        if exon.start <= current_exon.end:
            # Merge the current exon with the next overlapping exon
            current_exon.end = max(current_exon.end, exon.end)
        else:
            merged.append(current_exon)
            current_exon = exon
    merged.append(current_exon)
    return merged

def merge_putative_exons(dct):
    # Convert keys to intervals and sort
    intervals = [(int(k.split('_')[0]), int(k.split('_')[1])) for k in dct.keys()]
    intervals.sort(key=lambda x: x[0])
    
    # Merge overlapping intervals and sum their values
    merged_intervals = []
    current_interval = intervals[0]
    current_value = dct['_'.join(map(str, current_interval))]
    
    for i in range(1, len(intervals)):
        # If the current interval overlaps with the next one
        if current_interval[1] >= intervals[i][0]:
            current_value += dct['_'.join(map(str, intervals[i]))]
            current_interval = (current_interval[0], max(current_interval[1], intervals[i][1]))
        else:
            merged_intervals.append((current_interval, current_value))
            current_interval = intervals[i]
            current_value = dct['_'.join(map(str, current_interval))]
    merged_intervals.append((current_interval, current_value))
    
    # Convert merged intervals back to dictionary format
    merged_dct = {'_'.join(map(str, interval)): value for interval, value in merged_intervals}
    
    return merged_dct


#######################################################################
# process the acceptor data and get read distribution of all donors ###
#######################################################################

def process_drug_acceptor_data_raw(temp_acceptor_custom):
    def extracted_info(sample):
        sample_id, reads = sample.split(':')
        return str(sample_id), int(reads)

    # Get all the donors and the sample:read information
    donor_df = temp_acceptor_custom[['start', 'samples']]
    donor_df = donor_df.copy()
    donor_df['start'] = donor_df['start'].astype(str)
    donor_df['samples'] = donor_df['samples'].str.lstrip(',')
    split_samples = donor_df['samples'].str.split(',')
    extract_info = split_samples.apply(lambda x: [extracted_info(sample) for sample in x])

    # Convert the list of tuples into separate lists for sample IDs and read counts
    sample_ids, read_counts = zip(*[item for sublist in extract_info for item in sublist])

    # Create a set of unique sample IDs
    unique_sample_ids = set(sample_ids)

    # Initialize a dictionary to store read counts per sample per start
    reads_per_start_per_sample = {}

    # Iterate over each row in the DataFrame
    for index, row in donor_df.iterrows():
        start = row['start']
        samples = row['samples'].split(',')

        # Initialize a dictionary to store read counts per sample for this start
        reads_per_sample = {sample_id: 0 for sample_id in unique_sample_ids}  # Initialize with 0 instead of 1

        for sample in samples:
            sample_id, reads = extracted_info(sample)
            reads_per_sample[sample_id] = reads  # Store raw read counts

        reads_per_start_per_sample[start] = reads_per_sample

    reads_df = pd.DataFrame(reads_per_start_per_sample).fillna(0)
    reads_df.reset_index(inplace=True)
    reads_df.rename(columns={'index': 'Samples'}, inplace=True)

    return reads_df

def plot_distributions_data(data):
    # Determine the number of columns in the data
    num_columns = data.shape[1]
    
    # Create a figure with subplots
    fig, axes = plt.subplots(num_columns, 1, figsize=(8, num_columns * 4),dpi=200)
    
    # Flatten the axes array if only one column (i.e., if num_columns is 1)
    if num_columns == 1:
        axes = [axes]
        
    # Loop through each column and plot
    for i in range(num_columns):
        sns.histplot(data[:, i], kde=True, ax=axes[i])
        axes[i].set_title(f'Distribution of Column {i+1}')
        axes[i].set_xlabel('Value')
        axes[i].set_ylabel('Frequency')
    
    # Adjust the layout
    plt.tight_layout()
    plt.show()


def process_row(row):
    items = row.split(',')
    result = {}
    for item in items:
        if ':' in item:  # Check to make sure there is a colon to split on
            key, value = item.split(':')
            result[key] = int(value)
    return result


def get_gene_acceptor_data(selected_gene = 'SMN2'):
    # Grab custom data
    conn = sqlite3.connect('/Users/jaimetaitz/cci_internship/GeneAnalysis/junctions.sqlite')
    #conn = sqlite3.connect('/Users/paceramateos/projects/cryptic_exons/data/U1_brana_risdi/jc_custom_STARjunc.sqlite')

    # Create a cursor object
    cur = conn.cursor()

    ### Input of the app
    #selected_gene = 'SMN2'

    geneinfo = get_ensembl_gene_info(selected_gene)

    # get intron positions of canonical transcript
    enst = re.sub('[.].*$','',geneinfo['canonical_transcript'])
    goi_introns = get_enst_intron_df(enst)

    # get canonical transcript sequence
    trx_stat, seq = ensembl_trxseq(enst)

    # get snaptron data
    strand = {-1:'-', 1:'+'}[geneinfo['strand']]

    # Define your region and strand
    chromosome = trx_stat['chr']
    start_pos = trx_stat['start']
    end_pos = trx_stat['end']
    query = f"""
    SELECT snaptron_id as 'DataSource:Type', snaptron_id, chrom as 'chromosome', start, end,
    length, strand, annotated, donor as 'left_motif', acceptor as 'right_motif',
    left_annotated, right_annotated, samples, samples_count,
    coverage_sum, coverage_avg, coverage_median, source_dataset_id
    FROM intron 
    WHERE chrom = '{chromosome}' 
    AND strand = '{strand}' 
    AND start <= {end_pos} 
    AND end >= {start_pos};
    """
    # Execute the query and load the results into a DataFrame
    snap_custom = pd.read_sql_query(query, conn)

    # Close the connection
    conn.close()

    # only consider introns within the limit of the canonical transcript
    snap_custom = snap_custom.loc[(snap_custom['start'] >= goi_introns['start'].min()) & (snap_custom['end'] <= goi_introns['end'].max())]

    snap_custom = snap_custom.join(goi_introns[['start','end','transcript_name', 'annotation']].set_index(['start','end']), how='left', on=['start','end'])
    snap_custom = snap_custom.fillna('non-canonical')

    # sort snap0
    snap_custom = snap_custom.sort_values('start').reset_index(drop=True)
    snap_custom['genesymbol'] = geneinfo['display_name']
    snap_custom['ensembl'] = geneinfo['id']

    exon_seq = 4 # how many exonic nt to keep in reported motif
    intron_seq = 7 # how many intonic nt to keep in reported motif 

    #############################################
    # generate donor and acceptor data frames ###
    #############################################
    if geneinfo['strand'] == 1:
        acceptor_side = 'end'
        donor_side = 'start'
        left_seq = exon_seq
        right_seq = intron_seq
    else:
        acceptor_side = 'start'
        donor_side = 'end'
        left_seq = intron_seq-1
        right_seq = exon_seq+1

    canonical_acceptor = snap_custom.loc[snap_custom.transcript_name != 'non-canonical', acceptor_side].tolist()
    return snap_custom, canonical_acceptor, acceptor_side, donor_side
    #canonical_donor    = snap_custom.loc[snap_custom.transcript_name != 'non-canonical', donor_side].tolist()
    #canonical_lengths  = snap_custom.loc[snap_custom.transcript_name != 'non-canonical'].length.tolist()

#########################################################################
#########################################################################

### controls
### HTT k=47
### SMN2 k=7

## Start App
app = Dash(__name__)

# Generate dropdown options dynamically
def get_all_ensembl_genes():
    with open("./mart_gene_names.txt", "r") as file:
        lines = file.readlines()
        gene_names = []
        for x in lines[1:]:
            cols = x.split('\t')
            if len(cols) == 1:
                continue
            
            name = cols[1]
            if len(name) != 0 and name not in gene_names:
                name = name.rstrip('\n')
                gene_names.append(name)
            
    return gene_names

def get_acceptors_for_gene(gene_name):
    snap_custom, canonical_acceptors, acceptor_side, donor_side = get_gene_acceptor_data(gene_name)
    return canonical_acceptors

gene_options = [{"label": gene, "value": gene} for gene in get_all_ensembl_genes()]

# App layout
app.layout = html.Div([
    html.H1("Proportion of Counts in an Acceptor", style={'text-align': 'center'}),

    dcc.Dropdown(id="slct_gene",
                 options=gene_options,
                 multi=False,
                 value='SMN2', # default value
                 style={'width': "40%"}
                 ),

    dcc.Dropdown(id="slct_acceptor",
                 options=[],
                 multi=False,
                 value=None,
                 style={'width': "40%"}
                 ),
    
    html.Div(id='output_container', children=[]),

    dcc.Graph(id='my_acceptor_map', figure={})
])

@app.callback(
    Output(component_id='slct_acceptor', component_property='options'),
    [Input(component_id='slct_gene', component_property='value')])


def update_acceptor_options(selected_gene):
    canonical_acceptors = get_acceptors_for_gene(selected_gene)  # Get acceptors for the selected gene
    return [{"label": acceptor, "value": acceptor} for acceptor in canonical_acceptors] if canonical_acceptors else []


@app.callback(
    Output(component_id='slct_acceptor', component_property='value'),
    Input(component_id='slct_acceptor', component_property='options'))

def set_acceptor_value(available_options):
    return available_options[0]['value'] if available_options else None


# Connect the Plotly graphs with Dash Components
@app.callback(
    [Output(component_id='output_container', component_property='children'),
     Output(component_id='my_acceptor_map', component_property='figure')],
    [Input(component_id='slct_acceptor', component_property='value'),
     Input(component_id='slct_gene', component_property='value')]
)

def update_graph(slct_acceptor, slct_gene):
    print(f"Acceptor: {slct_acceptor}, Gene: {slct_gene}")

    container = f"The acceptor and gene chosen by user was: {slct_acceptor}, {slct_gene}"

    ''' 
    columns
    Index(['DataSource:Type', 'snaptron_id', 'chromosome', 'start', 'end',
            'length', 'strand', 'annotated', 'left_motif', 'right_motif',
            'left_annotated', 'right_annotated', 'samples', 'samples_count',
            'coverage_sum', 'coverage_avg', 'coverage_median', 'source_dataset_id',
            'transcript_name', 'annotation', 'genesymbol', 'ensembl'],
          dtype='object')
    '''

    # Filter based on user input 
    snap_custom, canonical_acceptor, acceptor_side, donor_side = get_gene_acceptor_data(slct_gene)
    temp_acceptor_custom = snap_custom[snap_custom[acceptor_side] == slct_acceptor] # changed this based on user input 
    
    temp_acceptor_custom_annotation = temp_acceptor_custom[[donor_side,'annotation']]
    temp_acceptor_custom_annotation.loc[:,'annotation'] = np.where(
    temp_acceptor_custom_annotation['annotation'] != 'non-canonical',
    'canonical',  # This will be set if condition is True
    'non-canonical'  # This will be set if condition is False
    )
    
    if temp_acceptor_custom.shape[0] < 2:
        #continue
        container = f"The acceptor and gene chosen by user was: {slct_acceptor}, {slct_gene}.\n Not enough data for this acceptor"
        sys.stdout.write('Not enough data')
        return container, None
    
    if temp_acceptor_custom.strand.tolist()[0] == '+':
        acceptor_name = temp_acceptor_custom.end.iloc[0]
    else:
        acceptor_name = temp_acceptor_custom.start.iloc[0]
    
    custom_data_raw = process_drug_acceptor_data_raw(temp_acceptor_custom)
    
    # Make a df with the proportion numbers
    numeric_columns = custom_data_raw.columns.drop('Samples')
    row_sums = custom_data_raw[numeric_columns].sum(axis=1)
    prop_data = custom_data_raw[numeric_columns].div(row_sums, axis=0)
    prop_data.insert(0, 'Samples', custom_data_raw['Samples'])
    
    ### discard columns where for all samples values are smaller than 1%
    valid_columns = prop_data.drop('Samples', axis=1).apply(lambda x: (x >= 0.01).any())

    # Filter the DataFrame to include only the valid columns plus 'Samples'
    prop_data = prop_data.loc[:, ['Samples'] + valid_columns[valid_columns].index.tolist()]    
    
    # Visualize data as a heatmap
    annotation_map = temp_acceptor_custom_annotation.set_index('start')['annotation'].to_dict()

    # Drop the 'Samples' column for the heatmap data but store it for use as y-axis labels
    heatmap_data = prop_data.set_index('Samples')

    custom_xtick_labels = [f"{col}\n({annotation_map.get(int(col), 'N/A')})" for col in heatmap_data.columns]
    
    # Plotly Express heatmap
    fig = px.imshow(heatmap_data, 
                    labels=dict(x="Acceptor sites", y="Sample ID", colour="Proportion"),
                    x=custom_xtick_labels,
                    y=heatmap_data.index,
                    text_auto=".2f", 
                    aspect="auto",
                    color_continuous_scale='Viridis')
    
    return container, fig


if __name__ == '__main__':
    app.run(debug=True) 



''' 
for k in range(len(canonical_acceptor)):
    
    ###### Make a boxplot (more convenient plot in case of lots of samples)
    
    boxplot_data = prop_data.drop('Samples', axis=1)
    custom_xtick_labels = [f"{col}\n({annotation_map.get(int(col), 'N/A')})" for col in boxplot_data.columns]
    
    # Set the style and increase font sizes using seaborn
    sns.set(style="whitegrid", font_scale=1.5)  # Adjust font_scale for bigger text
    
    # Create the plot
    plt.figure(figsize=(12, 9), dpi=200)  # Enhanced figure size and dpi for clarity
    sns.boxplot(data=boxplot_data)
    sns.stripplot(data=boxplot_data, color='grey', size=9, jitter=True, edgecolor='black', linewidth=1.5)
    
    # Setting the plot title and labels with increased font size directly
    plt.title(f"Proportion of counts in acceptor {canonical_acceptor[k]} {selected_gene}",
              fontsize=35,
              pad=20)
    plt.ylabel('Proportion', fontsize=30)
    plt.xlabel('Donor Sites', fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    
    # Setting custom x-tick labels with larger text
    plt.xticks(ticks=range(len(custom_xtick_labels)), labels=custom_xtick_labels, rotation=45)
    
    plt.show()
    
'''

    


    
    
        
        
        
   