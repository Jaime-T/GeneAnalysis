#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wednesday 18 Sep

@author: Analysis by Pablo Acera Mateos. GUI functionality by Jaime Taitz  

Description: Added GUI with Plotly Dash

Notes: Change the SQLite connection depending on where you store the dataset: 
        get_gene_acceptor_data connects to the database of choice
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
import dash_bio as dashbio
from scipy import special
import json
import math

description_data = """
 SRR15622469 - DMSO control_1
 SRR15622470 - DMSO control_2
 SRR15622463 - DMSO control_3
 SRR15622464 - Treated with 2 nM branaplam_1
 SRR15622465 - Treated with 2 nM branaplam_2 
 SRR15622466 - Treated with 2 nM branaplam_3
 SRR15622468 - Treated with 40 nM branaplam_1
 SRR15622467 - Treated with 40 nM branaplam_2
 SRR15622456 - Treated with 40 nM branaplam_3
 SRR15622457 - Treated with 50 nM risdiplam_1
 SRR15622458 - Treated with 50 nM risdiplam_2
 SRR15622459 - Treated with 50 nM risdiplam_3
 SRR15622461 - Treated with 1000 nM risdiplam_1
 SRR15622460 - Treated with 1000 nM risdiplam_2
 SRR15622462 - Treated with 1000 nM risdiplam_3
 SRR8697000 - HEK293T_U1Mut_1
 SRR8697001 - HEK293T_U1Mut_2
 SRR8697002 - HEK293T_U1WT_1
 SRR8697003 - HEK293T_U1WT_2
 SRR5206789 - Treated_SMA_404
 SRR5206788 - Treated_SMA_403
 SRR5206787 - Treated_SMA_402
 SRR5206786 - Treated_SMA_401
 SRR5206785 - Control_SMA_02
 SRR5206784 - Control_SMA_01
 SRR5206783 - Control_FB13c2
 SRR5206782 - Control_FB13c1
 SRR9674472 - MEC1_RNU1_3_WT_B
 SRR9674471 - MEC1_RNU1_3_AC_B
 SRR9674470 - MEC1_RNU1_3_WT_A
 SRR9674469 - MEC1_RNU1_3_AC_A
 SRR9674468 - HG3_RNU1_3_WT_B
 SRR9674467 - HG3_RNU1_3_AC_B
 SRR9674466 - HG3_RNU1_3_WT_A
 SRR9674465 - HG3_RNU1_3_AC_A
 SRR9674464 - JVM3_RNU1_3_WT_B
 SRR9674463 - JVM3_RNU1_3_AC_B
 SRR9674462 - JVM3_RNU1_3_WT_A
 SRR9674461 - JVM3_RNU1_3_AC_A
 SRR10485752 - hela_overexpression_1.5ug
 SRR10485751 - hela_overexpression_1ug
 SRR10485750 - hela_overexpression_control
 SRR10485749 - U1AMO_62.5pmol
 SRR10485748 - U1AMO_12.5pmol
 SRR10485747 - cAMO_for_U1AMO
 """

experiments = [
    {
        "control": ["SRR15622469", "SRR15622470", "SRR15622463"],
        "treatment": ["SRR15622464", "SRR15622465", "SRR15622466"],
        "label": "DMSO vs 2 nM branaplam"
    },
    {
        "control": ["SRR15622469", "SRR15622470", "SRR15622463"],
        "treatment": ["SRR15622457", "SRR15622458", "SRR15622459"],
        "label": "DMSO vs 50 nM risdiplam"
    },
    {
        "control": ["SRR15622469", "SRR15622470", "SRR15622463"],
        "treatment": ["SRR15622468", "SRR15622467", "SRR15622456"],
        "label": "DMSO vs 40 nM branaplam"
    },
    {
        "control": ["SRR15622469", "SRR15622470", "SRR15622463"],
        "treatment": ["SRR15622461", "SRR15622460", "SRR15622462"],
        "label": "DMSO vs 1000 nM risdiplam"
    },
    {
        "control": ["SRR8697000", "SRR8697001"],
        "treatment": ["SRR8697002", "SRR8697003"],
        "label": "HEK293T_U1Mut vs HEK293T_U1WT"
    },
    {
        "control": ["SRR5206783", "SRR5206782"],
        "treatment": ["SRR5206789", "SRR5206788", "SRR5206787", "SRR5206786"],
        "label": "Treated_SMA vs Control_FB13"
    },
    {
        "control": ["SRR9674472", "SRR9674470"],
        "treatment": ["SRR9674471", "SRR9674469"],
        "label": "MEC1_RNU1_3_WT vs MEC1_RNU1_3_AC"
    },
    {
        "control": ["SRR9674468", "SRR9674466"],
        "treatment": ["SRR9674467", "SRR9674465"],
        "label": "HG3_RNU1_3_WT vs HG3_RNU1_3_AC"
    },
    {
        "control": ["SRR9674464", "SRR9674462"],
        "treatment": ["SRR9674463", "SRR9674461"],
        "label": "JVM3_RNU1_3_WT vs JVM3_RNU1_3_AC"
    },
    {
        "control": ["SRR10485750"],
        "treatment": ["SRR10485752"],
        "label": "hela_overexpression_1.5ug vs Control"
    },
    {
        "control": ["SRR10485750"],
        "treatment": ["SRR10485751"],
        "label": "hela_overexpression_1ug vs Control"
    }
]

lines = description_data.strip().split('\n')
cleaned_lines = [item.strip() for item in lines]
# Split each line into two parts and create a list of tuples
data_tuples = [tuple(line.split(' - ')) for line in cleaned_lines]
# Create a DataFrame from the list of tuples
sample_condition = dict(data_tuples)
#print(sample_condition)

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
    conn = sqlite3.connect('/Users/jaimetaitz/cci_internship/GeneAnalysis/jc_custom_STARjunc.sqlite')
    #conn = sqlite3.connect('/Users/paceramateos/projects/cryptic_exons/data/U1_brana_risdi/jc_custom_STARjunc.sqlite')

    # Create a cursor object
    cur = conn.cursor()

    ### Input of the app
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

## Needle Values 
def calc_lnlr(r1, r2, alpha = None):

    if(r1.shape[0] != r2.shape[0]):
        raise Exception("Different K's implied by r1 and r2.")
        
    if(np.sum(alpha == None) == 0):
        if(alpha.shape[0] != r1.shape[0]):
            raise Exception("Different K's implied by alpha and r1.")
    else:
        alpha = np.ones(r1.shape[0])
    
    # calculate ln(LR) using counts and prior parameters
    lnlr = (np.sum(special.loggamma(alpha + r1 + r2)) - special.loggamma(np.sum(alpha + r1 + r2))
           - np.sum(special.loggamma(alpha + r1)) + special.loggamma(np.sum(alpha + r1))
           - np.sum(special.loggamma(alpha + r2)) + special.loggamma(np.sum(alpha + r2))
           + np.sum(special.loggamma(alpha)) - special.loggamma(np.sum(alpha))) 
    
    return lnlr


def calc_LR(alpha_estimates, raw_sample, prior_count=1):
    
    alpha_estimates_sample = raw_sample #.squeeze() 

    LR = np.exp(calc_lnlr(alpha = np.full((alpha_estimates.shape[0],), prior_count),
                          r1=alpha_estimates, # alphas calculated using MLE with ratios from GTEx
                          r2=alpha_estimates_sample))# alphas calculated using raw counts+1, to test against the background  
    return LR


def get_needle_value(mean, indiv_sample):
    bayes_factor = calc_LR(mean,indiv_sample,7)
    needle_value = 1/bayes_factor
    return needle_value

# Raw data for each acceptor of a gene
def raw_data(snap_custom,  acceptor_side, donor_side, slct_acceptor):
    temp_acceptor_custom = snap_custom[snap_custom[acceptor_side] == slct_acceptor] # changed this based on user input 
    
    if temp_acceptor_custom.shape[0] < 2:
        # not enough data
        return {} 
    
    custom_data_raw = process_drug_acceptor_data_raw(temp_acceptor_custom)
    raw_data = pd.DataFrame(custom_data_raw)
    return raw_data

def filtered_raw_data(raw_data):
    # Make a df with the proportion numbers
    numeric_columns = raw_data.columns.drop('Samples')
    row_sums = raw_data[numeric_columns].sum(axis=1)
    prop_data = raw_data[numeric_columns].div(row_sums, axis=0)
    prop_data.insert(0, 'Samples', raw_data['Samples'])
    
    ### discard columns where for all samples values are smaller than 1%
    valid_columns = prop_data.drop('Samples', axis=1).apply(lambda x: (x >= 0.01).any())

    # Filter the raw data to include only the valid columns plus 'Samples'
    filtered_data = raw_data.loc[:, ['Samples'] + valid_columns[valid_columns].index.tolist()]

    # Drop the 'Samples' column but store it for use as index
    filtered_data.set_index('Samples', inplace=True)

    # less than 2 valid acceptor sites
    if filtered_data.shape[1] < 2:
        return {}
    return filtered_data 


# HeatMap data for gui 
def heatmap(slct_gene, slct_acceptor):
    # Filter based on user input 
    if slct_acceptor is None:
            container = "Choose an acceptor. Click a point on the needleplot, or select from dropdown options"
            return {}, {}, container
    else:
        container = f"The acceptor and gene chosen by user was: {slct_acceptor}, {slct_gene}"
    
    snap_custom, canonical_acceptor, acceptor_side, donor_side = get_gene_acceptor_data(slct_gene)
    temp_acceptor_custom = snap_custom[snap_custom[acceptor_side] == slct_acceptor] # changed this based on user input 
    
    if temp_acceptor_custom.shape[0] < 2:
        container = f"The acceptor and gene chosen by user was: {slct_acceptor}, {slct_gene}.\n Not enough data for this acceptor"
        return {}, {}, container

    temp_acceptor_custom_annotation = temp_acceptor_custom[[donor_side,'annotation']]
    temp_acceptor_custom_annotation.loc[:,'annotation'] = np.where(
    temp_acceptor_custom_annotation['annotation'] != 'non-canonical',
    'canonical',  # This will be set if condition is True
    'non-canonical'  # This will be set if condition is False
    )
    
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

    # less than 2 valid acceptor sites
    if heatmap_data.shape[1] < 2:
        container = container = f"The acceptor and gene chosen by user was: {slct_acceptor}, {slct_gene}.\n Not enough donors for this acceptor"
        return {}, {}, container

    # change sampleid to sample name of condition, using the sample_condition dict 
    # case: if id doesnt exist in dict, then keep its SSR value
    new_index = [sample_condition.get(i,i) for i in heatmap_data.index] 
    heatmap_data.index = new_index
    
    # sort in alphabetical order
    heatmap_data = heatmap_data.sort_index()

    custom_xtick_labels = [f"{col}\n({annotation_map.get(int(col), 'N/A')})" for col in heatmap_data.columns]

    # Plotly Express heatmap
    fig = px.imshow(heatmap_data, 
                    labels=dict(x="Acceptor sites", y="Sample ID", color="Proportion"),
                    x=custom_xtick_labels,
                    y=heatmap_data.index,
                    text_auto=".2f", 
                    aspect="auto",
                    color_continuous_scale='Viridis')
    fig.update_layout(
        height=800,
        yaxis=dict(
            automargin=True,
            title='Sample ID',  
            tickmode='linear'  # Ensures ticks are shown at regular intervals
        ),
        margin=dict(l=60, r=20, t=40, b=40)  
    )
                    
    return heatmap_data, fig, container


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

def get_acceptor_coords(gene_name, acceptor):
    snap_custom, canonical_acceptors, acceptor_side, donor_side = get_gene_acceptor_data(gene_name)
    temp_acceptor_custom = snap_custom[snap_custom[acceptor_side] == acceptor]
    start = temp_acceptor_custom['start'].values
    end = temp_acceptor_custom['end'].values
    return start,end

gene_options = [{"label": gene, "value": gene} for gene in get_all_ensembl_genes()]

# Needle sample data:
needle_sample_data = {
    'x': ['271.0-279.0', '808.0-825.0', '661.0-672.0'], 
    'y': ['1', '1', '1'], 
    'mutationGroups': ['Helix', 'Helix', 'Beta strand'], 
    'domains': [
        {'name': 'Intron_p85B', 'coord': '32-107'},
        {'name': 'Intron_rbd', 'coord': '173-292'},
        {'name': 'Intron_C2', 'coord': '350-485'}
    ]  
}

## Start App
app = Dash(__name__, suppress_callback_exceptions=True, external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css'])

# App layout
app.layout = html.Div([
    html.H1("Proportion of Counts in an Acceptor", style={'text-align': 'center'}),

    # Gene View or Experiment View, choose a tab 
    dcc.Tabs(id='gene_or_exper_tab', value='gene_tab', children=[
        dcc.Tab(label='Gene View', value='gene_tab'),
        dcc.Tab(label='Experiment View', value='exper_tab'),
    ]),
    html.Div(id='gene_or_exper_tab_content'),



])

@app.callback(
    Output(component_id='gene_or_exper_tab_content', component_property='children'),
    Input(component_id='gene_or_exper_tab', component_property='value')
)

def render_content(tab):
    if tab == 'gene_tab': # Gene View 
        return html.Div([
            html.H3('Gene View!'),
            html.P("Select Gene:"),

            dcc.Dropdown(id="slct_gene",
                        options=gene_options,
                        multi=False,
                        value='H3-3A', # default value
                        style={'width': "40%", 'display':'inline-block'}
                        ),

            html.P("Show or hide range slider:"),

            dcc.Dropdown(
                id='needleplot-rangeslider',
                options=[
                    {'label': 'Show', 'value': 1},
                    {'label': 'Hide', 'value': 0}
                ],
                clearable=False,
                multi=False,
                value=1,
                style={'width': "40%"}
            ),

            html.H1("NeedlePlot"),
            
            html.P("Maximum needle value:"),
            dcc.Dropdown(
                id='needleplot_max',
                options=[{"label": num, "value": num} for num in range(100,1100,100)],
                value=200,
                multi=False,
                style={'width': "40%"}
            ),

            html.P("Needle value options:"),
            dcc.Dropdown(
                id='needleplot_option',
                options=[{"label": "Actual", "value": "actual"}, {"label": "Log", "value": "log"}],
                value="actual",
                multi=False,
                style={'width': "40%"}
            ),

            # Centered NeedlePlot
            html.Div(
                style={
                    'display': 'flex',
                    'justify-content': 'center',
                    'align-items': 'center',
                    'margin': '20px',
                    'width': '90%'
                },
                children=[
                    dashbio.NeedlePlot(
                        id='gene-needleplot',
                        mutationData=needle_sample_data,
                        xlabel='Gene',
                        ylabel='Score',
                        width=1200,
                        domainStyle={'displayMinorDomains': True}
                    )
                ]
            ),
            html.P("Click on an intron of interest", style={'color': 'blue', 'fontSize': 15}),

            dcc.Dropdown(id="slct_acceptor",
                        options=[],
                        multi=False,
                        value='70067293', # default 
                        style={'width': "40%"}
                        ),
            
            html.Div(id='output_container', children=[], style={'color': 'blue'}),

            # Heat map
            html.H1("Heat Map"),
            dcc.Graph(id='my_acceptor_map', figure={})

        ])
    
    elif tab == 'exper_tab': # Experiment view
        return html.Div([
            html.H3('Experiment Tab'),
            html.P('Still being developed! But.. you can select an experiment of interest :)'),

            dcc.Dropdown(id="slct_exper",
                         options=[exp["label"] for exp in experiments],
                         multi=False,
                         value='DMSO vs 40 nM branaplam',
                         style={'width': "40%"}
                         )
        ])


@app.callback(
    Output(component_id='gene-needleplot', component_property='rangeSlider'),
    Input(component_id='needleplot-rangeslider', component_property='value')
)
def update_needleplot(show_rangeslider):
    return True if show_rangeslider else False

@app.callback(
    Output(component_id='gene-needleplot', component_property='mutationData'),
    [Input(component_id='slct_gene', component_property='value'), 
     Input(component_id='needleplot_max', component_property='value'),
     Input(component_id='needleplot_option', component_property='value')])

def update_needleplot(selected_gene, needle_max, needle_option):
    snap_custom, canonical_acceptors, acceptor_side, donor_side = get_gene_acceptor_data(selected_gene)
    if len(canonical_acceptors) == 0:
        return needle_sample_data
    
    domains = []
    x = []
    y = []

    # find max needle value for each acceptor
    acceptor_needles = {}

    for acceptor in canonical_acceptors:
        
        samples_raw_data = raw_data(snap_custom, acceptor_side, donor_side, acceptor)

        # if not a dataframe, means no info, so skip this acceptor 
        if not isinstance(samples_raw_data, pd.DataFrame):
            continue

        valid_raw_data = filtered_raw_data(samples_raw_data)
        
        # if not a dataframe, means no info, so skip this acceptor 
        if not isinstance(valid_raw_data, pd.DataFrame):
            continue

        if valid_raw_data.empty:
            continue
        
        acceptor_needles[str(acceptor)] = -math.inf

        for experiment in experiments:
            
            # process the control samples. Samples column index is SRR id 
            controls = experiment["control"]
            control_samples = [valid_raw_data[valid_raw_data.index.str.contains(control)] for control in controls]
            control_samples_df = pd.concat(control_samples)

            # if theres no samples for this experiment 
            if control_samples_df.empty:
                continue
            else:
                control_avg = control_samples_df.mean(axis=0)

            # process the treatment samples
            treatments = experiment["treatment"]
            treatment_samples = [valid_raw_data[valid_raw_data.index.str.contains(treatment)] for treatment in treatments]
            treatment_samples_df = pd.concat(treatment_samples)

            # if theres no samples for this experiment 
            if treatment_samples_df.empty:
                continue
            
            # Calculate all needles at once
            needles = treatment_samples_df.apply(lambda row: get_needle_value(control_avg, row), axis=1)
            print('all needles:', needles)

            # Get the maximum needle value
            max_needle = needles.max()

            # Update the acceptor_needles dictionary
            acceptor_needles[str(acceptor)] = max(max_needle, acceptor_needles[str(acceptor)])
            print('new max is:', acceptor_needles[str(acceptor)] )

            label = experiment["label"]
            if label == 'Treated_SMA vs Control_FB13' and acceptor == 226071350:
                print('gene is: ' + str(selected_gene) + 'acceptor is: ' + str(acceptor) + 'label: ' + label)
                print('treatment samples: ')
                print(treatment_samples_df)
                print('control sample: ')
                print(control_samples_df)
                print('needle values would be... ', needles)

                

        coord = str(acceptor) + '-' + str(acceptor + 100)
        domains.append({"name": str(acceptor), "coord": coord})
        

    if needle_option == 'log':
        acceptor_needles = {key: math.log(value) for key, value in acceptor_needles.items()}
    
    else:
        acceptor_needles = {key: min(needle_max,value) for key, value in acceptor_needles.items()}
    
    print(acceptor_needles)


    plot_data = {"x": list(acceptor_needles.keys()), 
                "y": list(acceptor_needles.values()), 
                "domains": domains}
    return plot_data


@app.callback(
    Output(component_id='slct_acceptor', component_property='value'),
    Input(component_id='gene-needleplot', component_property='clickData')
)

def update_acceptor_value(clickData):
    # Check if clickData is None
    if clickData is None:
        return "No selection"  # or any default value you want to return

    # Assuming clickData is structured correctly, call click_acceptor
    x_value = click_acceptor(clickData)
    return x_value

def click_acceptor(clickData):

    if clickData is None:
        sys.stdout.write('Error: clickData is None\n')
        return None
    
    if not isinstance(clickData, dict):
        sys.stdout.write('Error: clickData must be a dictionary\n')
        return None

    try:
        x = clickData['points'][0]['x']
    except (IndexError, KeyError) as e:
        sys.stdout.write(f'Error accessing x value: {e}\n')
        return None

    return int(x)


@app.callback(
    Output(component_id='slct_acceptor', component_property='options'),
    Input(component_id='slct_gene', component_property='value'))


def update_acceptor_options(selected_gene):
    canonical_acceptors = get_acceptors_for_gene(selected_gene)  # Get acceptors for the selected gene
    return [{"label": acceptor, "value": acceptor} for acceptor in canonical_acceptors] if canonical_acceptors else []

# Heatmap 
@app.callback(
    [Output(component_id='output_container', component_property='children'),
     Output(component_id='my_acceptor_map', component_property='figure')],
    [Input(component_id='slct_acceptor', component_property='value'),
     Input(component_id='slct_gene', component_property='value')]
)

def update_graph(slct_acceptor, slct_gene):
    print(f"Acceptor: {slct_acceptor}, Gene: {slct_gene}")
    heatmap_data, fig, container = heatmap(slct_gene, slct_acceptor)
    
    return container, fig


if __name__ == '__main__':
    app.run(debug=True) 



''' 
columns
Index(['DataSource:Type', 'snaptron_id', 'chromosome', 'start', 'end',
        'length', 'strand', 'annotated', 'left_motif', 'right_motif',
        'left_annotated', 'right_annotated', 'samples', 'samples_count',
        'coverage_sum', 'coverage_avg', 'coverage_median', 'source_dataset_id',
        'transcript_name', 'annotation', 'genesymbol', 'ensembl'],
        dtype='object')
'''

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

    


    
    
        
        
        
   