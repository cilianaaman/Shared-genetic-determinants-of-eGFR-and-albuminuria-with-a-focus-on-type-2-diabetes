
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: cilianaaman
"""
import openpyxl
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn import preprocessing
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from statsmodels.stats.multitest import multipletests



#%% Load eGFR Data
#########################
#### Load eGFR Data #####
#########################


# Load Excel file
df = pd.read_excel('your/path')

df=df.drop(labels="LocusCoordinates", axis=1)

# Create a matrix

eGFR_matrix_overall = df.drop(range(424,453))

del df

#%% UKBB Data import (run one time)
#####################
####  UKBB Data #####
#####################

#########################
### Load the alb data ###
#########################

df_Genomic_alb = pd.read_excel('your/path', sheet_name='Ark1')

# Select columns starting with "rs"
rs_columns = [col for col in df_Genomic_alb.columns if col.startswith('rs')]

# round the values 
df_Genomic_alb[rs_columns] = df_Genomic_alb[rs_columns].apply(np.round)


##########################
### Load the eGFR data ###
##########################


df_Genomic_eGFR = pd.read_excel('your/path', sheet_name='Ark1')

# Select columns starting with "rs"
rs_columns = [col for col in df_Genomic_eGFR.columns if col.startswith('rs')]

# round the values 
df_Genomic_eGFR[rs_columns] = df_Genomic_eGFR[rs_columns].apply(np.round)




###############################
### Load the phenotype data ###
###############################


df_phenotype = pd.read_excel('your/path', sheet_name='phenotype_data2', decimal=',')


import pickle
with open('variables.pkl', 'wb') as f:
    pickle.dump((df_phenotype, df_Genomic_eGFR, df_Genomic_alb), f)

#%% Getting the UKBB imported data fast + creating 2 matresis with the phenotype and genotype data
import pickle


# Load variables
with open('your/path/variables.pkl', 'rb') as f:
    df_phenotype, df_Genomic_eGFR, df_Genomic_alb = pickle.load(f)


# change pickle file 
df_Genomic_eGFR.drop(columns=['rs1004441_A.1'], inplace=True)


ethnicity = pd.read_excel('your/path/ethnicity.xlsx', sheet_name='ethnicity', decimal=',')

################
### Withdraw ###
################

# Load the participants who have withdrawn
withdraw = pd.read_csv('your/path/withdraw.csv', header=None)

# Rename the column to 'Eid'
withdraw.columns = ['Eid']

###############
### PC 1-10 ###
###############
pc = pd.read_csv(
    'your/path/GPC_1_to_10_participant.csv',
    skiprows=1,
    header=None,
    low_memory=False
)

pc.columns = ['eid'] + [f'PC{i}' for i in range(1, 11)]



### BMI and HbA1c
bmi_hba1c = pd.read_csv(
    'your/path/BMI_HbA1c_participant.csv',
    skiprows=1,
    header=None,
    low_memory=False
)

bmi_hba1c.columns = ['eid'] + ['hba1c'] + ['BMI'] 

################################################################
### creating 2 matresis with the phenotype and genotype data ###
################################################################


#####################
### the eGFR data ###
#####################

df_eGFR_data = pd.merge(df_phenotype, df_Genomic_eGFR, left_on='eid', right_on='IID')

# Specify the columns to delete
columns_to_delete = ['BMI', 'micro albumin in urin', 'Unnamed: 9', 'Unnamed: 10', 'Unnamed: 11', 'Unnamed: 12', 'Unnamed: 13','Unnamed: 0','SEX', 'IID']

# Drop the specified columns from the DataFrame
df_eGFR_data = df_eGFR_data.drop(columns_to_delete, axis=1)

# Filter out rows where 'Eid' in Overall_matrix matches 'Eid' in withdraw
df_eGFR_data = df_eGFR_data[~df_eGFR_data['eid'].isin(withdraw['Eid'])]

# Remove rows with NaN values
df_eGFR_data = df_eGFR_data.dropna()



######################
### overall matrix ###
######################

genomic_eGFR_columns = [col for col in df_Genomic_eGFR.columns if col.startswith('rs')]
rename_mapping = {col: col + '_eGFR' for col in genomic_eGFR_columns}
df_Genomic_eGFR = df_Genomic_eGFR.rename(columns=rename_mapping)


Overall_matrix = pd.merge(df_phenotype, df_Genomic_eGFR, left_on='eid', right_on='IID')
# Specify the columns to delete
columns_to_delete = ['Unnamed: 0', 'BMI', 'eGFR_Creatinine_2021_0','Unnamed: 9' ,'Unnamed: 10', 'Unnamed: 11', 'Unnamed: 12', 'Unnamed: 13']

# Drop the specified columns from the DataFrame
Overall_matrix = Overall_matrix.drop(columns_to_delete, axis=1)

Overall_matrix = pd.merge(Overall_matrix, ethnicity, left_on='eid', right_on='eid')

# Filter out rows where 'Eid' in Overall_matrix matches 'Eid' in withdraw
Overall_matrix = Overall_matrix[~Overall_matrix['eid'].isin(withdraw['Eid'])]


# Deleting variables that are not needed     
variables_to_delete = ['columns_to_delete','df_Genomic_eGFR','f','genomic_eGFR_columns','rename_mapping']


for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete
del ethnicity
del df_phenotype
del df_Genomic_alb



#%% Ethnicity european - How many ?

Overall_matrix['22006-0.0'] = Overall_matrix['22006-0.0'].fillna(0)

value_counts = Overall_matrix['21000-0.0'].value_counts()

total_count = value_counts.sum()

desired_values = [1001, 1002, 1003]
desired_counts = value_counts.loc[desired_values].sum()

# Calculate the overall percentage
overall_percentage = (desired_counts / total_count) * 100

# Print the overall percentage
print("Overall percentage: {:.2f}%".format(overall_percentage))

value_counts = Overall_matrix['22006-0.0'].value_counts()

# Get the count of 1 values
count_1 = value_counts.get(1, 0)

# Get the total number of values in the column
total_count = Overall_matrix['22006-0.0'].count()

# Calculate the percentage of 1s compared to the entire column
percentage_1 = (count_1 / total_count) * 100

print(f"Number of 1 values: {count_1}")
print(f"Percentage of 1 values: {percentage_1}%")

Overall_matrix = Overall_matrix.drop(['22006-0.0', '21000-0.0'], axis=1)


#%% Calculating Allele Frequency UK Biobank

allele_frequency = pd.DataFrame(columns=['rsID', 'Allele Frequency'])

for column in Overall_matrix.columns:
    if column.startswith('rs'):
        allele_sum = Overall_matrix[column].sum()
        total_count = len(Overall_matrix[column]) * 2
        allele_frequency.loc[len(allele_frequency)] = [column, allele_sum / total_count]

# Split the rsID into base rsID and EA
allele_frequency[['rsID_clean', 'EA']] = allele_frequency['rsID'].str.extract(r'^(rs\d+)_([ACGT])_eGFR$')

# Rename and reorder for clarity
allele_frequency = allele_frequency[['rsID_clean', 'EA', 'Allele Frequency']]
allele_frequency.columns = ['rsID', 'EA', 'EAF']


# Select only the needed columns from eGFR_matrix_overall
eGFR_meta = eGFR_matrix_overall[['rsID', 'chr', 'Gene']].copy()
# Force chr to stay as string
eGFR_meta['chr'] = eGFR_meta['chr'].astype(str)

# Merge on 'rsID'
merged_df = pd.merge(allele_frequency, eGFR_meta, on='rsID', how='left')

# Reorder columns to your desired structure
df_allele_fq = merged_df[['rsID', 'chr', 'Gene', 'EA', 'EAF']]
df_allele_fq['chr'] = df_allele_fq['chr'].astype(int)


# Export the final table to an Excel file
output_file = 'your/path/EAF.xlsx'
df_allele_fq.to_excel(output_file, index=False)



#%% Preprocessing data 

##########################
### Preprocessing data ###
##########################


######################
### Calculate UACR ###
######################
creatinine_data = pd.read_csv('your/path/creatinine_data.csv', header=None)
creatinine_data.columns = ['eid', 'creatinine_urine']
creatinine_data = creatinine_data.iloc[1:].reset_index(drop=True)

overall_copy = Overall_matrix[['eid', 'micro albumin in urin']].copy()
overall_copy['eid'] = overall_copy['eid'].astype(int)
creatinine_data['eid'] = creatinine_data['eid'].astype(int)

merged_df = pd.merge(overall_copy, creatinine_data, on='eid', how='left')

# First make sure both columns are numeric
merged_df['micro albumin in urin'] = pd.to_numeric(merged_df['micro albumin in urin'], errors='coerce')
merged_df['creatinine_urine'] = pd.to_numeric(merged_df['creatinine_urine'], errors='coerce')

# Compute UACR
merged_df['UACR'] = (merged_df['micro albumin in urin'] / merged_df['creatinine_urine']) * 1000  # Converts to mg/mmol


# Create a DataFrame with only eid and UACR to avoid duplicate columns during merge
uacr_data = merged_df[['eid', 'UACR']]

# Merge UACR into the original Overall_matrix
Overall_matrix = pd.merge(Overall_matrix, uacr_data, on='eid', how='left')

del creatinine_data
del uacr_data
del merged_df
del overall_copy

###############################
### calculating eGFR values ###
###############################


# Convert creatinine from µmol/L to mg/dL
Overall_matrix['Creatinine'] = Overall_matrix['Creatinine'] * 0.0113123


sex_column = Overall_matrix['Sex']
creatinine_column = Overall_matrix['Creatinine']
age_column = Overall_matrix['Age at recruitment']

   
def calculate_eGFR_2021(sex, creatinine, age):
    if sex == 1:
        return 142*((min(creatinine/0.9, 1))**(-0.302))*((max(creatinine/0.9, 1))**(-1.200))*(0.9938**age)
    elif sex == 0:
        return 142*(min(creatinine/0.7, 1))**(-0.241)*(max(creatinine/0.7, 1))**(-1.200)*0.9938**age*1.012

egfr_2021 = []

for eid, sex, creatinine, age in zip(Overall_matrix['eid'], sex_column, creatinine_column, age_column):
    if not np.isnan(sex) and not np.isnan(creatinine) and not np.isnan(age):
        eGFR_2021 = calculate_eGFR_2021(sex, creatinine, age)
    else:
        eGFR_2021 = np.nan
    egfr_2021.append((eid, eGFR_2021))

egfr_2021 = pd.DataFrame(egfr_2021, columns=['eid', 'eGFR_2021'])

Overall_matrix = Overall_matrix.merge(egfr_2021, on='eid', how='left')

Overall_matrix_backup = Overall_matrix


# Deleting variables that are not needed     
variables_to_delete = ['age','age_column','creatinine','creatinine_column','egfr_2021', 'eGFR_2021','eid','sex','sex_column']


for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete

########################################
### Overall matrix removing outliers ###
########################################

eid_column = Overall_matrix['eid']

# Select the columns with phenotype data
phenotype_columns = Overall_matrix[['Sex', 'Systolic blood pressure', 'eGFR_2021', 'Age at recruitment']]

# Calculate the mean and standard deviation for each phenotype column excluding NaN values
mean_values = np.nanmean(phenotype_columns, axis=0)
std_values = np.nanstd(phenotype_columns, axis=0)

# Define the lower and upper bounds for outlier removal
lower_bounds = mean_values - 4 * std_values
upper_bounds = mean_values + 4 * std_values

# Identify the rows where any phenotype value is outside the bounds, excluding NaN values
outlier_rows = np.any((phenotype_columns < lower_bounds) | (phenotype_columns > upper_bounds), axis=1) & np.all(~np.isnan(phenotype_columns), axis=1)

# Remove the outlier rows from the Overall_matrix matrix
Overall_matrix_filtered = Overall_matrix[~outlier_rows].copy()

Overall_matrix_filtered.loc[:, 'eid'] = eid_column.loc[~outlier_rows].values


# summary statistics 


selected_columns = [
    'Sex',
    'Systolic blood pressure',
    'Creatinine',
    'Age at recruitment',
    'eGFR_2021'
]

subset = Overall_matrix_filtered[selected_columns]

subset = subset.dropna()
Summary_overall = subset.describe(include='all')


selected_columns = ['eid',
    'Sex',
    'Systolic blood pressure',
    'Creatinine',
    'Age at recruitment',
    'eGFR_2021'
]

eGFR_population = Overall_matrix_filtered[selected_columns]

eGFR_population = eGFR_population.dropna()


# Deleting variables that are not needed     
variables_to_delete = ['eid_column','lower_bounds','mean_values','outlier_rows','phenotype_columns', 'selected_columns','std_values','subset','upper_bounds']


for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete



# Convert eGFR to ln(eGFR)

Overall_matrix_filtered['eGFR_2021'] = np.log(Overall_matrix_filtered['eGFR_2021'])


#############################
### creating eGFR matrix  ###
#############################

columns_to_keep = ['eid', 'Sex', 'Systolic blood pressure', 'Age at recruitment','eGFR_2021']

# Columns to filter and rename
columns_to_filter = [col for col in Overall_matrix_filtered.columns if col.endswith('_eGFR')]
filtered_columns = [col.replace('_eGFR', '') for col in columns_to_filter]

# Create the new DataFrame with filtered columns
df_eGFR_data_filtered = Overall_matrix_filtered[columns_to_keep + columns_to_filter].copy()
df_eGFR_data_filtered.columns = columns_to_keep + filtered_columns


# Drop rows with NaN values
df_eGFR_data_filtered = df_eGFR_data_filtered.dropna()




# Deleting variables that are not needed     
variables_to_delete = ['columns_to_filter','columns_to_keep',
                       'filtered_columns']

for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete



#%% preprocessing scale


##########################
### Preprocessing data ###
##########################



columns_to_scale = ['Sex', 'Systolic blood pressure', 'Age at recruitment']

# Select all columns that start with 'rs'
rs_columns = [col for col in Overall_matrix_filtered.columns if col.startswith('rs')]


df_selected_columns = Overall_matrix_filtered[columns_to_scale]

# Scale the selected columns
scaler = preprocessing.StandardScaler()
scaled_columns = scaler.fit_transform(df_selected_columns)

Overall_matrix_stand = pd.DataFrame(index=Overall_matrix_filtered.index, data=scaled_columns, columns=columns_to_scale)

# Include the columns that were not scaled
non_scaled_columns = [col for col in Overall_matrix_filtered.columns if col not in columns_to_scale]
Overall_matrix_stand[non_scaled_columns] = Overall_matrix_filtered[non_scaled_columns]


columns_to_drop = ['IID', 'SEX']
Overall_matrix_stand.drop(columns_to_drop, axis=1, inplace=True)




# Deleting variables that are not needed     
variables_to_delete = ['columns_to_drop','columns_to_scale',
                       'non_scaled_columns','rs_columns','scaled_columns']

for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete



#%% MLR 1 2021 eGFR 

###########################
### MLR 1 for eGFR 2021 ###
###########################

# Select the columns with '_eGFR' at the end and include 'eid', 'Sex', and 'age' and SBP
columns_to_select = [col for col in Overall_matrix_stand.columns if col.endswith('_eGFR')] + ['eid', 'Sex', 'Age at recruitment','Systolic blood pressure','eGFR_2021']

# Create a new DataFrame with the selected columns
new_matrix = Overall_matrix_stand[columns_to_select].copy()

pc_cols = [f'PC{i}' for i in range(1, 11)]

new_matrix = pd.merge(
    new_matrix,
    pc[['eid'] + pc_cols],
    on='eid',
    how='left'
)

# Remove rows with NaN values
new_matrix.dropna(inplace=True)
new_matrix.reset_index(drop=True, inplace=True)

# Compute mean and std
mean_egfr_overall = new_matrix['eGFR_2021'].mean()
std_egfr_overall = new_matrix['eGFR_2021'].std()

# Base covariates
X_base = new_matrix[['Sex', 'Age at recruitment',
                     'Systolic blood pressure'] + pc_cols].values

# add intercept
X_base = sm.add_constant(X_base)

# Outcome
y_log = new_matrix['eGFR_2021'].values

# SNP columns
rs_columns = [col for col in new_matrix.columns if col.startswith('rs')]

results = []


for rsID in rs_columns:

    snp = new_matrix[rsID].values.reshape(-1, 1)

    X_full = np.hstack((X_base, snp))

    model = sm.OLS(y_log, X_full)
    result = model.fit()

    beta = result.params[-1]
    p_value = result.pvalues[-1]
    ci_lower, ci_upper = result.conf_int()[-1]

    results.append((rsID, beta, p_value, ci_lower, ci_upper))

MLR_eGFR_2021_all = pd.DataFrame(
    results,
    columns=['rsID', 'Effect Size', 'P-value', 'CI Lower', 'CI Upper']
)

MLR_eGFR_2021_significant_results_all = MLR_eGFR_2021_all[MLR_eGFR_2021_all['P-value'] <= 0.05]

MLR_eGFR_2021_all['FDR'] = multipletests(
    MLR_eGFR_2021_all['P-value'], 
    method='fdr_bh'
)[1]

MLR_eGFR_2021_FDR_significant = MLR_eGFR_2021_all[
    MLR_eGFR_2021_all['FDR'] <= 0.05
]


variables_to_delete = ['result','results','rs_columns','rsID', 'scaler','snp','X_base','X_full','y_log',
                       'p_value','new_matrix','model','columns_to_select','ci_upper','ci_lower', 'beta']



for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete



#%% MLR 2
###########################################
### eGFR 2021 rsID associated with UACR ### 
###########################################

# including micro albumin
alb = Overall_matrix_backup[['eid', 'UACR']].dropna()
eid_column = alb['eid']

# Select the columns with phenotype data
phenotype_columns = alb[['UACR']]

# Calculate the mean and standard deviation for each phenotype column excluding NaN values
mean_values = np.nanmean(phenotype_columns, axis=0)
std_values = np.nanstd(phenotype_columns, axis=0)

# Define the lower and upper bounds for outlier removal
lower_bounds = mean_values - 4 * std_values
upper_bounds = mean_values + 4 * std_values

# Identify the rows where any phenotype value is outside the bounds, excluding NaN values
outlier_rows = np.any((phenotype_columns < lower_bounds) | (phenotype_columns > upper_bounds), axis=1) & np.all(~np.isnan(phenotype_columns), axis=1)

# Remove the outlier rows from the Overall_matrix matrix
alb_filtered = alb[~outlier_rows].copy()

alb_filtered.loc[:, 'eid'] = eid_column.loc[~outlier_rows].values
rsIDs_to_select = MLR_eGFR_2021_significant_results_all['rsID'].tolist()

# add PC1 to 10
pc_cols = [f'PC{i}' for i in range(1, 11)]

base = Overall_matrix_stand[['eid', 'Sex', 'Age at recruitment', 'Systolic blood pressure']].copy()
base = pd.merge(base, alb_filtered[['eid', 'UACR']], on='eid', how='inner')

# merge PC 
base = pd.merge(base, pc[['eid'] + pc_cols], on='eid', how='left')

# drop NA 
base.dropna(inplace=True)
base.reset_index(drop=True, inplace=True)

# summary statistics 
selected_columns = ['UACR']
subset = base[selected_columns].dropna()

Summary_overall2 = subset.describe(include='all')
Summary_combined = pd.concat([Summary_overall, Summary_overall2], axis=1)

UACR_df = pd.merge(Overall_matrix_filtered, alb_filtered[['eid', 'UACR']], on='eid', how='inner')

selected_columns = ['eid',
    'Sex',
    'Systolic blood pressure',
    'Age at recruitment',
    'UACR_x'
]

UACR_population = UACR_df[selected_columns]
UACR_population = UACR_population.dropna()


base['UACR'] = np.log(base['UACR'])
mean_uae_overall = base['UACR'].mean()
std_uae_overall  = base['UACR'].std()

y_log = base['UACR']

X = base[['Sex', 'Age at recruitment', 'Systolic blood pressure'] + pc_cols]

X_base = sm.add_constant(X.to_numpy())
y_np   = y_log.to_numpy()

# Bonferroni
bonferroni_threshold_eGFR_2021_associated_with_UACR = 0.05 / len(rsIDs_to_select)

# Loop one SNP at a time
results = []

base_eids = base[['eid']] 

for rsID in rsIDs_to_select:
    snp_df = Overall_matrix_stand[['eid', rsID]]
    snp_aligned = pd.merge(base_eids, snp_df, on='eid', how='left')[rsID].to_numpy().reshape(-1, 1)
    mask = ~np.isnan(snp_aligned[:, 0])
    if mask.sum() < 3:
        continue

    X_full = np.hstack((X_base[mask], snp_aligned[mask]))
    y_use  = y_np[mask]

    res = sm.OLS(y_use, X_full).fit()

    beta = res.params[-1]
    p_value = res.pvalues[-1]
    ci_low, ci_high = res.conf_int()[-1]

    results.append((rsID, beta, p_value, ci_low, ci_high))


# Store results
eGFR_2021_rsid_associated_with_UACR_model2 = pd.DataFrame(
    results, columns=['rsID', 'Effect Size', 'P-value', 'CI Lower', 'CI Upper']
)

eGFR_2021_rsid_associated_with_UACR_bonferroni_model2 = (
    eGFR_2021_rsid_associated_with_UACR_model2[
        eGFR_2021_rsid_associated_with_UACR_model2['P-value'] <= bonferroni_threshold_eGFR_2021_associated_with_UACR
    ]
)

eGFR_2021_rsid_associated_with_UACR_significant_model2 = (
    eGFR_2021_rsid_associated_with_UACR_model2[
        eGFR_2021_rsid_associated_with_UACR_model2['P-value'] <= 0.05
    ]
)

eGFR_2021_rsid_associated_with_UACR_model2['FDR'] = multipletests(
    eGFR_2021_rsid_associated_with_UACR_model2['P-value'],
    method='fdr_bh'
)[1]

eGFR_2021_rsid_associated_with_UACR_FDR_model2 = (
    eGFR_2021_rsid_associated_with_UACR_model2[
        eGFR_2021_rsid_associated_with_UACR_model2['FDR'] <= 0.05
    ]
)



variables_to_delete = [ 'base','beta', 'ci_high','ci_low', 'df_selected_columns', 'base_eids',
                       'eid_column', 'lower_bounds','upper_bounds', 'mask', 'outlier_rows','p_value',
                       'res','results','rsID','rsIDs_to_select','selected_columns','snp_df', 'snp_aligned', 'subset', 'X',
                       'X_base', 'X_full', 'y_log', 'y_np','y_use']

for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete



#%% Preprocessing - creating matrix for scatterplots
##############################################################
### eGFR_2021_rsid_associated_with_UACR_significant_model2 ###
##############################################################

# Load Excel file
lookup = pd.read_excel('your/path')

# Extract base rsID and EA from rsID
split_rsID = eGFR_2021_rsid_associated_with_UACR_model2['rsID'].str.split('_', expand=True)
eGFR_2021_rsid_associated_with_UACR_model2['base_rsID'] = split_rsID[0]
eGFR_2021_rsid_associated_with_UACR_model2['EA'] = split_rsID[1]

# Merge on base_rsID
merged = pd.merge(
    eGFR_2021_rsid_associated_with_UACR_model2,
    lookup[['base_rsID', 'chr', 'Gene']],
    on='base_rsID',
    how='left'
)

# Rename and format final output
MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2 = merged.rename(columns={
    'base_rsID': 'rsID_eGFR_2021',
    'Effect Size': 'Beta_UACR',
    'P-value': 'p_val_UACR',
    'FDR': 'FDR_UACR',
    'CI Lower': 'UACR CI Lower',
    'CI Upper': 'UACR CI Upper'
})[[
    'rsID_eGFR_2021', 'chr', 'Gene', 'EA',
    'Beta_UACR', 'p_val_UACR','FDR_UACR', 'UACR CI Lower', 'UACR CI Upper'
]]

#########################
### UACR X eGFR 2021  ### eGFR
#########################
Scatter_plot_eGFR_2021_associated_with_UACR = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2.copy()

# Step 1: Extract base rsID from MLR1 results
MLR_eGFR_2021_all['rsID_prefix'] = MLR_eGFR_2021_all['rsID'].str.split('_').str[0]

# Step 2: Merge MLR1 results into scatter plot matrix
Scatter_plot_eGFR_2021_associated_with_UACR = Scatter_plot_eGFR_2021_associated_with_UACR.merge(
    MLR_eGFR_2021_all[['rsID_prefix', 'Effect Size', 'P-value','FDR', 'CI Lower', 'CI Upper']],
    left_on='rsID_eGFR_2021',
    right_on='rsID_prefix',
    how='left'
)

# Step 3: Rename merged columns
Scatter_plot_eGFR_2021_associated_with_UACR.rename(columns={
    'Effect Size': 'Beta_eGFR_2021',
    'P-value': 'p_val_eGFR_2021',
    'FDR': 'FDR',
    'CI Lower': 'eGFR CI Lower',
    'CI Upper': 'eGFR CI Upper'
}, inplace=True)

# Step 4: Drop the helper merge key
Scatter_plot_eGFR_2021_associated_with_UACR.drop(columns='rsID_prefix', inplace=True)


#%% Scatter plot 2021 
############################################
### filtereing only the significant data ###
############################################
Scatter_plot_eGFR_2021_associated_with_UACR_significant = Scatter_plot_eGFR_2021_associated_with_UACR[
    (Scatter_plot_eGFR_2021_associated_with_UACR['p_val_eGFR_2021'] <=  0.05) &
    (Scatter_plot_eGFR_2021_associated_with_UACR['p_val_UACR'] <=  0.05)
].reset_index(drop=True)



############################################################
####  Scatterplot eGFR 2021 rsID association with UACR #####
############################################################

Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'] = \
    pd.to_numeric(Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'], errors='coerce')

Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'] = \
    pd.to_numeric(Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'], errors='coerce')

# Extract base rsIDs from Bonferroni dataset
uacr_bonferroni_rsIDs = (
    eGFR_2021_rsid_associated_with_UACR_bonferroni_model2['rsID']
    .str.split('_')
    .str[0]
)

# Assign plotting group
Scatter_plot_eGFR_2021_associated_with_UACR_significant['plotting_groups'] = np.where(
    Scatter_plot_eGFR_2021_associated_with_UACR_significant['rsID_eGFR_2021']
    .isin(uacr_bonferroni_rsIDs),
    'Bonferroni',
    'P<0.05'
)


# Create plot function
def plot_scatter(data, annotate=False, save_path=None):

    g = sns.relplot(
        data=data,
        x='Beta_eGFR_2021',
        y='Beta_UACR',
        hue='plotting_groups',
        palette={'P<0.05': 'black', 'Bonferroni': '#4A90E2'},
        aspect=2,
        s=45,
        alpha=0.8,
        legend=False
    )

    g.ax.set_xlabel(r'$\ln(eGFR_{2021})$ beta', fontsize=14)
    g.ax.set_ylabel('ln(UACR) beta', fontsize=14)

    g.ax.axhline(0, color='black', lw=0.5)
    g.ax.axvline(0, color='black', lw=0.5)
    g.ax.grid(True, linestyle='--', linewidth=0.5)

    g.ax.set_xlim(-0.015, 0.01)
    g.ax.set_ylim(-0.075, 0.06)

    # Optional annotation
    if annotate:
        bonf_df = data[data['plotting_groups'] == 'Bonferroni']

        annotations = bonf_df.apply(
            lambda p: g.ax.annotate(
                p['Gene'],
                (p['Beta_eGFR_2021'], p['Beta_UACR']),
                fontsize=11
            ),
            axis=1
        ).to_list()

        adjust_text(
            annotations,
            arrowprops=dict(arrowstyle="-", color='k', lw=0.4),
            force_points=0.2,
            force_text=0.1
        )

    handles = [
        plt.Line2D([], [], color='black', marker='o', linestyle='None'),
        plt.Line2D([], [], color='#4A90E2', marker='o', linestyle='None')
    ]

    labels = [
    'Nominal significance (P < 0.05)',
    r'Bonferroni significance ($P_{b}$ < 1.18 × 10$^{-4}$)'
]

    g.ax.legend(handles, labels, loc='upper left', fontsize=12, frameon=True)

    if save_path:
        plt.savefig(save_path, dpi=500, bbox_inches='tight')

    plt.show()


# plots

save_path = 'your/path/scatterplot_eGFR_2021_associated_with_UACR'

plot_scatter(
    Scatter_plot_eGFR_2021_associated_with_UACR_significant,
    annotate=True,
    save_path='your/path/scatterplot_eGFR_2021_associated_with_UACR_with_labels.jpg'
)


plot_scatter(
    Scatter_plot_eGFR_2021_associated_with_UACR_significant,
    annotate=False,
    save_path='your/path/scatterplot_eGFR_2021_associated_with_UACR_no_labels.jpg'
)


#%% Sensitivity analysis 
############################################
### Sensitivity analysis: add eGFR_2021 ####
### Only for Bonferroni-significant SNPs ###
############################################

rsIDs_to_select = eGFR_2021_rsid_associated_with_UACR_bonferroni_model2['rsID'].tolist()

base = Overall_matrix_stand[['eid', 'Sex', 'Age at recruitment',
                              'Systolic blood pressure', 'eGFR_2021']].copy()

# Merge UACR
base = pd.merge(base, alb_filtered[['eid', 'UACR']], on='eid', how='inner')

# Merge PC (remove this line if PCs already in Overall_matrix_stand)
base = pd.merge(base, pc[['eid'] + pc_cols], on='eid', how='left')

# Drop NA once
base.dropna(inplace=True)
base.reset_index(drop=True, inplace=True)

base['UACR'] = np.log(base['UACR'])

# Scale eGFR
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
base['eGFR_2021_std'] = scaler.fit_transform(base[['eGFR_2021']])

# Define outcome
y_np = base['UACR'].to_numpy()

X_cov = base[['Sex', 'Age at recruitment',
              'Systolic blood pressure',
              'eGFR_2021_std'] + pc_cols].to_numpy()

X_base = sm.add_constant(X_cov)

base_eids = base[['eid']]

results = []

for rsID in rsIDs_to_select:

    # align SNP to base via eid
    snp_df = Overall_matrix_stand[['eid', rsID]]
    snp_aligned = pd.merge(base_eids, snp_df, on='eid', how='left')[rsID].to_numpy().reshape(-1,1)

    mask = ~np.isnan(snp_aligned[:,0])

    if mask.sum() < 3:
        continue

    X_full = np.hstack((X_base[mask], snp_aligned[mask]))
    y_use  = y_np[mask]

    res = sm.OLS(y_use, X_full).fit()

    beta = res.params[-1]
    pval = res.pvalues[-1]
    ci_low, ci_high = res.conf_int()[-1]

    results.append((rsID, beta, pval, ci_low, ci_high))

# -----------------------------------
# 7) Store results
# -----------------------------------

UACR_sensitivity_model1 = pd.DataFrame(
    results, columns=['rsID', 'Effect Size', 'P-value', 'CI Lower', 'CI Upper']
)

UACR_sensitivity_model1_significant = \
    UACR_sensitivity_model1[UACR_sensitivity_model1['P-value'] <= 0.05]
           
#################################
### UACR missingness analysis ###
#################################

from scipy import stats

has_uacr = Overall_matrix_filtered['UACR'].notna()
has_egfr = Overall_matrix_filtered['eGFR_2021'].notna()

print("Total:", len(Overall_matrix_filtered))
print("With UACR:", has_uacr.sum())
print("Without UACR:", (~has_uacr).sum())
print("With eGFR:", has_egfr.sum())

summary_vars = [
    'Sex',
    'Systolic blood pressure',
    'Creatinine',
    'Age at recruitment',
    'eGFR_2021',
    'UACR'
]

Summary_with_UACR = Overall_matrix_filtered.loc[
    has_uacr, summary_vars
].describe(include='all')


Summary_with_eGFR = Overall_matrix_filtered.loc[
    has_egfr, summary_vars
].describe(include='all')


#################################
### Descriptive comparison ######
#################################

def standardized_mean_difference(x1, x2):
    x1 = pd.Series(x1).dropna()
    x2 = pd.Series(x2).dropna()

    mean1, mean2 = x1.mean(), x2.mean()
    sd1, sd2 = x1.std(), x2.std()

    pooled_sd = np.sqrt((sd1**2 + sd2**2) / 2)

    if pooled_sd == 0:
        return np.nan

    return (mean1 - mean2) / pooled_sd


def standardized_difference_binary(g1, g2):
    g1 = pd.Series(g1).dropna()
    g2 = pd.Series(g2).dropna()

    p1 = (g1 == 1).mean()
    p2 = (g2 == 1).mean()

    denom = np.sqrt((p1*(1-p1) + p2*(1-p2)) / 2)

    if denom == 0:
        return np.nan

    return (p1 - p2) / denom


results = []

continuous_vars = [
    'Systolic blood pressure',
    'Creatinine',
    'Age at recruitment',
    'eGFR_2021',
    'UACR'
]

for var in continuous_vars:

    valid_mask = Overall_matrix_filtered[var].notna()

    g1 = Overall_matrix_filtered.loc[has_uacr & valid_mask, var]
    g2 = Overall_matrix_filtered.loc[has_egfr & valid_mask, var]

    smd = standardized_mean_difference(g1, g2)

    results.append({
        'Variable': var,
        'UACR group': f"{g1.mean():.2f} ± {g1.std():.2f}",
        'eGFR group': f"{g2.mean():.2f} ± {g2.std():.2f}",
        'SMD': smd,
        'N UACR group': g1.count(),
        'N eGFR group': g2.count()
    })
    
    
sex_mask = Overall_matrix_filtered['Sex'].notna()

sex_g1 = Overall_matrix_filtered.loc[has_uacr & sex_mask, 'Sex']
sex_g2 = Overall_matrix_filtered.loc[has_egfr & sex_mask, 'Sex']

n_female_uacr = (sex_g1 == 0).sum()
n_male_uacr   = (sex_g1 == 1).sum()

n_female_egfr = (sex_g2 == 0).sum()
n_male_egfr   = (sex_g2 == 1).sum()

smd_sex = standardized_difference_binary(sex_g1, sex_g2)

results.append({
    'Variable': 'Sex',
    'UACR group': f"Female: {n_female_uacr}, Male: {n_male_uacr}",
    'eGFR group': f"Female: {n_female_egfr}, Male: {n_male_egfr}",
    'SMD': smd_sex,
    'N UACR group': sex_g1.count(),
    'N eGFR group': sex_g2.count()
})

comparison_table = pd.DataFrame(results)

n_row = pd.DataFrame([{
    'Variable': 'Total N',
    'UACR group': has_uacr.sum(),
    'eGFR group': has_egfr.sum(),
    'SMD': np.nan,
    'N UACR group': has_uacr.sum(),
    'N eGFR group': has_egfr.sum()
}])
comparison_table = pd.concat([n_row, comparison_table], ignore_index=True)

print(comparison_table)

)


output_file = 'your/path/UACR_missingness_analysis.xlsx'

with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
    
    comparison_table.to_excel(
        writer,
        sheet_name='Descriptive_UACR_vs_eGFR',
        index=False
    )
    


#%% Doing it only for diabetes 

#####################################
### creating diabetes matrix data ###
#####################################

# Load the diabetes csv data using the Pandas library
df_T2D = pd.read_csv('your/path')

# Select only the 'eid' and 'T2D' columns
df_T2D = df_T2D[['eid', 'T2D']]
t2d_eids = df_T2D.loc[df_T2D['T2D'] == 1, 'eid']
t2d_eids = t2d_eids[~t2d_eids.isin(withdraw['Eid'])]

# Filter Overall_matrix to T2D 
Overall_matrix_T2D = Overall_matrix.loc[Overall_matrix['eid'].isin(t2d_eids)].copy()

del df_T2D

# add BMI and Hba1c
Overall_matrix_T2D = Overall_matrix_T2D.merge(bmi_hba1c, on='eid', how='inner')

##########################
### Preprocessing data ###
##########################

eid_column = Overall_matrix_T2D['eid']

# Select the columns with phenotype data
phenotype_columns = Overall_matrix_T2D[['Sex', 'Systolic blood pressure', 'eGFR_2021', 'Age at recruitment', 'hba1c','BMI']]

# Calculate the mean and standard deviation for each phenotype column excluding NaN values
mean_values = np.nanmean(phenotype_columns, axis=0)
std_values = np.nanstd(phenotype_columns, axis=0)

# Define the lower and upper bounds for outlier removal
lower_bounds = mean_values - 4 * std_values
upper_bounds = mean_values + 4 * std_values

# Identify the rows where any phenotype value is outside the bounds, excluding NaN values
outlier_rows = np.any((phenotype_columns < lower_bounds) | (phenotype_columns > upper_bounds), axis=1) & np.all(~np.isnan(phenotype_columns), axis=1)

# Remove the outlier rows from the Overall_matrix matrix
Overall_matrix_T2D_filtered = Overall_matrix_T2D[~outlier_rows].copy()

Overall_matrix_T2D_filtered.loc[:, 'eid'] = eid_column.loc[~outlier_rows].values

Overall_matrix_T2D_backup = Overall_matrix_T2D

# summary statistics 
selected_columns = [
    'Sex',
    'Systolic blood pressure',
    'Creatinine',
    'Age at recruitment',
    'eGFR_2021',
    'BMI',
    'hba1c'
]

subset = Overall_matrix_T2D_filtered[selected_columns]

subset = subset.dropna()

Summary_T2D = subset.describe(include='all')



variables_to_delete = ['eid_column', 'phenotype_columns', 'mean_values','std_values', 'lower_bounds', 'upper_bounds',
                       'outlier_rows','selected_columns', 'subset']

for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete




# Convert  eGFR to ln(eGFR)
Overall_matrix_T2D_filtered['eGFR_2021'] = np.log(Overall_matrix_T2D_filtered['eGFR_2021'])

columns_to_scale = ['Sex', 'Systolic blood pressure', 'Age at recruitment', 'BMI', 'hba1c']

df_selected_columns = Overall_matrix_T2D_filtered[columns_to_scale]

# Scale the selected columns
scaled_columns = preprocessing.scale(df_selected_columns)

Overall_matrix_stand_T2D = pd.DataFrame(index=Overall_matrix_T2D_filtered.index, data=scaled_columns, columns=columns_to_scale)

# Include the columns that were not scaled
non_scaled_columns = [col for col in Overall_matrix_T2D_filtered.columns if col not in columns_to_scale]
Overall_matrix_stand_T2D[non_scaled_columns] = Overall_matrix_T2D_filtered[non_scaled_columns]


columns_to_drop = ['IID', 'SEX']
Overall_matrix_stand_T2D.drop(columns_to_drop, axis=1, inplace=True)

del columns_to_scale
del df_selected_columns
del scaled_columns
del non_scaled_columns
del columns_to_drop


###############################
### MLR 1 for eGFR 2021 T2D ###
###############################

base = Overall_matrix_stand_T2D[['eid','Sex','Age at recruitment','Systolic blood pressure','eGFR_2021','BMI','hba1c']].copy()
base = base.merge(pc[['eid'] + pc_cols], on='eid', how='left')
base = base.dropna().reset_index(drop=True)

# summery stats
mean_egfr_T2D = base['eGFR_2021'].mean()
std_egfr_T2D  = base['eGFR_2021'].std()


y_np = base['eGFR_2021'].to_numpy()

# X_base (BMI+hba1c + PC)
X_cov = base[['Sex', 'Age at recruitment', 'Systolic blood pressure', 'BMI', 'hba1c'] + pc_cols].to_numpy()

X_base = sm.add_constant(X_cov)

base_eids = base[['eid']]

# SNP-list
rs_columns = [col for col in Overall_matrix_stand_T2D.columns if col.endswith('_eGFR') and col.startswith('rs')]

# Loop
results = []

for rsID in rs_columns:
    snp_df = Overall_matrix_stand_T2D[['eid', rsID]]

    snp = pd.merge(base_eids, snp_df, on='eid', how='left')[rsID].to_numpy().reshape(-1, 1)

    mask = ~np.isnan(snp[:, 0])
    if mask.sum() < 3:
        continue

    X_full = np.hstack((X_base[mask], snp[mask]))
    y_use  = y_np[mask]

    res = sm.OLS(y_use, X_full).fit()

    beta = res.params[-1]
    p_value = res.pvalues[-1]
    ci_low, ci_high = res.conf_int()[-1]

    results.append((rsID, beta, p_value, ci_low, ci_high))

# Results
MLR_eGFR_2021_all_T2D = pd.DataFrame(
    results, columns=['rsID', 'Effect Size', 'P-value', 'CI Lower', 'CI Upper']
)

MLR_eGFR_2021_significant_results_all_T2D = MLR_eGFR_2021_all_T2D[
    MLR_eGFR_2021_all_T2D['P-value'] <= 0.05
]
MLR_eGFR_2021_all_T2D['FDR'] = multipletests(
    MLR_eGFR_2021_all_T2D['P-value'],
    method='fdr_bh'
)[1]

MLR_eGFR_2021_significant_FDR_results_all_T2D = (
    MLR_eGFR_2021_all_T2D[
        MLR_eGFR_2021_all_T2D['FDR'] <= 0.05
    ]
)


###########################################
### eGFR 2021 rsID associated with UACR ### 
###########################################
# including micro albumin
alb_T2D = Overall_matrix_T2D_backup[['eid', 'UACR']]
alb_T2D=alb_T2D.dropna()

eid_column = alb_T2D['eid']

# Select the columns with phenotype data
phenotype_columns = alb_T2D[['UACR']]

# Calculate the mean and standard deviation for each phenotype column excluding NaN values
mean_values = np.nanmean(phenotype_columns, axis=0)
std_values = np.nanstd(phenotype_columns, axis=0)

# Define the lower and upper bounds for outlier removal
lower_bounds = mean_values - 4 * std_values
upper_bounds = mean_values + 4 * std_values

# Identify the rows where any phenotype value is outside the bounds, excluding NaN values
outlier_rows = np.any((phenotype_columns < lower_bounds) | (phenotype_columns > upper_bounds), axis=1) & np.all(~np.isnan(phenotype_columns), axis=1)

# Remove the outlier rows from the Overall_matrix matrix
alb_filtered_T2D = alb_T2D[~outlier_rows].copy()

alb_filtered_T2D.loc[:, 'eid'] = eid_column.loc[~outlier_rows].values



variables_to_delete = ['eid_column','phenotype_columns', 'mean_values', 'std_values', 'lower_bounds', 'upper_bounds',
                       'outlier_rows',]

for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete



rsIDs_to_select = MLR_eGFR_2021_significant_results_all_T2D['rsID'].tolist()

base = Overall_matrix_stand_T2D[['eid', 'Sex', 'Age at recruitment', 'Systolic blood pressure', 'BMI', 'hba1c']].copy()

# merge UACR in
base = base.merge(alb_filtered_T2D, on='eid', how='inner')

# include PCs:
base = base.merge(pc[['eid'] + pc_cols], on='eid', how='left')

# drop NA once
base = base.dropna().reset_index(drop=True)


# Summary stats 
Summary_T2D_2 = base[['UACR']].describe(include='all')
Summary_combined_T2D = pd.concat([Summary_T2D, Summary_T2D_2], axis=1)

# Excel export (same)
output_file = r'your/path\Summary.xlsx'
with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    Summary_combined.to_excel(writer, sheet_name='Overall')
    Summary_combined_T2D.to_excel(writer, sheet_name='T2D')


# log(UACR) + mean/std 

base['UACR'] = np.log(base['UACR'].to_numpy())
mean_uae_T2D = base['UACR'].mean()
std_uae_T2D  = base['UACR'].std()


# Build X_base once 
y_np = base['UACR'].to_numpy()
X_cols = ['Sex', 'Age at recruitment', 'Systolic blood pressure', 'BMI', 'hba1c'] + pc_cols
X_cov = base[X_cols].to_numpy()
X_base = sm.add_constant(X_cov)

base_eids = base[['eid']]

# Bonferroni threshold based on number of SNPs tested
bonferroni_threshold_eGFR_2021_associated_with_UACR_T2D = 0.05 / len(rsIDs_to_select)


# SNP loop
results = []  

for rsID in rsIDs_to_select:
    snp_df = Overall_matrix_stand_T2D[['eid', rsID]]

    snp = base_eids.merge(snp_df, on='eid', how='left')[rsID].to_numpy().reshape(-1, 1)

    mask = ~np.isnan(snp[:, 0])
    if mask.sum() < 3:
        continue

    X_full = np.hstack((X_base[mask], snp[mask]))
    y_use  = y_np[mask]

    res = sm.OLS(y_use, X_full).fit()

    beta = res.params[-1]
    pval = res.pvalues[-1]
    ci_low, ci_high = res.conf_int()[-1]

    results.append((rsID, beta, pval, ci_low, ci_high))


# results

eGFR_2021_rsid_associated_with_UACR_model2_T2D = pd.DataFrame(
    results, columns=['rsID', 'Effect Size', 'P-value', 'CI Lower', 'CI Upper']
)

eGFR_2021_rsid_associated_with_UACR_bonferroni_model2_T2D = (
    eGFR_2021_rsid_associated_with_UACR_model2_T2D[
        eGFR_2021_rsid_associated_with_UACR_model2_T2D['P-value'] <= bonferroni_threshold_eGFR_2021_associated_with_UACR_T2D
    ]
)

eGFR_2021_rsid_associated_with_UACR_significant_model2_T2D = (
    eGFR_2021_rsid_associated_with_UACR_model2_T2D[
        eGFR_2021_rsid_associated_with_UACR_model2_T2D['P-value'] <= 0.05
    ]
)
# Add FDR (Benjamini–Hochberg)
eGFR_2021_rsid_associated_with_UACR_model2_T2D['FDR'] = multipletests(
    eGFR_2021_rsid_associated_with_UACR_model2_T2D['P-value'],
    method='fdr_bh'
)[1]

# FDR-significant
eGFR_2021_rsid_associated_with_UACR_FDR_model2_T2D = (
    eGFR_2021_rsid_associated_with_UACR_model2_T2D[
        eGFR_2021_rsid_associated_with_UACR_model2_T2D['FDR'] <= 0.05
    ]
)


##################################################################
### eGFR_2021_rsid_associated_with_UACR_significant_model2_T2D ###
##################################################################
split_rsID = eGFR_2021_rsid_associated_with_UACR_model2_T2D['rsID'].str.split('_', expand=True)
eGFR_2021_rsid_associated_with_UACR_model2_T2D['base_rsID'] = split_rsID[0]
eGFR_2021_rsid_associated_with_UACR_model2_T2D['EA'] = split_rsID[1]

lookup = pd.read_excel('your/path')

merged = pd.merge(
    eGFR_2021_rsid_associated_with_UACR_model2_T2D,
    lookup[['base_rsID', 'chr', 'Gene']],
    on='base_rsID',
    how='left'
)
MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2_T2D = merged.rename(columns={
    'base_rsID': 'rsID_eGFR_2021',
    'Effect Size': 'Beta_UACR',
    'P-value': 'p_val_UACR',
    'FDR': 'FDR_UACR',
    'CI Lower': 'UACR CI Lower',
    'CI Upper': 'UACR CI Upper'
})[[
    'rsID_eGFR_2021', 'chr', 'Gene', 'EA',
    'Beta_UACR', 'p_val_UACR','FDR_UACR', 'UACR CI Lower', 'UACR CI Upper'
]]

#########################
### UACR X eGFR 2021  ### 
#########################

Scatter_plot_eGFR_2021_associated_with_UACR_T2D = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2_T2D.copy()

MLR_eGFR_2021_all_T2D['rsID_prefix'] = MLR_eGFR_2021_all_T2D['rsID'].str.split('_').str[0]

Scatter_plot_eGFR_2021_associated_with_UACR_T2D = Scatter_plot_eGFR_2021_associated_with_UACR_T2D.merge(
    MLR_eGFR_2021_all_T2D[['rsID_prefix', 'Effect Size', 'P-value', 'CI Lower', 'CI Upper']],
    left_on='rsID_eGFR_2021',
    right_on='rsID_prefix',
    how='left'
)

Scatter_plot_eGFR_2021_associated_with_UACR_T2D.rename(columns={
    'Effect Size': 'Beta_eGFR_2021',
    'P-value': 'p_val_eGFR_2021',
    'CI Lower': 'eGFR CI Lower',
    'CI Upper': 'eGFR CI Upper'
}, inplace=True)

Scatter_plot_eGFR_2021_associated_with_UACR_T2D.drop(columns='rsID_prefix', inplace=True)


############################################
### filtereing only the significant data ###
############################################



Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D = Scatter_plot_eGFR_2021_associated_with_UACR_T2D[
    (Scatter_plot_eGFR_2021_associated_with_UACR_T2D['p_val_eGFR_2021'] <=  0.05) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_T2D['p_val_UACR'] <=  0.05)
].reset_index(drop=True)


#%%

############################################################
####  Scatterplot eGFR 2021 rsID association with UACR ##### 
############################################################

Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_UACR'] = pd.to_numeric(Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_UACR'], errors='coerce')
Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_eGFR_2021'] = pd.to_numeric(Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_eGFR_2021'], errors='coerce')

Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['plotting_groups'] = 'A'
# Extract base rsIDs from the Bonferroni dataset
uacr_bonferroni_rsIDs_T2D = eGFR_2021_rsid_associated_with_UACR_bonferroni_model2_T2D['rsID'].str.split('_').str[0]

# Assign 'B' to rows where rsID_eGFR matches a Bonferroni rsID
Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['plotting_groups'] = np.where(
    Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['rsID_eGFR_2021'].isin(uacr_bonferroni_rsIDs_T2D),
    'B',
    Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D.get('plotting_groups', 'A')  # Default to 'A'
)


# Create the scatter plot 
g = sns.relplot(data=Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D,
                x='Beta_eGFR_2021',
                y='Beta_UACR',
                aspect=2,
                hue='plotting_groups',
                palette=['black', '#4A90E2'],
                linewidth=0,
                s=45,
                alpha=0.8,
                legend=None)


# Add labels and title to the plot
g.ax.set_xlabel(r'$\ln(eGFR_{2021})$ beta', fontsize=14)
g.ax.set_ylabel('ln(UACR) beta', fontsize=14)
g.ax.axhline(0, color='black', lw=0.5)
g.ax.axvline(0, color='black', lw=0.5)

# Set x and y axis limits
g.ax.set_xlim(-0.01, 0.01)
g.ax.set_ylim(-0.25, 0.25)

# Create the annotations
#Anno = Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D.apply(
#    lambda p: g.ax.annotate(p['Gene'], (p['Beta_eGFR_2021'], p['Beta_UACR']), fontsize=11),
#    axis=1
#).to_list()

#adjust_text(Anno, arrowprops=dict(arrowstyle="-", color='k', lw=0.4), force_points=0.1, force_text=0.2)


# Define legend labels and colors
legend_labels = [
    'Nominal significance (P < 0.05)',
    r'Bonferroni significance ($P_{b}$ < 6.58 × 10$^{-4}$)'
]
legend_markers = ['o','o']  # Specify the markers for each label
legend_colors = ['black','#26538d']

handles = [plt.Line2D([], [], color=color, marker=marker, linestyle='None') for color, marker in zip(legend_colors, legend_markers)]

g.ax.grid(True, which='both', linestyle='--', linewidth=0.5)
g.ax.legend(handles, legend_labels, loc='upper right', markerscale=1, fontsize=12, frameon=True, edgecolor='black')

# save the plot
save_path = 'your/path/scatterplot_eGFR_2021_associated_with_UACR_T2D_text.jpg'
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.show()

#%% Sensitivity analysis T2D
############################################
### Sensitivity analysis: add eGFR_2021 ####
### Only for Bonferroni-significant SNPs ###
############################################

rsIDs_to_select = eGFR_2021_rsid_associated_with_UACR_bonferroni_model2_T2D['rsID'].tolist()

base_cols = [
    'eid',
    'Sex',
    'Age at recruitment',
    'Systolic blood pressure',
    'eGFR_2021',
    'BMI',
    'hba1c'
]

base = Overall_matrix_stand_T2D[base_cols].copy()
base = base.merge(alb_filtered_T2D[['eid', 'UACR']], on='eid', how='inner')
base = base.merge(pc[['eid'] + pc_cols], on='eid', how='left')
base = base.dropna().reset_index(drop=True)

base['UACR'] = np.log(base['UACR'].to_numpy())
 
# Scale eGFR
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
base['eGFR_2021_std'] = scaler.fit_transform(base[['eGFR_2021']])

# Define outcome
egfr = base['eGFR_2021'].to_numpy()
base['eGFR_2021_std'] = (egfr - egfr.mean()) / egfr.std()

y_np = base['UACR'].to_numpy()


X_cols = ['Sex', 'Age at recruitment', 'Systolic blood pressure',
          'eGFR_2021_std', 'BMI', 'hba1c'] + pc_cols

X_base = sm.add_constant(base[X_cols].to_numpy())

base_eids = base[['eid']]
results = []  # vigtigt i Spyder: reset hver gang

for rsID in rsIDs_to_select:

    # hent SNP fra T2D SNP-matrixen (brug den samme som du brugte til model2_T2D)
    snp_df = Overall_matrix_stand_T2D[['eid', rsID]]

    # align til base via eid (samme rækkefølge som X_base/y_np)
    snp = base_eids.merge(snp_df, on='eid', how='left')[rsID].to_numpy().reshape(-1, 1)

    mask = ~np.isnan(snp[:, 0])
    if mask.sum() < 3:
        continue

    X_full = np.hstack((X_base[mask], snp[mask]))
    y_use  = y_np[mask]

    res = sm.OLS(y_use, X_full).fit()

    beta = res.params[-1]
    pval = res.pvalues[-1]
    ci_low, ci_high = res.conf_int()[-1]

    results.append((rsID, beta, pval, ci_low, ci_high))

UACR_sensitivity_model1_T2D = pd.DataFrame(
    results, columns=['rsID', 'Effect Size', 'P-value', 'CI Lower', 'CI Upper']
)

UACR_sensitivity_model1_significant_T2D = UACR_sensitivity_model1_T2D[
    UACR_sensitivity_model1_T2D['P-value'] <= 0.05
]



# Export the final table to an Excel file
output_file = r'your/path/sensitivity_Results.xlsx'

with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    UACR_sensitivity_model1.to_excel(writer, sheet_name='Overall', index=False)
    UACR_sensitivity_model1_significant.to_excel(writer, sheet_name='Overall significant', index=False)
    UACR_sensitivity_model1_T2D.to_excel(writer, sheet_name='T2D', index=False)
    UACR_sensitivity_model1_significant_T2D.to_excel(writer, sheet_name='T2D significant', index=False)



#################################
### Descriptive comparison ######
#################################

def standardized_mean_difference(x1, x2):
    x1 = pd.Series(x1).dropna()
    x2 = pd.Series(x2).dropna()

    mean1, mean2 = x1.mean(), x2.mean()
    sd1, sd2 = x1.std(), x2.std()

    pooled_sd = np.sqrt((sd1**2 + sd2**2) / 2)

    if pooled_sd == 0:
        return np.nan

    return (mean1 - mean2) / pooled_sd


def standardized_difference_binary(g1, g2):
    g1 = pd.Series(g1).dropna()
    g2 = pd.Series(g2).dropna()

    p1 = (g1 == 1).mean()
    p2 = (g2 == 1).mean()

    denom = np.sqrt((p1*(1-p1) + p2*(1-p2)) / 2)

    if denom == 0:
        return np.nan

    return (p1 - p2) / denom



has_uacr = Overall_matrix_T2D_backup['UACR'].notna()
has_egfr = Overall_matrix_T2D_backup['eGFR_2021'].notna()


results = []

continuous_vars = [
    'Age at recruitment',
    'Systolic blood pressure',
    'eGFR_2021',
    'BMI',
    'hba1c',
    'UACR',
    'Creatinine',
]

for var in continuous_vars:

    valid_mask = Overall_matrix_T2D_backup[var].notna()

    g1 = Overall_matrix_T2D_backup.loc[has_uacr & valid_mask, var]
    g2 = Overall_matrix_T2D_backup.loc[has_egfr & valid_mask, var]

    smd = standardized_mean_difference(g1, g2)

    results.append({
        'Variable': var,
        'UACR group': f"{g1.mean():.2f} ± {g1.std():.2f}",
        'eGFR group': f"{g2.mean():.2f} ± {g2.std():.2f}",
        'SMD': smd,
        'N UACR group': g1.count(),
        'N eGFR group': g2.count()
    })
    
    

sex_mask = Overall_matrix_T2D_backup['Sex'].notna()

sex_g1 = Overall_matrix_T2D_backup.loc[has_uacr & sex_mask, 'Sex']
sex_g2 = Overall_matrix_T2D_backup.loc[has_egfr & sex_mask, 'Sex']

n_female_uacr = (sex_g1 == 0).sum()
n_male_uacr   = (sex_g1 == 1).sum()

n_female_egfr = (sex_g2 == 0).sum()
n_male_egfr   = (sex_g2 == 1).sum()

smd_sex = standardized_difference_binary(sex_g1, sex_g2)

results.append({
    'Variable': 'Sex',
    'UACR group': f"Female: {n_female_uacr}, Male: {n_male_uacr}",
    'eGFR group': f"Female: {n_female_egfr}, Male: {n_male_egfr}",
    'SMD': smd_sex,
    'N UACR group': sex_g1.count(),
    'N eGFR group': sex_g2.count()
})

comparison_table = pd.DataFrame(results)

n_row = pd.DataFrame([{
    'Variable': 'Total N',
    'UACR group': has_uacr.sum(),
    'eGFR group': has_egfr.sum(),
    'SMD': np.nan,
    'N UACR group': has_uacr.sum(),
    'N eGFR group': has_egfr.sum()
}])
comparison_table = pd.concat([n_row, comparison_table], ignore_index=True)

print(comparison_table)


output_file = 'your/path/UACR_missingness_analysis_T2D.xlsx'

with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
    
    comparison_table.to_excel(
        writer,
        sheet_name='Descriptive_UACR_vs_eGFR',
        index=False
    )


#%% count how many points are in each quadrent 

###############
### For T2D ###
###############

# Counting points in the 1. quadrant 
count = Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_UACR'] > 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_eGFR_2021'] > 0)
].shape[0]

print(f"Number of points in the 1. quadrant {count}")

#  Counting points in the 2. quadrant 
count = Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_UACR'] > 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_eGFR_2021'] < 0)
].shape[0]

print(f"Number of points in the 2. quadrant: {count}")

# Counting points in the 3. quadrant 
count = Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_UACR'] < 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_eGFR_2021'] < 0)
].shape[0]

print(f"Number of points in the 3. quadrant: {count}")

# Counting points in the 4. quadrant 
count = Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_UACR'] < 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_eGFR_2021'] > 0)
].shape[0]

print(f"Number of points in the 4. quadrant: {count}")


###################
### For overall ###
###################

# Counting points in the 1st quadrant and their occurrences in plotting_groups 'A' and 'B'
count_1st_quadrant = Scatter_plot_eGFR_2021_associated_with_UACR_significant[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'] > 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'] > 0)
]
count_1st_quadrant_A = count_1st_quadrant[
    count_1st_quadrant['plotting_groups'] == 'A'
].shape[0]
count_1st_quadrant_B = count_1st_quadrant[
    count_1st_quadrant['plotting_groups'] == 'B'
].shape[0]
print(f"Number of points in the 1st quadrant: {count_1st_quadrant.shape[0]}")
print(f"Number of points in the 1st quadrant with 'A' in plotting_groups: {count_1st_quadrant_A}")
print(f"Number of points in the 1st quadrant with 'B' in plotting_groups: {count_1st_quadrant_B}")

# Counting points in the 2nd quadrant and their occurrences in plotting_groups 'A' and 'B'
count_2nd_quadrant = Scatter_plot_eGFR_2021_associated_with_UACR_significant[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'] > 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'] < 0)
]
count_2nd_quadrant_A = count_2nd_quadrant[
    count_2nd_quadrant['plotting_groups'] == 'A'
].shape[0]
count_2nd_quadrant_B = count_2nd_quadrant[
    count_2nd_quadrant['plotting_groups'] == 'B'
].shape[0]
print(f"Number of points in the 2nd quadrant: {count_2nd_quadrant.shape[0]}")
print(f"Number of points in the 2nd quadrant with 'A' in plotting_groups: {count_2nd_quadrant_A}")
print(f"Number of points in the 2nd quadrant with 'B' in plotting_groups: {count_2nd_quadrant_B}")

# Counting points in the 3rd quadrant and their occurrences in plotting_groups 'A' and 'B'
count_3rd_quadrant = Scatter_plot_eGFR_2021_associated_with_UACR_significant[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'] < 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'] < 0)
]
count_3rd_quadrant_A = count_3rd_quadrant[
    count_3rd_quadrant['plotting_groups'] == 'A'
].shape[0]
count_3rd_quadrant_B = count_3rd_quadrant[
    count_3rd_quadrant['plotting_groups'] == 'B'
].shape[0]
print(f"Number of points in the 3rd quadrant: {count_3rd_quadrant.shape[0]}")
print(f"Number of points in the 3rd quadrant with 'A' in plotting_groups: {count_3rd_quadrant_A}")
print(f"Number of points in the 3rd quadrant with 'B' in plotting_groups: {count_3rd_quadrant_B}")

# Counting points in the 4th quadrant and their occurrences in plotting_groups 'A' and 'B'
count_4th_quadrant = Scatter_plot_eGFR_2021_associated_with_UACR_significant[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'] < 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'] > 0)
]
count_4th_quadrant_A = count_4th_quadrant[
    count_4th_quadrant['plotting_groups'] == 'A'
].shape[0]
count_4th_quadrant_B = count_4th_quadrant[
    count_4th_quadrant['plotting_groups'] == 'B'
].shape[0]
print(f"Number of points in the 4th quadrant: {count_4th_quadrant.shape[0]}")
print(f"Number of points in the 4th quadrant with 'A' in plotting_groups: {count_4th_quadrant_A}")
print(f"Number of points in the 4th quadrant with 'B' in plotting_groups: {count_4th_quadrant_B}")

#%%
# Sort the DataFrame by both p_val_UACR and p_val_eGFR_2021 in ascending order
sorted_df = Scatter_plot_eGFR_2021_associated_with_UACR_significant.sort_values(by=['p_val_UACR', 'p_val_eGFR_2021'], ascending=True)

# Select the top 5 rows
top_5_rows = sorted_df.head(20)

# Print the values in the rsID_eGFR_2021, p_val_UACR, and p_val_eGFR_2021 columns for the top 5 rows
print(top_5_rows[['rsID_eGFR_2021', 'p_val_UACR', 'p_val_eGFR_2021']])

# Sort the DataFrame by both p_val_UACR and p_val_eGFR_2021 in ascending order
sorted_df_T2D = Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D.sort_values(by=['p_val_UACR', 'p_val_eGFR_2021'], ascending=True)

# Select the top 5 rows
top_5_rows = sorted_df_T2D.head(10)

# Print the values in the rsID_eGFR_2021, p_val_UACR, and p_val_eGFR_2021 columns for the top 5 rows
print(top_5_rows[['rsID_eGFR_2021', 'p_val_UACR', 'p_val_eGFR_2021']])

#%%


# Save to a specific folder 
output_path = 'your/path/Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D.xlsx'

# Export to Excel
Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D.to_excel(output_path, index=False, header=True)

# Save to a specific folder (replace the path with your desired folder)
output_path2 = 'your/path/Scatter_plot_eGFR_2021_associated_with_UACR_significant.xlsx'

# Export to Excel
Scatter_plot_eGFR_2021_associated_with_UACR_significant.to_excel(output_path2, index=False, header=True)

#%% Result table export 

###################################################
### result table significant overall population ###
###################################################

# Merge dataframes on 'rsID' (or the equivalent shared identifier)
result_table = pd.merge(
    lookup,
    Scatter_plot_eGFR_2021_associated_with_UACR_significant,
    left_on='base_rsID',
    right_on='rsID_eGFR_2021',
    how='inner'  # Use 'inner' to keep only matching rows
)

# Select and rename columns for the result table
result_table3 = result_table[[
    'base_rsID',            
    'chr_x',
    'Gene_x',
    'EA_x',             
    'OA',               
    'EAF',             
    'Beta_eGFR_2021',  
    'p_val_eGFR_2021', 
    'eGFR CI Lower',
    'eGFR CI Upper',
    'Beta_UACR',       
    'p_val_UACR',
    'FDR_UACR',        
    'UACR CI Lower',
    'UACR CI Upper'
]].rename(columns={
    'base_rsID': 'rsID',
    'chr_x': 'Chr',
    'Gene_x': 'Gene',
    'EA_x': 'Effect Allele',
    'OA': 'Other Allele',
    'EAF': 'EAF',
    'Beta_eGFR_2021': 'Beta ln(eGFR)',
    'p_val_eGFR_2021': 'P-value eGFR',
    'Beta_UACR': 'Beta ln(UACR)',
    'p_val_UACR': 'P-value UACR',
    'FDR_UACR': 'FDR UACR',
    'eGFR CI Lower': 'CI Lower eGFR',
    'eGFR CI Upper': 'CI Upper eGFR',
    'UACR CI Lower': 'CI Lower UACR',
    'UACR CI Upper': 'CI Upper UACR'

})

result_table3 = result_table3.drop_duplicates(subset='rsID', keep='first')

#############################################
### result table MLR 1 overall population ###
#############################################
# Extract the clean rsID from the MLR_eGFR_2021_all DataFrame
MLR_eGFR_2021_all['clean_rsID'] = MLR_eGFR_2021_all['rsID'].str.split('_').str[0]

# Merge dataframes on 'rsID' (or the equivalent shared identifier)
result_table = pd.merge(
    lookup,
    MLR_eGFR_2021_all,
    left_on='base_rsID',
    right_on='clean_rsID',
    how='inner'  # Use 'inner' to keep only matching rows
)

# Select and rename columns for the result table
result_table1 = result_table[[
    'base_rsID',            
    'chr',
    'Gene',
    'EA',               
    'OA',             
    'EAF',             
    'Effect Size',   
    'P-value', 
    'CI Lower',
    'CI Upper'
]].rename(columns={
    'base_rsID': 'rsID',
    'chr': 'Chr',
    'EA': 'Effect Allele',
    'OA': 'Other Allele',
    'EAF': 'EAF',
    'Effect Size': 'Beta ln(eGFR)',
    'P-value': 'P-value eGFR',
    'CI Upper': 'CI Upper eGFR',
    'CI Lower': 'CI Lower eGFR'

})
result_table1 = result_table1.drop_duplicates(subset='rsID', keep='first')


#############################################
### result table MLR 2 overall population ###
#############################################

# Merge dataframes on 'rsID' (or the equivalent shared identifier)
result_table = pd.merge(
    lookup,
    Scatter_plot_eGFR_2021_associated_with_UACR,
    left_on='base_rsID',
    right_on='rsID_eGFR_2021',
    how='inner'  # Use 'inner' to keep only matching rows
)


# Select and rename columns for the result table
result_table2 = result_table[[
    'base_rsID',            
    'chr_x',
    'Gene_x',
    'EA_x',            
    'OA',             
    'EAF',             
    'Beta_eGFR_2021',   
    'p_val_eGFR_2021', 
    'eGFR CI Lower',
    'eGFR CI Upper',
    'Beta_UACR',       
    'p_val_UACR',
    'FDR_UACR',       
    'UACR CI Lower',
    'UACR CI Upper'
]].rename(columns={
    'base_rsID': 'rsID',
    'chr_x': 'Chr',
    'Gene_x': 'Gene',
    'EA_x': 'Effect Allele',
    'OA': 'Other Allele',
    'EAF': 'EAF',
    'Beta_eGFR_2021': 'Beta ln(eGFR)',
    'p_val_eGFR_2021': 'P-value eGFR',
    'Beta_UACR': 'Beta ln(UACR)',
    'p_val_UACR': 'P-value UACR',
    'FDR_UACR': 'FDR UACR',
    'eGFR CI Lower': 'CI Lower eGFR',
    'eGFR CI Upper': 'CI Upper eGFR',
    'UACR CI Lower': 'CI Lower UACR',
    'UACR CI Upper': 'CI Upper UACR'
    
})
    
result_table2 = result_table2.drop_duplicates(subset='rsID', keep='first')


# Export the final table to an Excel file
output_file = 'your/path/Overall Population Results.xlsx'

with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    result_table1.to_excel(writer, sheet_name='MLR 1', index=False)
    result_table2.to_excel(writer, sheet_name='MLR 2', index=False)
    result_table3.to_excel(writer, sheet_name='MLR 2 Significant', index=False)


###############################################
### result table significant T2D population ###
###############################################
result_table = pd.merge(
    lookup,
    Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D,
    left_on='base_rsID',
    right_on='rsID_eGFR_2021',
    how='inner'  # Use 'inner' to keep only matching rows
)

# Select and rename columns for the result table
result_table3 = result_table[[
    'base_rsID',          
    'chr_x',
    'Gene_x',
    'EA_x',             
    'OA',              
    'EAF',              
    'Beta_eGFR_2021',  
    'p_val_eGFR_2021',  
    'eGFR CI Lower',
    'eGFR CI Upper',
    'Beta_UACR',       
    'p_val_UACR', 
    'FDR_UACR',      
    'UACR CI Lower',
    'UACR CI Upper'
]].rename(columns={
    'base_rsID': 'rsID',
    'chr_x': 'Chr',
    'Gene_x': 'Gene',
    'EA_x': 'Effect Allele',
    'OA': 'Other Allele',
    'EAF': 'EAF',
    'Beta_eGFR_2021': 'Beta ln(eGFR)',
    'p_val_eGFR_2021': 'P-value eGFR',
    'Beta_UACR': 'Beta ln(UACR)',
    'p_val_UACR': 'P-value UACR',
    'FDR_UACR': 'FDR UACR',
    'eGFR CI Lower': 'CI Lower eGFR',
    'eGFR CI Upper': 'CI Upper eGFR',
    'UACR CI Lower': 'CI Lower UACR',
    'UACR CI Upper': 'CI Upper UACR'
    
})

result_table3 = result_table3.drop_duplicates(subset='rsID', keep='first')


#############################################
### result table MLR 1 overall population ###
#############################################
MLR_eGFR_2021_all_T2D['clean_rsID'] = MLR_eGFR_2021_all_T2D['rsID'].str.split('_').str[0]

result_table = pd.merge(
    lookup,
    MLR_eGFR_2021_all_T2D,
    left_on='base_rsID',
    right_on='clean_rsID',
    how='inner'  
)

# Select and rename columns for the result table
result_table1 = result_table[[
    'base_rsID',          
    'chr',
    'Gene',
    'EA',              
    'OA',             
    'EAF',            
    'Effect Size',  
    'P-value', 
    'CI Lower',
    'CI Upper'
]].rename(columns={
    'base_rsID': 'rsID',
    'chr': 'Chr',
    'EA': 'Effect Allele',
    'OA': 'Other Allele',
    'EAF': 'EAF',
    'Effect Size': 'Beta ln(eGFR)',
    'P-value': 'P-value eGFR',
    'CI Lower': 'CI Lower eGFR',
    'CI Upper': 'CI Upper eGFR'

})

result_table1 = result_table1.drop_duplicates(subset='rsID', keep='first')


#############################################
### result table MLR 2 overall population ###
#############################################

# Merge dataframes on 'rsID' (or the equivalent shared identifier)
result_table = pd.merge(
    lookup,
    Scatter_plot_eGFR_2021_associated_with_UACR_T2D,
    left_on='base_rsID',
    right_on='rsID_eGFR_2021',
    how='inner'  # Use 'inner' to keep only matching rows
)

# Select and rename columns for the result table
result_table2 = result_table[[
    'base_rsID',             
    'chr_x',
    'Gene_x',
    'EA_x',              
    'OA',            
    'EAF',              
    'Beta_eGFR_2021',   
    'p_val_eGFR_2021',  
    'eGFR CI Lower',
    'eGFR CI Upper',
    'Beta_UACR',        
    'p_val_UACR',
    'FDR_UACR',        
    'UACR CI Lower',
    'UACR CI Upper'
]].rename(columns={
    'base_rsID': 'rsID',
    'chr_x': 'Chr',
    'Gene_x': 'Gene',
    'EA_x': 'Effect Allele',
    'OA': 'Other Allele',
    'EAF': 'EAF',
    'Beta_eGFR_2021': 'Beta ln(eGFR)',
    'p_val_eGFR_2021': 'P-value eGFR',
    'Beta_UACR': 'Beta ln(UACR)',
    'p_val_UACR': 'P-value UACR',
    'FDR_UACR': 'FDR UACR',
    'eGFR CI Lower': 'CI Lower eGFR',
    'eGFR CI Upper': 'CI Upper eGFR',
    'UACR CI Lower': 'CI Lower UACR',
    'UACR CI Upper': 'CI Upper UACR'
    
})

# Export the final table to an Excel file
output_file = 'your/path/T2D Results.xlsx'

with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    result_table1.to_excel(writer, sheet_name='MLR 1', index=False)
    result_table2.to_excel(writer, sheet_name='MLR 2', index=False)
    result_table3.to_excel(writer, sheet_name='MLR 2 Significant', index=False)







