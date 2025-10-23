
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




#%% Load eGFR Data
#########################
#### Load eGFR Data #####
#########################


# Load Excel file
df = pd.read_excel('file/path/eGFR.xlsx')

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

df_Genomic_alb = pd.read_excel('file/path/Data_UKBB_alb.xlsx', sheet_name='Ark1')

# Select columns starting with "rs"
rs_columns = [col for col in df_Genomic_alb.columns if col.startswith('rs')]

# round the values 
df_Genomic_alb[rs_columns] = df_Genomic_alb[rs_columns].apply(np.round)


##########################
### Load the eGFR data ###
##########################


df_Genomic_eGFR = pd.read_excel('file/path/Data_ukbb_eGFR.xlsx', sheet_name='Ark1')

# Select columns starting with "rs"
rs_columns = [col for col in df_Genomic_eGFR.columns if col.startswith('rs')]

# round the values 
df_Genomic_eGFR[rs_columns] = df_Genomic_eGFR[rs_columns].apply(np.round)




###############################
### Load the phenotype data ###
###############################


df_phenotype = pd.read_excel('file/path/phenotype_data_alb_eGFR.xlsx', sheet_name='phenotype_data2', decimal=',')


import pickle
with open('variables.pkl', 'wb') as f:
    pickle.dump((df_phenotype, df_Genomic_eGFR, df_Genomic_alb), f)

#%% Getting the UKBB imported data fast + creating 2 matresis with the phenotype and genotype data
from openpyxl import load_workbook
import pickle


# Load variables
with open('file/path/variables.pkl', 'rb') as f:
    df_phenotype, df_Genomic_eGFR, df_Genomic_alb = pickle.load(f)


# change pickle file 
df_Genomic_eGFR.drop(columns=['rs1004441_A.1'], inplace=True)


ethnicity = pd.read_excel('file/path/ethnicity.xlsx', sheet_name='ethnicity', decimal=',')

################
### Withdraw ###
################


# Load the participants who have withdrawn
withdraw = pd.read_csv('file/path/withdraw.csv', header=None)

# Rename the column to 'Eid'
withdraw.columns = ['Eid']



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
output_file = 'file/path/EAF.xlsx'
df_allele_fq.to_excel(output_file, index=False)



#%% Preprocessing data 

##########################
### Preprocessing data ###
##########################


######################
### Calculate UACR ###
######################
creatinine_data = pd.read_csv('file/path/creatinine_data.csv', header=None)
creatinine_data.columns = ['eid', 'creatinine_urine']
creatinine_data = creatinine_data.iloc[1:].reset_index(drop=True)

overall_copy = Overall_matrix.copy()
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

# Remove rows with NaN values
new_matrix.dropna(inplace=True)

new_matrix.reset_index(drop=True, inplace=True)

# Compute mean and std
mean_egfr_overall = new_matrix['eGFR_2021'].mean()
std_egfr_overall = new_matrix['eGFR_2021'].std()

y_log = new_matrix['eGFR_2021']

# Define the independent variables (X) as the phenotype data and the rsIDs
X = new_matrix[['Sex', 'Age at recruitment','Systolic blood pressure']]  # Selecting Age and Sex SBP

# Perform individual regressions for each rsID
rs_columns = [col for col in new_matrix.columns if col.startswith('rs')]
results = []
#bonferroni_threshold = 0.05 / len(rs_columns)  # Bonferroni corrected threshold
#import scipy.stats as stats
#import matplotlib.pyplot as plt

for rsID in rs_columns:
    # Add the current rsID column to X
    X_with_rsID = X.join(new_matrix[rsID])

    # Add a constant column to X_with_rsID (for the intercept)
    X_with_rsID = sm.add_constant(X_with_rsID)

    y_log.reset_index(drop=True, inplace=True)

    # Fit the regression model
    model = sm.OLS(y_log, X_with_rsID)
    result = model.fit()
    # Only plot Q-Q plot for specific SNPs (so you don't get 1000s of plots)
    #if rsID in ['rs2433601_C_eGFR', 'rs2412608_T_eGFR']:
    #    residuals = result.resid
        
    #    stats.probplot(residuals, dist="norm", plot=plt)
    #    plt.title(f"Q-Q plot for residuals of {rsID}")
    #    plt.xlabel("Theoretical Quantiles")
    #    plt.ylabel("Sample Quantiles")
    #    plt.grid(True, linestyle='--', linewidth=0.5)
    #    plt.tight_layout()
    #    plt.show()
    # Standard p-value from statsmodels (may round to 0 if very small)
    p_value = result.pvalues[rsID]
    
    # If it's 0, recompute with high precision from t-statistic
    if p_value == 0.0:
        from scipy.stats import t
        t_stat = result.tvalues[rsID]
        df_resid = result.df_resid
        p_value = 2 * t.sf(np.abs(t_stat), df=df_resid)
        print(f"[INFO] Corrected p-value for {rsID}")
        print(f"  t-statistic: {t_stat:.10f}")
        print(f"  Degrees of freedom: {df_resid}")
        print(f"  Recomputed p-value: {p_value:.10e}")

    conf_int = result.conf_int().loc[rsID]  # Confidence interval for this rsID


    results.append((rsID, result.params[rsID], p_value, conf_int[0], conf_int[1]))

# Create a DataFrame to store the results
MLR_eGFR_2021_all = pd.DataFrame(results, columns=['rsID', 'Effect Size', 'P-value', 'CI Lower', 'CI Upper'])

# Filter non-significant rsIDs based on Bonferroni correction
#MLR_eGFR_2021_significant_results_bonferroni_all = MLR_eGFR_2021_all[MLR_eGFR_2021_all['P-value'] <= bonferroni_threshold]

# Filter non-significant rsIDs based on a significance level of 0.05
MLR_eGFR_2021_significant_results_all = MLR_eGFR_2021_all[MLR_eGFR_2021_all['P-value'] <= 0.05]



variables_to_delete = ['model','p_value','result','results','rs_columns','new_matrix','df_selected_columns', 't_stat',
                       'rsID','X','X_with_rsID','y_log','columns_to_select','conf_int', 'df_resid']

for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete






#%% MLR 2


###########################################
### eGFR 2021 rsID associated with UACR ### 
###########################################

# including micro albumin

alb = Overall_matrix_backup[['eid', 'UACR']]
alb=alb.dropna()


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

new_matrix = Overall_matrix_stand[rsIDs_to_select + ['eid', 'Sex', 'Age at recruitment','Systolic blood pressure']].copy()

new_matrix = pd.merge(new_matrix, alb_filtered[['eid','UACR']], on='eid', how='inner')



# Remove rows with NaN values
new_matrix.dropna(inplace=True)

new_matrix.reset_index(drop=True, inplace=True)


# summary statistics

selected_columns = [
    'UACR'

]


subset = new_matrix[selected_columns]

subset = subset.dropna()


Summary_overall2 = subset.describe(include='all')


Summary_combined = pd.concat([Summary_overall, Summary_overall2], axis=1)

new_matrix['UACR'] = np.log(new_matrix['UACR'])

# Compute mean and std
mean_uae_overall = new_matrix['UACR'].mean()
std_uae_overall = new_matrix['UACR'].std()



y_log = new_matrix['UACR']

# Define the independent variables (X) as the phenotype data and the rsIDs
X = new_matrix[['Sex', 'Age at recruitment','Systolic blood pressure']]  # Selecting Age and Sex


# Perform individual regressions for each rsID
rs_columns = [col for col in new_matrix.columns if col.startswith('rs')]
results = []
bonferroni_threshold_eGFR_2021_associated_with_UACR = 0.05 / len(rs_columns)  # Bonferroni corrected threshold

for rsID in rs_columns:
    # Add the current rsID column to X
    X_with_rsID = X.join(new_matrix[rsID])

    # Add a constant column to X_with_rsID (for the intercept)
    X_with_rsID = sm.add_constant(X_with_rsID)

    y_log.reset_index(drop=True, inplace=True)

    # Fit the regression model
    model = sm.OLS(y_log, X_with_rsID)
    result = model.fit()
    p_value = result.pvalues[rsID]
    conf_int = result.conf_int().loc[rsID]  # Confidence interval for this rsID


    results.append((rsID, result.params[rsID], p_value, conf_int[0], conf_int[1]))

# Create a DataFrame to store the results
eGFR_2021_rsid_associated_with_UACR_model2 = pd.DataFrame(results, columns=['rsID', 'Effect Size', 'P-value','CI Lower', 'CI Upper'])

# Filter non-significant rsIDs based on Bonferroni correction
eGFR_2021_rsid_associated_with_UACR_bonferroni_model2 = eGFR_2021_rsid_associated_with_UACR_model2[eGFR_2021_rsid_associated_with_UACR_model2['P-value'] <= bonferroni_threshold_eGFR_2021_associated_with_UACR]

# Filter non-significant rsIDs based on a significance level of 0.05
eGFR_2021_rsid_associated_with_UACR_significant_model2 = eGFR_2021_rsid_associated_with_UACR_model2[eGFR_2021_rsid_associated_with_UACR_model2['P-value'] <= 0.05]




variables_to_delete = ['conf_int','eid_column', 'lower_bounds','mean_values','model', 'new_matrix', 'upper_bounds', 'X',
                       'outlier_rows', 'p_value', 'phenotype_columns', 'result', 'rsID','y_log', 'X_with_rsID', 'Summary_overall2',
                       'subset','std_values','selected_columns','rsIDs_to_select','results','rs_columns', 'Summary_overall']

for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete



#%% Preprocessing - creating matrix for scatterplots


##############################################################
### eGFR_2021_rsid_associated_with_UACR_significant_model2 ###
##############################################################



# Load Excel file
lookup = pd.read_excel('file/path/lookup.xlsx')


# Step 1: Extract base rsID and EA from rsID
split_rsID = eGFR_2021_rsid_associated_with_UACR_model2['rsID'].str.split('_', expand=True)
eGFR_2021_rsid_associated_with_UACR_model2['base_rsID'] = split_rsID[0]
eGFR_2021_rsid_associated_with_UACR_model2['EA'] = split_rsID[1]

# Step 3: Merge on base_rsID
merged = pd.merge(
    eGFR_2021_rsid_associated_with_UACR_model2,
    lookup[['base_rsID', 'chr', 'Gene']],
    on='base_rsID',
    how='left'
)

# Step 4: Rename and format final output
MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2 = merged.rename(columns={
    'base_rsID': 'rsID_eGFR_2021',
    'Effect Size': 'Beta_UACR',
    'P-value': 'p_val_UACR',
    'CI Lower': 'UACR CI Lower',
    'CI Upper': 'UACR CI Upper'
})[[
    'rsID_eGFR_2021', 'chr', 'Gene', 'EA',
    'Beta_UACR', 'p_val_UACR', 'UACR CI Lower', 'UACR CI Upper'
]]



#########################
### UACR X eGFR 2021  ### eGFR
#########################

Scatter_plot_eGFR_2021_associated_with_UACR = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2.copy()

# Step 1: Extract base rsID from MLR1 results
MLR_eGFR_2021_all['rsID_prefix'] = MLR_eGFR_2021_all['rsID'].str.split('_').str[0]

# Step 2: Merge MLR1 results into scatter plot matrix
Scatter_plot_eGFR_2021_associated_with_UACR = Scatter_plot_eGFR_2021_associated_with_UACR.merge(
    MLR_eGFR_2021_all[['rsID_prefix', 'Effect Size', 'P-value', 'CI Lower', 'CI Upper']],
    left_on='rsID_eGFR_2021',
    right_on='rsID_prefix',
    how='left'
)

# Step 3: Rename merged columns
Scatter_plot_eGFR_2021_associated_with_UACR.rename(columns={
    'Effect Size': 'Beta_eGFR_2021',
    'P-value': 'p_val_eGFR_2021',
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
####  Scatterplot eGFR 2021 rsID association with UACR ##### all without names 
############################################################


Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'] = pd.to_numeric(Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'], errors='coerce')
Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'] = pd.to_numeric(Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'], errors='coerce')

Scatter_plot_eGFR_2021_associated_with_UACR_significant['plotting_groups'] = 'A'
# Extract base rsIDs from the Bonferroni dataset
uacr_bonferroni_rsIDs = eGFR_2021_rsid_associated_with_UACR_bonferroni_model2['rsID'].str.split('_').str[0]

# Check if rsID_eGFR is in the UACR Bonferroni list and assign 'B' to plotting_groups
Scatter_plot_eGFR_2021_associated_with_UACR_significant['plotting_groups'] = np.where(
    Scatter_plot_eGFR_2021_associated_with_UACR_significant['rsID_eGFR_2021'].isin(uacr_bonferroni_rsIDs),
    'B',
    Scatter_plot_eGFR_2021_associated_with_UACR_significant.get('plotting_groups', np.nan)  # Preserve existing values if any
)



# Create the scatter plot 
g = sns.relplot(data=Scatter_plot_eGFR_2021_associated_with_UACR_significant,
                x='Beta_eGFR_2021',
                y='Beta_UACR',
                aspect=2,
                hue='plotting_groups',
                palette=['black','#26538d'],
                linewidth=0,
                s=45,
                alpha=0.7,
                legend=None)



# Add labels and title to the plot

g.ax.set_xlabel(r'$\ln(eGFR_{2021})$ beta', fontsize=14)
g.ax.set_ylabel('ln(UACR) beta', fontsize=14)
g.ax.axhline(0, color='black', lw=0.5)
g.ax.axvline(0, color='black', lw=0.5)

# Set x and y axis limits
g.ax.set_xlim(-0.013, 0.013)
g.ax.set_ylim(-0.08, 0.05)
g.ax.grid(True, which='both', linestyle='--', linewidth=0.5)

# Define legend labels and colors
legend_labels = ['Significance level: 0.05', 'Bonferroni corrected: 0.000137']
legend_markers = ['o', 'o']  # Specify the markers for each label
legend_colors = ['black', '#26538d']

handles = [plt.Line2D([], [], color=color, marker=marker, linestyle='None') for color, marker in zip(legend_colors, legend_markers)]


g.ax.legend(handles, legend_labels, loc='upper right', markerscale=1, fontsize=12, frameon=True, edgecolor='black')

# save the plot
save_path = 'file/path/scatterplot_eGFR_2021_associated_with_UACR_a.jpg'
plt.savefig(save_path, dpi=500, bbox_inches='tight')

plt.show()


############################################################ 
####  Scatterplot eGFR 2021 rsID association with UACR ##### 
############################################################

Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'] = pd.to_numeric(Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'], errors='coerce')
Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'] = pd.to_numeric(Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'], errors='coerce')


Scatter_plot_eGFR_2021_associated_with_UACR_significant['plotting_groups'] = 'A'

# Iterate over the rows and check for matching rsIDs
# Extract base rsIDs from the Bonferroni dataset
uacr_bonferroni_rsIDs = eGFR_2021_rsid_associated_with_UACR_bonferroni_model2['rsID'].str.split('_').str[0]

# Assign 'B' to rows where rsID_eGFR matches a Bonferroni rsID
Scatter_plot_eGFR_2021_associated_with_UACR_significant['plotting_groups'] = np.where(
    Scatter_plot_eGFR_2021_associated_with_UACR_significant['rsID_eGFR_2021'].isin(uacr_bonferroni_rsIDs),
    'B',
    Scatter_plot_eGFR_2021_associated_with_UACR_significant.get('plotting_groups', 'A')  # Default to 'A'
)
filtered_df = Scatter_plot_eGFR_2021_associated_with_UACR_significant[Scatter_plot_eGFR_2021_associated_with_UACR_significant['plotting_groups'] == 'B'].copy()



# Create the scatter plot
g = sns.relplot(data=filtered_df,
                x='Beta_eGFR_2021',
                y='Beta_UACR',
                aspect=2,
                hue='plotting_groups',
                palette=['#26538d'],
                linewidth=0,
                s=45,
                alpha=0.8,
                legend=None)



# Add labels and title to the plot
g.ax.set_xlabel(r'$\ln(eGFR_{2021})$ beta', fontsize=14)
g.ax.set_ylabel('ln(UACR) beta', fontsize=14)
g.ax.axhline(0, color='black', lw=0.5)
g.ax.axvline(0, color='black', lw=0.5)
g.ax.grid(True, which='both', linestyle='--', linewidth=0.5)
# Set x and y axis limits
g.ax.set_xlim(-0.01, 0.01)
g.ax.set_ylim(-0.08, 0.05)


# Create the annotations
Anno = filtered_df.apply(
    lambda p: g.ax.annotate(p['Gene'], (p['Beta_eGFR_2021'], p['Beta_UACR']), fontsize=11),
    axis=1
).to_list()

adjust_text(Anno, arrowprops=dict(arrowstyle="-", color='k', lw=0.4), force_points=0.1, force_text=0.1)




# Define legend labels and colors
legend_labels = ['Bonferroni corrected: 0.000136']
legend_markers = ['o']  # Specify the markers for each label
legend_colors = ['#26538d']

handles = [plt.Line2D([], [], color=color, marker=marker, linestyle='None') for color, marker in zip(legend_colors, legend_markers)]


g.ax.legend(handles, legend_labels, loc='upper right', markerscale=1, fontsize=12, frameon=True, edgecolor='black')

# save the plot
save_path = 'file/path/scatterplot_eGFR_2021_associated_with_UACR_b.jpg'
plt.savefig(save_path, dpi=500, bbox_inches='tight')


plt.show()





variables_to_delete = ['uacr_bonferroni_rsIDs','split_rsID','save_path','merged','legend_markers','legend_labels', 'legend_colors',
                       'handles', 'g', 'Anno']

for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete


#%% Sensitivity analysis 
############################################
### Sensitivity analysis: add eGFR_2021 ####
### Only for Bonferroni-significant SNPs ###
############################################


rsIDs_to_select = eGFR_2021_rsid_associated_with_UACR_bonferroni_model2['rsID'].tolist()

new_matrix = Overall_matrix_stand[rsIDs_to_select + ['eid', 'Sex', 'Age at recruitment','Systolic blood pressure', 'eGFR_2021']].copy()

new_matrix = pd.merge(new_matrix, alb_filtered[['eid','UACR']], on='eid', how='inner')


# Remove rows with NaN values
new_matrix.dropna(inplace=True)
new_matrix.reset_index(drop=True, inplace=True)


new_matrix['UACR'] = np.log(new_matrix['UACR'])

# Scale eGFR
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
new_matrix['eGFR_2021_std'] = scaler.fit_transform(new_matrix[['eGFR_2021']])

# Define outcome
y_log = new_matrix['UACR']

# Define covariates INCLUDING standardized eGFR
X = new_matrix[['Sex', 'Age at recruitment','Systolic blood pressure','eGFR_2021_std']]

# Run regressions for each Bonferroni-significant SNP
rs_columns = [col for col in new_matrix.columns if col.startswith('rs')]
results = []

for rsID in rs_columns:
    X_with_rsID = sm.add_constant(X.join(new_matrix[rsID]))
    model = sm.OLS(y_log, X_with_rsID)
    result = model.fit()

    beta = result.params[rsID]
    p_value = result.pvalues[rsID]
    ci_lower, ci_upper = result.conf_int().loc[rsID]

    results.append((rsID, beta, p_value, ci_lower, ci_upper))

# Store results
UACR_sensitivity_model1 = pd.DataFrame(results, columns=['rsID', 'Effect Size', 'P-value', 'CI Lower', 'CI Upper'])



UACR_sensitivity_model1_significant = UACR_sensitivity_model1[UACR_sensitivity_model1['P-value'] <= 0.05]




#%% Doing it only for diabetes 

#####################################
### creating diabetes matrix data ###
#####################################

# Load the diabetes csv data using the Pandas library
df_T2D = pd.read_csv('file/path/df_preprocessed.csv')

# Select only the 'eid' and 'T2D' columns
df_T2D = df_T2D[['eid', 'T2D']]

Overall_matrix_T2D = pd.merge(Overall_matrix, df_T2D, on='eid', how='inner')

# Filter out rows where 'Eid' in Overall_matrix matches 'Eid' in withdraw
Overall_matrix_T2D = Overall_matrix_T2D[~Overall_matrix_T2D['eid'].isin(withdraw['Eid'])]

# only the people with T2D = 1

Overall_matrix_T2D = Overall_matrix_T2D[Overall_matrix_T2D['T2D'] == 1]

del df_T2D


##########################
### Preprocessing data ###
##########################


eid_column = Overall_matrix_T2D['eid']

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
Overall_matrix_T2D_filtered = Overall_matrix_T2D[~outlier_rows].copy()

Overall_matrix_T2D_filtered.loc[:, 'eid'] = eid_column.loc[~outlier_rows].values

Overall_matrix_T2D_backup = Overall_matrix_T2D


# summary statistics 


selected_columns = [
    'Sex',
    'Systolic blood pressure',
    'Creatinine',
    'Age at recruitment',
    'eGFR_2021'
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




# Convert UACR to ln(UACR) and eGFR to ln(eGFR)
Overall_matrix_T2D_filtered['eGFR_2021'] = np.log(Overall_matrix_T2D_filtered['eGFR_2021'])




columns_to_scale = ['Sex', 'Systolic blood pressure', 'Age at recruitment']

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


#####################
### for all pheno ###
#####################

# Select the columns with '_eGFR' at the end and include 'eid', 'Sex', 'age' and SBP
columns_to_select = [col for col in Overall_matrix_stand_T2D.columns if col.endswith('_eGFR')] + ['eid', 'Sex', 'Age at recruitment','Systolic blood pressure','eGFR_2021']


# Create a new DataFrame with the selected columns
new_matrix = Overall_matrix_stand_T2D[columns_to_select].copy()

# Remove rows with NaN values
new_matrix.dropna(inplace=True)

new_matrix.reset_index(drop=True, inplace=True)

# Compute mean and std
mean_egfr_T2D = new_matrix['eGFR_2021'].mean()
std_egfr_T2D = new_matrix['eGFR_2021'].std()



y_log = new_matrix['eGFR_2021']

# Define the independent variables (X) as the phenotype data and the rsIDs
X = new_matrix[['Sex', 'Age at recruitment','Systolic blood pressure']]  # Selecting Age and Sex ans SBP

# Perform individual regressions for each rsID
rs_columns = [col for col in new_matrix.columns if col.startswith('rs')]
results = []
#bonferroni_threshold_T2D = 0.05 / len(rs_columns)  # Bonferroni corrected threshold

for rsID in rs_columns:
    # Add the current rsID column to X
    X_with_rsID = X.join(new_matrix[rsID])

    # Add a constant column to X_with_rsID (for the intercept)
    X_with_rsID = sm.add_constant(X_with_rsID)

    y_log.reset_index(drop=True, inplace=True)

    # Fit the regression model
    model = sm.OLS(y_log, X_with_rsID)
    result = model.fit()
    p_value = result.pvalues[rsID]
    conf_int = result.conf_int().loc[rsID]  # Confidence interval for this rsID


    results.append((rsID, result.params[rsID], p_value, conf_int[0], conf_int[1]))

# Create a DataFrame to store the results
MLR_eGFR_2021_all_T2D = pd.DataFrame(results, columns=['rsID', 'Effect Size', 'P-value','CI Lower', 'CI Upper'])

# Filter non-significant rsIDs based on Bonferroni correction
#MLR_eGFR_2021_significant_results_bonferroni_all_T2D = MLR_eGFR_2021_all_T2D[MLR_eGFR_2021_all_T2D['P-value'] <= bonferroni_threshold_T2D]

# Filter non-significant rsIDs based on a significance level of 0.05
MLR_eGFR_2021_significant_results_all_T2D = MLR_eGFR_2021_all_T2D[MLR_eGFR_2021_all_T2D['P-value'] <= 0.05]



variables_to_delete = ['model','p_value','result','results','rs_columns', 'columns_to_select', 'conf_int',
                       'rsID','X','X_with_rsID','y_log']

for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete

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

new_matrix = Overall_matrix_stand_T2D[rsIDs_to_select + ['eid', 'Sex', 'Age at recruitment','Systolic blood pressure']].copy()


new_matrix = pd.merge(new_matrix, alb_filtered_T2D[['eid','UACR']], on='eid', how='inner')




# Remove rows with NaN values
new_matrix.dropna(inplace=True)

new_matrix.reset_index(drop=True, inplace=True)


# summary statistics


selected_columns = [
    'UACR'

]

subset = new_matrix[selected_columns]

subset = subset.dropna()


Summary_T2D_2 = subset.describe(include='all')


Summary_combined_T2D = pd.concat([Summary_T2D, Summary_T2D_2], axis=1)

del Summary_T2D
del Summary_T2D_2




output_file = 'file/path/Summary.xlsx'

with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    Summary_combined.to_excel(writer, sheet_name='Overall')
    Summary_combined_T2D.to_excel(writer, sheet_name='T2D')
    
del output_file

new_matrix['UACR'] = np.log(new_matrix['UACR'])

# Compute mean and std
mean_uae_T2D = new_matrix['UACR'].mean()
std_uae_T2D = new_matrix['UACR'].std()


y_log = new_matrix['UACR']

# Define the independent variables (X) as the phenotype data and the rsIDs
X = new_matrix[['Sex', 'Age at recruitment','Systolic blood pressure']]  # Selecting Age and Sex


# Perform individual regressions for each rsID
rs_columns = [col for col in new_matrix.columns if col.startswith('rs')]
results = []
bonferroni_threshold_eGFR_2021_associated_with_UACR_T2D = 0.05 / len(rs_columns)  # Bonferroni corrected threshold

for rsID in rs_columns:
    # Add the current rsID column to X
    X_with_rsID = X.join(new_matrix[rsID])

    # Add a constant column to X_with_rsID (for the intercept)
    X_with_rsID = sm.add_constant(X_with_rsID)

    y_log.reset_index(drop=True, inplace=True)

    # Fit the regression model
    model = sm.OLS(y_log, X_with_rsID)
    result = model.fit()
    p_value = result.pvalues[rsID]
    conf_int = result.conf_int().loc[rsID]  # Confidence interval for this rsID


    results.append((rsID, result.params[rsID], p_value, conf_int[0], conf_int[1]))

# Create a DataFrame to store the results
eGFR_2021_rsid_associated_with_UACR_model2_T2D = pd.DataFrame(results, columns=['rsID', 'Effect Size', 'P-value', 'CI Lower', 'CI Upper'])

# Filter non-significant rsIDs based on Bonferroni correction
eGFR_2021_rsid_associated_with_UACR_bonferroni_model2_T2D = eGFR_2021_rsid_associated_with_UACR_model2_T2D[eGFR_2021_rsid_associated_with_UACR_model2_T2D['P-value'] <= bonferroni_threshold_eGFR_2021_associated_with_UACR_T2D]

# Filter non-significant rsIDs based on a significance level of 0.05
eGFR_2021_rsid_associated_with_UACR_significant_model2_T2D = eGFR_2021_rsid_associated_with_UACR_model2_T2D[eGFR_2021_rsid_associated_with_UACR_model2_T2D['P-value'] <= 0.05]





variables_to_delete = ['model','p_value','result','results','rs_columns', 'conf_int',
                       'rsID','X','X_with_rsID','y_log', 'rsIDs_to_select', 'subset', 'selected_columns', 'writer']

for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete



##############################################################
### eGFR_2021_rsid_associated_with_UACR_significant_model2_T2D ###
##############################################################
# Step 1: Split rsID into base ID and EA
split_rsID = eGFR_2021_rsid_associated_with_UACR_model2_T2D['rsID'].str.split('_', expand=True)
eGFR_2021_rsid_associated_with_UACR_model2_T2D['base_rsID'] = split_rsID[0]
eGFR_2021_rsid_associated_with_UACR_model2_T2D['EA'] = split_rsID[1]



# Load Excel file
lookup = pd.read_excel('file/path/lookup.xlsx')



# Step 3: Merge on base_rsID
merged = pd.merge(
    eGFR_2021_rsid_associated_with_UACR_model2_T2D,
    lookup[['base_rsID', 'chr', 'Gene']],
    on='base_rsID',
    how='left'
)

# Step 4: Clean final output
MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2_T2D = merged.rename(columns={
    'base_rsID': 'rsID_eGFR_2021',
    'Effect Size': 'Beta_UACR',
    'P-value': 'p_val_UACR',
    'CI Lower': 'UACR CI Lower',
    'CI Upper': 'UACR CI Upper'
})[[
    'rsID_eGFR_2021', 'chr', 'Gene', 'EA',
    'Beta_UACR', 'p_val_UACR', 'UACR CI Lower', 'UACR CI Upper'
]]

#########################
### UACR X eGFR 2021  ### eGFR
#########################

# Step 1: Start from your UACR-associated matrix
Scatter_plot_eGFR_2021_associated_with_UACR_T2D = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2_T2D.copy()

# Step 2: Extract base rsID from eGFR MLR results
MLR_eGFR_2021_all_T2D['rsID_prefix'] = MLR_eGFR_2021_all_T2D['rsID'].str.split('_').str[0]

# Step 3: Merge on rsID prefix
Scatter_plot_eGFR_2021_associated_with_UACR_T2D = Scatter_plot_eGFR_2021_associated_with_UACR_T2D.merge(
    MLR_eGFR_2021_all_T2D[['rsID_prefix', 'Effect Size', 'P-value', 'CI Lower', 'CI Upper']],
    left_on='rsID_eGFR_2021',
    right_on='rsID_prefix',
    how='left'
)

# Step 4: Rename columns
Scatter_plot_eGFR_2021_associated_with_UACR_T2D.rename(columns={
    'Effect Size': 'Beta_eGFR_2021',
    'P-value': 'p_val_eGFR_2021',
    'CI Lower': 'eGFR CI Lower',
    'CI Upper': 'eGFR CI Upper'
}, inplace=True)

# Step 5: Drop the helper column
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
                palette=['black', '#26538d'],
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
g.ax.set_xlim(-0.014, 0.016)
g.ax.set_ylim(-0.075, 0.075)

# Create the annotations
#Anno = Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D.apply(
#    lambda p: g.ax.annotate(p['Gene'], (p['Beta_eGFR_2021'], p['Beta_UACR']), fontsize=11),
#    axis=1
#).to_list()

#adjust_text(Anno, arrowprops=dict(arrowstyle="-", color='k', lw=0.4), force_points=0.1, force_text=0.2)





# Define legend labels and colors
legend_labels = ['Significance level: 0.05',
                 'Bonferroni corrected: 0.000331']
legend_markers = ['o','o']  # Specify the markers for each label
legend_colors = ['black','#26538d']

handles = [plt.Line2D([], [], color=color, marker=marker, linestyle='None') for color, marker in zip(legend_colors, legend_markers)]

g.ax.grid(True, which='both', linestyle='--', linewidth=0.5)
g.ax.legend(handles, legend_labels, loc='upper right', markerscale=1, fontsize=12, frameon=True, edgecolor='black')

# save the plot
save_path = 'file/path/scatterplot_eGFR_2021_associated_with_UACR_T2D_no_text.jpg'
plt.savefig(save_path, dpi=300, bbox_inches='tight')


plt.show()




#%% Sensitivity analysis T2D
############################################
### Sensitivity analysis: add eGFR_2021 ####
### Only for Bonferroni-significant SNPs ###
############################################


rsIDs_to_select = eGFR_2021_rsid_associated_with_UACR_bonferroni_model2_T2D['rsID'].tolist()

new_matrix = Overall_matrix_stand[rsIDs_to_select + ['eid', 'Sex', 'Age at recruitment','Systolic blood pressure', 'eGFR_2021']].copy()

new_matrix = pd.merge(new_matrix, alb_filtered[['eid','UACR']], on='eid', how='inner')


# Remove rows with NaN values
new_matrix.dropna(inplace=True)
new_matrix.reset_index(drop=True, inplace=True)


new_matrix['UACR'] = np.log(new_matrix['UACR'])

# Scale eGFR
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
new_matrix['eGFR_2021_std'] = scaler.fit_transform(new_matrix[['eGFR_2021']])

# Define outcome
y_log = new_matrix['UACR']

# Define covariates INCLUDING standardized eGFR
X = new_matrix[['Sex', 'Age at recruitment','Systolic blood pressure','eGFR_2021_std']]

# Run regressions for each Bonferroni-significant SNP
rs_columns = [col for col in new_matrix.columns if col.startswith('rs')]
results = []

for rsID in rs_columns:
    X_with_rsID = sm.add_constant(X.join(new_matrix[rsID]))
    model = sm.OLS(y_log, X_with_rsID)
    result = model.fit()

    beta = result.params[rsID]
    p_value = result.pvalues[rsID]
    ci_lower, ci_upper = result.conf_int().loc[rsID]

    results.append((rsID, beta, p_value, ci_lower, ci_upper))

# Store results
UACR_sensitivity_model1_T2D = pd.DataFrame(results, columns=['rsID', 'Effect Size', 'P-value', 'CI Lower', 'CI Upper'])



UACR_sensitivity_model1_significant_T2D = UACR_sensitivity_model1_T2D[UACR_sensitivity_model1_T2D['P-value'] <= 0.05]





# Export the final table to an Excel file
output_file = 'file/path/sensitivity_Results.xlsx'

with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
    UACR_sensitivity_model1.to_excel(writer, sheet_name='Overall', index=False)
    UACR_sensitivity_model1_significant.to_excel(writer, sheet_name='Overall significant', index=False)
    UACR_sensitivity_model1_T2D.to_excel(writer, sheet_name='T2D', index=False)
    UACR_sensitivity_model1_significant_T2D.to_excel(writer, sheet_name='T2D significant', index=False)






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

top_5_rows = sorted_df.head(20)

print(top_5_rows[['rsID_eGFR_2021', 'p_val_UACR', 'p_val_eGFR_2021']])

# Sort the DataFrame by both p_val_UACR and p_val_eGFR_2021 in ascending order
sorted_df_T2D = Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D.sort_values(by=['p_val_UACR', 'p_val_eGFR_2021'], ascending=True)


top_5_rows = sorted_df_T2D.head(10)

print(top_5_rows[['rsID_eGFR_2021', 'p_val_UACR', 'p_val_eGFR_2021']])

#%%


# Save to a specific folder 
output_path = 'file/path/Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D.xlsx'

# Export to Excel
Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D.to_excel(output_path, index=False, header=True)

# Save to a specific folder (replace the path with your desired folder)
output_path2 = 'file/path/Scatter_plot_eGFR_2021_associated_with_UACR_significant.xlsx'

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
    how='inner' 
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
    'Beta_UACR': 'Beta ln(UAE)',
    'p_val_UACR': 'P-value UAE',
    'eGFR CI Lower': 'CI Lower eGFR',
    'eGFR CI Upper': 'CI Upper eGFR',
    'UACR CI Lower': 'CI Lower UAE',
    'UACR CI Upper': 'CI Upper UAE'

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
    'CI Upper': 'CI Upper eGFR',
    'CI Lower': 'CI Lower eGFR'

})
result_table1 = result_table1.drop_duplicates(subset='rsID', keep='first')


#############################################
### result table MLR 2 overall population ###
#############################################

# Merge dataframes on 'rsID' 
result_table = pd.merge(
    lookup,
    Scatter_plot_eGFR_2021_associated_with_UACR,
    left_on='base_rsID',
    right_on='rsID_eGFR_2021',
    how='inner'  
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
    'Beta_UACR': 'Beta ln(UAE)',
    'p_val_UACR': 'P-value UAE',
    'eGFR CI Lower': 'CI Lower eGFR',
    'eGFR CI Upper': 'CI Upper eGFR',
    'UACR CI Lower': 'CI Lower UAE',
    'UACR CI Upper': 'CI Upper UAE'
    
})
    
result_table2 = result_table2.drop_duplicates(subset='rsID', keep='first')


# Export the final table to an Excel file
output_file = 'file/path/Overall Population Results.xlsx'

with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
    result_table1.to_excel(writer, sheet_name='MLR 1', index=False)
    result_table2.to_excel(writer, sheet_name='MLR 2', index=False)
    result_table3.to_excel(writer, sheet_name='MLR 2 Significant', index=False)


###############################################
### result table significant T2D population ###
###############################################

# Merge dataframes on 'rsID'
result_table = pd.merge(
    lookup,
    Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D,
    left_on='base_rsID',
    right_on='rsID_eGFR_2021',
    how='inner' 
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
    'Beta_UACR': 'Beta ln(UAE)',
    'p_val_UACR': 'P-value UAE',
    'eGFR CI Lower': 'CI Lower eGFR',
    'eGFR CI Upper': 'CI Upper eGFR',
    'UACR CI Lower': 'CI Lower UAE',
    'UACR CI Upper': 'CI Upper UAE'
    
})

result_table3 = result_table3.drop_duplicates(subset='rsID', keep='first')




#############################################
### result table MLR 1 overall population ###
#############################################
MLR_eGFR_2021_all_T2D['clean_rsID'] = MLR_eGFR_2021_all_T2D['rsID'].str.split('_').str[0]

# Merge dataframes on 'rsID' (or the equivalent shared identifier)
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
    'Beta_UACR': 'Beta ln(UAE)',
    'p_val_UACR': 'P-value UAE',
    'eGFR CI Lower': 'CI Lower eGFR',
    'eGFR CI Upper': 'CI Upper eGFR',
    'UACR CI Lower': 'CI Lower UAE',
    'UACR CI Upper': 'CI Upper UAE'
    
})

# Export the final table to an Excel file
output_file = 'file/path/T2D Results.xlsx'

with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
    result_table1.to_excel(writer, sheet_name='MLR 1', index=False)
    result_table2.to_excel(writer, sheet_name='MLR 2', index=False)
    result_table3.to_excel(writer, sheet_name='MLR 2 Significant', index=False)







