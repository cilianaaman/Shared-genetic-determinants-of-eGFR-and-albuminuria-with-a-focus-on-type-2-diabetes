#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: cilianaaman
"""

import openpyxl
import numpy as np
import pandas as pd
#%% overall population 

# Load the Excel file
file_path = 'file/path/Overall Population Results.xlsx'
output_path = 'file/path/Overall Population Results With Percent Change.xlsx'

# -------- MLR 1 --------
df_mlr1 = pd.read_excel(file_path, sheet_name='MLR 1')
idx_egfr = df_mlr1.columns.get_loc('Beta ln(eGFR)') + 1
df_mlr1.insert(idx_egfr, '% change mean eGFR', (np.exp(df_mlr1['Beta ln(eGFR)']) - 1) * 100)

# -------- MLR 2 --------
df_mlr2 = pd.read_excel(file_path, sheet_name='MLR 2')
idx_egfr = df_mlr2.columns.get_loc('Beta ln(eGFR)') + 1
idx_uae = df_mlr2.columns.get_loc('Beta ln(UAE)') + 1
df_mlr2.insert(idx_egfr, '% change mean eGFR', (np.exp(df_mlr2['Beta ln(eGFR)']) - 1) * 100)
df_mlr2.insert(idx_uae + 1, '% change mean UAE', (np.exp(df_mlr2['Beta ln(UAE)']) - 1) * 100)

# -------- MLR 2 Significant --------
df_mlr2_sig = pd.read_excel(file_path, sheet_name='MLR 2 Significant')
idx_egfr = df_mlr2_sig.columns.get_loc('Beta ln(eGFR)') + 1
idx_uae = df_mlr2_sig.columns.get_loc('Beta ln(UAE)') + 1
df_mlr2_sig.insert(idx_egfr, '% change mean eGFR', (np.exp(df_mlr2_sig['Beta ln(eGFR)']) - 1) * 100)
df_mlr2_sig.insert(idx_uae + 1, '% change mean UAE', (np.exp(df_mlr2_sig['Beta ln(UAE)'] )-1) * 100)

# -------- Save to new Excel file --------
with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
    df_mlr1.to_excel(writer, sheet_name='MLR 1', index=False)
    df_mlr2.to_excel(writer, sheet_name='MLR 2', index=False)
    df_mlr2_sig.to_excel(writer, sheet_name='MLR 2 Significant', index=False)

print("done :)")

#%% T2D population

# Load the Excel file
file_path = 'file/path/T2D Results.xlsx'
output_path = 'file/path/T2D Results With Percent Change.xlsx'

# -------- MLR 1 --------
df_mlr1 = pd.read_excel(file_path, sheet_name='MLR 1')
idx_egfr = df_mlr1.columns.get_loc('Beta ln(eGFR)') + 1
df_mlr1.insert(idx_egfr, '% change mean eGFR', (np.exp(df_mlr1['Beta ln(eGFR)']) - 1) * 100)

# -------- MLR 2 --------
df_mlr2 = pd.read_excel(file_path, sheet_name='MLR 2')
idx_egfr = df_mlr2.columns.get_loc('Beta ln(eGFR)') + 1
idx_uae = df_mlr2.columns.get_loc('Beta ln(UAE)') + 1
df_mlr2.insert(idx_egfr, '% change mean eGFR', (np.exp(df_mlr2['Beta ln(eGFR)']) - 1) * 100)
df_mlr2.insert(idx_uae + 1, '% change mean UAE', (np.exp(df_mlr2['Beta ln(UAE)']) - 1) * 100)

# -------- MLR 2 Significant --------
df_mlr2_sig = pd.read_excel(file_path, sheet_name='MLR 2 Significant')
idx_egfr = df_mlr2_sig.columns.get_loc('Beta ln(eGFR)') + 1
idx_uae = df_mlr2_sig.columns.get_loc('Beta ln(UAE)') + 1
df_mlr2_sig.insert(idx_egfr, '% change mean eGFR', (np.exp(df_mlr2_sig['Beta ln(eGFR)']) - 1) * 100)
df_mlr2_sig.insert(idx_uae + 1, '% change mean UAE', (np.exp(df_mlr2_sig['Beta ln(UAE)'] )-1) * 100)

# -------- Save to new Excel file --------
with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
    df_mlr1.to_excel(writer, sheet_name='MLR 1', index=False)
    df_mlr2.to_excel(writer, sheet_name='MLR 2', index=False)
    df_mlr2_sig.to_excel(writer, sheet_name='MLR 2 Significant', index=False)

print("done :)")
