# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 10:32:56 2024
Health risk assessment. Using the RfC database from USEPA. Molecular weights database was used too. 
@author: lb945465
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

# Define file paths and directories
concentration_dir = r'C:\Users\LB945465\OneDrive - University at Albany - SUNY\State University of New York\Spyder\1_Preparation of the Database\PMF_Dispersion'
constrained_dir = r'C:\Users\LB945465\Desktop\DN PMF Optimal'
source_profiles_file = os.path.join(constrained_dir, 'Source_Profiles_named.xlsx')
rfc_file = r'C:\Users\LB945465\OneDrive - University at Albany - SUNY\State University of New York\Spyder\7_Health risk assessment\RfC_Toxicity_Values.xlsx'
iur_file = r'C:\Users\LB945465\OneDrive - University at Albany - SUNY\State University of New York\Spyder\7_Health risk assessment\RfC_Toxicity_Values.xlsx'
molecular_weight_file = r'C:\Users\LB945465\OneDrive - University at Albany - SUNY\State University of New York\Spyder\7_Health risk assessment\VOC unit conversion_detected species.xlsx'

# Load RfC values from the Excel file
rfc_data = pd.read_excel(rfc_file, sheet_name='Constants', index_col='Specie')

# Remove special characters and reassign to the index
rfc_data.index = rfc_data.index.str.replace(',', '', regex=True)

rfc_values = rfc_data['RfC'].to_dict()  # Convert RfC values to dictionary

# Load RfC values from the Excel file
iur_data = pd.read_excel(iur_file, sheet_name='Constants', index_col='Specie')

# Remove special characters and reassign to the index
iur_data.index = iur_data.index.str.replace(',', '', regex=True)

iur_values = iur_data['IUR'].to_dict()  # Convert RfC values to dictionary

# Load molecular weights from the Excel file
mw_data = pd.read_excel(molecular_weight_file, index_col='Specie')

# Remove special characters and reassign to the index
mw_data.index = mw_data.index.str.replace(',', '', regex=True)

# Convert molecular weights to dictionary
molecular_weights = mw_data['MW'].to_dict()

# List of sites
sites = ['Bronx', 'Queens', 'Kings', 'Richmond', 'Elizabeth', 'Chester']

# Load concentration data
concentration_data = {}
for site in sites:
    conc_file = os.path.join(concentration_dir, f'{site}_conc.csv')
    df = pd.read_csv(conc_file, parse_dates=['Date'], index_col='Date')
    df = df.drop(columns=['TVOC', 'TVOC_recons'], errors='ignore')
    
    # remove special character
    df.columns = df.columns.str.replace(',', '')
    
    concentration_data[site] = df

# Load source profiles
source_profiles_data = pd.read_excel(source_profiles_file, sheet_name=None, index_col='Species')
for sheet_name, df in source_profiles_data.items():
    # remove special character
    df.index = df.index.str.replace(',', '')
    
    source_profiles_data[sheet_name] = df.drop(index=['TVOC', 'TVOC_recons'], errors='ignore')
    
# Load VC_ratio data
vc_ratio_data = {}
for site in sites:
    vc_file = os.path.join(concentration_dir, f'VC_{site}.xlsx')
    vc_df = pd.read_excel(vc_file, parse_dates=['Date'], index_col='Date')
    vc_ratio_data[site] = vc_df['VC_ratio']

# Create VOC_modelled dataframe for each site
voc_modelled_data = {}
for site in sites:
    VOC_conc = concentration_data[site]
    
    # Create VOC_modelled by subtracting Residuals from VOC_conc
    VOC_modelled = VOC_conc
    
    # Merge VOC_modelled with VC_ratio
    vc_ratio = vc_ratio_data[site]
    VOC_unnorm = VOC_modelled.merge(vc_ratio, left_index=True, right_index=True, how='left')
    
    # Drop rows where VC_ratio is zero or NaN
    VOC_unnorm = VOC_unnorm[VOC_unnorm['VC_ratio'] > 0]
    
    # Normalize VOC_modelled by VC_ratio
    for col in VOC_modelled.columns:
        VOC_unnorm[col] = VOC_unnorm[col] / VOC_unnorm['VC_ratio']
    
    # Drop the VC_ratio column
    VOC_unnorm.drop(columns=['VC_ratio'], inplace=True)
    
    voc_modelled_data[site] = VOC_unnorm

# Function to convert ppbv to µg/m3 using molecular weight
def ppbv_to_ugm3(concentration_ppbv, molecular_weight):
    # Convert concentration from ppbv to µg/m³
    return concentration_ppbv * (molecular_weight / 24.45)

# Convert VOC_modelled data to µg/m3 for selected species using the molecular weights from the loaded database
for site, df in voc_modelled_data.items():
    # Create a new dataframe to store the converted values
    df_ugm3 = pd.DataFrame(index=df.index)
    
    for voc in df.columns:
        if voc in molecular_weights:
            # Convert each VOC to µg/m³ using the molecular weight from the database
            df_ugm3[voc] = ppbv_to_ugm3(df[voc], molecular_weights[voc])
        else:
            print(f"Skipping {voc} - not in molecular weight database.")
    
    # Replace the original ppbv dataframe with the converted µg/m³ dataframe
    voc_modelled_data[site] = df_ugm3

# Define constants for EC calculation
ET = 24     # Exposure time in hours/day
EF = 350   # Exposure frequency in days/year
ED = 24    # Exposure duration in years
AT_general = 25 * 365 * 24  # Averaging time in hours for non-cancer risk (25 years)
AT_lifetime = 70 * 365 * 24 # Averaging time in hours for lifetime cancer risk (70 years)

# Function to calculate Exposure Concentration (EC)
def calculate_ec(concentration, et, ef, ed, at):
    # EC = (C * ET * EF * ED) / AT
    return (concentration * et * ef * ed) / at

# Function to calculate HQ based on selected scenario with debug statements
def calculate_hq(voc_modelled_data, scenario='general'):
    # Set AT based on scenario
    if scenario == 'general':
        at = AT_general
    elif scenario == 'lifetime':
        at = AT_lifetime
    else:
        raise ValueError("Invalid scenario. Choose 'general' or 'lifetime'.")

    # Create a dictionary to store the EC and HQ values separately for each site
    ec_values = {}  # Store EC values separately
    hq_values = {}  # Store HQ values separately

    # Calculate EC and HQ for each site and VOC species
    for site, df in voc_modelled_data.items():
        # Create new dataframes to store EC and HQ values
        ec_df = pd.DataFrame(index=df.index)
        hq_df = pd.DataFrame(index=df.index)
        
        print(f"\nCalculating EC and HQ for site: {site}")
        
        for voc in df.columns:
            if voc in rfc_values:
                # Calculate EC based on the selected scenario
                ec = calculate_ec(df[voc], ET, EF, ED, at)
                
                # Calculate HQ
                hq = ec / rfc_values[voc]
                
                # Store the calculated EC and HQ values in separate dataframes
                ec_df[f'{voc}'] = ec
                hq_df[f'{voc}'] = hq
                
                print(f"VOC: {voc}, EC: {ec.mean()}, HQ: {hq.mean()}")  # Debugging: Print mean values for verification
            else:
                print(f"Skipping {voc} - RfC value not available.")
        
        # Store the calculated EC and HQ values for the site
        ec_values[site] = ec_df
        hq_values[site] = hq_df

    return ec_values, hq_values

# Calculate EC and HQ values based on the selected scenario
selected_scenario = 'lifetime'  # Change to 'general' or 'lifetime' for lifetime risk assessment

# Calculate EC and HQ values separately
ec_results, hq_results = calculate_hq(voc_modelled_data, scenario=selected_scenario)

# Calculate HQ per factor for each site
hq_factors = {}

for site in sites:
    print(f"Processing site: {site}")
    
    site_hq_data = hq_results[site]  # Use hq_results instead of hq_data as per our updated structure

    # Ensure alignment of columns between HQ data and source profiles
    site_source_profiles = source_profiles_data[site] / 100  # Convert percentages to fractions
    
    # Filter only common species between site_hq_data and site_source_profiles
    common_species = site_hq_data.columns.intersection(site_source_profiles.index)
    site_hq_data_filtered = site_hq_data[common_species]
    site_source_profiles_filtered = site_source_profiles.loc[common_species]

    # Create an empty dataframe to store the HQ values per factor
    site_hq_factors = pd.DataFrame(index=site_hq_data.index)
    
    # Iterate over each factor (column) in source profiles
    for factor in site_source_profiles_filtered.columns:
        # Multiply each VOC HQ value by its respective source profile contribution for this factor
        hq_factor = site_hq_data_filtered.mul(site_source_profiles_filtered[factor], axis=1)
        
        # Sum across all VOC species to get total HQ per factor
        site_hq_factors[factor] = hq_factor.sum(axis=1)
    
    # Rename the index to 'Date'
    site_hq_factors.index.name = 'Date'
    
    # Store the HQ factors for the current site
    hq_factors[site] = site_hq_factors

# Output for verification
for site, df in hq_factors.items():
    print(f"\nHQ Factors for {site} (first 5 rows):")
    print(df.head())

# Create a dictionary to store LCR values for each site
lcr_results = {}

# Calculate LCR for each site and VOC species
for site, ec_df in ec_results.items():
    # Create a new dataframe to store LCR values
    lcr_df = pd.DataFrame(index=ec_df.index)
    
    for voc in ec_df.columns:
        if voc in iur_values and iur_values[voc] is not None:
            # Calculate LCR by multiplying EC with IUR
            lcr = ec_df[voc] * iur_values[voc]
            lcr_df[voc] = lcr
    
    # Store the calculated LCR values for the site
    lcr_results[site] = lcr_df

# Output for verification
lcr_results_summary = {site: df.mean().to_dict() for site, df in lcr_results.items()}

# Calculate LCR per factor for each site
lcr_factors = {}

for site in sites:
    print(f"Processing site: {site}")
    
    site_lcr_data = lcr_results[site]  # Use lcr_results as the base for LCR data
    
    # Ensure alignment of columns between LCR data and source profiles
    site_source_profiles = source_profiles_data[site] / 100  # Convert percentages to fractions
    
    # Filter only common species between site_lcr_data and site_source_profiles
    common_species = site_lcr_data.columns.intersection(site_source_profiles.index)
    site_lcr_data_filtered = site_lcr_data[common_species]
    site_source_profiles_filtered = site_source_profiles.loc[common_species]

    # Create an empty dataframe to store the LCR values per factor
    site_lcr_factors = pd.DataFrame(index=site_lcr_data.index)
    
    # Iterate over each factor (column) in source profiles
    for factor in site_source_profiles_filtered.columns:
        # Multiply each VOC LCR value by its respective source profile contribution for this factor
        lcr_factor = site_lcr_data_filtered.mul(site_source_profiles_filtered[factor], axis=1)
        
        # Sum across all VOC species to get total LCR per factor
        site_lcr_factors[factor] = lcr_factor.sum(axis=1)
    
    # Rename the index to 'Date'
    site_lcr_factors.index.name = 'Date'
    
    # Store the LCR factors for the current site
    lcr_factors[site] = site_lcr_factors

# Output for verification
for site, df in lcr_factors.items():
    print(f"\nLCR Factors for {site} (first 5 rows):")
    print(df.head())

# Create lists to keep track of missing RfC and IUR values
missing_rfc = []
missing_iur = []

# Calculate LCR for each site and VOC species
for site, ec_df in ec_results.items():
    # Create a new dataframe to store LCR values
    lcr_df = pd.DataFrame(index=ec_df.index)
    
    for voc in ec_df.columns:
        if voc in iur_values and iur_values[voc] is not None:
            # Calculate LCR by multiplying EC with IUR
            lcr = ec_df[voc] * iur_values[voc]
            lcr_df[voc] = lcr
        elif voc not in iur_values or iur_values[voc] is None:
            print(f"Skipping {voc} - IUR value not available.")
            missing_iur.append(voc)  # Track missing IUR values
    
    # Store the calculated LCR values for the site
    lcr_results[site] = lcr_df

# Output the list of missing RfC and IUR values
print("\nMissing RfC values:")
print(set(missing_rfc))

print("\nMissing IUR values:")
print(set(missing_iur))
    
# Function to visualize results (HQ or LCR) with 95% confidence intervals
def visualize_risk_assessment_with_ci(factors_data, assessment_type='HQ'):
    # Combine factors from all sites into a single DataFrame
    combined_factors = pd.DataFrame()

    for site, df in factors_data.items():
        combined_factors = pd.concat([combined_factors, df], axis=0)

    # Calculate the average value per factor across all data
    average_value_by_factor = combined_factors.mean()

    # Calculate the standard error of the mean (SEM)
    sem_by_factor = combined_factors.sem()

    # Calculate the 95% confidence interval (mean ± 1.96 * SEM)
    ci_lower = average_value_by_factor - 1.96 * sem_by_factor
    ci_upper = average_value_by_factor + 1.96 * sem_by_factor

    # Sort the average values and confidence intervals in ascending order
    sorted_factors = average_value_by_factor.sort_values(ascending=True).index
    average_value_by_factor_sorted = average_value_by_factor[sorted_factors]
    ci_lower_sorted = ci_lower[sorted_factors]
    ci_upper_sorted = ci_upper[sorted_factors]

    # Create a bar plot of the sorted average value by factor with error bars
    plt.figure(figsize=(15, 5), dpi=300)
    plt.bar(average_value_by_factor_sorted.index, 
            average_value_by_factor_sorted, 
            yerr=[average_value_by_factor_sorted - ci_lower_sorted, 
                  ci_upper_sorted - average_value_by_factor_sorted],
            capsize=5, 
            color='lightcoral',
            alpha=0.7)

    plt.xlabel('Factors')
    plt.ylabel(f'Average {assessment_type}')
    plt.title(f'Average {assessment_type} by Factor Across All Sites (with 95% CI)')
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()
    
    # Calculate total value by summing all factors in the provided data for each site
    total_values = {}
    
    for site, df in factors_data.items():
        total_values[site] = df.sum(axis=1)
    
    # Create a time series plot for total values (HQ/LCR) in each site
    plt.figure(figsize=(15, 5), dpi=300)
    
    for site, value_series in total_values.items():
        plt.plot(value_series.index, value_series, label=site, alpha=0.7)
    
    plt.xlabel('Date')
    plt.ylabel(f'Total {assessment_type}')
    plt.title(f'Total {assessment_type} Over Time for Each Site')
    plt.legend(title='Sites')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Visualize HQ with confidence intervals
visualize_risk_assessment_with_ci(hq_factors, assessment_type='HQ')

# Visualize LCR with confidence intervals
visualize_risk_assessment_with_ci(lcr_factors, assessment_type='LCR')
