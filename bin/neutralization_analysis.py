#!/usr/bin/env python

import pandas as pd
import openpyxl
import numpy as np
import warnings
from scipy.optimize import curve_fit, OptimizeWarning
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import math
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import os

parser = argparse.ArgumentParser(description="Analyzes raw neutralization data starting with Excel file output from BioTek Neo2. Normalizes data, determines IC50s, and plots data and curve for each sample.")
parser.add_argument('--excel', '-e', required=True, help='Path to Excel file produced by plate reader')
parser.add_argument('--layout', '-l', required=True, help='Plate to layout csv file')
parser.add_argument('--metadata', '-m', required=True, help='Path to metadata csv file')
parser.add_argument('--output', '-o', required=True, help='Prefix for output files')
parser.add_argument('--no-plot', '-np', action='store_true', help='Do not generate individual dose-response plots')
parser.add_argument('--fit-mean', '-fm', action='store_true', help='Fit mean values instead of individual points')
args = parser.parse_args()

plate_reader_file_path = os.path.abspath(args.excel)
plate_layout_file_path = os.path.abspath(args.layout)
medata_file_path = os.path.abspath(args.metadata)
output_prefix = args.output
parent_dir = os.path.dirname(plate_reader_file_path)
output_directory = os.path.join(parent_dir, output_prefix)
os.makedirs(output_directory, exist_ok=True)

#Loading Excel, metadata, and plate layout files

input_excel_sheet = openpyxl.load_workbook(plate_reader_file_path)

sheet_i=1

for sheet in input_excel_sheet:
    new_sheet_title = 'Plate_' + str(sheet_i)
    sheet.title = new_sheet_title
    sheet_i = sheet_i+1

renamed_sheet_file = output_prefix + '_renamed.xlsx'
renamed_sheet_file_path = os.path.join(output_directory, renamed_sheet_file)

input_excel_sheet.save(renamed_sheet_file_path)

renamed_excel_sheet = pd.read_excel(renamed_sheet_file_path, sheet_name=None)
plate_layout = pd.read_csv(plate_layout_file_path , index_col=[0])
metadata = pd.read_csv(medata_file_path)

#Normalizing RLU values in each well

normalized_RLU_d = {}

for key, value in renamed_excel_sheet.items():
    results_row = value['Unnamed: 0'] == 'Results'
    rlu_start = results_row.shift(4)
    rlu_start_index = rlu_start[rlu_start == True].index[0]
    clean_excel_sheet = value.drop(range(0,rlu_start_index)).drop(['Unnamed: 0', 'Unnamed: 1', 'Unnamed: 14'], axis = 1)
    clean_excel_sheet.columns = [1,2,3,4,5,6,7,8,9,10,11,12]
    clean_excel_sheet.reset_index(drop=True, inplace=True)
    clean_excel_sheet = clean_excel_sheet.rename(index={0: 'A',1: 'B',2: 'C',3: 'D',4: 'E',5: 'F',6: 'G',7:'H'})
    renamed_excel_sheet[key] = clean_excel_sheet

    clean_excel_sheet_list = clean_excel_sheet.values.flatten().tolist()
    plate_layout_list = plate_layout.values.flatten().tolist()
    layout_value_tuple_list = list(zip(plate_layout_list, clean_excel_sheet_list))

    virus_only_total_RLU = 0
    virus_only_counter = 0

    for tuple in layout_value_tuple_list:
        if tuple[0] == 'virus_only':
            virus_only_total_RLU = virus_only_total_RLU + tuple[1]
            virus_only_counter = virus_only_counter + 1

    virus_only_avg_RLU = virus_only_total_RLU/virus_only_counter

    cells_only_total_RLU = 0
    cells_only_counter = 0

    for tuple in layout_value_tuple_list:
        if tuple[0] == 'cells_only':
            cells_only_total_RLU = cells_only_total_RLU + tuple[1]
            cells_only_counter = cells_only_counter + 1

    cells_only_avg_RLU = cells_only_total_RLU/cells_only_counter

    clean_excel_sheet.replace('OVRFLW', -999999, inplace=True)
    clean_excel_sheet.fillna(-999999, inplace=True)

    def normalize(x):
        if x == -999999:
            return -999999
        if x != -999999:
            return (x - cells_only_avg_RLU) / (virus_only_avg_RLU - cells_only_avg_RLU)*100

    normalized_RLU_df = clean_excel_sheet.applymap(normalize)
    normalized_RLU_d.update({key:normalized_RLU_df})

#Determining location (plate and wells) of each sample

samples_per_plate = []
unique_samples_per_plate = []

plate_layout_list = plate_layout.values.flatten().tolist()

for item in plate_layout_list:
    if item == 'blank' or item == 'cells_only' or item == 'virus_only':
        continue
    else:
        item_list = item.split('_')
        samples_per_plate.append(item_list[0])

for item in samples_per_plate:
    if item not in unique_samples_per_plate:
        unique_samples_per_plate.append(item)

number_unique_samples_per_plate = len(unique_samples_per_plate)

#Creating a tidy dataframe, averaging normalized RLU values, and assigning averaged values to each sample & dilution. Writes out the tidy dataframe as a csv file.

sample_i = 1
plate_number = 1

tidy_df = pd.DataFrame(columns=['Sample_ID','Bleed', 'Treatment', 'Group', 'Virus', 'Sample_Type', 'Dilution','Plate_Number', 'Sample_Code'])

for row in metadata.itertuples():

    sample_ID = getattr(row, 'Sample_ID')
    bleed = getattr(row, 'Bleed_or_Dose')
    treatment = getattr(row, 'Experimental_Manipulation')
    group = getattr(row, 'Group')
    virus = getattr(row, 'Virus')
    dilutions = getattr(row, 'Number_of_Dilutions')
    starting_dilution = getattr(row, 'Starting_Dilution_or_Concentration')
    dilution_factor = getattr(row, 'Dilution_Factor')
    sample_type = getattr(row, 'Sample_Type')

    sample_number = 'sample' + str(sample_i)

    for i in range(dilutions):
        if sample_type.lower() == 'sera':
            total_dilution = starting_dilution / ((1/dilution_factor) ** i)
        if sample_type.lower() == 'antibody':
            total_dilution = starting_dilution / (dilution_factor ** i)
        sample_code = sample_number + '_' + str(i+1)
        new_tidy_data_d = {'Sample_ID': sample_ID, 'Bleed' : bleed, 'Treatment' : treatment, 'Virus' : virus, 'Group' : group, 'Sample_Type' : sample_type, 'Dilution': total_dilution, 'Plate_Number': plate_number, 'Sample_Code': sample_code}
        new_tidy_data_df = pd.DataFrame([new_tidy_data_d])
        tidy_df = pd.concat([tidy_df, new_tidy_data_df], ignore_index = True)

    sample_i = sample_i + 1
    if sample_i > number_unique_samples_per_plate:
        plate_number = plate_number + 1
        sample_i = 1

normalized_avg_list = []
normalized_std_list = []
normalized_n_list = []

for row in tidy_df.itertuples():

    key_number = getattr(row, 'Plate_Number')
    key = 'Plate_' + str(key_number)
    sample_code = getattr(row, 'Sample_Code')

    normalized_RLU_list = normalized_RLU_d[key].values.flatten().tolist()
    plate_layout_list = plate_layout.values.flatten().tolist()
    layout_normalized_tuple_list = list(zip(plate_layout_list, normalized_RLU_list))

    normalized_values = []

    for tuple in layout_normalized_tuple_list:
        if tuple[0] == sample_code:
            if tuple[1] != -999999:
                normalized_values.append(tuple[1])

    normalized_avg = np.mean(normalized_values)
    normalized_std = np.std(normalized_values)
    normalized_n = len(normalized_values)

    normalized_avg_list.append(normalized_avg)
    normalized_std_list.append(normalized_std)
    normalized_n_list.append(normalized_n)

tidy_df['Averaged_Normalized_Entry'] = normalized_avg_list
tidy_df['Standard_Deviation_Normalized_Entry'] = normalized_std_list
tidy_df['Technical_replicates'] = normalized_n_list

tidy_file_name = output_prefix + '_tidy_data.csv'
tidy_file_path = os.path.join(output_directory, tidy_file_name)

tidy_df.to_csv(tidy_file_path, header = True, index = False)

#Creates a new tidy dataframe where the normalized RLU values aren't averaged across technical replicates and kept separate instead.

tidy_unaveraged_df = pd.DataFrame(columns=['Sample_ID','Bleed', 'Treatment', 'Group', 'Virus', 'Sample_Type', 'Dilution','Plate_Number', 'Sample_Code', 'Technical_Replicate', 'Normalized_Entry'])

for row in tidy_df.itertuples():

    key_number = getattr(row, 'Plate_Number')
    key = 'Plate_' + str(key_number)
    sample_code = getattr(row, 'Sample_Code')

    normalized_RLU_list = normalized_RLU_d[key].values.flatten().tolist()
    plate_layout_list = plate_layout.values.flatten().tolist()
    layout_normalized_tuple_list = list(zip(plate_layout_list, normalized_RLU_list))

    normalized_values = []

    for tuple in layout_normalized_tuple_list:
        if tuple[0] == sample_code:
            if tuple[1] != -999999:
                normalized_values.append(tuple[1])

    tech_rep = 0

    for i in range(len(normalized_values)):
        pos = tech_rep
        tech_rep = tech_rep + 1
        new_row = [getattr(row,'Sample_ID'),getattr(row,'Bleed'), getattr(row,'Treatment'), getattr(row,'Group'), getattr(row,'Virus'), getattr(row,'Sample_Type'), getattr(row,'Dilution'),getattr(row,'Plate_Number'), getattr(row,'Sample_Code'), tech_rep, normalized_values[pos]]
        tidy_unaveraged_df.loc[len(tidy_unaveraged_df)] = new_row


#Creates and writes a csv file with the averaged normalzied RLU values, standard deviation, and number of technical replicates for each sample. Can be directly pasted into Prism for plotting or IC50 determination.

i = 0

for row in metadata.itertuples():

    sample_id = getattr(row, 'Sample_ID')
    bleed = str(getattr(row, 'Bleed_or_Dose'))
    treatment = str(getattr(row, 'Experimental_Manipulation'))
    virus = str(getattr(row, 'Virus'))
    sample_type = str(getattr(row, 'Sample_Type'))


    tidy_df['Bleed'].fillna('nan', inplace=True)
    tidy_df['Treatment'].fillna('nan', inplace=True)
    tidy_df['Virus'].fillna('nan', inplace=True)

    temp_df = tidy_df[(tidy_df['Sample_ID'] == sample_id) & (tidy_df['Bleed'] == bleed) & (tidy_df['Treatment'] == treatment) & (tidy_df['Virus'] == virus)]

    if bleed == 'nan':
        bleed = ''

    if treatment == 'nan':
        treatment = ''

    if virus == 'nan':
        virus = ''

    name = [sample_id, bleed, treatment, virus]
    filtered_name = list(filter(None, name))
    sample_name = ' '.join(filtered_name)

    if i == 0:
        normalized_results_df = temp_df.filter(['Dilution', 'Averaged_Normalized_Entry', 'Standard_Deviation_Normalized_Entry', 'Technical_replicates'])
        normalized_results_df.rename(columns={'Averaged_Normalized_Entry': sample_name + ' Mean', 'Standard_Deviation_Normalized_Entry': sample_name + ' STD', 'Technical_replicates' : sample_name + ' N'}, inplace=True)

    if i != 0:
        temp_df2 = temp_df.filter(['Dilution', 'Averaged_Normalized_Entry', 'Standard_Deviation_Normalized_Entry', 'Technical_replicates'])
        temp_df2.rename(columns={'Averaged_Normalized_Entry': sample_name + ' Mean', 'Standard_Deviation_Normalized_Entry': sample_name + ' STD', 'Technical_replicates' : sample_name + ' N'}, inplace=True)
        normalized_results_df = pd.merge(normalized_results_df, temp_df2, on='Dilution', how='outer')

    i = i + 1

if sample_type.lower() == 'sera':
    normalized_results_df.sort_values(by='Dilution', ascending=True, inplace=True)

if sample_type.lower() == 'antibody':
    normalized_results_df.sort_values(by='Dilution', ascending=False, inplace=True)
    normalized_results_df.rename(columns={'Dilution' : 'Concentration'}, inplace = True)

normalized_file_name = output_prefix + '_normalized_values.csv'
normalized_file_path = os.path.join(output_directory, normalized_file_name)

normalized_results_df.to_csv(normalized_file_path, header = True, index = False)

#Defining Hill Equation

def hill_fit(x, IC50, s):
    return 100/(1+((IC50/x)**s))

#Calculates IC50 and Hill slope for each sample. If --fit-mean, calculates IC50 and Hill slope from averaged normalized RLU, otherwise calculates from individual data points. Writes out IC50 results in a csv file.

ic50_list = []
hill_slope_list = []
r_squared_list = []
rmse_list = []
warning_list = []

for row in metadata.itertuples():
    with np.errstate(invalid='ignore'), warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=OptimizeWarning)

        sample_id = getattr(row, 'Sample_ID')
        bleed = str(getattr(row, 'Bleed_or_Dose'))
        treatment = str(getattr(row, 'Experimental_Manipulation'))
        virus = str(getattr(row, 'Virus'))
        sample_type = str(getattr(row, 'Sample_Type'))

        if args.fit_mean:
            print('fit mean')
            tidy_df['Bleed'].fillna('nan', inplace=True)
            tidy_df['Treatment'].fillna('nan', inplace=True)
            tidy_df['Virus'].fillna('nan', inplace=True)

            temp_df = tidy_df[(tidy_df['Sample_ID'] == sample_id) & (tidy_df['Bleed'] == bleed) & (tidy_df['Treatment'] == treatment) & (tidy_df['Virus'] == virus)]
            dilution_column = temp_df['Dilution']
            x = dilution_column.to_numpy()
            p0_dilution = np.median(x)
            xmin = np.min(x)
            xmax = np.max(x)
            entry_column = temp_df['Averaged_Normalized_Entry']
            y = entry_column.to_numpy()
        else:
            tidy_unaveraged_df['Bleed'].fillna('nan', inplace=True)
            tidy_unaveraged_df['Treatment'].fillna('nan', inplace=True)
            tidy_unaveraged_df['Virus'].fillna('nan', inplace=True)

            temp_df = tidy_unaveraged_df[(tidy_unaveraged_df['Sample_ID'] == sample_id) & (tidy_unaveraged_df['Bleed'] == bleed) & (tidy_unaveraged_df['Treatment'] == treatment) & (tidy_unaveraged_df['Virus'] == virus)]
            dilution_column = temp_df['Dilution']
            x = dilution_column.to_numpy()
            p0_dilution = np.median(x)
            if sample_type.lower() == 'sera':
                s0 = 1
            if sample_type.lower() == 'antibody':
                s0 = -1
            xmin = np.min(x)
            xmax = np.max(x)
            entry_column = temp_df['Normalized_Entry']
            y = entry_column.to_numpy()

        param, param_cov = curve_fit(hill_fit, x, y, p0 = [p0_dilution, 1],  maxfev=10000)

        ic50 = param[0]
        if ic50 < xmin:
            warning = 'IC50 below LOD'
        elif ic50 > xmax:
            warning = 'IC50 above LOD'
        else:
            warning = ''

        ic50_list.append(ic50)
        hill_slope_list.append(param[1])

        y_pred = hill_fit(x, *param)
        r_squared = r2_score(y, y_pred)
        mse = np.mean((y-y_pred)**2)
        rmse = np.sqrt(mse)
        r_squared_list.append(r_squared)
        rmse_list.append(rmse)
        warning_list.append(warning)

ic50_df = metadata.filter(['Sample_ID', 'Bleed_or_Dose', 'Experimental_Manipulation', 'Virus', 'Group', 'Sample_Type'])
ic50_df['IC50'] = ic50_list
ic50_df['Hill_Slope'] = hill_slope_list
ic50_df['r_squared'] = r_squared_list
ic50_df['rmse'] = rmse_list
ic50_df['Warning'] = warning_list

ic50_results_file = output_prefix + '_ic50_results.csv'
ic50_results_path = os.path.join(output_directory, ic50_results_file)

ic50_df.to_csv(ic50_results_path, header = True, index = False)

#Plots data creating an individual plot for each sample. Writes out a PDF with all plots. If --no-plot, does not plot data.

if args.no_plot:
    pass

else:
    sample_name_list = []
    for row in ic50_df.itertuples():

        sample_id = getattr(row, 'Sample_ID')
        bleed = str(getattr(row, 'Bleed_or_Dose'))
        treatment = str(getattr(row, 'Experimental_Manipulation'))
        virus = str(getattr(row, 'Virus'))

        if bleed == 'nan':
            bleed = ''

        if treatment == 'nan':
            treatment = ''

        if virus == 'nan':
            virus = ''

        name = [sample_id, bleed, treatment, virus]
        filtered_name = list(filter(None, name))
        sample_name = ' '.join(filtered_name)

        sample_name_list.append(sample_name)


    internal_ic50_df = ic50_df.filter(['IC50', 'Hill_Slope'])
    internal_ic50_df['Sample_Name'] = sample_name_list


    hill_slope_df = internal_ic50_df.filter(['Sample_Name','Hill_Slope'])
    hill_slope_d = dict(hill_slope_df.values)
    ic50_only_df = internal_ic50_df.filter(['Sample_Name','IC50'])
    ic50_only_d = dict(ic50_only_df.values)

    fig_size = len(internal_ic50_df['Sample_Name'].tolist())
    fig_len = math.ceil(fig_size/4)

    fig = plt.figure(figsize = (12,fig_len*3))

    fig_i = 1

    for row in metadata.itertuples():

        sample_id = getattr(row, 'Sample_ID')
        bleed = str(getattr(row, 'Bleed_or_Dose'))
        treatment = str(getattr(row, 'Experimental_Manipulation'))
        virus = str(getattr(row, 'Virus'))
        sample_type = str(getattr(row, 'Sample_Type'))
        unit = str(getattr(row, 'Unit'))


        tidy_df['Bleed'].fillna('nan', inplace=True)
        tidy_df['Treatment'].fillna('nan', inplace=True)
        tidy_df['Virus'].fillna('nan', inplace=True)

        temp_plotting_df = tidy_df[(tidy_df['Sample_ID'] == sample_id) & (tidy_df['Bleed'] == bleed) & (tidy_df['Treatment'] == treatment) & (tidy_df['Virus'] == virus)]

        if bleed == 'nan':
            bleed = ''

        if treatment == 'nan':
            treatment = ''

        if virus == 'nan':
            virus = ''

        name = [sample_id, bleed, treatment, virus]
        filtered_name = list(filter(None, name))
        sample_name = ' '.join(filtered_name)

        dilution_column = temp_plotting_df['Dilution']
        x = dilution_column.to_numpy()
        entry_column = temp_plotting_df['Averaged_Normalized_Entry']
        y = entry_column.to_numpy()
        std_column = temp_plotting_df['Standard_Deviation_Normalized_Entry']
        e = std_column.to_numpy()

        hill_slope = hill_slope_d[sample_name]
        ic50 = ic50_only_d[sample_name]

        new_x = np.logspace(np.log10(min(x)), np.log10(max(x)), 1000)
        y_smooth = hill_fit(new_x, ic50, hill_slope)

        ax = fig.add_subplot(fig_len,4,fig_i)
        ax.plot(new_x, y_smooth, 'k-', color = 'black')
        ax.errorbar(x, y, e, linestyle = 'None', marker = 'o', mec = 'black', mfc = 'white', ecolor = 'black', ms = 8, capsize = 5)
        ax.set_xscale('log')
        y_range = ax.get_ylim()
        if y_range[0] > 0:
            ax.set_ylim(bottom = 0)
        if y_range[1] < 100:
            ax.set_ylim(top = 100)
        if sample_type.lower() == 'sera':
            ax.set_xlabel('Reciprocal Dilution')
        if sample_type.lower() == 'antibody':
            if unit == None:
                ax.set_xlabel('Concentration')
            else:
                ax.set_xlabel('Concentration' + ' [' + unit + ']')
        ax.set_ylabel('Percent Entry (%)')
        ax.set_title(sample_name)
        fig_i = fig_i + 1

    plots_file_name = output_prefix + '_individual_plots.pdf'
    plots_file_path = os.path.join(output_directory, plots_file_name)

    with PdfPages(plots_file_path) as pdf:
        fig.tight_layout()
        pdf.savefig()

os.remove(renamed_sheet_file_path)
