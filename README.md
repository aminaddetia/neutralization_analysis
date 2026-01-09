# Neutralization Analysis

This program takes Excel workbooks generated from the 96-well neutralization protocol on the BioTek Neo2 plate reader and normalizes, fits, and plots the data. By default, all technical replicates are kept as individual data points and the data are fit using the Hill equation for a normalized response using non-linear least squares curve fitting. The program generates a single PDF file with a plot for each of the samples along with a CSV file containing IC50 values, the Hill slope, R-squared, and RMSE values and a CSV file with the normalized entry values. 

## Initial Set-up

1. Install Python 3 using Homebrew.
```bash
brew install python
```
2. Install pandas, openpyxl, numpy, scipy, scikit-learn, matplotlib.
```bash
pip install numpy pandas openpyxl scipy scikit-learn matplotlib
```
3. Clone Git Repo locally.
```bash
cd ~
git clone https://github.com/aminaddetia/neutralization_analysis
```
4. Make neutralization_analysis.py executable from the command line. 
```bash
cd ~/neutralization_analysis/bin
chmod +x neutralization_analysis.py
echo 'export PATH="<absolute/path/to>/neutralization_analysis/bin:$PATH"' >> ~/.zshrc 
source ~/.zshrc
```
5. Confirm everything is working properly.
```bash
mkdir ~/Documents/neutralization_example
cd ~/Documents/neutralization_example
cp ~/neutralization_analysis/example/* .
neutralization_analysis.py -e example_Data.xlsx -m example_metadata.csv -l example_plate_layout.csv -o setup_test
```
Check if a new directory titled 'setup_test' was created inside of the neutralization_example directory and contains four ouput files: setup_test_ic50_results.csv, setup_test_tidy_data.csv, setup_test_individual_plots.pdf, setup_test_normalized_values.csv. If so, you ready to run the program with your data.

## Running neutralization_analysis.py

### Experimental set-up

Follow the standards below when setting up your neutralization experiment:

1. All plates must have the exact same layout - if you are planning on running some samples with more dilutions, different numbers of technical replicates, etc. either read the plates separately so the results end up in separate excel files and provide a different plate layout file for each excel file or analyze with Prism.
   
2. The 0% entry (cells only) wells and 100% entry (virus + cells) wells must be in the same location on every plate.

3. If columns/rows were skipped on a plate (beside the final plate), this can be indicated in the metadata file. Importantly, for example, if samples are run down columns in technical duplicate, two columns should be skipped, not shifted by one so the plate still retains the same layout.

4. If the final plate is not filled, start filling samples in the same order as the other plates (ie. if the first sample on the plate is in column 1, start filling the final plate at column 1 and leave the higher columns empty).

### File set-up

1. **Results excel file:** No changes should be necessary for the results file if it was directly output from the plate reader.

2. **Metada csv file:**
  - Each sample in the metadata file have a unique name when the values for the "Sample_ID", "Bleed_or_Dose", "Experimental_Manipulation", and "Virus" columns are combined.
  - Values are required for "Sample_ID", "Sample_Type", "Starting_Dilution_or_Concentration", "Dilution_Factor", and "Number_of_Dilutions" columns. The other columns are optional, but can be useful for sorting data or grouping plots.
  - Samples must be in the same order as they are on the plate(s). If columns/rows were skipped in any plate besides the final plate, an entry with a "Sample_ID" of "skip" must be included in the metadata file. For example, if samples 1-4 were supposed to be on plate 1, but due to a set-up error, sample 4 was moved to the next plate, you must create an entry named "skip" between samples 3 and 4. If the final plate has empty columns/rows, this can be ignored and no additional entries need to be created.
  - The options for "Sample_Type" are either "sera" or "antibody". Use "sera" if you are plotting reciprocal dilution on the x-axis and "antibody" if you are plotting concentration on the a-axis.
  - The "Dilution_Factor" should be a positive integer. For example, if you do a 1:3 serial dilution, put 3 in this column not 1/3.

3. **Plate Layout csv file:**
- To indicate which wells reflect 0% entry (cells only or mock infected), use "cells_only". To indicate which wells reflect 100% entry (virus + cells), use "virus_only". To indicate unused wells across all plates, use "blank".
- To indicate where the samples are on the plate, use sample1, sample2, sample3, etc. To indicate the dilutions, add _1, _2, _3, etc., to the end of each sample#, so what you are filling in on the plate layout is sample1_1, sample1_2, sample2_1, sample2_2, etc.  If technical replicates are conducted, there will mutliple wells with the same name (eg. if two technical replicates, two wells will have the name sample1_1).

Example results, metadata, and plate layout files are included in the "example" directory for reference. Empty metadata and plate layout files are included in the "input_files" directory.

### Running the program

The suggested directory structure is recommended:
  1. Create a folder for the entire project (eg. January 2026 Nanoparticle Study).
  2. Within the project folder, create a new folder for the specific experiment ran (eg. 260109 Wuhan Neutralization Rep1)
  3. Transfer the results excel file and completed metadata and plate layout files to this folder.
  4. After running neutralization_analysis.py, the results will show up in a new folder in the experiment folder.

To run the program, use the following command:
```bash
cd ~/<absolute-path-to-experiment-folder>
neutralization_analysis.py --excel <name-of-results-excel-file> --metadata <name-of-metadata-csv-file> --layout <name-of-plate-layout> --output <name-for-output-directory-and-files>
```
Alternatively, the program can be run using shorthand notion:
```bash
neutralization_analysis.py -e <name-of-results-excel-file> -m <name-of-metadata-csv-file> -l <name-of-plate-layout> -o <name-for-output-directory-and-files>
```

### Advanced options

To run neutralization_analysis.py without producing plots, use the --no-plots or -np flag:
```bash
neutralization_analysis.py -e <name-of-results-excel-file> -m <name-of-metadata-csv-file> -l <name-of-plate-layout> -o <name-for-output-directory-and-files> -np
```

To run the curve fitting using the mean value for technical replicates, rather than the individual data points, use the --fit-mean or -fm flag:
```bash
neutralization_analysis.py -e <name-of-results-excel-file> -m <name-of-metadata-csv-file> -l <name-of-plate-layout> -o <name-for-output-directory-and-files> -fm
```

## Understanding the outputs

All results are output into the newly created directory in the experiment folder. Both the directory and files are named for the value given for the --output/-o flag. Note: if the suggested directory structure wasn't followed, the new directory will be created in the same directory containing the results excel file.

1. **_normalized_values.csv:** this file contains the averaged normalized entry, standard deviation, and number of technical replicates for each dilution for each sample. This can be directly pasted into Prism (XY data with mean, std, n calculated elsewhere).

2. **_ic50_results.csv:** this file contains the IC50 and Hill slope values calculated for each sample. It additionally contains metrics about how well the curve fit: R-squared (ideally 1) and RMSE (ideally 0). There is an additional warning column that indicates if the calculated IC50 is outside the limit of detection (either IC50 > highest dilution or IC50 < lowest dilution).

3. **_individual_plots.csv:** this is a single PDF containing an individual plot for each of the samples analyzed.

4. **_tidy_data.csv:** this file is for subsequent plotting programs.

## Additional Information

This workflow was set up using the following software versions:
- Python 3.7.7
- Pandas 1.3.5
- Numpy 1.17.4
- Openpyxl 3.1.3
- Scipy 1.6.2
- Scikit-learn 1.0.2
- Matplotlib 3.2.2
