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
echo 'export PATH="<absolute/path/to>/neutralization_analysis/bin:$PATH"' >> ~/.zshrc #Replace <absolute/path/to> based on your local naming (will likely be /Users/Username/)
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


