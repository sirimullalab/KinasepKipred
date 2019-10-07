These codes are taken from https://github.com/PatWalters/metk.git
I modified files and used for different evaluation metrics, specially for confidence interval calculation.

Usage:
python metk.py --in input_file_name.csv --prefix results_filename

(--in and --prefix should be provided as shown above)
(two output files will be generated with .txt and .pdf extension)
(so we don't have to provide the extension for the output file, just prefix (filename) should be provided instead)

Note: 
If an error is occurred with matplotlib (in python2.7)

Create a file ~/.matplotlib/matplotlibrc and add backend:TkAgg
