#!/usr/bin/env python3

# Usage: a.py [OPTIONS] [barcode_files] > barcode_combinations.out

import sys
import argparse
import itertools

def parse_args():
    parser = argparse.ArgumentParser(description='Enumerate combinatorial barcodes from input files.')
    parser.add_argument('barcode_files', nargs='+', help='Input barcode files')
    parser.add_argument('--input-translate', action='store_true', help='Flag to indicate if input files are for barcode translation, where the barcode is in the second column')
    parser.add_argument('--input-translate-delimiter', help='Delimiter for input translator files (default: tab,space)', default=None)
    parser.add_argument('--translate-output', help='Output file for barcode translation (default: no translation)', default=None)
    parser.add_argument('--translate-delimiter', default='-', help='Delimiter for translation output file')
    return parser.parse_args()

params = parse_args()
enumerated_barcodes = []
enumerated_barcodes_translate = {}
barcodes = []
barcode_translation = [] # A list of dictionaries

# When output translate file, the input has to be specificed as a translate file
if (params.translate_output is not None
    and params.input_translate is False):
    sys.exit("Error: When output translate file is specified, the input has to be specified as a translate file.")

# Iterate through the barcode files
for i, barcode_file in enumerate(params.barcode_files):
    with open(barcode_file, 'r') as f:
        if (params.input_translate):
            # If the input files are translators, read the second column as barcodes
            barcodes.append([line.rstrip().split( params.input_translate_delimiter )[1] for line in f])
        else:
            barcodes.append([line.rstrip() for line in f])

# Enumerate all combinations of barcodes
enumerated_barcodes = list(itertools.product(*barcodes))
for barcode in enumerated_barcodes:
    print(''.join(barcode))

if (params.translate_output is not None):
    # Read in the translation files and create a list of dictionaries for each barcode file
    for i, barcode_file in enumerate(params.barcode_files):
        with open(barcode_file, 'r') as f:
            if (params.translate_output is not None):
                if (params.input_translate_delimiter == None):
                    barcode_translation.append({line.rstrip().split()[1]: line.rstrip().split()[0] for line in f})
                else:
                    barcode_translation.append({line.rstrip().split( params.input_translate_delimiter )[1]: line.rstrip().split( params.input_translate_delimiter )[0] for line in f})
    
    # Create a translation dictionary for the enumerated barcodes
    for barcode in enumerated_barcodes:
        translated_barcode = []
        for i, b in enumerate(barcode):
            translated_barcode.append(barcode_translation[i][b])
        enumerated_barcodes_translate[''.join(barcode)] = params.translate_delimiter.join(translated_barcode)

    # Write the translation to the output file
    with open(params.translate_output, 'w') as f:
        for barcode, translation in enumerated_barcodes_translate.items():
            f.write(f"{translation}\t{barcode}\n")
