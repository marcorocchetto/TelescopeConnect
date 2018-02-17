'''

Download frames

The input to the script is a json file with the following format. See input_examples

Example:

python /home/telescope/TelescopeConnect/download.py --json=/home/telescope/TelescopeConnect/input_examples/download_files.json

'''

import argparse
from astropy.io import fits
import sys
import numpy as np
from time import strftime
import time
import os
import zipfile, random, string
import shutil

import library.general
from library.general import *

#loading parameter file parser
parser = argparse.ArgumentParser()
parser.add_argument('--json',
                    dest='json_filename',
                    type=str,
                    default=False,
                    )

options = parser.parse_args()

if not options.json_filename:
    RaiseError('You need an inpunt JSON file')

json_data = json.loads(open(options.json_filename).read())

# read all files from list
# images = []
# for image_fname in json_data['input_fits']:
#     images.append(get_ImageData(image_fname))

# TEMPORARY. THIS DOES NOT DO ANYTHING, JUST COMPRESS ALL FILES INTO SINGLE ZIP
random_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))

# create directory for files to be compressed
zip_directory = os.path.abspath(os.path.join(json_data['output_folder'], random_string))
if not os.path.exists(zip_directory):
    os.makedirs(zip_directory)

for image_fname in json_data['input_fits']:

    # set filename
    hdulist = fits.open(image_fname)
    header = hdulist[0].header

    # ccd temperature
    ccdtemp = header['CCD-TEMP']

    # Date and time of observation: keyword DATE-OBS
    dateobs_utc_str, dateobs_utc_datetime = get_DateObs(header['DATE-OBS'])
    filename = dateobs_utc_datetime.strftime('%Y-%m-%dT%H-%M-%S')

    # determine Image type from header IMAGETYP
    imagetype = get_ImageType(header['IMAGETYP'])
    if imagetype == 'LIGHT':
        if 'OBJECT' in header:
            filename += '_' + header['OBJECT']
        filename += header['FILTER']
        filename += '_T'
        filename += str(ccdtemp)
    if imagetype == 'BIAS':
        filename += '_Bias'
    if imagetype == 'MASTER BIAS':
        filename += '_MasterBias'
    if imagetype == 'DARK':
        filename += '_Dark'
        filename += '_T'
        filename += str(ccdtemp)
    if imagetype == 'MASTER DARK':
        filename += '_MasterDark'
        filename += '_T'
        filename += str(ccdtemp)
    if imagetype == 'FLAT':
        filename += '_Flat_'
        filename += header['FILTER']
        filename += '_T'
        filename += str(ccdtemp)
    if imagetype == 'MASTER FLAT':
        filename += '_MasterFlat_'
        filename += header['FILTER']
        filename += '_T'
        filename += str(ccdtemp)

    filename = filename.replace(' ', '_')
    filename_fits = filename + '.fits'

    # copy file to archive directory
    shutil.copy(image_fname, os.path.join(zip_directory, filename_fits))


output_filename_zip = 'TelescopeConnect_' + random_string + '.zip'
output_filename = 'TelescopeConnect_' + random_string

shutil.make_archive(os.path.join(json_data['output_folder'], output_filename), 'zip', zip_directory)

shutil.rmtree(zip_directory)

print(json.dumps({'result': 'SUCCESS',
                  'output_file': os.path.abspath(os.path.join(json_data['output_folder'], output_filename_zip))},
                 separators=(',', ':'), indent=4))
sys.exit()
