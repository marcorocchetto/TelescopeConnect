'''

Processing of raw Light frames

The input to the script is a json file with the following format. See input_examples

Example:

python /home/telescope/TelescopeConnect/process_fits.py --json=/home/telescope/TelescopeConnect/input_examples/process_fits.json

'''

import argparse
from astropy.io import fits
import sys
import numpy as np
from time import strftime
import os

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
images = []
for image_fname in json_data['input_fits']:
    images.append(get_ImageData(image_fname))

# todo: check size of frame

output_fits = []
output_jpg_large = []
output_jpg_thumb = []

result = 'SUCCESS'

for idx, image in enumerate(images):

    # get output filenames for fits and jpg previews
    filename, file_extension = os.path.splitext(os.path.basename(json_data['input_fits'][idx]))

    filename_fits = '%s_processed.fits' % filename
    path_fits = os.path.abspath(os.path.join(json_data['output_folder'], filename_fits))

    filename_jpg_large = '%s_processed.jpg' % filename
    path_jpg_large  = os.path.abspath(os.path.join(json_data['output_folder'], filename_jpg_large))

    filename_jpg_thumb = '%s_processed_thumb.jpg' % filename
    path_jpg_thumb = os.path.abspath(os.path.join(json_data['output_folder'], filename_jpg_thumb))


    # check if output fits file already exists
    if os.path.isfile(filename_fits):
        RaiseError('File %s already exists' % filename_fits)

    try:

        # subtract MasterBias, if present
        if "master_bias" in json_data:
            master_bias = get_ImageData(json_data['master_bias'])[0]
            images[idx][0][:, :] = images[idx][0] - master_bias

        # subtract MasterDark, if present
        if "master_dark" in json_data:
            master_dark = get_ImageData(json_data['master_dark'])[0]
            images[idx][0][:, :] = images[idx][0] - master_dark
            exptime = images[idx][1]['EXPTIME']
            expfactor = exptime / 60.
            images[idx][0][:, :] = images[idx][0] - master_dark * expfactor

        # divide by MasterFlat
        if "master_flat" in json_data:
            master_flat = get_ImageData(json_data['master_flat'])[0]
            images[idx][0][:, :] = images[idx][0] / master_flat

        # use header of first file and add some comments
        header_out = images[0][1]
        header_out['PROC'] = 'True'
        header_out['COMMENT'] = 'Processed %s on %s' % (get_Version(), strftime("%Y-%m-%dT%H-%M-%S"))

        hdu = fits.PrimaryHDU(images[idx][0], header=header_out)
        hdu.writeto(path_fits)

        # generate previews
        os.system("/usr/bin/convert '" + path_fits + "' -linear-stretch 600x1500 -resize 1024x '" + path_jpg_large + "'")
        os.system("/usr/bin/convert '" + path_fits + "' -linear-stretch 600x1500 -resize 100x '" + path_jpg_thumb + "'")

    except:

        result = 'WARNING'
        path_fits = ''
        path_jpg_large = ''
        path_jpg_thumb = ''
        warning_description = 'One or more images have not been processed'

    # append filename to list of fits images
    output_fits.append(path_fits)
    output_jpg_large.append(path_jpg_large)
    output_jpg_thumb.append(path_jpg_thumb)

output = {
    'result': result,
    'output_fits': output_fits,
    'output_jpg_large': output_jpg_large,
    'output_jpg_thumb': output_jpg_thumb,
}

if result == 'WARNING':
    output['warning_description'] = warning_description


print(json.dumps(output, separators=(',',':'), sort_keys=True, indent=4))