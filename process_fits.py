'''

FPS.1b Part 2 – Processing of raw Light frames

The input to the script is a json file with the following format. See input_examples

Example:

python process_fits.py --json=input_examples/process_fits.json

'''

import argparse
from astropy.io import fits
import sys
import numpy as np
from time import strftime
import os

sys.path.append(os.path.join(os.path.realpath(__file__), 'library'))

import general
from general import *

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

result = 'SUCCESS'

for idx, image in enumerate(images):

    # get output filename
    filename, file_extension = os.path.splitext(os.path.basename(json_data['input_fits'][idx]))
    fname_out = '%s_processed.fits' % filename
    fname_out_fullpath = os.path.join(json_data['output_folder'], fname_out)

    # check if file already exists
    if os.path.isfile(fname_out_fullpath):
        RaiseError('File %s already exists' % fname_out_fullpath)

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
        header_out['COMMENT'] = 'Processed with CORTEX %s on %s' % (get_Version(), strftime("%Y-%m-%dT%H-%M-%S"))
        hdu = fits.PrimaryHDU(images[idx][0], header=header_out)
        hdu.writeto(fname_out_fullpath)

    except:

        result = 'WARNING'
        fname_out_fullpath = ''

        warning_description = 'One or more images have not been processed'

    # append filename to list of fits images
    output_fits.append(fname_out_fullpath)

output = {
    'result': result,
    'output_fits': output_fits,
}

if result == 'WARNING':
    output['warning_description'] = warning_description


print(json.dumps(output, separators=(',',':'), sort_keys=True, indent=4))