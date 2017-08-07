'''

FPS.1a Part 1 – Creation of Master Calibration Files: Master Bias, Master Dark and Master Flat

The input to the script is a json file with the following format. See input_examples

For Master Bias, Dark, Flat creation:

python make_master.py --json=input_examples/make_bias.json
python make_master.py --json=input_examples/make_dark.json
python make_master.py --json=input_examples/make_flat.json

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
for image_path in json_data['input_fits']:
    images.append(get_ImageData(image_path))

# todo: check size of frame

if json_data['make_type'].upper() == 'BIAS'\
        or json_data['make_type'].upper() == 'DARK'\
        or json_data['make_type'].upper() == 'FLAT':

    # load header of first image, will be used as header for output file
    header_out = images[0][1]

    # create filename
    dateobs_utc_str, dateobs_utc_datetime = get_DateObs(header_out['DATE-OBS'])
    filename = dateobs_utc_datetime.strftime('%Y-%m-%dT%H-%M-%S')
    if json_data['make_type'].upper() == 'BIAS':
            filename += '_MasterBias'
    if json_data['make_type'].upper() == 'DARK':
            filename += '_T%sC' % json_data['ccd_temp_avg']
            filename += '_60s'
            filename += '_MasterDark'
    if json_data['make_type'].upper() == 'FLAT':
            filename += '_T%sC' % json_data['ccd_temp_avg']
            filename += '_%s' % json_data['filter']
            filename += '_MasterFlat'

    filename = filename.replace(' ', '-')

    filename_fits = filename + '.fits'
    filename_jpg_large = filename + '.jpg'
    filename_jpg_thumb = filename + '_thumb.jpg'

    path_fits = os.path.abspath(os.path.join(json_data['output_folder'], filename_fits))
    path_jpg_large = os.path.abspath(os.path.join(json_data['output_folder'], filename_jpg_large))
    path_jpg_thumb = os.path.abspath(os.path.join(json_data['output_folder'], filename_jpg_thumb))

    # check ifoutput fits file already exists
    if os.path.isfile(path_fits):
        RaiseError('File %s already exists' % path_fits)

    if json_data['make_type'].upper() == 'DARK':

        # if dark, apply master bias to all raw darks
        master_bias = get_ImageData(json_data['master_bias'])[0]
        for idx, image in enumerate(images):
            images[idx][0][:,:] = images[idx][0]-master_bias
            # todo after bias subtraction we get darks with negative values...

            # normalise to exptime of 60 sec
            exptime = images[idx][1]['EXPTIME']
            expfactor = (60./exptime)
            images[idx][0][:, :] = images[idx][0]*expfactor

    if json_data['make_type'].upper() == 'FLAT':
        # if flat, apply master bias and master dark to all raw flats
        master_bias = get_ImageData(json_data['master_bias'])[0]
        master_dark = get_ImageData(json_data['master_dark'])[0]
        for idx, image in enumerate(images):
            images[idx][0][:, :] = images[idx][0] - master_bias

            # subtract dark; scale dark to exptime of flat
            exptime = images[idx][1]['EXPTIME']
            expfactor = exptime/60.
            images[idx][0][:, :] = images[idx][0] - master_dark*expfactor



    # combine frames to create Master frame

    # create 3-axis array (all images are paced in a cube)
    master_image_all = np.zeros((len(images), np.shape(images[0][0])[0], np.shape(images[0][0])[1]), dtype=float)
    for idx, image in enumerate(images):
        master_image_all[idx, :, :] = image[0] # 0 is imagedata, 1 is header. See get_ImageData

    # get average or median of all corrected input frames to create master
    if json_data['combine_method'] == 'average':
        master_image = np.average(master_image_all, axis=0)
    elif json_data['combine_method'] == 'median':
        master_image = np.median(master_image_all, axis=0)

    # if flat, normalise
    if json_data['make_type'].upper() == 'FLAT':
        master_image = master_image / np.average(master_image)

    # create output fits frame

    # use header of first file and add some comments
    header_out['NCOMBINE'] = len(images)
    header_out['COMMENT'] = 'Processed with CORTEX %s on %s' % (get_Version(), strftime("%Y-%m-%dT%H-%M-%S"))
    if json_data['make_type'].upper() == 'DARK':
        header_out['EXPTIME'] = 60. # exp time is 60 if Master Dark frame

    hdu = fits.PrimaryHDU(master_image, header=header_out)

    hdu.writeto(path_fits)

    os.system("/usr/bin/convert '" + path_fits + "' -linear-stretch 600x1500 -resize 1024x '" + path_jpg_large + "'")
    os.system("/usr/bin/convert '" + path_fits + "' -linear-stretch 600x1500 -resize 100x '" + path_jpg_thumb + "'")

    output = {
        'result': 'SUCCESS',
        'output_fits': path_fits,
        'output_jpg_large': path_jpg_large,
        'output_jpg_thumb': path_jpg_thumb,
        'make_type': json_data['make_type'].upper(),
    }

print(json.dumps(output, separators=(',',':'), sort_keys=True, indent=4))