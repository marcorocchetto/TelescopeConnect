'''

Creation of Master Calibration Files: Master Bias, Master Dark and Master Flat

The input to the script is a json file with the following format. See input_examples

For Master Bias, Dark, Flat creation:

python /home/telescope/TelescopeConnect/make_master.py --json=/home/telescope/TelescopeConnect/input_examples/make_bias.json
python /home/telescope/TelescopeConnect/make_master.py --json=/home/telescope/TelescopeConnect/input_examples/make_dark.json
python /home/telescope/TelescopeConnect/make_master.py --json=/home/telescope/TelescopeConnect/input_examples/make_flat.json

'''

import time

begin = time.time()

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

# read all darks from list if type is FLATS
if json_data['make_type'].upper() == 'FLAT':
    darks = []
    for image_path in json_data['master_dark']:
        darks.append(get_ImageData(image_path))


if json_data['make_type'].upper() == 'BIAS'\
        or json_data['make_type'].upper() == 'DARK'\
        or json_data['make_type'].upper() == 'FLAT':

    if json_data['make_type'].upper() == 'DARK':

        # check that we have uniform exposure times
        for image_idx, image_val in enumerate(images):
            if image_idx == 0:
                continue
            if image_val[1]['EXPTIME'] != images[0][1]['EXPTIME']:
                RaiseError('Dark frames do not have the same exposure time, cannot continue')

    if json_data['make_type'].upper() == 'FLAT':

        # exp times of flats and darks
        exptimes_flats = [image[1]['EXPTIME'] for image in images]
        exptimes_darks = [image[1]['EXPTIME'] for image in darks]

        darks_exp = {}
        for image in darks:
            darks_exp[image[1]['EXPTIME']] = image


    # load header of first image, will be used as header for output file
    header_out = images[0][1]

    # create filename
    dateobs_utc_str, dateobs_utc_datetime = get_DateObs(header_out['DATE-OBS'])
    filename = dateobs_utc_datetime.strftime('%Y-%m-%dT%H-%M-%S')

    if json_data['make_type'].upper() == 'BIAS':
            filename += '_MasterBias'

    if json_data['make_type'].upper() == 'DARK':
            filename += '_T%sC' % json_data['ccd_temp_avg']
            filename += '_%is' % int(round(images[0][1]['EXPTIME']))
            filename += '_MasterDark'

    if json_data['make_type'].upper() == 'FLAT':
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
    # if os.path.isfile(path_fits):
    #     RaiseError('File %s already exists' % path_fits)


    if json_data['make_type'].upper() == 'DARK':

        # if dark, apply master bias to all raw darks
        master_bias = get_ImageData(json_data['master_bias'])[0]
        for idx, image in enumerate(images):

            # remove bias
            images[idx][0][:,:] = images[idx][0]-master_bias

            # normalise to exptime of 60 sec
            # exptime = images[idx][1]['EXPTIME']
            # expfactor = (60./exptime)
            # images[idx][0][:, :] = images[idx][0]*expfactor

    if json_data['make_type'].upper() == 'FLAT':

        # if flat, apply master bias and master dark to all raw flats
        master_bias = get_ImageData(json_data['master_bias'])[0]

        for idx, image in enumerate(images):
            #remove bias
            images[idx][0][:, :] = images[idx][0] - master_bias


            if images[idx][1]['EXPTIME'] in darks_exp.keys():
                # master dark same exptime as frame
                master_dark = darks_exp[images[idx][1]['EXPTIME']][0]
            else:
                # scale master dark
                # take first master dark
                # todo we should actually take master dark with closesest exptime
                dark_exptime = next (iter (darks_exp.values()))
                image_exptime = images[idx][1]['EXPTIME']
                master_dark_unscaled = darks_exp[dark_exptime][0]
                master_dark = master_dark_unscaled * (image_exptime/dark_exptime)

            # scale flat to 60sec exp; subtract dark;
            images[idx][0][:, :] = images[idx][0][:, :] - master_dark
            images[idx][0][:, :][images[idx][0][:, :] < 0] = 0.

            # normalise
            images[idx][0][:, :] = images[idx][0][:, :]/np.average(images[idx][0][:, :])

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

    if json_data['make_type'].upper() == 'FLAT':
        # rescale to ~ 30k counts max
        master_image = 3e4 * master_image / np.average(master_image)

    # create output fits frame

    # use header of first file and add some comments
    header_out['NCOMBINE'] = len(images)
    header_out['COMMENT'] = 'Processed with CORTEX %s on %s' % (get_Version(), strftime("%Y-%m-%dT%H-%M-%S"))
    # if json_data['make_type'].upper() == 'DARK' or json_data['make_type'].upper() == 'FLAT':
    #     header_out['EXPTIME'] = 60. # exp time is 60 if Master Dark frame or Master Flat

    hdu = fits.PrimaryHDU(master_image.astype(np.uint16), header=header_out)

    hdu.writeto(path_fits, overwrite=True)

    os.system("/usr/bin/convert '" + path_fits + "' -contrast-stretch 3%  -resize 1024x '" + path_jpg_large + "'")
    os.system("/usr/bin/convert '" + path_fits + "' -contrast-stretch 3% -resize 100x '" + path_jpg_thumb + "'")

    output = {
        'result': 'SUCCESS',
        'output_fits': path_fits,
        'output_jpg_large': path_jpg_large,
        'output_jpg_thumb': path_jpg_thumb,
        'make_type': json_data['make_type'].upper(),
    }

print(json.dumps(output, separators=(',',':'), sort_keys=True, indent=4))