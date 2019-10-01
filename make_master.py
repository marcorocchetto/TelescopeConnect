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
            filename += '_T%sC' % int(json_data['ccd_temp_avg'])
            filename += '_%is' % int(round(images[0][1]['EXPTIME']))
            filename += '_MasterDark'

    if json_data['make_type'].upper() == 'FLAT':
            filename += '_%s' % json_data['filter'].replace("'", "prime")
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

    master_image_all = np.zeros((len(images), np.shape(images[0][0])[0], np.shape(images[0][0])[1]), dtype=float)

    for idx, image in enumerate(images):
        master_image_all[idx, :, :] = image[0] # 0 is imagedata, 1 is header. See get_ImageData

    if json_data['make_type'].upper() == 'DARK':
        master_bias = np.zeros((np.shape(images[0][0])[0], np.shape(images[0][0])[1]), dtype=float)

        # if dark, apply master bias to all raw darks
        master_bias[:,:] = get_ImageData(json_data['master_bias'])[0]

        # pedestal
        pedestal = 100

        for idx, image in enumerate(images):

            # remove bias
            master_image_all[idx, :, :] = images[idx][0] + pedestal - master_bias

    if json_data['make_type'].upper() == 'FLAT':

        master_dark = np.zeros((np.shape(images[0][0])[0], np.shape(images[0][0])[1]), dtype=float)
        master_bias = np.zeros((np.shape(images[0][0])[0], np.shape(images[0][0])[1]), dtype=float)

        # if flat, apply master bias and master dark to all raw flats
        master_bias[:,:] = get_ImageData(json_data['master_bias'])[0]
        master_bias[:,:] = master_bias.astype(np.float)


        for idx, image in enumerate(images):
            #remove bias
            master_image_all[idx, :, :] = master_image_all[idx, :, :] - master_bias

            # get master dark
            if images[idx][1]['EXPTIME'] in darks_exp.keys():

                # master dark same exptime as frame
                master_dark[:,:] = darks_exp[images[idx][1]['EXPTIME']][0]

            else:
                # scale master dark
                # take first master dark
                # todo we should actually take master dark with closesest exptime
                dark_exptime = darks[0][1]['EXPTIME']
                image_exptime = images[idx][1]['EXPTIME']
                master_dark_unscaled = np.zeros((np.shape(images[0][0])[0], np.shape(images[0][0])[1]), dtype=float)
                master_dark_unscaled[:,:] = darks[0][0]
                master_dark[:,:] = master_dark_unscaled * (image_exptime/dark_exptime)

            # subtract dark;
            master_image_all[idx, :, :] = master_image_all[idx, :, :] - master_dark

            # normalise
            master_image_all[idx, :, :] = master_image_all[idx, :, :]  / np.average(master_image_all[idx, :, :])


    # combine frames to create Master frame

    # get average or median of all corrected input frames to create master
    if json_data['combine_method'] == 'average':
        master_image = np.average(master_image_all, axis=0)

    elif json_data['combine_method'] == 'median':
        master_image = np.median(master_image_all, axis=0)

    if json_data['make_type'].upper() == 'FLAT':
        # rescale to 30k counts
        master_image = 3e4 * master_image / np.average(master_image)

    # create output fits frame

    # use header of first file and add some comments
    # todo this is not ideal, we should modify DATE-OBS and EXPTIME to match master frame

    header_out['NCOMBINE'] = len(images)
    header_out['MASTER'] = True

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