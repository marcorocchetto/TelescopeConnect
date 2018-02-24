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
from astropy.stats import sigma_clip
import shutil
import glob

import psutil


process = psutil.Process(os.getpid())


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

if not 'preprocess' in json_data:
    json_data['preprocess'] = []

random_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))

# create directory for files to be compressed
zip_directory = os.path.abspath(os.path.join(json_data['output_folder'], random_string))
if not os.path.exists(zip_directory):
    os.makedirs(zip_directory)

# create working directory
working_directory = os.path.abspath(os.path.join(json_data['output_folder'], 'wdir'))
if not os.path.exists(working_directory):
    os.makedirs(working_directory)

# create directory for copy of original files
orig_directory = os.path.abspath(os.path.join(json_data['output_folder'], 'orig'))
if not os.path.exists(orig_directory):
    os.makedirs(orig_directory)

# logging to text file, in zip folder
log_filename = os.path.join(zip_directory, 'log.txt')
sys.stdout = open(log_filename, "w")

# preprocessing

# first copy files to working directory
for image_fname in json_data['input_fits']:

    # set filename
    hdulist = fits.open(image_fname)
    header = hdulist[0].header

    # ccd temperature
    ccdtemp = round(header['CCD-TEMP'])

    # Date and time of observation: keyword DATE-OBS
    dateobs_utc_str, dateobs_utc_datetime = get_DateObs(header['DATE-OBS'])
    filename = dateobs_utc_datetime.strftime('%Y-%m-%dT%H-%M-%S')

    # determine Image type from header IMAGETYP
    imagetype = get_ImageType(header['IMAGETYP'])
    if imagetype == 'LIGHT':
        if 'OBJECT' in header:
            filename += '_' + header['OBJECT']
        filename += '_' + header['FILTER']
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

    # copy file to orig directory
    shutil.copy(image_fname, os.path.join(orig_directory, filename_fits))

print('****** MEMORY MB %.1f' % (float(process.memory_info().rss) / 1024 / 1024))

# no align, no stack. Deliver original files
if not 'align' in json_data['preprocess'] and not 'stack' in json_data['preprocess']:
    if not 'stack' in json_data['preprocess']:
        # if you are not stacking but only alignining, then output will be aligned images
        images = glob.glob(os.path.join(orig_directory, '*.fits'))
        for image_fname in images:
            if json_data['output_format'] == 'fits':
               shutil.copy(image_fname, zip_directory)
            elif json_data['output_format'] == 'tiff':
                output_filename_tiff = '%s.tiff' % os.path.splitext(os.path.basename(image_fname))[0]
                os.system("/usr/bin/convert '" + os.path.join(zip_directory, image_fname) + \
                          "' -contrast-stretch 50% '" + os.path.join(zip_directory, output_filename_tiff) + "'")
            elif json_data['output_format'] == 'jpg':
                output_filename_jpg = '%s.jpg' % os.path.splitext(os.path.basename(image_fname))[0]
                os.system("/usr/bin/convert '" + os.path.join(zip_directory, image_fname) + \
                          "' -contrast-stretch 20% '" + os.path.join(zip_directory, output_filename_jpg) + "'")
            print('****** MEMORY MB %.1f' % (float(process.memory_info().rss) / 1024 / 1024))

# align
if 'align' in json_data['preprocess'] or 'stack' in json_data['preprocess']:
    print("Align")

    import alipy

    aligned_directory = os.path.join(working_directory, 'aligned')
    if not os.path.isdir(aligned_directory):
        os.mkdir(aligned_directory)

    # all images and reference image
    images = glob.glob(os.path.join(orig_directory, '*.fits'))
    ref_image = images[0]

    # reference image
    identifications = alipy.ident.run(ref_image, images, visu=False, stdout_file=log_filename)

    for id in identifications:
        if id.ok == True:
            print("%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio))
        else:
            print("%20s : no transformation found !" % (id.ukn.name))
    print('****** MEMORY MB %.1f' % (float(process.memory_info().rss) / 1024 / 1024))

    outputshape = alipy.align.shape(images[0])
    for id in identifications:
        if id.ok == True:
            alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, outdir=aligned_directory, makepng=False)
    print('****** MEMORY MB %.1f' % (float(process.memory_info().rss) / 1024 / 1024))

    print("Prepare for copy")

    if not 'stack' in json_data['preprocess']:
        # if you are not stacking but only alignining, then output will be aligned images
        images_aligned = glob.glob(os.path.join(aligned_directory, '*.fits'))
        for image_fname in images_aligned:
            if json_data['output_format'] == 'fits':
               shutil.copy(image_fname, zip_directory)

            elif json_data['output_format'] == 'tiff':
                output_filename_tiff = '%s.tiff' % os.path.splitext(os.path.basename(image_fname))[0]
                os.system("/usr/bin/convert '" + os.path.join(zip_directory, image_fname) + \
                          "' -contrast-stretch 50% '" + os.path.join(zip_directory, output_filename_tiff) + "'")
            elif json_data['output_format'] == 'jpg':
                output_filename_jpg = '%s.jpg' % os.path.splitext(os.path.basename(image_fname))[0]
                os.system("/usr/bin/convert '" + os.path.join(zip_directory, image_fname) + \
                          "' -contrast-stretch 50% '" + os.path.join(zip_directory, output_filename_jpg) + "'")

    print('****** MEMORY MB %.1f' % (float(process.memory_info().rss) / 1024 / 1024))

# align and stack (images already aligned in step above)
if 'stack' in json_data['preprocess']:

    print("Stack")

    stack_directory = os.path.join(working_directory, 'stacked')
    if not os.path.isdir(stack_directory):
        os.mkdir(stack_directory)

    images_aligned = glob.glob(os.path.join(aligned_directory, '*.fits'))

    ref_image = get_ImageData(images_aligned[0])

    filters = {}

    # stacking method
    stack_method = json_data['stack_type']

    shape = np.shape(ref_image[0])

    # sort images by filter
    for image in images_aligned:
        load_image = get_ImageData(image)
        header = load_image[1]
        filter = header['FILTER']
        if not filter in filters:
            filters[filter] = []
        filters[filter].append(image)
    print('****** MEMORY MB %.1f' % (float(process.memory_info().rss) / 1024 / 1024))


    # loop through each filter
    for filter in filters:

        print("Stack filter %s" % filter)

        nimages = len(filters[filter])
        stack_array = np.zeros((shape[0], shape[1], nimages), dtype=np.uint16)

        exp_total = 0

        for image_idx, image_val in enumerate(filters[filter]):
            load_image = get_ImageData(image_val)
            exp_total += np.float(load_image[1]['EXPTIME'])
            stack_array[:, :, image_idx] = load_image[0]

        exp_total_round = round(exp_total / 60.)
        print('****** MEMORY MB %.1f' % (float(process.memory_info().rss) / 1024 / 1024))

        stack_final = np.zeros((shape[0], shape[1]), dtype=np.uint16)

        if stack_method == 'median':
            stack_final[:, :] = np.median(stack_array, axis=2)
        elif stack_method == 'mean':
            stack_final[:, :] = np.average(stack_array, axis=2)
        elif stack_method == 'sum':
            stack_final[:, :] = np.sum(stack_array, axis=2)
        print('****** MEMORY MB %.1f' % (float(process.memory_info().rss) / 1024 / 1024))

        output_filename = ref_image[1]['OBJECT'].replace(' ', '_')
        output_filename += '_%s' % filter
        output_filename += '_stack-%s_%imin' % (stack_method, exp_total_round)
        print(output_filename)

        output_filename_fits = '%s.fits' % output_filename
        output_filename_tiff = '%s.tiff' % output_filename
        output_filename_jpg = '%s.jpg' % output_filename

        hdu = fits.PrimaryHDU(stack_final)
        hdul = fits.HDUList([hdu])
        print('****** MEMORY MB %.1f' % (float(process.memory_info().rss) / 1024 / 1024))

        if not os.path.isfile(os.path.join(stack_directory, output_filename_fits)):
            hdul.writeto(os.path.join(stack_directory, output_filename_fits))

        # convert to tiff or jpg if needed, save output to zip folder
        if json_data['output_format'] == 'tiff':
            os.system("/usr/bin/convert '" + os.path.join(stack_directory, output_filename_fits) + \
                     "' -contrast-stretch 50% '" + os.path.join(zip_directory, output_filename_tiff) + "'")
        elif json_data['output_format'] == 'jpeg' or json_data['output_format'] == 'jpg':
            os.system("/usr/bin/convert '" + os.path.join(stack_directory, output_filename_fits) + \
                      "' -contrast-stretch 50% '" + os.path.join(zip_directory, output_filename_jpg) + "'")
        else:
            shutil.copy(os.path.join(stack_directory, output_filename_fits),
                    os.path.join(zip_directory, output_filename_fits))

print('****** MEMORY MB %.1f' % (float(process.memory_info().rss) / 1024 / 1024))

sys.stdout = sys.__stdout__

output_filename_zip = 'TelescopeConnect_' + random_string + '.zip'
output_filename = 'TelescopeConnect_' + random_string

shutil.make_archive(os.path.join(json_data['output_folder'], output_filename), 'zip', zip_directory)
shutil.rmtree(zip_directory)

print(json.dumps({'result': 'SUCCESS',
                  'output_file': os.path.abspath(os.path.join(json_data['output_folder'], output_filename_zip))},
                 separators=(',', ':'), indent=4))

sys.exit()