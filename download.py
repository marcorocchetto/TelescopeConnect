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

if not 'preprocess' in json_data or not json_data['preprocess']:
    json_data['preprocess'] = []
if not 'stack_type' in json_data or not json_data['stack_type']:
    json_data['stack_type'] = []

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


# stats_filename = os.path.join(zip_directory, 'stats.txt')
# stats_file = open(stats_filename, 'w')

# preprocessing
#
#

print('****** MEMORY MB %.1f' % (float(process.memory_info().rss) / 1024 / 1024))

# sort all files
all_images = []
if 'input_fits' in json_data:
    all_images += json_data['input_fits']
if 'input_calib' in json_data:
    all_images += json_data['input_calib']

all_images_dict = {
    'Flat': [],
    'Dark': [],
    'Bias': [],
    'Master': [],
    'Processed': [],
    'Raw': []
}
for image_fname in all_images:

    # get header
    hdulist = fits.open(image_fname)
    header = hdulist[0].header

    # determine Image type from header IMAGETYP
    imagetype = get_ImageType(header['IMAGETYP'])
    if imagetype == 'LIGHT':
        if 'PROC' in header:
            if header['PROC']:
                all_images_dict['Processed'].append(image_fname)
            else:
                all_images_dict['Raw'].append(image_fname)
        else:
            all_images_dict['Raw'].append(image_fname)
    elif 'NCOMBINE' in header: # ncombine is only inserted for Master frames
        all_images_dict['Master'].append(image_fname)
    elif imagetype == 'BIAS':
        all_images_dict['Bias'].append(image_fname)
    elif imagetype == 'DARK':
        all_images_dict['Dark'].append(image_fname)
    elif imagetype == 'FLAT':
        all_images_dict['Flat'].append(image_fname)

# if we specify Master or Raw calibs, Raw frames, then we need to deliver them
for imgtype in ['Flat', 'Dark', 'Bias', 'Master', 'Raw']:

    if len(all_images_dict[imgtype]) > 0:

        os.makedirs(os.path.join(zip_directory, imgtype))

        for image_fname in all_images_dict[imgtype]:

            hdulist = fits.open(image_fname)
            header = hdulist[0].header
            filename = os.path.splitext(os.path.basename(image_fname))[0]
            filename = filename.replace(' ', '_')
            output_filename_fits = filename + '.fits'
            output_filename_tiff = filename + '.tiff'
            output_filename_jpg = filename + '.jpg'

            if json_data['output_format'] == 'fits':
                shutil.copy(image_fname, os.path.join(zip_directory, imgtype, output_filename_fits))
            elif json_data['output_format'] == 'tiff':
                os.system("/usr/bin/convert '" + image_fname + \
                          "' -contrast-stretch 50% '" + os.path.join(zip_directory, imgtype, output_filename_tiff) + "'")
            elif json_data['output_format'] == 'jpg':
                os.system("/usr/bin/convert '" + image_fname + \
                          "' -contrast-stretch 20% '" + os.path.join(zip_directory, imgtype, output_filename_jpg) + "'")

print('****** MEMORY MB %.1f' % (float(process.memory_info().rss) / 1024 / 1024))

# no align, no stack. Deliver original Processed files if present
if not 'align' in json_data['preprocess'] and not 'stack' in json_data['preprocess']:

    if len(all_images_dict['Processed']) > 0:

        imgtype = 'Processed'

        os.makedirs(os.path.join(zip_directory, imgtype))

        for image_fname in all_images_dict[imgtype]:

            hdulist = fits.open(image_fname)
            header = hdulist[0].header
            filename = os.path.splitext(os.path.basename(image_fname))[0]
            filename = filename.replace(' ', '_')
            output_filename_fits = filename + '.fits'
            output_filename_tiff = filename + '.tiff'
            output_filename_jpg = filename + '.jpg'

            if json_data['output_format'] == 'fits':
                shutil.copy(image_fname, os.path.join(zip_directory, imgtype, output_filename_fits))
            elif json_data['output_format'] == 'tiff':
                os.system("/usr/bin/convert '" + image_fname + \
                          "' -contrast-stretch 50% '" + os.path.join(zip_directory, imgtype, output_filename_tiff) + "'")
            elif json_data['output_format'] == 'jpg':
                os.system("/usr/bin/convert '" + image_fname + \
                          "' -contrast-stretch 20% '" + os.path.join(zip_directory, imgtype, output_filename_jpg) + "'")

    print('****** MEMORY MB %.1f' % (float(process.memory_info().rss) / 1024 / 1024))


# align
aligned_success = False
if 'align' in json_data['preprocess'] or 'stack' in json_data['preprocess']:

    if not len(all_images_dict['Processed']) > 1:

        # cannot align images if we have less than 1 processed image
        exit()

    else:

        print("Align")

        import alipy

        aligned_directory = os.path.join(working_directory, 'aligned')
        if not os.path.isdir(aligned_directory):
            os.mkdir(aligned_directory)

        # all images and reference image
        images = all_images_dict['Processed']
        ref_image = images[0] # ref image is first one

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

            os.makedirs(os.path.join(zip_directory, 'Aligned'))

            for image_fname in images_aligned:
                if json_data['output_format'] == 'fits':
                   shutil.copy(image_fname, os.path.join(zip_directory, 'Aligned'))

                elif json_data['output_format'] == 'tiff':
                    output_filename_tiff = '%s.tiff' % os.path.splitext(os.path.basename(image_fname))[0]

                    os.system("/usr/bin/convert '" + image_fname + \
                              "' -contrast-stretch 50% '" + os.path.join(zip_directory, 'Aligned', output_filename_tiff) + "'")

                elif json_data['output_format'] == 'jpg':
                    output_filename_jpg = '%s.jpg' % os.path.splitext(os.path.basename(image_fname))[0]
                    os.system("/usr/bin/convert '" + image_fname + \
                              "' -contrast-stretch 50% '" + os.path.join(zip_directory, 'Aligned', output_filename_jpg) + "'")

        aligned_success = True
        print('****** MEMORY MB %.1f' % (float(process.memory_info().rss) / 1024 / 1024))



# align and stack (images already aligned in step above)
if 'stack' in json_data['preprocess']:


    if not len(all_images_dict['Processed']) > 1:

        # cannot stack images if we have less than 1 processed image
        exit()

    elif not aligned_success:
        # cannot stack images if aligned failed
        exit()
    else:

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

            # get header of first image
            header_out = get_ImageData(filters[filter][0])[1]

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

            elif stack_method == 'mean' or stack_method == 'average':
                stack_final[:, :] = np.mean(stack_array, axis=2)
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

            header_out['NCOMBINE'] = nimages
            header_out['EXPTIME'] = exp_total

            # todo need to fix OBSTIME for stacked image.
            # todo we should rerun plate solution for this field

            hdu = fits.PrimaryHDU(stack_final, header=header_out)
            hdul = fits.HDUList([hdu])
            print('****** MEMORY MB %.1f' % (float(process.memory_info().rss) / 1024 / 1024))

            if not os.path.isfile(os.path.join(stack_directory, output_filename_fits)):
                hdul.writeto(os.path.join(stack_directory, output_filename_fits))

            os.makedirs(os.path.join(zip_directory, 'Stacked'))

            # convert to tiff or jpg if needed, save output to zip folder
            if json_data['output_format'] == 'tiff':
                os.system("/usr/bin/convert '" + os.path.join(stack_directory, output_filename_fits) + \
                         "' -contrast-stretch 50% '" + os.path.join(zip_directory, 'Stacked', output_filename_tiff) + "'")
            elif json_data['output_format'] == 'jpeg' or json_data['output_format'] == 'jpg':
                os.system("/usr/bin/convert '" + os.path.join(stack_directory, output_filename_fits) + \
                          "' -contrast-stretch 50% '" + os.path.join(zip_directory, 'Stacked', output_filename_jpg) + "'")
            else:
                shutil.copy(os.path.join(stack_directory, output_filename_fits),
                        os.path.join(zip_directory, 'Stacked', output_filename_fits))

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