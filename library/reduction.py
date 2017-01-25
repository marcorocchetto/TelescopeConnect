
import numpy as np
import logging
import sys
import os
from astropy.io import fits
from time import strftime

from common import *

def create_master(self, frames, bias=None, dark=None, flat=None, tag=None):

    ''' Function create_master_bias

    Create master bias frames. Input frames can be a dictionary of grouped frames (by date/temp)
    or a simple list of bias. If not specified, frames are automatically set to self.bias_frames

    '''

    master = combine(frames)

    if bias:
        pass
    if dark:
        pass
    if flat:
        pass

    # create master frame
    timeproc = strftime("%Y-%m-%dT%H-%M-%S")
    header = get_header(frames[0])
    #header['NCOMBINE'] = len(frames)
    #header['COMMENT'] = 'Processed with robopipe %s on %s' % (get_version(), timeproc)
    hdu = fits.PrimaryHDU(master, header=header)
    fname = os.path.join(self.params.data_working_dir, self.data.rename(frames[0]).replace('.fits', '-master.fits'))
    hdu.writeto(fname)

    print
    frames
    sys.exit()

    return os.path.join(self.params.data_working_dir, fname)


def combine(self, frames, tag=None):

    zerocorr = np.zeros((len(frames), self.params.red_img_y, self.params.red_img_x))

    for idx, fname in enumerate(frames):
        zerocorr[idx, :, :] = get_Image(fname)

    if self.params.red_bias_combine == 'average':
        out = np.average(zerocorr, axis=0)
    elif self.params.red_bias_combine == 'median':
        out = np.median(zerocorr, axis=0)

    return out