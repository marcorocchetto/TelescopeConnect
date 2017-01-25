#from library_general import *

import numpy as np
import logging
import sys
import os
from astropy.io import fits
from time import strftime


class reduction(object):

    def __init__(self, data):

        self.data = data

        # divide frames by date, temp etc
        self.light_frames = self.group_frames_by_date_temp(self.data.all_light)
        self.bias_frames = self.group_frames_by_date_temp(self.data.all_bias)
        self.dark_frames = self.group_frames_by_date_temp(self.data.all_dark)
        self.flat_frames = self.group_frames_by_date_temp(self.data.all_flat)


    def create_master(self, frames, bias=None, dark=None, flat=None, tag=None):

        ''' Function create_master_bias

        Create master bias frames. Input frames can be a dictionary of grouped frames (by date/temp)
        or a simple list of bias. If not specified, frames are automatically set to self.bias_frames

        '''

        master = self.combine(frames)

        if bias:
            pass
        if dark:
            pass
        if flat:
            pass

        # create master frame
        timeproc = strftime("%Y-%m-%dT%H-%M-%S")
        header = get_header(frames[0])
        header['NCOMBINE'] = len(frames)
        header['COMMENT'] = 'Processed with robopipe %s on %s' % (get_version(), timeproc)
        hdu = fits.PrimaryHDU(master, header=header)
        fname = os.path.join(self.params.data_working_dir, self.data.rename(frames[0]).replace('.fits', '-master.fits'))
        hdu.writeto(fname)

        print frames
        sys.exit()

        return os.path.join(self.params.data_working_dir, fname)

    def combine(self, frames, tag=None):

        zerocorr = np.zeros((len(frames), self.params.red_img_y, self.params.red_img_x))
        for idx, fname in enumerate(frames):
            zerocorr[idx,:,:] = get_image(fname)

        if self.params.red_bias_combine == 'average':
            out = np.average(zerocorr, axis=0)
        elif self.params.red_bias_combine == 'median':
            out = np.median(zerocorr, axis=0)

        return out

    def group_frames_by_date_temp(self, frames):

        # todo need to reject bias/dark/flat if min_nbias < param etc

        out_frames = {}

        for fname in frames:

            frame_dict = {}
            frame_dict['path'] = fname

            header = get_header(fname)
            imgtype = self.data.get_frame_imgtyp(header)

            # get date of observation
            timestamp = self.data.get_timestamp_from_dateobs(header)
            daten = datenight(timestamp)  #convert obs-date in night of observation (see function datenight in library)
            if not daten in out_frames:
                out_frames[daten] = {}

            # get ccd temperature
            ccdtemp = get_headerval_from_keywords(header, self.params.header_ccdtemp)
            ccdtemp_round = round_base(ccdtemp, self.params.red_ccdtemp_round) #round ccdtemp
            if not ccdtemp_round in out_frames[daten]:
                out_frames[daten][ccdtemp_round] = []

            # get filter fot light and flat frames
            if imgtype == 'light' or imgtype == 'flat':
                filter = get_headerval_from_keywords(header, self.params.header_filter)
                frame_dict['filter'] = filter

            # create a dictionary for the frame
            frame_dict['timestamp'] = timestamp
            frame_dict['datenight'] = daten
            frame_dict['ccdtemp'] = ccdtemp
            out_frames[daten][ccdtemp_round].append(frame_dict)

        return out_frames

