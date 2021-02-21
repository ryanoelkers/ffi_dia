""" This set of functions will function for the MAST release. It will prepare both the RAW and CLEAN light curves
in the standard .fits format for MAST."""

from config import Configuration
from libraries.utils import Utils
from astropy.io import fits
import os
import pandas as pd
import numpy as np


class MastRelease:

    @staticmethod
    def make_mast_lc(star_list):
        """ This function will convert all of the .lc files in the give directory to the MAST compliant release form.

        :parameter star_list - The star list for the current sector CCD combination

        :return Nothing is returned, but the new fits files are output to the output directory
        """

        # get the file list of the differenced image flux information
        image_list = Utils.get_file_list(Configuration.DIFFERENCED_DIRECTORY, '-' + Configuration.SECT + '-' +
                                         Configuration.CAMERA + '-' + Configuration.CCD + '-' +
                                         Configuration.SECT_NUM + '-s_ffic_sd.fits.gz', "N")

        # get the image headers for the first and the last image
        img1 = fits.getheader(Configuration.DIFFERENCED_DIRECTORY + image_list[0])
        img2 = fits.getheader(Configuration.DIFFERENCED_DIRECTORY + image_list[-1])

        # the prefix and suffix names for the files
        file_pre = 'hlsp_filtg-imsub_ffi_tic000'
        file_suf = '-s' + Configuration.SECTOR_NUMBER + '_tess_v01_lc.fits'

        # read in the filtergraph file
        tic_file = pd.read_csv(Configuration.RELEASE_SECTOR_DIRECTORY + Configuration.SECTOR + "_" +
                               Configuration.CAMERA + "_" + Configuration.CCD + "_tic7_data.csv",
                               header=0, delimiter=',')
        mast_file = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" +
                                Configuration.CAMERA + "_" + Configuration.CCD + "_master.ap",
                                header=0, delimiter=',')
        filterg = pd.merge(tic_file, mast_file, left_on='ticid', right_on='TICID', suffixes=['', '_master'])

        for idx, row in filterg.iterrows():

            if os.path.exists(Configuration.DETREND_LC_DIRECTORY + str(row.ticid) + "_" + Configuration.SECTOR + "_" +
                              Configuration.CAMERA + "_" + Configuration.CCD + ".lc"):
                
                if idx % 1000 == 0:
                    Utils.log(str(idx) + " files written. " + str(len(filterg) - idx) + " files remain.",
                              "info", Configuration.LOG_SCREEN)

                # read in the raw light curve
                raw_lc = pd.read_csv(Configuration.RAW_LC_DIRECTORY + str(row.ticid) + "_" +
                                     Configuration.SECTOR + "_" + Configuration.CAMERA + "_" +
                                     Configuration.CCD + ".lc", sep=" ", names=['JD', 'raw_mag', 'mag_err'])
                raw_lc.raw_mag.replace('*********', 99.999999, inplace=True)

                # get the flux and the normalized flux
                raw_lc['raw_flux'] = raw_lc.apply(lambda x: 10. ** (-1 * (float(x.raw_mag) - 25.) / 2.5), axis=1)

                # read in the cleaned light curve
                cln_lc = pd.read_csv(Configuration.DETREND_LC_DIRECTORY + str(row.ticid) + "_" +
                                     Configuration.SECTOR + "_" + Configuration.CAMERA + "_" +
                                     Configuration.CCD + ".lc",
                                     sep=" ", names=['JD', 'cln_mag', 'cln_mag_err'])
                cln_lc.cln_mag.replace('*********', 99.999999, inplace=True)

                # get the flux and the normalized flux
                cln_lc['cln_flux'] = cln_lc.apply(lambda x: 10. ** (-1 * (float(x.cln_mag) - 25.) / 2.5), axis=1)

                # merge the two data frames on the time stamps
                lc = pd.merge(raw_lc, cln_lc, on='JD', how='left')

                # get the trend in magnitude
                lc['trend'] = lc.apply(lambda x: float(x.raw_mag) - float(x.cln_mag), axis=1)

                # fill in the null values
                lc = lc.fillna(-99.99999)

                # we need to create a fits file with all of the data columns
                c1 = fits.Column(name='BTJD', array=lc['JD'].to_numpy(), format='F')
                c2 = fits.Column(name='raw_mag', array=lc['raw_mag'].astype(float).to_numpy(), format='F')
                c3 = fits.Column(name='raw_flux', array=lc['raw_flux'].to_numpy(), format='F')
                c4 = fits.Column(name='cln_mag', array=lc['cln_mag'].astype(float).to_numpy(), format='F')
                c5 = fits.Column(name='cln_flux', array=lc['cln_flux'].to_numpy(), format='F')
                c6 = fits.Column(name='mag_err', array=lc['mag_err'].to_numpy(), format='F')
                c7 = fits.Column(name='trend', array=lc['trend'].to_numpy(), format='F')

                # add the date specific information from the images
                mjd_beg = np.around(img1['TSTART'] + img1['BJDREFI'] + img1['BJDREFF'] - 2400000.5, decimals=6)
                mjd_end = np.around(img2['TSTOP'] + img2['BJDREFI'] + img2['BJDREFF'] - 2400000.5, decimals=6)

                # add the necessary MAST information
                hdr = fits.Header()
                hdr.append(('TELESCOP', 'TESS', 'telescope'))
                hdr.append(('INSTRUME', 'TESS Photometer', 'detector type'))
                hdr.append(('FILTER', 'TESS    ', 'filter type'))
                hdr.append(('DATE-OBS', img1['DATE-OBS'], 'TSTART as UTC calendar date'))
                hdr.append(('BJDREFI', img1['BJDREFI'], 'integer part of the BTJD reference date'))
                hdr.append(('BJDREFF', img1['BJDREFF'], 'fraction of the day in the BTJD reference date'))
                hdr.append(('TSTART', img1['TSTART'], 'TSTART of the first frame in BTJD'))
                hdr.append(('TSTOP', img2['TSTOP'], 'TSTOP of the last frame in BTJD'))
                hdr.append(('MJD-BEG', mjd_beg, 'MJD of the first frame'))
                hdr.append(('MJD-END', mjd_end, 'MJD of the last frame'))
                hdr.append(('OBJECT', 'TIC ' + str(int(row.ticid)), 'object'))
                hdr.append(('TICID', int(row.ticid), 'TIC ID of the object'))
                hdr.append(('CAMERA', int(Configuration.CAMERA), 'Camera number object was observed'))
                hdr.append(('CCD', int(Configuration.CCD), 'CCD number object was observed'))
                hdr.append(('SECTOR', int(Configuration.SECTOR_NUMBER), 'Sector number of observation'))
                hdr.append(('XPOSURE', 1800.0, 'exposure time in seconds'))
                hdr.append(('RADESYS', 'ICRS    ', 'Coordinate system'))
                hdr.append(('EQUINOX', 2000.0, 'Equinox of observation'))
                hdr.append(('RA_OBJ', filterg.ra[0], 'Right Ascension of object in degrees'))
                hdr.append(('DEC_OBJ', filterg.dec[0], 'Declination of object in degrees'))
                hdr.append(('DETREND', 'YES', 'Was the light curve detrended'))
                hdr.append(('TICVER', 'v7', 'TIC version of parameters'))
                hdr.append(('X_Pixel', filterg['x'][0], 'X pixel of star centroid'))
                hdr.append(('Y_Pixel', filterg['y'][0], 'Y pixel of star centroid'))
                hdr.append(('InstMag', filterg['mag'][0], 'Instrumental magnitude'))
                hdr.append(('APERErr', filterg['mag_err'][0], 'Photometric error from APER routine'))

                # add the TIC information
                for colum in filterg.columns[13:99]:

                    bd = 1
                    if isinstance(filterg[colum][0], str) == 0:
                        if np.isnan(filterg[colum][0]) == 1:
                            bd = 0

                    if bd == 1:
                        hdr.append((colum, filterg[colum][0], 'see TICv7 for details'))

                # add the descriptions to the header information
                hdr.append(('TDESC1', 'BTJD', 'The mean BTJD time of exposure, in days.'))
                hdr.append(('TDESC2', 'raw_mag', 'Raw instrumental magnitude'))
                hdr.append(('TDESC3', 'raw_flux', 'Raw instrumental flux'))
                hdr.append(('TDESC4', 'cln_mag', 'Detrended instrumental magnitude'))
                hdr.append(('TDESC5', 'cln_flux', 'Detrended instrumental flux'))
                hdr.append(('TDESC6', 'mag_err', 'Photometric error in magnitude'))
                hdr.append(('TDESC7', 'trend', 'Photometric trend in magnitude'))

                # create the FITS file binary table, then write out the light curve
                lc_table = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7], header=hdr)

                # check to make sure the directory is made (if so, then write; if not, make directory then write)
                # make the 16 integer TICID if it doesn't exist
                tic_ids = str(int(row.ticid))
                while len(tic_ids) < 16:
                    tic_ids = '0' + tic_ids

                # split up the IDs
                tic_ids1 = tic_ids[0] + '000'
                tic_ids2 = tic_ids[4] + '000'
                tic_ids3 = tic_ids[8] + '000'
                tic_ids4 = tic_ids[12] + '000'

                # make the directories
                sects = Configuration.RELEASE_SECTOR_DIRECTORY + Configuration.SECT
                tid1 = sects + '/' + tic_ids1
                tid2 = tid1 + '/' + tic_ids2
                tid3 = tid2 + '/' + tic_ids3
                tid4 = tid3 + '/' + tic_ids4

                # make the directories if they do not exist
                mast_dirs = [sects, tid1, tid2, tid3, tid4]
                Utils.create_directories(mast_dirs)

                # write out the file
                lc_table.writeto(tid4 + '/' + file_pre + str(int(row.ticid)) + file_suf, overwrite=True)

        return
