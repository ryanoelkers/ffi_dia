""" This script will work on the TESS FFIs from year 1, which were created using IDL. This script will need to be
updated slightly for future versions. Please confirm directory location of files etc. """

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

# set this flag to be SPEARMAN or ZEROPOINT depending on how you want to detrend
METHOD = 'SPEARMAN'
SPEARMAN_CUT = 0.85

# change this for the sector camera ccd you are looking for
SECTOR = 'sector01'
CAMERA = '1'
CCD = '1'

# these can be updated for whatever you want
TREND_STARS = 500
MG_CLIP = 13

# change these for the specific directories
# SECTOR_DIRECTORY = '/media/oelkerrj/Yavin/TESS/data/'
SECTOR_DIRECTORY = '/net/jovan/export/jovan/oelkerrj/FFI/'
# MASTER_DIRECTORY = SECTOR_DIRECTORY + SECTOR + '/master/'
MASTER_DIRECTORY = SECTOR_DIRECTORY + SECTOR + '/master/star_list/'
RAW_DIRECTORY = SECTOR_DIRECTORY + SECTOR + '/lc/raw/'
DETREND_DIRECTORY = SECTOR_DIRECTORY + SECTOR + '/lc/detrend_test/'
OUTPUT_DIRECTORY = SECTOR_DIRECTORY + SECTOR + '/output/'  # UPDATE THIS DIRECTORY TO BE WHERE YOU WANT TO WRITE TO

# read in the master magnitude file from the master frame directory
master = pd.read_csv(MASTER_DIRECTORY + SECTOR + '_' + CAMERA + '_' + CCD + '_master.ap',
                     delimiter=' ', names=['TICID', 'x', 'y', 'tessmag', 'mag', 'er', 'flux', 'flux_er'])
# master = pd.read_csv(MASTER_DIRECTORY + SECTOR + '_' + CAMERA + '_' + CCD + '_master.data',
#                     delimiter=',', index_col=0)

# push TICIDs to a list
ticid = master['TICID'].to_list()
master_mags = master['mag'].to_numpy()
zpt_ticid = master[master['tessmag'] < MG_CLIP].TICID.to_list()
zpt_mags = master[master['tessmag'] < MG_CLIP]['mag'].to_numpy()

# set up the light curve vectors based on the method you are using
if METHOD == 'SPEARMAN':
    # read in the necessary light curves
    filenames = [RAW_DIRECTORY + str(x) + "_" + SECTOR + "_" + CAMERA + "_" + CCD + ".lc" for x in ticid]

    # read in the raw lightcurves to find the zeropoint
    # lc_data = np.array([np.genfromtxt(f, delimiter=' ')[:, 2] for f in filenames[0:TREND_STARS]])
    lc_data = np.array([np.genfromtxt(f, delimiter=' ')[:, 1] for f in filenames[0:TREND_STARS]])
    # subtract the master frame magnitude
    median_df = pd.DataFrame(lc_data).replace('*********', 99.999999).apply(lambda x: x-x.median(), axis=1)


if METHOD == 'ZEROPOINT':
    filenames = [DETREND_DIRECTORY + str(x) + "_" + SECTOR + "_" + CAMERA + "_" + CCD + ".lc" for x in zpt_ticid]

    # read in the raw lightcurves to find the zeropoint
    lc_data = np.array([np.genfromtxt(f, delimiter=' ')[:, 2] for f in filenames])

    # subtract the master frame magnitude
    median_df = pd.DataFrame(lc_data).replace('*********', 99.999999) - zpt_mags.reshape(-1, 1)

    # calculat the zeropoint offset
    zpt_offset = median_df.median(axis=0)

iter_num = 0
for idx, tic in enumerate(ticid):

    if METHOD == 'ZEROPOINT':
        # read in the raw light curve
        lc = pd.read_csv(DETREND_DIRECTORY + str(tic) + "_" + SECTOR + "_" + CAMERA + "_" + CCD + ".lc",
                         names=['JD', 'detrend', 'raw', 'err'], sep=" ", na_values='*********')

        # this is the current version of 'clean' in the pipeline
        lc['clean'] = lc['detrend'] - zpt_offset

        # write the cleaned light curve
        lc[['JD', 'clean', 'detrend', 'raw', 'err']].to_csv(OUTPUT_DIRECTORY +
                                                            str(tic) + "_" + SECTOR + "_" +
                                                            CAMERA + "_" + CCD + ".lc", sep=" ",
                                                            float_format='%.6f', header=False, index=False)

        if (idx % 100 == 0) & (idx > 0):
            print("100 stars have been zeropoint corrected in sector: " + SECTOR + " camera: " +
                  CAMERA + " ccd: " + CCD + ". Now working on the next 100. " +
                  str(len(ticid) - idx - 1) + " stars remain to be zeropoint corrected.")

    if METHOD == 'SPEARMAN':
        # read in the raw light curve
        # lc = pd.read_csv(RAW_DIRECTORY + str(tic) + "_" + SECTOR + "_" + CAMERA + "_" + CCD + ".lc",
        #                 names=['JD', 'raw', 'err'], sep=" ", na_values='*********', usecols=[0, 2, 3])
        lc = pd.read_csv(RAW_DIRECTORY + str(tic) + "_" + SECTOR + "_" + CAMERA + "_" + CCD + ".lc",
                         names=['JD', 'raw', 'err'], sep=" ", na_values='*********')
        
        # get the light curve offset from median
        spear_lc = lc['raw'] - lc['raw'].median()

        # get the spearman rank relationship for each light curve
        corr_lc = median_df.apply(lambda x: x.corr(pd.Series(spear_lc), 'spearman'), axis=1)

        # check for like stars based on the spearman rank
        if len(corr_lc[(corr_lc > SPEARMAN_CUT) & (corr_lc < .99999)]) >= 1:
            # be sure to remove any likely variables based on the standard deviation
            std_df = median_df[(corr_lc > SPEARMAN_CUT) & (corr_lc < .99999)].std(axis=1)

            # pull out stars with high spearman ranks correlation
            trend = np.average(median_df[(corr_lc > SPEARMAN_CUT) &
                                         (corr_lc < .99999) &
                                         (std_df < std_df.min() * 3)],
                               weights=corr_lc[(corr_lc > SPEARMAN_CUT) &
                                               (corr_lc < .99999) &
                                               (std_df < std_df.min() * 3)], axis=0)
            lc['detrend'] = lc['raw'] - trend
        else:
            if len(corr_lc[(corr_lc > SPEARMAN_CUT - 0.05) & (corr_lc < .99999)]) >= 1:
                # be sure to remove any likely variables based on the standard deviation
                std_df = median_df[(corr_lc > SPEARMAN_CUT - 0.05) & (corr_lc < .99999)].std(axis=1)

                # pull out stars with high spearman ranks correlation
                trend = np.average(median_df[(corr_lc > SPEARMAN_CUT - 0.05) &
                                             (corr_lc < .99999) &
                                             (std_df < std_df.min() * 3)],
                                   weights=corr_lc[(corr_lc > SPEARMAN_CUT - 0.05) &
                                                   (corr_lc < .99999) &
                                                   (std_df < std_df.min() * 3)], axis=0)
                lc['detrend'] = lc['raw'] - trend
            else:
                lc['detrend'] = lc['raw']

        # remove NaN values
        lc = lc.fillna(-9.999999)

        # write the detrended light curve
        lc[['JD', 'detrend', 'raw', 'err']].to_csv(DETREND_DIRECTORY +
                                                   str(tic) + "_" + SECTOR + "_" +
                                                   CAMERA + "_" + CCD + ".lc", sep=" ",
                                                   float_format='%.6f', header=False, index=False)

        # we are looking for stars within +/- 200 nearest magnitudes
        if (idx % 100 == 0) & (idx > 300):
            print("Moving the trend star filter by 100 stars.")

            # set the bounds on the light curves to read in
            lw_bnd = TREND_STARS + (iter_num * 100)
            up_bnd = TREND_STARS + (iter_num * 100) + 100

            # get the light curves for the initial number of light curves
            fname_hold = [RAW_DIRECTORY + str(x) + "_" + SECTOR + "_" + CAMERA + "_" + CCD + ".lc"
                          for x in ticid[lw_bnd:up_bnd]]
            # raw_hold = np.array([np.genfromtxt(f, delimiter=' ')[:, 2] for f in fname_hold])
            raw_hold = np.array([np.genfromtxt(f, delimiter=' ')[:, 1] for f in fname_hold])
            median_hold = pd.DataFrame(raw_hold).replace('*********', 99.999999) - \
                          np.median(raw_hold, axis=1).reshape(-1, 1)

            # add the next 100 stars, clip the top 100 and reorganize the data frame
            median_df = median_df.append(median_hold)
            median_df = median_df.drop(np.arange(100)).reset_index(drop=True)

            iter_num += 1

        if (idx % 100 == 0) & (idx > 0):
            print("100 stars have been de-trended in sector: " + SECTOR + " camera: " +
                  CAMERA + " ccd: " + CCD + ". Now working on the next 100. " +
                  str(len(ticid) - idx - 1) + " stars remain to be cleaned.")
