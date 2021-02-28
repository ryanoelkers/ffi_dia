""" This script will work on the TESS FFIs from year 1, which were created using IDL. This script will need to be
updated slightly for future versions. Please confirm directory location of files etc. """

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

# set this flag to 1 if you want to experiment with the new SPEARMAN RANK method
EXPERIMENT = 1

# change this for the sector camera ccd you are looking for
SECTOR = 'sector01'
CAMERA = '1'
CCD = '1'

# these can be updated for whatever you want
TREND_STARS = 500
MG_CLIP = 9

# change these for the specific directories
SECTOR_DIRECTORY = '/net/jovan/export/jovan/oelkerrj/FFI/'
MASTER_DIRECTORY = SECTOR_DIRECTORY + SECTOR + '/master/star_list/'
RAW_DIRECTORY = SECTOR_DIRECTORY + SECTOR + '/lc/raw/'
OUTPUT_DIRECTORY = SECTOR_DIRECTORY + SECTOR + '/output/'  # UPDATE THIS DIRECTORY TO BE WHERE YOU WANT TO WRITE TO

# read in the master magnitude file from the master frame directory
master = pd.read_csv(MASTER_DIRECTORY + SECTOR + '_' + CAMERA + '_' + CCD + '_master.ap',
                     delimiter=' ', names=['TICID', 'x', 'y', 'tessmag', 'mag', 'er', 'flux', 'flux_er'])

# push TICIDs to a list
ticid = master['TICID'].to_list()
master_mags = master['mag'].to_numpy()
zpt_ticid = master[master['tessmag'] < MG_CLIP].TICID.to_list()
zpt_mags = master[master['tessmag'] < MG_CLIP]['mag'].to_numpy()

filenames = [RAW_DIRECTORY + str(x) + "_" + SECTOR + "_" + CAMERA + "_" + CCD + ".lc" for x in zpt_ticid]

# read in the raw lightcurves to find the zeropoint
lc_data = np.array([np.genfromtxt(f, delimiter=' ')[:, 1] for f in filenames])

# subtract the master frame magnitude, and then median combine to determine the zeropoint trend
median_df = pd.DataFrame(lc_data).replace('*********', 99.999999) - zpt_mags.reshape(-1, 1)
zpt_offset = median_df.median(axis=0)

if EXPERIMENT == 1:
    # make the experimental dataframe
    exp_array = median_df[0:TREND_STARS]
    exp_df = pd.DataFrame(exp_array)

iter_num = 0
for idx, tic in enumerate(ticid):
    # read in the raw light curve
    lc = pd.read_csv(RAW_DIRECTORY + str(tic) + "_" + SECTOR + "_" + CAMERA + "_" + CCD + ".lc",
                     names=['JD', 'mag', 'err'], sep=" ", na_values='*********')

    # this is the current version of 'clean' in the pipeline
    lc['clean'] = lc['mag'] - zpt_offset

    if (idx % 100 == 0) & (idx > 0):
        print("100 stars have been zeropoint corrected in sector: " + SECTOR + " camera: " +
              CAMERA + " ccd: " + CCD + ". Now working on the next 100. " +
              str(len(ticid) - idx - 1) + " stars remain to be zeropoint corrected.")

    # if you'd like to experiment with the SPEARMAN RANK gradient method it will be done here
    if EXPERIMENT == 0:
        # write the cleaned light curve
        lc[['JD', 'clean', 'mag', 'err']].to_csv(OUTPUT_DIRECTORY +
                                                 str(tic) + "_" + SECTOR + "_" +
                                                 CAMERA + "_" + CCD + ".lc", sep=" ",
                                                 float_format='%.6f', header=False, index=False)
    else:

        # get the light curve offset from median
        exp_lc = lc['clean'] - lc['clean'].median()

        # get the spearman rank relationship for each light curve
        corr_lc = exp_df.apply(lambda x: x.corr(pd.Series(exp_lc), 'spearman'), axis=1)

        # check for like stars
        if len(corr_lc[(corr_lc > 0.85) & (corr_lc < .99999)]) > 1:
            # pull out stars with high spearman ranks correlation
            trend = np.average(exp_df[(corr_lc > 0.85) & (corr_lc < .99999)],
                               weights=corr_lc[(corr_lc > 0.85) & (corr_lc < .99999)], axis=0)
            lc['detrend'] = lc['clean'] - trend
        else:
            lc['detrend'] = -99.999999

            # remove NaN values
            lc = lc.fillna(-9.999999)

            # write the detrended light curve
            lc[['JD', 'detrend', 'clean', 'mag', 'err']].to_csv(OUTPUT_DIRECTORY +
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
            raw_hold = np.array([np.genfromtxt(f, delimiter=' ')[:, 1] for f in filenames])
            median_hold = pd.DataFrame(raw_hold).replace('*********', 99.999999) - \
                          np.median(raw_hold, axis=1).reshape(-1, 1)

            # add the next 100 stars, clip the top 100 and reorganize the data frame
            exp_df = exp_df.append(median_hold)
            exp_df = exp_df.drop(np.arange(100)).reset_index(drop=True)

            iter_num += 1

        if (idx % 100 == 0) & (idx > 0):
            print("100 stars have been de-trended in sector: " + SECTOR + " camera: " +
                  CAMERA + " ccd: " + CCD + ". Now working on the next 100. " +
                  str(len(ticid) - idx - 1) + " stars remain to be cleaned.")
