#%%
import glob
import mne
import pandas as pd
import numpy as np
from fooof import FOOOF
from scipy.signal import welch

workingpath = 'D:\\aperiod\\'

segmentcode = [56, 63, 78]
segment = ['baseline', 'training', 'stress']
segmentlength = [5, 5, 5]
roi_channels = ['Fp1', 'Fp2', 'Fz', 'F4', 'F3', 'F8', 'F7', 'Cz', 'C4', 'C3', 'T4', 'T3', 'Pz', 'P4', 'P3', 'P8', 'P7', 'O2', 'O1']

subjectlist = []
no_subject = len(glob.glob(workingpath + "eeg\\*.vhdr"))
for i in range(0, no_subject):
    subjectlist.append(glob.glob(workingpath + "eeg\\*.vhdr")[i][-9:-5])
subjectlist.sort()

with open('fooof_results.txt', 'w') as f:
    f.write("Subject\tsegment\tROI\tOffset\tExponent\tR2\tMAE\tResidual_Delta\tResidual_Theta\tResidual_Alpha\tResidual_Beta\tDelta\tTheta\tAlpha\tBeta\t\n")

for subject in subjectlist:

    filepathevent = workingpath + 'eeg\\' + str(subject) + '.vhdr'
    filepath = workingpath + 'eeg_clean\\' + str(subject) + '.edf'
    raw = mne.io.read_raw_brainvision(filepathevent, misc='auto', scale=1.0, preload=True, verbose=None)
    raw_clean = mne.io.read_raw_edf(filepath, preload=True, verbose=None, stim_channel='Status')
    raw_clean._data *= 1e6

    eventcoding = pd.read_csv((workingpath + 'master_eventsII.txt'), delimiter='\t', header=0)
    events = mne.find_events(raw, stim_channel="TRIGGER")
    tempeventdf = pd.DataFrame(events)
    tempeventdf = tempeventdf.loc[:, (tempeventdf != 0).any(axis=0)]
    tempeventdf = tempeventdf.loc[:, (tempeventdf != 47831).any(axis=0)]
    tempeventdf.columns = range(tempeventdf.shape[1])
    tempeventdf = tempeventdf.rename(columns={0: 'trigger_time', 1: 'Code'})
    tempeventdf = tempeventdf[tempeventdf['Code'] >= 50]

    eventdf = pd.merge(eventcoding, tempeventdf, on='Code', how='inner')

    for s in range(0, len(segment)):
        starttime = eventdf[eventdf['Code'] == segmentcode[s]]['trigger_time'].values[0]

        if segment[s] == 'baseline':
            starttime = starttime + (15 * 60 * 500)

        for roi in roi_channels:
            rawdata = raw_clean.get_data(roi)
            filterdata = rawdata[0]
            data = filterdata[starttime:starttime + (segmentlength[s] * 60 * 500)]

            freqs, psd = welch(
                data,
                fs=500,
                window='hann',
                nperseg=1000,
                noverlap=500,
                scaling='density',
                average='mean',
                detrend='constant'
            )

            df = pd.DataFrame({'Frequency_Hz': freqs, 'PSD': psd})
            df.to_csv('D:\\aperiod\\export_eeg_psd\\' + segment[s] + '\\' + subject + '_' + roi + '.csv', index=False)

            fm = FOOOF(
                peak_width_limits=[1, 6],
                max_n_peaks=6,
                min_peak_height=0.05,
                peak_threshold=1.5,
                aperiodic_mode='fixed'
            )
            fm.fit(freqs, psd, [2, 40])

            offset, exponent = fm.aperiodic_params_
            aperiodic_freqs = fm.freqs
            aperiodic_fit = fm._ap_fit

            df_aperiod = pd.DataFrame({'Frequency_Hz': aperiodic_freqs, 'PSD': aperiodic_fit})
            df_aperiod.to_csv('D:\\aperiod\\export_eeg_aperiod\\' + segment[s] + '\\' + subject + '_' + roi + '.csv', index=False)

            mask = freqs > 0
            freqs = freqs[mask]
            psd = psd[mask]
            aperiodic_fit_full = offset - exponent * np.log10(freqs)
            peak_only_psd = np.log10(psd) - aperiodic_fit_full

            with open('fooof_results.txt', 'a') as f:
                f.write(str(subject) + '\t')
                f.write(str(segment[s]) + '\t')
                f.write(str(roi) + '\t')
                f.write(str(offset) + '\t')
                f.write(str(exponent) + '\t')
                f.write(str(fm.r_squared_) + '\t')
                f.write(str(fm.error_) + '\t')

            bands = {
                'delta': (1, 4),
                'theta': (4, 8),
                'alpha': (8, 12),
                'beta': (12, 30)
            }

            for band_name, (fmin, fmax) in bands.items():
                mask = (freqs >= fmin) & (freqs < fmax)
                band_power = np.mean(peak_only_psd[mask])
                with open('fooof_results.txt', 'a') as f:
                    f.write(str(band_power) + '\t')

            for band_name, (fmin, fmax) in bands.items():
                mask = (freqs >= fmin) & (freqs < fmax)
                band_power = np.mean(psd[mask])
                with open('fooof_results.txt', 'a') as f:
                    f.write(str(band_power) + '\t')

            with open('fooof_results.txt', 'a') as f:
                f.write('\n')
