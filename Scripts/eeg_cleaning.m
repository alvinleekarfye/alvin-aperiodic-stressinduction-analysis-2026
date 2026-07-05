rawPath = 'D:\aperiod\eeg\';
savePath = 'D:\aperiod\eeg_clean\';
fileList = dir(fullfile(rawPath, '*.vhdr'));

if isempty(gcp('nocreate'))
    parpool(16);
end

parfor k = 1:length(fileList)
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui');
    
    filePath = fullfile(rawPath, fileList(k).name);
    subject = filePath(16:19);
    
    EEG = pop_loadbv(rawPath, [subject '.vhdr'], [], [1:19]);
    EEG = pop_chanedit(EEG, {'lookup','standard_1005.elc'});
    EEG = pop_eegfiltnew(EEG, 'locutoff', 0.5, 'hicutoff', 45, 'plotfreqz',0);
    EEG = pop_reref(EEG, []);
    EEG = GEDAI(EEG, 'auto', 12,  0.5, 'precomputed', false, false);
    
    pop_writeeeg(EEG, [savePath subject '.edf'], 'TYPE','EDF');
 
    disp(subject);
end