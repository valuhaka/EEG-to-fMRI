clc; clear;

addpath('/data/tu_steinfath_cloud/PhD/Tools/eeglab2021.1/')
addpath('/data/tu_steinfath_cloud/PhD/Code/HEP/MPI/Supplement/')

pathEEG = '/data/pt_02035/Data/eeg_data/novelty/';
pathICA ='/data/pt_02584/Oddball/Preproc/ICA/';
pathPostICA = '/data/pt_02584/Oddball/Preproc/postICA/';
savePath = '/data/pt_02584/Oddball/Preproc/postICA_1_45Hz/data/';
QA_save_path='/data/pt_02584/Oddball/Preproc/postICA_1_45Hz/QA/';
save_error = '/data/pt_02584/Oddball/Preproc/postICA_1_45Hz/Logfiles/';

eeglab

fs=500; % downsample the data to 250

files = dir(fullfile(pathEEG, '*.vhdr'));
prepFiles = dir(fullfile(savePath, '*.set'));
prepFileList = {prepFiles.name};

rawFiles = {files.name}
ICAFiles = dir(fullfile(pathICA, '*.set'));
ICAFilesList = {ICAFiles.name};

postICAFiles = dir(fullfile(pathPostICA, '*.set'));
postICAFilesList = {postICAFiles.name};

% assigning standard electrode positions
elecfile='/data/tu_steinfath_cloud/PhD/Tools/eeglab2021.1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp';

error_log = struct('ID', {}, 'str', {});

removedChannels = {};

%for debugging
% s = find(count({postICAFiles.name}, 'LI02255015_Novelty'));


% randSub_ids = {ICAFiles(randi(length(ICAFiles),100,1)).name}
% save( '/data/tu_steinfath_cloud/PhD/Code/HEP/MPI/randSub_100.mat', 'randSub_ids')
load('/data/tu_steinfath_cloud/PhD/Code/HEP/MPI/randSub_100.mat')

%% downsample save as set
parfor s=1:length(randSub_ids)
    
    try
        
        
        subjid = randSub_ids{s}(1:end-4);
        
        %create QA folders
        if not(isfolder([QA_save_path, subjid]))
            mkdir([QA_save_path, subjid])
        end
        
        
        %         subjid = postICAFiles(s).name(1:end-4);
        %
        %         %Check if subject was processed already
        newID = [subjid '.set']
        rawID = [subjid '.vhdr']
        %
        %         if any(strcmp(prepFileList, newID))
        %             continue
        %         end
        
        %match post ICA and Raw data
        post_ica_ID =  find(strcmp(postICAFilesList, newID))
        ica_ID =  find(strcmp(ICAFilesList, newID))
        raw_ID = find(strcmp(rawFiles, rawID))
        
        %       If subjids dont mepochatch - give error
        %         if ~(strcmp(postICAFilesList{post_ica_ID}(1:end-4), files(raw_ID).name(1:end-5)) && strcmp(files(raw_ID).name(1:end-5), ICAFilesList{ica_ID}(1:end-4)))
        %             fileID = fopen([save_error, subjid, 'no_match_IDs_error_log.txt'],'w');
        %             fprintf(fileID,'%6s\n',subjid, postICAFilesList{post_ica_ID}(1:end-4), ICAFilesList{ica_ID}(1:end-4));
        %             fprintf(fileID,'%6s\n',getReport(ME));
        %             fclose(fileID);
        %             continue
        %         end
        
        %load ICA data
        EEG_ica = pop_loadset(ICAFilesList{ica_ID}, pathICA);
        
        %load post ICA data
        EEG_post_ica = pop_loadset(postICAFilesList{post_ica_ID}, pathPostICA);
        
        %load Raw
        [EEG, com] = pop_loadbv(pathEEG, files(raw_ID).name);
        
        % resample to fs
        EEG = pop_resample(EEG,fs);
        
        EEG.setname = files(s).name(1:end-5)
        
        %Add montage
        EEG=pop_chanedit(EEG, 'lookup',elecfile);
        
        %%% For testing only
        %EEG = pop_select(EEG, 'point', [1 5000])
        
        %select channels
        EEG.ECG =   pop_select(EEG, 'channel',{'EKG'});
        EEG.VEOG =  pop_select(EEG, 'channel',{'VEOG'});
        EEG.HEOG =  pop_select(EEG, 'channel',{'HEOG'});
        EEG =       pop_select(EEG, 'channel',[1:31]);
        
        
        % Band-pass filter combined with the notch
        %         [b1,a1]=butter(2,[0.3 45]/(EEG.srate/2));
        %         [b2,a2]=butter(2,[49 51]/(EEG.srate/2),'stop');
        %
        %         EEG.data = filtfilt(conv(b1,b2),conv(a1,a2),double(EEG.data)')';
        
        
        % %       low pass than high pass filter the data
        [b,a]=butter(2,1/(EEG.srate/2),'high'); % high pass filter
        EEG.data=filtfilt(b,a,double(EEG.data)')'; % this function doubles the filter order -> what does it mean?
        [c,d]=butter(2,45/(EEG.srate/2), 'low'); % low pass filter
        EEG.data=filtfilt(c,d,double(EEG.data)')';
        
        EEG = pop_eegfiltnew(EEG, 'locutoff',49,'hicutoff',51,'revfilt',1,'plotfreqz',0);
        
        %keep backup for channel interpolation
        originalEEG=EEG;
        
        %number of interpolated electrodes
        %rejnum(s)=originalEEG.nbchan-EEG.nbchan;
        
        
        %run trim_outlier to detect bad epochs
        [EEG, rejData] = PS_trimOutlier_adjust_new(EEG, 80,40, 500);
        EEG.reject.trimOutlier_reject = rejData;
        
        %Add ICA results to EEG struct
        EEG.icawinv = EEG_ica.icawinv;
        EEG.icasphere = EEG_ica.icasphere;
        EEG.icaweights = EEG_ica.icaweights;
        EEG.icachansind = EEG_ica.icachansind;
        eeg_checkset(EEG);
        
        
        reject = unique([EEG_post_ica.etc.ic_remove.heart', EEG_post_ica.etc.ic_remove.muscle', EEG_post_ica.etc.ic_remove.eye', EEG_post_ica.etc.ic_remove.LineNoise', EEG_post_ica.etc.ic_remove.ChannelNoise', EEG_post_ica.etc.ic_remove.other']);
        EEG = pop_subcomp(EEG, reject);
        EEG.etc = EEG_post_ica.etc
        EEG = eeg_checkset(EEG);
        
        
        % remove flat channels & channels with low correlation
        %         EEG = clean_artifacts(EEG,'ChannelCriterion','off', 'FlatlineCriterion', 5, 'BurstCriterion', 'off', 'WindowCriterion','off');
        
        EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',4,'ChannelCriterion',0.8,'LineNoiseCriterion',5,'Highpass','off','BurstCriterion','off','WindowCriterion','off','BurstRejection','off','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
        
        allchan = {originalEEG.chanlocs.labels}
        EEG.reject.removed_channels = allchan(~ismember({originalEEG.chanlocs.labels}, {EEG.chanlocs.labels}))
        
        %         removedChannels{s,1} = subjid;
        %         removedChannels{s,2} = numel(EEG.reject.removed_channels);
        %         removedChannels{s,3} = EEG.reject.removed_channels;
        
        %interpolate
        EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');
        eeg_checkset(EEG);
        
        % Re-reference the data to common average
        EEG = pop_reref(EEG,[]);
        
        
        %Plot spectogram
        figure; plot_spec(EEG.data', 500, 'f_max',70);
        saveas(gcf, [QA_save_path, subjid, '/', subjid, '_PSD_post_05_45Hz_ICA.png'],'png');
        
        disp([subjid, ' done'])
        
        %save results
        pop_saveset(EEG,'filename',subjid,'filepath',savePath);
        
        
        %if Loop completed and previously a error log was created, delete
        %the error log
        if isfile([save_error, subjid, '_error_log.txt'])
            delete([save_error, subjid, '_error_log.txt'])
        end
        
        
    catch ME
        
        fileID = fopen([save_error, subjid, '_error_log.txt'],'w');
        fprintf(fileID,'%6s\n',subjid);
        fprintf(fileID,'%6s\n',getReport(ME));
        fclose(fileID);
        
        
    end
end
