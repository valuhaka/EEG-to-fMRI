%%Preprocess Wormcry task data
%script to filter EEG data at 1 hz and calculate ICA

%% set up
clc; clear; %housekeeping
addpath('/data/pt_02584/Tools/eeglab2022.0')                          %add path for eeglab software
addpath('/data/p_02191/Analysis/Nadine/EEG/task/preprocessing/functions/') %add path for functions called in this script

%define paths
basedir =  '/data/p_02191/';    % define base directory 
datadir = [basedir,'Data/'];     % path to EEG data
preprocdir = [basedir,'Analysis/Nadine/EEG/task/preprocessing/'];
save_data = [preprocdir, 'data/']; %path to store outcome of this script

%save error logs here
error_log = struct('ID', {}, 'str', {});
save_error = [preprocdir,'error_log/']; %path to save error_log

%get sublist
sublist = dir(fullfile(datadir, 'S0*'));
%sublist.name(17) 

%% start preprocessing loop

%start eeglab
eeglab

% set downsample rate for data to 500
fs=500; 

% assigning standard electrode positions
elecfile='/data/u_naherzog_software/eeglab/eeglab2022.0/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp'; 
  
parfor s=33:length(sublist)

    try 
        %get subject ID and task data path
        subjid = sublist(s).name;
        sub_data_path = [datadir, subjid, '/EEG/colourWM/']        
        
        %create saveDir folders
        if not(isfolder([save_data, subjid]))
            mkdir([save_data, subjid])
        end
        
        %load EEG data
        subData =  dir(fullfile(sub_data_path, '*.vhdr'));
        [EEG, com] = pop_loadbv(sub_data_path, subData.name);

        %select data, remove irrelevant data befor and after task
        count = length(EEG.event);
        endlatency = EEG.event(count).latency/2500;  %latency of last event in seconds (therefore devided by 2500 (sampling rate)
        EEG = pop_rmdat( EEG, {'S  1'},[0 endlatency] ,0);
            
        % resample to fs
        EEG = pop_resample(EEG,fs); 

        %Add montage 
        EEG=pop_chanedit(EEG, 'lookup',elecfile);
                 
        %select channels
        EEG.ECG =   pop_select(EEG, 'channel',{'ECG'});
        EEG.VEOG =  pop_select(EEG, 'channel',{'VEOG'}); 
        EEG =       pop_select(EEG, 'nochannel', {'ECG', 'VEOG'});

        % Band-pass filter combined with the notch
        [b1,a1]=butter(2,[0.5 45]/(EEG.srate/2));
        [b2,a2]=butter(2,[49 51]/(EEG.srate/2),'stop');

        EEG.data = filtfilt(conv(b1,b2),conv(a1,a2),double(EEG.data)')';

        %keep backup for channel interpolation
        originalEEG=EEG;

        % remove flat channels & channels with low correlation 
        EEG = clean_artifacts(EEG,'ChannelCriterion','off', 'FlatlineCriterion', 5, 'BurstCriterion', 'off', 'WindowCriterion','off');

        %number of interpolated electrodes
        %rejnum(s)=originalEEG.nbchan-EEG.nbchan; 

        allchan = {originalEEG.chanlocs.labels}
        EEG.reject.removed_channels = allchan(~ismember({originalEEG.chanlocs.labels}, {EEG.chanlocs.labels}))

        %interpolate
        EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');
        eeg_checkset(EEG);        
 
        %run trim_outlier to detect bad epochs
        [EEG, rejData] = PS_trimOutlier_adjust_new(EEG, 150,40, 500);
        EEG.reject.trimOutlier_reject = rejData;
        
        if isempty(rejData)
            
            EEG_rej = EEG;
            eeg_checkset(EEG_rej);
            
            %rank
            dataRank = rank(double(EEG_rej.data'));
            
            %run ICA
            EEG_ica = pop_runica(EEG_rej, 'icatype', 'runica', 'pca', dataRank , 'options', {'extended' 1});
        
        elseif rejData(2) > 10*EEG.srate*60
            err_msg = 'More than 10 minutes of Data were rejected by Trim Outlier. ICA calculated on RAW data. Check EEG Timecourse and ICA components carefully'
            fileID = fopen([save_error, subjid, '_error_log.txt'],'w');
            fprintf(fileID,'%6s\n',subjid);
            fprintf(fileID,'%6s\n',err_msg);
            fclose(fileID);
            
            %rank
            dataRank = rank(double(EEG.data'));
            
            %run ICA
            EEG_ica = pop_runica(EEG, 'icatype', 'runica', 'pca', dataRank , 'options', {'extended' 1});
       
        else
            
            EEG_rej = pop_select(EEG, 'nopoint', [rejData(:,1) rejData(:,2)]); %reject data for ICA. TO DO: for manual artifact detection
            eeg_checkset(EEG_rej)
            
            %rank
            dataRank = rank(double(EEG_rej.data'));
            
            %run ICA
            EEG_ica = pop_runica(EEG_rej, 'icatype', 'runica', 'pca', dataRank , 'options', {'extended' 1});
        end
    
        %Add ICA results to EEG struct
        EEG.icawinv = EEG_ica.icawinv;
        EEG.icasphere = EEG_ica.icasphere;
        EEG.icaweights = EEG_ica.icaweights;
        EEG.icachansind = EEG_ica.icachansind;

        disp([subjid, ' ICA done'])
        
        %save results
        pop_saveset(EEG,'filename',subjid,'filepath',[save_data, subjid]);

    catch ME

        fileID = fopen([save_error, subjid, '_error_log.txt'],'w');
        fprintf(fileID,'%6s\n',subjid);  
        fprintf(fileID,'%6s\n',getReport(ME));
        fclose(fileID);
        
        
    end
end
