%%script to select noise components from ICA calculated in a_1
%we will remove bad timepoints only after epoching the data, since
%data_scroll with then be much shorter and we only need to inspect the
%relevant timepoints related to the events/triggers

%% set up
clc; clear; %housekeeping
addpath('/data/pt_02584/Tools/eeglab2022.0')                          %add path for eeglab software
addpath('/data/p_02191/Analysis/Nadine/EEG/rest/preprocessing/code/functions/') %add path for functions called in this script

%define paths
basedir =  '/data/p_02191/';    % define base directory 
datadir = [basedir,'Data/'];     % path to EEG data
preprocdir = [basedir,'Analysis/Nadine/EEG/task/preprocessing/data/'];
save_data = '/data/p_02191/Nadine/'; %path to store outcome of this script
QA_save_path= [basedir,'Analysis/Nadine/EEG/task/preprocessing/', 'QA/']; %path to store quality assesmt images

%save error logs here
error_log = struct('ID', {}, 'str', {});
save_error = [basedir,'Analysis/Nadine/EEG/task/preprocessing/','error_log/']; %path to save error_log

%get sublist
sublist = dir(fullfile(datadir, 'S0*'));
%sublist.name(17) 

eeglab
%% start processing

remChan = [];

parfor s=13:length(sublist)
    
    try

        %get subject ID and task path
        subid = sublist(s).name;
                
        %create QA folders
        if not(isfolder([QA_save_path, subid]))
            mkdir([QA_save_path, subid])
        end

        EEG = pop_loadset('filename',[subid, '.set'], 'filepath', [preprocdir, subid]);
        %eeglab redraw

        %run IClable for classification of components
        %looks at components from step 1 and applies label to them (noise vs. eye vs. movement)
        EEGorig = pop_iclabel(EEG, 'default');

        %%Add component idx to original set, saves index of components that
        %%we want to reject later
        EEGorig.etc.ic_remove.muscle = find(EEGorig.etc.ic_classification.ICLabel.classifications(:,2) >= 0.5);
        EEGorig.etc.ic_remove.eye = find(EEGorig.etc.ic_classification.ICLabel.classifications(:,3) >= 0.60);
        EEGorig.etc.ic_remove.LineNoise = find(EEGorig.etc.ic_classification.ICLabel.classifications(:,5) >= 0.5);
        EEGorig.etc.ic_remove.ChannelNoise = find(EEGorig.etc.ic_classification.ICLabel.classifications(:,6) >= 0.4);      
       %EEGorig.etc.ic_remove.other = find(EEGorig.etc.ic_classification.ICLabel.classifications(:,7) >= 0.5)
        EEGorig.etc.ic_remove.heart = find(EEGorig.etc.ic_classification.ICLabel.classifications(:,4) >= 0.8);

        reject = unique([EEGorig.etc.ic_remove.heart', EEGorig.etc.ic_remove.muscle', EEGorig.etc.ic_remove.eye', EEGorig.etc.ic_remove.LineNoise', EEGorig.etc.ic_remove.ChannelNoise']);
        %reject selected components
        EEGclean = pop_subcomp(EEGorig, reject);
        EEGclean = eeg_checkset(EEGclean);

        %re-reference data
        EEGclean= pop_reref(EEGclean, []);

        %save data
        EEGclean = pop_saveset(EEGclean,'filename',[subid, 'post_ICA'],'filepath',[preprocdir, subid]);
   
    catch ME

        fileID = fopen([save_error, subid, '_error_log_ICA_comp_select.txt'],'w');
        fprintf(fileID,'%6s\n',subid);  
        fprintf(fileID,'%6s\n',getReport(ME));
        fclose(fileID);
        
        
    end
        
end
