% 1st level analysis for BEDOB/GREADT WM task 
% 
% the script:
% - reads in functional files for WM task
% - creates conditions (ignore vs. updaten + respective control, depending on batch also other regressors + motion etc.)
% - reads in the onsets and durations of these conditions
% - all other setting are default in the batch
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Written by:      Nadine Herzog (naherzog@cbs.mpg.de)
%  Created on:      04.05.2020
%  Edited by:       Hendrik Hartmann (hehartmann@cbs.mpg.de)
%  
%  TO DO:  check if the timing variables are right!?
%          substract startT_i1_T from TestSessionStart? gives values between 5.08 and 1.012+03 (=1.01*1000 = 1012)
%          is this seconds? 
%  SOLVED: because 5.08 seconds is for the first to-be-ignored stimulus (first trial) 
%          so: stimulus 1 (2secs) - delay 2-6 secs - stimulus2 (at ~ sec 5)
%          indeed: startT.delay1_on - startT.delay1_off - first trial = 3 secs... so that would fit
%
%  TO DO:  how to align stimulus onsets with scan?
%  SOLVED: talked to Maria W. and she said because of the trigger the stimuli timings are aligned
%
%  TO DO:  incorporate 'check variable' to check for if motion should be included or not
%  DONE
%
%  TO DO:  check with people (lab and methods club) & read about whether to include: 
%          (1) regressor for error trials? 
%          (2) regressors for ALL events, or only events of interest?
%  SOLVED: in methods club Karsten said you should always use the least complex model possible, 
%          hence I will run the analysis including regressors such as probe display, response etc., motion... 
%          and then compare whether my results are the same when I exclude all these 
%                
%  TO DO: recommendation methods club: incorporate a 'junk' variable and store ALL incorrect trials in it, for the conditions of interest, take only the timings and durations
%         of the correct trials  
%  DONE
%
%  re-run analyses 3 times:
%       (1) using batch with only conditions of interst
%       (2) using batch including motion
%       (3) using batch including all regressors (encoding_T, encoding_fix, probe, response, junk, motion)
%
%  Updates by Hendrik Hartmann:
%         - include 24 motion regressors (calculated by realign_regressors.m)
%         - code m2 fixation cross for correct and incorrect trials
%         - remove response regressor and add parametric modulator for RTs to probe event
%         - add regressor for feedback screen
%         
%         ToDo:
%         - remove batch-for-loop - DONE 09.02.2021
%         - split probe  events for ignore and update condition into probe
%         type - DONE 11.02.2021
%
%     Changes after discussion with Sean (11.02.2021)
%         - too few trials to model probe type events separated per task condition
%             - split probe type collapsed for all task condition (novel, dist, target)
%         - don't separate correct and incorrect trials - DONE 12.02.2021
% 
%      16.02.2021
%     Run several models (on random sample) and compare the output.
%     1) correct and incorrect trials combined
%     2) only correct trials
%     3) separate model for probe
% 
%     18.03.2021
%     - replace negative RTs with mean
%     - exclude trials (separate regressor), where one of the timings
%     (presentation duration of encoding, interference or probe) was off.
%     Or maybe just remove the "excess time" of those trials from the implicit baseline and put them into separate regressor 

%% HOUSEKEEPING
clear all
%clearvars -except typetrial_cond  ... this was because 116 did not have trialtype_cond saved

%initialze directories
%first_level_dir = '/data/tu_hehartmann/GREADT/analyses/workingmemory/fMRI/1st_level'; %for HENDRIK
first_level_dir = '/data/p_nro188/Nadine/Analysis_scripts/MRI/1st_level/';    %for Nadine
%behavdata = '/data/p_01923/included/';                                          %directory where the behavioral data and also onsets and durations of stimuli are stored
behavdata = '/data/p_nro188/Data/';
taskdir = 'WM';
%datadir =  '/data/pt_01923/TmpOut/final_script/Derivatives/';
datadir =  '/data/pt_nro188/Data/Derivatives/';
preprocdata = [datadir,'Preproc/sub-'];                                     %directory where the pre-processed data is stored

%data stored for GREADT
% preprocdata = '/data/gh_gr_adipositas_share/Daten/NRO-189_GREADT/analysis/pt_01923/TmpOut/final_script/Derivatives/Preproc/'

%initialize datasets to convert
sublist = [102, 104, 105, 106, 107, 108, 109, 110, 111, 113, 114, 116, 120, 121, 122, 124, 125, 127, 128, 129, 130, 131, 132, 133, 135, 136, 137, 138, 139, 140, 144, 145, 148, 149, 151, 152, 153, 155, 157, 158, 162, 501, 504, 506, 507, 508, 509, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531 ,532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 546, 548]';
%sublist = [116]'; need to do 116 separately 
%sublist = 
subprefix = ''; %for greadt = S, for bedob it is '';

%required jobfiles/batches
batch_1st_level = fullfile(first_level_dir, 'batch_1st_level.m');  %load predefined batch

%%
for isub = 1:length(sublist)
      %prepare variables that need to be read in to the batch
      clearvars -except first_level_dir behavdata taskdir datadir preprocdata sublist subprefix batch_1st_level isub
  try
    %READ ONSETS & DURATIONS
    %initialize subject-specific folders and read onset file
        onsetdir = sprintf([behavdata,subprefix,'%s/','Behavioral/', taskdir], sprintf('%03d',sublist(isub)));
    %cd into folder to read data file & save timings into variable    
        cd(onsetdir)
        expfile = dir('*exp*WS.mat');
        load(expfile.name, 'startT');
        load(expfile.name, 'typetrial_cond');  
        load(expfile.name, 'correct');
        correct(correct(:,:) == 666) = 0;
        t0 = startT.ScanStart;
        
        % create RT vector
        clear resp_time vect_rt mc_rt
        resp_time = startT.rtimekeeper;
        vect_rt = reshape(resp_time,128,1);
        
        % check for negative RTs and replace with mean RT
        negatives = find(vect_rt < 0);
        positives = find(vect_rt > 0);
        vect_rt(negatives) = mean(vect_rt(positives));
        
        % mean correct RTs for later analysis
        meancorr_rt = vect_rt - mean(vect_rt);
        
        % calculate and save event onsets
%         0 = ignore
%         1 = m1
%         2 = update
%         3 = m2
        
        clear onsets_update durations_update onsets_ignore durations_ignore onsets_m1_stim durations_m1_stim onsets_ctrl_update         clear durations_ctrl_update

                for ii = 1:length(typetrial_cond(:,1))                      % update
                    for kk = 1:length(typetrial_cond(1,:))    
                        if typetrial_cond(ii,kk) == 2                                                         
                            onsets_update(ii,kk) = startT.i2_T_on(ii,kk);              % saves the onsets for each interfence stimulus in a 4x32 matrix
                            durations_update(ii,kk) = startT.i2_T_off(ii,kk) - startT.i2_T_on(ii,kk);
                        else
                            onsets_update(ii,kk) = NaN;
                            durations_update(ii,kk) = NaN;
                        end
                    end
                end
                
                for ii = 1:length(typetrial_cond(:,1))                      % ignore
                    for kk = 1:length(typetrial_cond(1,:)) 
                        if typetrial_cond(ii,kk) == 0                                                 
                            onsets_ignore(ii,kk) = startT.i2_N_on(ii,kk);
                            durations_ignore(ii,kk) = startT.i2_N_off(ii,kk)- startT.i2_N_on(ii,kk);
                        else
                            onsets_ignore(ii,kk) = NaN;
                            durations_ignore(ii,kk)= NaN;
                        end
                    end
                end
                 
                for ii = 1:length(typetrial_cond(:,1))                      % m1 fix cross
                    for kk = 1:length(typetrial_cond(1,:))                
                        if typetrial_cond(ii,kk) == 1            
                            onsets_m1_fix(ii,kk) = startT.i2_fix_on(ii,kk);
                            durations_m1_fix(ii,kk) = startT.i2_fix_off(ii,kk)- startT.i2_fix_on(ii,kk);
                        else
                            onsets_m1_fix(ii,kk) = NaN;
                            durations_m1_fix(ii,kk)= NaN;
                        end
                    end
                end
                                      
                for ii = 1:length(typetrial_cond(:,1))                      % m2 stimulus
                    for kk = 1:length(typetrial_cond(1,:))
                        if typetrial_cond(ii,kk) == 3                                             
                            onsets_m2_stim(ii,kk) = startT.i2_T_on(ii,kk);
                            durations_m2_stim(ii,kk) = startT.i2_T_off(ii,kk)- startT.i2_T_on(ii,kk);
                        else                        
                            onsets_m2_stim(ii,kk)= NaN;
                            durations_m2_stim(ii,kk)= NaN;
                        end
                    end
                end
                
                for ii = 1:length(typetrial_cond(:,1))                      % m2 fix cross
                    for kk = 1:length(typetrial_cond(1,:))
                        if typetrial_cond(ii,kk) == 3                                 
                            onsets_m2_fix(ii,kk) = startT.encoding_fix_on(ii,kk);
                            durations_m2_fix(ii,kk) = startT.encoding_off(ii,kk)- startT.encoding_fix_on(ii,kk);
                        else                        
                            onsets_m2_fix(ii,kk)= NaN;
                            durations_m2_fix(ii,kk)= NaN;
                        end
                    end
                end
                
        onsets_encoding_T = startT.encoding_T_on - t0;                      % 1st T (in ignore, update and m1)
        onsets_encoding_T = onsets_encoding_T(~isnan(onsets_encoding_T));          
        durations_encoding_T = startT.encoding_off - startT.encoding_T_on;   
        durations_encoding_T = durations_encoding_T(~isnan(durations_encoding_T));
                
        onsets_update = onsets_update(~isnan(onsets_update));               % remove NaN
        onsets_update = onsets_update - t0;                                 % substract onsets from start_task to get real timings
        durations_update = durations_update(~isnan(durations_update));      % remove NaN 
                  
        onsets_ignore = onsets_ignore(~isnan(onsets_ignore));               % remove NaN 
        onsets_ignore = onsets_ignore - t0;                                 % substract onsets from start_task to get real timings
        durations_ignore = durations_ignore(~isnan(durations_ignore));      % remove NaN 
         
        onsets_m1_fix = onsets_m1_fix(~isnan(onsets_m1_fix));               % remove NaN 
        onsets_m1_fix = onsets_m1_fix - t0;                                 % substract onsets from start_task to get real timings
        durations_m1_fix = durations_m1_fix(~isnan(durations_m1_fix));      % remove NaN 
        
        onsets_m2_fix = onsets_m2_fix(~isnan(onsets_m2_fix));               % remove NaN 
        onsets_m2_fix = onsets_m2_fix - t0;                                 % substract onsets from start_task to get real timings
        durations_m2_fix = durations_m2_fix(~isnan(durations_m2_fix));      % remove NaN 
        
        onsets_m2_stim = onsets_m2_stim(~isnan(onsets_m2_stim));            % remove NaN 
        onsets_m2_stim = onsets_m2_stim - t0;                               % substract onsets from start_task to get real timings
        durations_m2_stim = durations_m2_stim(~isnan(durations_m2_stim));   % remove NaN 
                         
        % probes
        clear onsets_probe durations_probe
            
        onsets_probe = startT.probe_on - startT.TestSessionStart;        
        onsets_probe = reshape(onsets_probe,1,4*(length(onsets_probe)));
            
        durations_probe = startT.probe_off - startT.probe_on;
        durations_probe = reshape(durations_probe,4*(length(durations_probe)),1);
        
        % feedback screen
        onsets_feedback = startT.show_accuracy_on - t0;
        
        % responses
        onsets_response = startT.response_on - t0;
        onsets_response = reshape(onsets_response,4*(length(onsets_response)),1);
        
        % check for longer stimulus presentations
        long_encoding_T = find(durations_encoding_T > 2.04);
        normal_encoding_T = find(durations_encoding_T < 2.04);
        durations_excess_encoding_T = durations_encoding_T(long_encoding_T) - mean(durations_encoding_T(normal_encoding_T));
        
        long_ignore = find(durations_ignore > 2.04);
        normal_ignore = find(durations_ignore < 2.04);
        durations_excess_ignore = durations_ignore(long_ignore) - mean(durations_ignore(normal_ignore));
        
        long_update = find(durations_update > 2.04);
        normal_update = find(durations_update < 2.04);
        durations_excess_update = durations_update(long_update) - mean(durations_update(normal_update));
        
        long_m1_fix = find(durations_m1_fix > 2.04);
        normal_m1_fix = find(durations_m1_fix < 2.04);
        durations_excess_m1_fix = durations_m1_fix(long_m1_fix) - mean(durations_m1_fix(normal_m1_fix));
        
        long_m2_fix = find(durations_m2_fix > 2.04);
        normal_m2_fix = find(durations_m2_fix < 2.04);
        durations_excess_m2_fix = durations_m2_fix(long_m2_fix) - mean(durations_m2_fix(normal_m2_fix));
        
        long_m2_stim = find(durations_m2_stim > 2.04);
        normal_m2_stim = find(durations_m2_stim < 2.04);
        durations_excess_m2_stim = durations_m2_stim(long_m2_stim) - mean(durations_m2_stim(normal_m2_stim));
        
        % adjust wanted stimulus durations
        durations_encoding_T(long_encoding_T) = mean(durations_encoding_T(normal_encoding_T));
        durations_ignore(long_ignore) = mean(durations_ignore(normal_ignore));
        durations_update(long_update) = mean(durations_update(normal_update));
        durations_m1_fix(long_m1_fix) = mean(durations_m1_fix(normal_m1_fix));
        durations_m2_fix(long_m2_fix) = mean(durations_m2_fix(normal_m2_fix));
        durations_m2_stim(long_m2_stim) = mean(durations_m2_stim(normal_m2_stim));
          
        %RUN THE ACTUAL BATCH AND READ IN ALL VARIABLES 
        run(batch_1st_level)
        
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; %defines the units of timing in batch as seconds
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;         %defines the repetition time in batch as 2 secs
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 60;    %number of slices
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 30;   %because multiband - middle slice
        
        %define output dir
        outputdir = [datadir,'1st_level/WM/',sprintf('%03d',sublist(isub))];       
        mkdir(outputdir)
        matlabbatch{1}.spm.stats.fmri_spec.dir = {outputdir};    %reads in the above defined outputdir in the batch
        
        %read in pre-processed functional scans
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList', fullfile([preprocdata,sprintf('%03d',sublist(isub)),'/ses-1/workingmemory']),['swuasub-',sprintf('%03d',sublist(isub)),'_ses-1_task_workingmemory_bold.nii*'],1:1000));
        
        %IGNORE 
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'ignore';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = [onsets_ignore];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = [durations_ignore];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod.name = 'accuracy';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod.param = '<UNDEFINED>';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod.poly = 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        %UPDATE
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'update';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = [onsets_update];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = [durations_update];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        %m1_fix
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'm1_fix';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset = [onsets_m1_fix];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = [durations_m1_fix];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).orth = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        %m2_stim
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'm2_stim';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset = [onsets_m2_stim];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = [durations_m2_stim];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
         
        %m2_fix
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).name = 'm2_fix';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).onset = [onsets_m2_fix];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).duration = [durations_m2_fix];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).orth = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).pmod = struct('name', {}, 'param', {}, 'poly', {});

        %probe
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).name = 'probe';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).onset = [onsets_probe];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).orth = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).pmod = struct('name', {'reaction time'}, 'param', {meancorr_rt}, 'poly', {1});
        
        %encoding_T
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).name = 'm1_1stT';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).onset = [onsets_encoding_T];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).duration = [durations_encoding_T];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).orth = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).pmod = struct('name', {}, 'param', {}, 'poly', {});
            
        % feedback screen
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).name = 'feedback';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).onset = [onsets_feedback];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).duration = 5;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).orth = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        % responses
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).name = 'response';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).onset = [onsets_response];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).duration = 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).orth = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        %read in motion parameters
        motion_pars = [preprocdata,sprintf('%03d',sublist(isub)), '/ses-1/workingmemory/rp_asub-',sprintf('%03d',sublist(isub)),'_ses-1_task_workingmemory_bold.txt'];
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {motion_pars}; 
            
        
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

        spm_jobman('interactive', matlabbatch); %if you put to 'interactive', batch manager will open, if you set it to 'run', batch will run, w0 opening
    catch ME
        errors(1,isub) = strcat("Error in subject ",num2str(sublist(isub)))
        errors(2,isub) = ME.message
   end
 end %end subject loop


%%
%  ---- DEFINE CONTRASTS ----
  
clear all
sublist = [102, 104, 105, 106, 107, 108, 109, 110, 111, 113, 114, 116, 120, 121, 122, 124, 125, 127, 128, 129, 130, 131, 132, 133, 135, 136, 137, 138, 139, 140, 144, 145, 148, 149, 151, 152, 153, 155, 157, 158, 162, 501, 504, 506, 507, 508, 509, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531 ,532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 546, 548]';

%specify directories
datadir =  '/data/pt_nro188/Data/Derivatives/1st_level/WM/';
first_level_dir = '/data/p_nro188/Nadine/Analysis_scripts/MRI/1st_level/';    
%required jobfiles/batches
jobfile = fullfile(first_level_dir, 'contrasts_job.m');

% coding indices: [ignore   update    m1_fix    m2_stim    m2_fix    probe    m1_1stT    feedback   response    upt > m2_stim    ign > m1_fix   EOI     upt > ign   ign > upt]
% coding indices: [   1        2        3         4          5         6        7           8           9            10               11        12          13          14]

for isub = 1:length(sublist)
    
    try
    %initialize subject-specific folders and read onset file
        spmdir = sprintf([datadir,sprintf('%03d',sublist(isub)),'/SPM.mat']);
        matlabbatch{1}.spm.stats.con.spmmat = {spmdir};
        
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'ignore';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1];
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'update';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [0 1];
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'm1_fix';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0 0 1];
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'm2_stim';
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 1];
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'm2_fix';
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 0 1];
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'probe';
        matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 0 0 1];
        matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'm1_1stT';
        matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 0 0 0 1];
        matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'feedback';
        matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = [0 0 0 0 0 0 0 1];
        matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = 'response';
        matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = [0 0 0 0 0 0 0 0 1];
        matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = 'upt > m2_stim';
        matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = [0 1 0 -1];
        matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

        matlabbatch{1}.spm.stats.con.consess{11}.tcon.name = 'ign > m1_fix';
        matlabbatch{1}.spm.stats.con.consess{11}.tcon.weights = [1 0 -1];
        matlabbatch{1}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

        matlabbatch{1}.spm.stats.con.consess{12}.fcon.name = 'EOI';
        matlabbatch{1}.spm.stats.con.consess{12}.fcon.weights = eye(10);
        matlabbatch{1}.spm.stats.con.consess{12}.fcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{13}.tcon.name = 'upt > ign';
        matlabbatch{1}.spm.stats.con.consess{13}.tcon.weights = [-1 1];
        matlabbatch{1}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

        matlabbatch{1}.spm.stats.con.consess{14}.tcon.name = 'ign > upt';
        matlabbatch{1}.spm.stats.con.consess{14}.tcon.weights = [1 -1];
        matlabbatch{1}.spm.stats.con.consess{14}.tcon.sessrep = 'none';
        
        
        matlabbatch{1}.spm.stats.con.delete = 1;
        spm_jobman('run', matlabbatch);
    %end
    catch ME
        errors(1,isub) = strcat("Error in subject ",num2str(sublist(isub)))
        errors(2,isub) = ME.message
    end
end
