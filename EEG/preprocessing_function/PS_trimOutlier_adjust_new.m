%% Modified from trimOutlier() %%

function [EEG rejectDataIntervals] = trimOutlier_adjust(EEG, amplitudeThreshold, amplitudeThreshold_hf, pointSpreadWidth)

if ~(nargin==4)
    error('trimOutlier() requires 4 input arguments.')
end

% return if 3-D
if length(size(EEG.data))==3
    disp('Epoched data detected: datapoint rejection will not be performed.')
    return
end
EEG_mfront = pop_select(EEG, 'nochannel', {'Fp1', 'Fp2'}); %reject frontal channels not to detect eye blinks

%% set the threshold for lower freq bands %%
[a,b] = butter(2,[1 15]/(EEG.srate/2));
EEG_lowfreq.data = filtfilt(a,b,double(EEG_mfront.data)')';
meanAllChan = mean(EEG_lowfreq.data(:,:));
stdAllChan  = std( EEG_lowfreq.data(:,:),0,1);
posi2SDChan = meanAllChan + 3*stdAllChan;
nega2SDChan = meanAllChan - 3*stdAllChan;
newThre = max(max(abs(nega2SDChan), posi2SDChan));

if newThre > 150
    amplitudeThreshold = 150;
else
    amplitudeThreshold = newThre;
end

%% remove bad datapoints

% obtain the window size
windowSize = pointSpreadWidth; % millisecond
windowSizeInFrame = round(windowSize/(1000/EEG.srate)); % frame

% compute bad datapoints
absMinMaxAllChan = max([abs(min(EEG_mfront.data)); abs(max(EEG_mfront.data))],[],1);
badPoints  = absMinMaxAllChan > amplitudeThreshold;

% Adjust parameters for higher frequency activity %
% high-pass filter
[d,c] = butter(2,[15 45]/(EEG.srate/2));
EEG_hf.data = filtfilt(d,c,double(EEG_mfront.data)')';

% compute bad datapoints

absMinMaxAllChan_hf = max([abs(min(EEG_hf.data)); abs(max(EEG_hf.data))],[],1);
badPoints_hf  = absMinMaxAllChan_hf > amplitudeThreshold_hf;
    

if any(badPoints)
    % expand badPoints
    badPointsExpanded = logical(conv(single(badPoints), ones(1,windowSizeInFrame), 'same'));
    
    % start Christian's impressive code
    rejectDataIntervals = reshape(find(diff([false badPointsExpanded false])),2,[])';
    rejectDataIntervals(:,2) = rejectDataIntervals(:,2)-1;
    
    % expand badPoints for high frequencies
    badPointsExpanded_hf = logical(conv(single(badPoints_hf), ones(1,windowSizeInFrame), 'same'));
    
    % start Christian's impressive code
    rejectDataIntervals_hf = reshape(find(diff([false badPointsExpanded_hf false])),2,[])';
    rejectDataIntervals_hf(:,2) = rejectDataIntervals_hf(:,2)-1;
    
    for n = 1:size(rejectDataIntervals_hf,1)
        l = size(rejectDataIntervals,1) +1;
        rejectDataIntervals(l,:)=rejectDataIntervals_hf(n,:);
    end
    
     [~,idx] = sort(rejectDataIntervals(:,1)); % sort just the first column
     rejectDataIntervals = rejectDataIntervals(idx,:);
    
      % merge overlapping segments and segments that are less than 3sec appart
        isb = 2;
    while size(rejectDataIntervals,1) > 1 && isb < size(rejectDataIntervals,1)
        prev_minus_next = rejectDataIntervals(isb+1,1) - rejectDataIntervals(isb,2);
        if (prev_minus_next) < (3 * EEG.srate)
            rejectDataIntervals(isb,2) = rejectDataIntervals(isb+1,2);
            rejectDataIntervals(isb+1,:) = [];
            isb = isb - 1;
        end
        isb = isb + 1 ;
    end
     rejectDataIntervals = rejectDataIntervals; %to account for datapoints taken away at the beginning.

     
%      if  rejectDataIntervals(1,1)~=1
%          sn = size(rejectDataIntervals,1)+1;
%          rejectDataIntervals(sn,1) = 1;
%          rejectDataIntervals(sn,2) = 2*EEG.srate;
%          [~,idx] = sort(rejectDataIntervals(:,1)); % sort just the first column
%          rejectDataIntervals = rejectDataIntervals(idx,:);
%      else
%          return
%          
%      end
          
    %    add one second before and after the artifact
    for kk = 1:size(rejectDataIntervals,1) 
        
        if rejectDataIntervals(kk,1) - EEG.srate > EEG.srate*2 % if minusing a second will not verlap with the first 2 seconds
            rejectDataIntervals(kk,1) = rejectDataIntervals(kk,1) - EEG.srate; %start minus a second
        end
        
        if rejectDataIntervals(kk,2) <= length(EEG.data)- EEG.srate
            rejectDataIntervals(kk,2) = rejectDataIntervals(kk,2) + EEG.srate; %end plus 1 second
        else
            rejectDataIntervals(kk,2) = length(EEG.data);
        end
    end
    
    % merge overlapping segments and segments that are less than 3sec appart
        isb = 1;
    while size(rejectDataIntervals,1) > 1 && isb < size(rejectDataIntervals,1)
        prev_minus_next = rejectDataIntervals(isb+1,1) - rejectDataIntervals(isb,2);
        if (prev_minus_next) < (3 * EEG.srate)
            rejectDataIntervals(isb,2) = rejectDataIntervals(isb+1,2);
            rejectDataIntervals(isb+1,:) = [];
            isb = isb - 1;
        end
        isb = isb + 1 ;
    end

      
    % mark them
    for ii = 1:size(rejectDataIntervals,1)
        n = length(EEG.event);
    if ii > 1
            EEG.event(n+1).latency = rejectDataIntervals(ii,1);
            EEG.event(n+1).duration = rejectDataIntervals(ii,2)-rejectDataIntervals(ii,1); % + half a second to compensate for the start extension and + 1 second to compensate for the end extension
            EEG.event(n+1).type = 'auto_start';
        
            EEG.event(n+2).latency = rejectDataIntervals(ii,2);
            EEG.event(n+2).duration = 0;
            EEG.event(n+2).type = 'auto_end';
    else
       EEG.event(n+1).latency = rejectDataIntervals(ii,1);
       EEG.event(n+1).duration = rejectDataIntervals(ii,2);
       EEG.event(n+1).type = 'auto_start';
       
       EEG.event(n+2).latency = rejectDataIntervals(ii,2);
       EEG.event(n+2).duration = 0;
       EEG.event(n+2).type = 'auto_end';
    end
            
    end
%     EEG = eeg_checkset(EEG,'eventconsistency');   

    % display log
    badPointsInSec = length(find(badPointsExpanded))*1000/EEG.srate/1000; %#ok<*NASGU>
    disp(sprintf('\n%2.0fuV threshold with %2.0fms spreading rejected %2.1fsec data, added %1.0f boundaries.', amplitudeThreshold, windowSize, badPointsInSec, size(rejectDataIntervals,1)));
else
    disp('No datapoint rejected.');
%     rejectDataIntervals(1,1) = 1;
%     rejectDataIntervals(1,2) = 2*EEG.srate;
%     n = length(EEG.event);
%      EEG.event(n+1).latency = rejectDataIntervals(1,1);
%      EEG.event(n+1).duration = rejectDataIntervals(1,2);
%      EEG.event(n+1).type = 'auto_start';
%        
%      EEG.event(n+2).latency = rejectDataIntervals(1,2);
%      EEG.event(n+2).duration = 0;
%      EEG.event(n+2).type = 'auto_end';
       rejectDataIntervals = []
end
     
     
     %(EEG_31ch.data, 'srate' , EEG_31ch.srate, 'winlength', 60, 'command');
%      V_Command = 'V_Rejected_Sample_Range = TMPREJ(:,1:2);';   % REJECT button action
%      eegplot(EEG_31ch.data, 'srate' , EEG_31ch.srate, 'winlength', 60, 'command',V_Command,'eloc_file',EEG_31ch.chanlocs);
%      waitfor(gcf);
    
         
% EEG_test = pop_eegthresh(EEG_31ch,1,[1:31] ,-60,60,-0.2,0.495,2,1, 0);
%   pop_eegplot_adjust(EEG, 1,1,1)
