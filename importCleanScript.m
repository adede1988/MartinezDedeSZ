
%adjust these for paths on your local machine: 
codePre = 'G:\My Drive\GitHub\';
datPre = 'G:\My Drive\Milne\SZproject';
addpath('C:\Users\dtf8829\Documents\MATLAB\eeglab2023.0')

%setting code paths: 
addpath([codePre 'SheffieldAutismBiomarkers'])
addpath([codePre 'MartinezDedeSZ'])
cd([codePre 'SheffieldAutismBiomarkers'])
standardTrodes

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


files = dir(datPre);
%select the Filt.set files: 
files = files(cellfun(@(x) contains(x, 'Filt.set'), {files.name}));

%loop subject files: 
for ii = 1:length(files)

    %% read in the data
    EEG = pop_loadset('filename',files(ii).name,'filepath',files(ii).folder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
   

    %% cut the EEG data to the rest periods: 

    %get latency of 1001 and 1002 event codes (eyes open and eyes closed)
    openIDX = find([EEG.event.type] == 1001);
    closedIDX = find([EEG.event.type] == 1002); 
    %deal with multiple stamps: 
    if length(openIDX) > 1
        openIDX = openIDX(length(openIDX)); %take the last one
    end
    if length(closedIDX) > 1
        closedIDX = closedIDX(length(closedIDX)); %take the last one
    end


    %check if the open or closed index is first!
    if openIDX < closedIDX
        EEG.data = EEG.data(:, EEG.event(openIDX).latency:end);
        openStart = 1; 
        closedStart = EEG.event(closedIDX).latency - ...
                      EEG.event(openIDX).latency;
    else
        EEG.data = EEG.data(:,EEG.event(closedIDX).latency:end);
        closedStart = 1; 
        openStart = EEG.event(openIDX).latency - ...
                      EEG.event(closedIDX).latency;
    end

    %how long is the total rest recording? 
    recLen = size(EEG.data,2)*(1000/EEG.srate);

    EEG.times = [1:1000/EEG.srate:recLen];
    %check taht time is aligned to the data: 
    if length(EEG.times) ~= size(EEG.data,2)
        disp(['Time and data are different lengths! file: ' num2str(ii)])
    end

    
    %% epoch the data into 2 second periods
    %how many epochs can we fit into the recording? 
    epochs = floor(size(EEG.data,2)/(2*EEG.srate));
    %trim the data and the time to a whole number epoch length:
    EEG.data = EEG.data(:,1:epochs*(2*EEG.srate)); 
    time = EEG.times(1:epochs*(2*EEG.srate)); 
    %reshape the data and the time stamps into epoched matrices
    EEG.data = reshape(EEG.data(:, 1:(2*EEG.srate)*epochs),...
        [size(EEG.data,1),(2*EEG.srate),epochs ]);
    time = reshape(time(1:(2*EEG.srate)*epochs), [(2*EEG.srate),epochs ]);
    %record that there are no problems so far: 
    EEG.probFlag = 0; 
    %hard code the time stamps for a single epoch
    EEG.times = [1:1000/EEG.srate:2000];
    %record the times across individual epochs and the open / closed start
    %times so that we can easily separate open and closed data later: 
    EEG.epochStamps = time(1,:)<closedStart; % boolean indicator where 1 = open 
    EEG.openStart = openStart; 
    EEG.closedStart = closedStart; 

    %% filter the data: 

    %design notch filter
%     wo = 50/(EEG.srate/2); 
%     bw = wo/35; 
%     [b,a] = iirnotch(wo,bw); 
%   
%     for seg = 1:epochs
%         nextSeg = EEG.data(:,:,seg);
%         nextSeg = nextSeg - mean(nextSeg,2); 
%         padDat = flip(nextSeg,2);
%         padDat = [padDat, nextSeg, padDat]; 
%         test =  highpass(double(padDat)', .5, EEG.srate); 
%         test = lowpass(test, 200, EEG.srate); 
%         test = filtfilt(b,a,test); 
%         EEG.data(:,:, seg) = test(EEG.srate*2+1:EEG.srate*4,:)';       
%     end



    %% clean the data: 
    [data.nbChanOrig, data.nbChanFinal, data.nbTrialOrig, data.nbTrialFinal,EEG]...
            = removeNoiseChansVoltSZ(EEG, "no auto save");

    %% interpolation to standard 32-channel montage 
%     standardTrodes; 
%     EEG = convertCoordinatesSZ(EEG, standardEEGlocs); 

    %% do connectivity analysis since this is easier to do before splitting the participant
    %frequency params
%     frex = logspace(log10(2),log10(80),100);
%     numfrex = length(frex); 
%     stds = linspace(2,5,numfrex);
%     %calculating the connectivity, but this is flawed right now because it
%     %will need to be split for eyes open versus closed data
%     [data.ispc, EEG] = getISPCSZ(EEG, frex, numfrex, stds);

    %% things to do next: 
    %the data need to be split into individual channel files for
    %compatibility with the code used for the discovery dataset
    %Also, we need to create data structs that hold the relevant metadata
    %for each subject. See getBioConsortDat.m for example of metadata
    %creation

    %read through fileSplitter.m to see how file splitting is done. We will
    %probably want to alter the file splitting so that it separates and
    %saves out two different files for each channel: eyes open and eyes
    %closed. 

    %once the files are split, we'll need to 

end
















