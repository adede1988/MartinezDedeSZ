
%adjust these for paths on your local machine: 
codePre = 'G:\My Drive\GitHub\';
datPre = 'H:\SZ_anton_data\EEG Resting State Data';
addpath('C:\Users\dtf8829\Documents\MATLAB\eeglab2023.0')

%setting code paths: 
addpath([codePre 'MartinezDedeSZ'])
addpath([codePre 'SheffieldAutismBiomarkers'])
cd([codePre 'MartinezDedeSZ'])
standardTrodes

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


files = dir(datPre);
%select the Filt.set files: 
%NOTE: I think we should really start with raw data, needs editing to
%select raw data
files = files(cellfun(@(x) contains(x, '.cnt'), {files.name}));

%loop subject files: 
for ii = 1:length(files)
ii
try
    %% read in the data
    EEG = pop_loadcnt([datPre '\' files(ii).name],...
        'dataformat', 'auto', 'memmapfile', '');
    chanDat = pop_loadset('filename','C2_ChannelLoc_Filt.set',...
        'filepath','G:\My Drive\Milne\SZproject');
    EEG.chanlocs = chanDat.chanlocs; 
%     [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
   

   



    %% make meta data struct
    data = struct; 
    %col1: subject key
    subName = split(files(ii).name, '_'); 
    data.key = subName{1};
    %col2: dir (directory of the data)
    data.dir = files(ii).folder; 
    %col3: fn (file name of data)
    data.fn = files(ii).name;
    %col4: sex
    data.sex = "M"; % HARD CODE NEEDS FIX
    %col5: age
    data.age = 25; % HARD CODE NEEDS FIX
    %col6: task
    data.task = "rest";
    %col7: diag (group)
    if contains(data.key, "C")
        data.group = "CON";
    else 
        data.group = "SZ"; 
    end
    %col8: paranoia measure 1, value? 
    data.para1 = 100; % HARD CODE NEEDS FIX
    %col9: paranoia measure 2, value? 
    data.para2 = 100; % HARD CODE NEEDS FIX
    %col10: para1_measure
    data.para1_measure = "trait variable"; %WHAT'S THE NAME? 
    %col11: para2_measure
    data.para2_measure = "diagnostic measure"; %WHAT'S THE NAME? 
    %col12: ANY OTHER MEASURE WE WANT? 
    data.BLANK = 100; %HARD CODE NEEDS FIX, AND CHANGE NAME
    %col13: EXTRA FIELD TO DESCRIBE OTHER MEASURE
    data.BLANK_version = "VARIABLE TYPE"; %DESCRIPTIVE TAG
    %col14: eyes
    data.eyes = "open";
    %col15: dataSet
    data.dataSet = "AntonDat"; %overall tag for the dataset in general


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
        flipFlag = 0; %flip to closed first if needed when splitting later
    else
        EEG.data = EEG.data(:,EEG.event(closedIDX).latency:end);
        closedStart = 1; 
        openStart = EEG.event(openIDX).latency - ...
                      EEG.event(closedIDX).latency;
        flipFlag = 1; 
    end

    %how long is the total rest recording? 
    recLen = size(EEG.data,2)*(1000/EEG.srate);

    EEG.times = [1:1000/EEG.srate:recLen];
    %check taht time is aligned to the data: 
    if length(EEG.times) ~= size(EEG.data,2)
        disp(['Time and data are different lengths! file: ' num2str(ii)])
    end

    
    %% filter data and epoch the data into 2 second periods
    %how many epochs can we fit into the recording? 
    epochs = floor(size(EEG.data,2)/(2*EEG.srate));
    %trim the data and the time to a whole number epoch length:
    EEG.data = EEG.data(:,1:epochs*(2*EEG.srate)); 
    time = EEG.times(1:epochs*(2*EEG.srate)); 


     %design notch filter
    wo = 50/(EEG.srate/2); 
    bw = wo/10; 
    [b,a] = iirnotch(wo,bw); 

%     wo = 100/(EEG.srate/2); 
%     bw = wo/10; 
%     [b2,a2] = iirnotch(wo,bw); 
% 
%     wo = 150/(EEG.srate/2); 
%     bw = wo/10; 
%     [b3,a3] = iirnotch(wo,bw); 
% 
%     wo = 200/(EEG.srate/2); 
%     bw = wo/10; 
%     [b4,a4] = iirnotch(wo,bw); 
  
    %apply notch, high pass, and low pass
    L = size(EEG.data,2); 
    nextSeg = EEG.data;
    nextSeg = nextSeg - mean(nextSeg,2); 
    padDat = flip(nextSeg,2);
    padDat = [padDat, nextSeg, padDat]; 
    test =  highpass(double(padDat)', .5, EEG.srate); 
    test = lowpass(test, 200, EEG.srate); 
    test = filtfilt(b,a,test); 
%     test = filtfilt(b2,a2,test);
%     test = filtfilt(b3,a3,test);
%     test = filtfilt(b4,a4,test);
    EEG.data(:,:) = test(L+1:L*2, :)';       
   



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
    EEG.openStart = openStart; 
    EEG.closedStart = closedStart; 
    
    if flipFlag == 0
        EEG.epochStamps = time(1,:)<(closedStart*(1000/EEG.srate)); % boolean indicator where 1 = open 
    else
        EEG.epochStamps = time(1,:)>(openStart*(1000/EEG.srate)); % boolean indicator where 1 = open 
    end

    %% filter the data: 
% 
%     %design notch filter
%     wo = 50/(EEG.srate/2); 
%     bw = wo/10; 
%     [b,a] = iirnotch(wo,bw); 
%   
%     %apply notch, high pass, and low pass
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
            = removeNoiseChansVoltSZ(EEG);

    %% reality check alpha plot
%     frex = logspace(log10(2),log10(80),100);
%     transforms = [1,2,3,4]; 
%     numfrex = length(frex); 
%     stds = linspace(2,5,numfrex);
%     powSpect = zeros(size(EEG.data,1), 1); 
%     for chan = 1:size(EEG.data,1)
%        
%     
%     test = getChanTrialTF(squeeze(EEG.data(chan,:,:)),...
%                         frex(43), 1, stds(43),500);
%     
%     powSpect(chan, :) = squeeze(mean(abs(test), [1,2])); 
%     
%     
%     
%     
%     
%     end
% %     
% %     figure
% %     imagesc(powSpect)
% %     yticks([1:size(EEG.data,1)])
% %     yticklabels({EEG.chanlocs.labels})
% %     title(data.key)
% %    
%     figure
%     subplot 211
%     topoplot(powSpect, EEG.chanlocs,'electrodes','on','plotrad',.53); 
%     
%     title(data.key)

    %% interpolation to standard 32-channel montage 
    standardTrodes; 

   
    EEG = convertCoordinatesSZ(EEG, standardEEGlocs, pi/2); 
    
%     powSpect = zeros(size(EEG.dataSTD,1), 1); 
%     for chan = 1:size(EEG.dataSTD,1)
%        
%     
%     test = getChanTrialTF(squeeze(EEG.dataSTD(chan,:,:)),...
%                         frex(43), 1, stds(43),500);
%     
%     powSpect(chan) = squeeze(mean(abs(test), [1,2])); 
%     
%     
%     
%     
%     
%     end
%     plot(powSpect)
%     topoplot(powSpect, EEG.chanlocsSTD,'electrodes','on','plotrad',1); 
    

    %% 

    chanLocs = EEG.chanlocsSTD; 
    save([datPre '/' 'exampleChanLocs.mat'], 'chanLocs'); 
    %% split to open and closed Dat
    EEGclosed = EEG; 
    EEGopen = EEG; 
    
   
    EEGopen.dataSTD = EEG.dataSTD(:,:, EEG.epochStamps==1); 
    EEGclosed.dataSTD = EEG.dataSTD(:,:,EEG.epochStamps==0); 
    EEGopen.dataLap = EEG.dataLap(:,:, EEG.epochStamps==1); 
    EEGclosed.dataLap = EEG.dataLap(:,:,EEG.epochStamps==0);

    %trim 4 seconds on either end
    EEGopen.dataSTD = EEGopen.dataSTD(:,:, 3:end-3); 
    EEGclosed.dataSTD = EEGclosed.dataSTD(:,:,3:end-3); 
    EEGopen.dataLap = EEGopen.dataLap(:,:, 3:end-3); 
    EEGclosed.dataLap = EEGclosed.dataLap(:,:,3:end-3);
   
    %make copies of meta data: 
    dataOpen = data; 
    dataClosed = data; 
    dataOpen.eyes = 'open'; 
    dataClosed.eyes = 'closed'; 
    %% do connectivity analysis since this is easier to do before splitting the participant
    %frequency params
    frex = logspace(log10(2),log10(80),100);
    numfrex = length(frex); 
    stds = linspace(2,5,numfrex);
    %calculating the connectivity, but this is flawed right now because it
    %will need to be split for eyes open versus closed data
    [dataOpen.ispc, EEGopen] = getISPCSZ(EEGopen, frex, numfrex, stds);
    [dataClosed.ispc, EEGclosed] = getISPCSZ(EEGclosed, frex, numfrex, stds);

    %% split files
 
    %NOTE: you need to make a CHANFILES folder inside the data directory
    fileSplitterSZ(EEGopen, join([data.dir '/CHANFILES/OPEN_' data.key], ''))
    fileSplitterSZ(EEGclosed, join([data.dir '/CHANFILES/CLOSED_' data.key], ''))


    %% save out the subject level files for open and closed data: 

    %NOTE: you need to make a SUMDAT folder inside the data directory
    save(join([data.dir '/SUMDAT/OPEN_' data.key '.mat']), 'dataOpen')
    save(join([data.dir '/SUMDAT/CLOSED_' data.key '.mat']), 'dataClosed')
catch
    files(ii).name
end
end
















