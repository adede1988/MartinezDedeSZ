%data combine


codePre = 'G:/My Drive/GitHub/';
datPre = 'H:\SZ_anton_data\EEG Resting State Data\';

addpath([codePre 'MartinezDedeSZ'])
summaryDatSave = [datPre 'SUMDAT/'];
chanDatSave = [datPre 'CHANFILES/'];
sumDatSave = [datPre 'SUMDAT/'];


filenames = dir(summaryDatSave);
filenames = filenames(contains({filenames.name}, '.mat'), :); 

%frequency params
frex = logspace(log10(2),log10(80),100);
transforms = [1,2,3,4]; 
numfrex = length(frex); 
stds = linspace(2,5,numfrex);

demographics = readtable([datPre 'demographics.csv']);

for ii = 1:length(filenames)
    tic
    filenames(ii).errorFlag = 0; 
    data = load([filenames(ii).folder '/' filenames(ii).name]);
   
    if isfield(data, 'data')
        data = data.data; 
    elseif isfield(data, 'dataClosed')
        data = data.dataClosed; 
    else
        data = data.dataOpen; 
    end
     %split the key 
    key = strsplit(data.key, '.cnt'); 
    data.key = key{1}; 
     if strcmp(data.eyes, 'open')
        out = load([chanDatSave 'OPEN_' data.key '_processed_temp' num2str(1) '.mat']).out; 
    else
        out = load([chanDatSave 'CLOSED_' data.key '_processed_temp' num2str(1) '.mat']).out; 
    end
    

    %update the demographic info
    if contains(data.key, "S")
        subNum = strsplit(data.key, 'SCH');
        subNum = ['s' subNum{2}];
        data.group = 'SZ'; 
    else
        subNum = strsplit(data.key, 'CTOL');
        subNum = ['c' subNum{2}];
    end
    if strcmp('c14a', subNum)
        subNum = 'c14'; 
    end

    IDs = table2cell(demographics(:,1)); 
    tmpDemo = demographics(cellfun(@(x) strcmp(subNum, x),... 
                            IDs), :);
    data.age = tmpDemo.Age; 
    data.sex = tmpDemo.Gender; 
    data.para1 = tmpDemo.Ideas_Persecution; 
    data.para1_measure = "Ideas_Persecution"; 
    data.para2 = tmpDemo.Paranoia_Traits; 
    data.para2_measure = "Paranoia_Traits"; 


    trialCount = size(out.data,3);
    missingData = false; 
    %preallocate power variable
    if ~isfield(data, 'aapower')
        
        data.power = zeros(32,numfrex,trialCount); 
        missingData = true; 
    end
    
    %preallocate PAC variables
    if ~isfield(data, 'aaMINull')
        data.PAC = zeros(10,21,18,32); %low freq X high freq X phase bins X channels
        data.PACMI = cell(32,1); 
        data.PACsig = cell(32,1); 
        data.phasePref = cell(32,1); 
        data.MINull = cell(32,1); 
        missingData = true; 
      
    end

    %preallocate entropy variable
    if ~isfield(data, 'aasampEnt')
        data.sampEnt = zeros(32,trialCount, 20); 
        missingData = true; 
    end

    if ~isfield(data, 'aafuzEnt')
        data.fuzEnt = zeros(32, trialCount, 20); 
        data.sampEnt = zeros(32, trialCount, 20); 
        missingData = true; 
    end

    %preallocate 1/f slope and peak alpha fits

    if ~isfield(data, 'aaslopeValsLog')
        data.slopeValsLog = zeros(32,2); 
        data.slopeValsRel = zeros(32,2); 
        data.alphaPeakLog = zeros(32,1); 
        data.alphaPeakRel = zeros(32,1); 
        missingData = true; 
    end





    if missingData %no point in opening all the channel files if there's no missing data
    for chan = 1:32
        try %try for the channel file in general
            if strcmp(data.eyes, 'open')
                out = load([chanDatSave 'OPEN_' data.key '_processed_temp' num2str(chan) '.mat']).out; 
            else
                out = load([chanDatSave 'CLOSED_' data.key '_processed_temp' num2str(chan) '.mat']).out; 
            end
            
    %% try for power values
            try
                data.power(chan,:,:) = out.power; 
            catch
                filenames(ii).errorFlag = [filenames(ii).errorFlag, chan]; %+0 give errors for power
            end
    %% try for PAC values
            try 
                if data.PAC(1,1,1,chan) == 0 %do we need to fill in the PAC values?
                    data.PAC(:,:,:,chan) = out.PAC;
                    data.PACMI{chan} = out.PACMI;
                    data.PACsig{chan} = out.PACsig;
                    data.phasePref{chan} = out.phasePref;
                    data.MINull{chan} = out.MINull;
                end
            catch
                filenames(ii).errorFlag = [filenames(ii).errorFlag, chan+100]; %+100 give errors for PAC
            end
            
    %% try for sample entropy values
            try 
                if data.sampEnt(chan,1,1) == 0 
                    data.sampEnt(chan,:,:) = out.sampEnt; 
                    data.fuzEnt(chan,:,:) = out.fuzEnt; 
                end
            catch
                filenames(ii).errorFlag = [filenames(ii).errorFlag, chan+200]; %+200 give errors for sampEnt
            end

    %% try for 1/f slopes and peak alpha values
            try
               if data.slopeValsLog(chan,1)==0
                   data.slopeValsLog(chan,:) = out.slopeValsLog; 
                   data.slopeValsRel(chan,:) = out.slopeValsRel; 
                   data.alphaPeakLog(chan) = out.alphaPeakLog; 
                   data.alphaPeakRel(chan) = out.alphaPeakRel; 
               end
            catch
                filenames(ii).errorFlag = [filenames(ii).errorFlag, chan+300]; %+300 give errors for slope/alpha
            end

 
        catch
        filenames(ii).errorFlag = [filenames(ii).errorFlag, chan+400]; %+400 gives global error
        
        end


    end
    end

  

    %% do a check for ISPC
    if ~isfield(data, 'ispc')
        filenames(ii).errorFlag = [filenames(ii).errorFlag, 500]; %500 means ispc missing
    end

    


% save the summary data file back out
% parsave([filenames(ii).folder '/' filenames(ii).name], data); 
save([filenames(ii).folder '/' filenames(ii).name], 'data');


disp(['subject: ' num2str(ii) ' errors: ' num2str(length(filenames(ii).errorFlag)-1) ' time: ' num2str(round(toc))])

end






