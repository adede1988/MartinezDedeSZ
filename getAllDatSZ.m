
codePre = 'G:/My Drive/GitHub/';
datPre = 'H:\SZ_anton_data\EEG Resting State Data/';

addpath([codePre 'MartinezDedeSZ'])
summaryDatSave = [datPre 'SUMDAT/'];
chanDatSave = [datPre 'CHANFILES/'];

filenames = dir(chanDatSave);
filenames = filenames(contains({filenames.name}, '.mat'),:);


parfor ii = 1:length(filenames)
tic
singleChanPipelineSZ(filenames(ii));
toc   
end

