
codePre = 'G:/My Drive/GitHub/';
datPre = 'G:/My Drive/Milne/SZproject/';

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

