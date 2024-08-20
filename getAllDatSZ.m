%local: 
% prefix = 'Z:/User/pc1aod/'; subNum = 1; chan = 1; 
%HPC: 
codePre = 'G:/My Drive/GitHub/';
datPre = 'G:/My Drive/Milne/SZproject/';

addpath([codePre 'MartinezDedeSZ'])
summaryDatSave = [datPre 'SUMDAT/'];
chanDatSave = [datPre 'CHANFILES/'];

filenames = dir(chanDatSave);
filenames = filenames(contains({filenames.name}, '.mat'),:);


for ii = 1:length(filenames)
tic
singleChanPipelineSZ(filenames(ii));
toc   
end

