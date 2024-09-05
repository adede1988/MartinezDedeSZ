

codePre = 'G:/My Drive/GitHub/';
datPre = 'H:\SZ_anton_data\EEG Resting State Data/';



files = dir('H:\SZ_anton_data\EEG Resting State Data\SUMDAT');


files = files(cellfun(@(x) contains(x, 'CLOSED'), {files.name}));

exampleDat = load([datPre 'exampleChanLocs.mat']).chanLocs; 
labels = {exampleDat.labels}; 
Th = [exampleDat.theta]; 
Rd = [exampleDat.radius]; 



meanClosed = zeros(length(files), 32,100); 

for ii = 1:length(files)

    data = load([files(ii).folder '\' files(ii).name]).data; 
    meanClosed(ii, :,:) = squeeze(mean(data.power(:,:,:), 3)); 

end

figure
quickTopo(squeeze(mean(meanClosed(:,:,40), 1)),Th, Rd, 1)