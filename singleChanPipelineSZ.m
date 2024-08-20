function [] = singleChanPipelineSZ(fileInfo)

%overall subject summary data
out = load([fileInfo.folder '/' fileInfo.name]).out;


%frequency params
frex = logspace(log10(2),log10(80),100);
numfrex = length(frex); 
stds = linspace(2,5,numfrex);



%% power
if ~isfield(out, 'fg') %check that power has not been calculated for this channel
    out.power = chanPowerSZ(out, frex, numfrex, stds);
end

%% PAC
if ~isfield(out, 'gfh') %check that PAC has not been calculated for this channel
    
     [out.PAC, out.PACMI, out.PACsig, out.phasePref, out.MINull] = chanPACSZ(out);
    
end


%% entropy
out.fuzEnt = chanFuzEntSZ(out); 
out.sampEnt = chanEntropyValsSZ(out);

%% 1/f and peak alpha
[out.slopeValsLog, out.slopeValsRel,out.alphaPeakLog, out.alphaPeakRel] = chanSlopeAlphaSZ(out,frex);


%% 
save([fileInfo.folder '/' fileInfo.name], 'out')


  


end