function [] = singleChanPipelineSZ(fileInfo)

%overall subject summary data
out = load([fileInfo.folder '/' fileInfo.name]).out;


%frequency params
frex = logspace(log10(2),log10(80),100);
numfrex = length(frex); 
stds = linspace(2,5,numfrex);
change = false; 


%% power
if ~isfield(out, 'power') %check that power has not been calculated for this channel
    out.power = chanPowerSZ(out, frex, numfrex, stds);
    change = true; 
end

%% PAC
if ~isfield(out, 'PAC') %check that PAC has not been calculated for this channel
    
     [out.PAC, out.PACMI, out.PACsig, out.phasePref, out.MINull] = chanPACSZ(out);
    change = true; 
    
end


%% entropy
if ~isfield(out, 'fuzEnt')
    out.fuzEnt = chanFuzEntSZ(out); 
    out.sampEnt = chanEntropyValsSZ(out);
    change = true; 
end
%% 1/f and peak alpha
if ~isfield(out, 'alphaPeakLog')
    [out.slopeValsLog, out.slopeValsRel,out.alphaPeakLog, out.alphaPeakRel] = chanSlopeAlphaSZ(out,frex);
    change = true; 
end

%% save out if anything has been added
if change
    save([fileInfo.folder '/' fileInfo.name], 'out')
end

  


end