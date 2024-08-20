function [sampEnt] = chanEntropyValsSZ(out)
%input: 
%       out     a struct with the fields data and srate (not used)
%               assumes that data is a channels X time X epochs array
%               and that each epoch is 2000ms long regardless of sampling
%               rate. This is critical, because all data are resampled to
%               1000 Hz sampling rate, but the way the code is written,
%               this only works with the assumption of 2s long epochs. 

%output: 
%       sampEnt Channels X epochs X scale; scale is hard coded to range
%       from 1 to 20. 

%function will first check to see if sampleEnt has already been calcualted
%for this data file. If yes, then don't waste time recalcualting it.
%Otherwise, do the calculation. 

%Note, the function sampEnt.m relies on functions coursegraining.m and
%msentropy.m. All three of these functions were obtained from 
% https://physionet.org/content/sampen/1.0.0/
% based on work in: 
% Lake, D. E., J. S. Richman, M. P. Griffin, and J. R. Moorman. 
% Sample entropy analysis of neonatal heart rate variability. 
% Am J Physiol 2002; 283(3):R789-R797

%this function was written by Adam Dede (adam.osman.dede@gmail.com)
%fall, 2022


if isfield(out, 'sampEnt')
   sampEnt = out.sampEnt; 
%    fuzEnt = out.fuzEnt; 
%    compSampEnt = out.compSampEnt; 
%    compFuzEnt = out.compFuzEnt; 
else
  
   sampEnt = nan([size(out.data,1), size(out.data,3), 20 ]);

   
   for snip = 1:size(out.data,3) %loop epochs
        curSnip = out.data(:,:,snip); 
        %resample to 1000 Hz sampling rate
        times = round(linspace(1,2000, size(out.data,2))); 
        curSnip = timeseries(curSnip', times ); 
        curSnip = resample(curSnip, [1:2000]);
        curSnip = curSnip.Data'; 
      
          
               
       sampEnt(1,snip, :) = msentropy(curSnip,2, .3, 20);
               
  

   end

  
 

end
end