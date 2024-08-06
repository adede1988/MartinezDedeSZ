% EEGLAB history file generated on the 06-Aug-2024
% ------------------------------------------------
EEG = pop_loadset('filename','C2_ChannelLoc_Filt.set','filepath','G:\\My Drive\\Milne\\SZproject\\');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
pop_saveh( EEG.history, 'eeglabhist.m', 'G:\My Drive\GitHub\MartinezDedeSZ\');
eeglab redraw;
