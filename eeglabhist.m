% EEGLAB history file generated on the 01-Sep-2024
% ------------------------------------------------
eeglab('redraw');
EEG = pop_loadcnt('H:\SZ_anton_data\EEG Resting State Data\CTOL1.cnt' , 'dataformat', 'auto', 'memmapfile', '');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
eeglab redraw;
