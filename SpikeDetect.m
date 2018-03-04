function [SpP_offset, dt_offset] = SpikeDetect( Vm, fs_ephys, thres, c )
%SPIKEDETECT selects the timepoints of spikes based on the given threshold

%Vm         2d matrix   voltage sample points
%thres      scalar      threshold value to get rid of noise
%fs_ephys   scalar      sampling frequency
%c          scalar      number of triggers

SpPVal = zeros(1,c);
SpP = zeros(1,c);
SpSTime = zeros(1,c);
SpETime = zeros(1,c);
SpPTime = zeros(1,c);

%Separate Channel Analysis
for i = 1:c
  
    Sp = find(Vm(:,i)>thres);               %Specify indices of threshold crossing
    SpS(i) = Sp(find(diff([0; Sp])>10));    %Start points of threshold crossing, ignore multiple threshold crossings due to noise
    SpE(i) = Sp(find(diff([Sp; inf])>10));  %End points of threshold crossing
    %SpP=round((SpE-SpS)/2+SpS);             %find index of threshold crossing peaks

    [SpPVal(i),SpP(i)]= max(Vm(:,i));
    SpSTime(i) = SpS(i)/fs_ephys; % Time of start
    SpETime(i) = SpE(i)/fs_ephys; % Time of end
    SpPTime(i) = SpP(i)/fs_ephys; % Time of peak
end

for i = 1:c-1
    SpP_offset(i) = SpP(i+1) - SpP(1);   % Find offset between first spike and subsequent spikes
end

dt_offset = SpP_offset/fs_ephys; % offset wrt time

clear SpPVal SpP SpS SpE SpP SpSTime SpETime SpPTime