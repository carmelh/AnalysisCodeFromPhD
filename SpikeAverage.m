function [Vm_avg, timeAP] = SpikeAverage (Vm, SpP_offset, c, timeEphys, numTriggers)
%SPIKEAVERAGE - average the spikes to give a standard voltage response

%Vm         2d matrix   voltage sample points
%SpP_offset      scalar      offset between first spike and subsequent spikes
%c               scalar      number of sample points
%timeEphys       vector      time points
%numTriggers     scalar      number of triggers


Vm_aligned = zeros(numTriggers, c);   % assign space
Vm_aligned = [Vm(:,1);];    % add first column of Vm to aligned matrix

%aliging spikes(n+1) to the first spike
for i = 1:c-1
   if SpP_offset(i) > 0 % if ith spike comes after the first
       Vm_temp = Vm((SpP_offset(i):end),i+1);       % deleting the first elements to bring back spike to meet the first
       Vm_temp = padarray(Vm_temp,[(SpP_offset(i)-1) 0],'replicate','post'); % padding the end of the array to make same size as original
   else        % if ith spike comes before the first
       Vm_temp = circshift(Vm(:,i+1),(length(Vm)-SpP_offset(i)));  % shifting spikes before the first spike to be aligned with the first
   end
   Vm_aligned(:,i+1) = Vm_temp;
end

clear Vm_temp

Vm_avg = mean(Vm_aligned.'); % Average spike

% Time of AP
[peakVoltage, indexAP] = max(Vm_avg);
timeAP = timeEphys(indexAP,:);


figure;
subplot(2,1,1);
p=plot(timeEphys*1000,Vm_aligned);

title(['Aligned Spikes. n = ' num2str(numTriggers)])
ylabel('Voltage (mV)');
axis([0 max(timeEphys*1000) min(Vm_aligned(:))*1.1 max(Vm_aligned(:))*1.1]);
grid on
clrs = winter(numel(p)); % Nx3 array of RGB values
for ii=1:numel(p)
    set(p(ii),'color',clrs(ii,:));
end

subplot(2,1,2)
plot(timeEphys*1000,Vm_avg,'Color',[160/255 14/255 138/255],'LineWidth',1.5);
title('Average Spike');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
axis([0 max(timeEphys*1000) min(Vm_aligned(:))*1.1 max(Vm_aligned(:))*1.1]);
grid on

clear Vm_aligned