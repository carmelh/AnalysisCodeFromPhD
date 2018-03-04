function[] = furtherAnalysis(fs_spr,fs_ephys,Vm_avg,rowAvg,newTimeSPR,nbins)
%furtherAnalysis - further statistical analysis
%1. Crops both time series to be +/- time from AP 
%2. Comparing histogram distributions from row average to noise
%3. Cross correlation
%4. FFT and PSD analysis
%5. Spectral coherence


%fs_spr          scalar      sampling frequency of imaging
%fs_ephys        scalar      sampling frequency of ephys
%Vm_avg         vector      average voltage wavefrom from ephys (spikeAverage)
%rowAvg         vector      average optical wavefrom from SPR (FFT_analysis)
%newTimeSPR     vector      time points
%nbins       scalar      number of bins for the histograms


% Resample Vm_avg for cross correlation
[p,qS] = rat(fs_spr/fs_ephys,0.0001);   % Rational fraction approximation
vmAvg8kHz = resample(Vm_avg, p,qS); 

%plot new vm Avg
figure;
plot(newTimeSPR,vmAvg8kHz,'Color',[13/255 165/255 138/255],'LineWidth',1)
xlabel('Time (ms)','FontName','Times New Roman');
ylabel('Voltage (mV)','FontName','Times New Roman');
grid on

%%%%%%%%%%%%
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'YColor',[0.0509803921568627 0.647058823529412 0.541176470588235]);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(vmAvg8kHz,'Parent',axes1,'LineWidth',2,...
    'Color',[0.0509803921568627 0.647058823529412 0.541176470588235]);
ylabel({'Voltage (mV)'},'FontName','Times New Roman');

% Create axes
axes2 = axes('Parent',figure1,...
    'YAxisLocation','right',...
    'YColor',[0.0549019607843137 0.254901960784314 0.627450980392157],...
    'ColorOrder',[0 0.5 0;1 0 0;0 0.75 0.75;0.75 0 0.75;0.75 0.75 0;0.25 0.25 0.25;0 0 1],...
    'Color','none');
hold(axes2,'all');

% Create plot
plot(rowAvg,'Parent',axes2,...
    'Color',[0.0549019607843137 0.254901960784314 0.627450980392157]);
ylabel({'Average Intensity (a.u.)'},'FontName','Times New Roman');
grid on

% crop time series to be around +/- 10 ms from time of AP
prompt = 'Crop time series? Lower... ';
cropValues(1) = input(prompt);
prompt = 'Crop time series? Higher... ';
cropValues(2) = input(prompt);

cropVm = vmAvg8kHz(1:cropValues(2));
cropVm = cropVm(cropValues(1):end);
cropTime = newTimeSPR(1:cropValues(2));
cropTime = cropTime(cropValues(1):end);
cropSPR = rowAvg(1:cropValues(2));
cropSPR = cropSPR(cropValues(1):end);

[meanCropRowAvg,stdCropRowAvg,meanCropRowAvgCI,stdCropRowAvgCI] = normfit(cropSPR);
varCropRowAvg = var(cropSPR);

%save mean, stddev and variance
fileID = fopen('statistics.txt','a');
formatSpec = 'cropSPR: mean = %8.4e (confidence intervals %8.4e and %8.4e), standard deviation = %8.4e (confidence intervals %8.4e and %8.4e) and variance = %8.4e \n';
fprintf(fileID,formatSpec,meanCropRowAvg,meanCropRowAvgCI(1),meanCropRowAvgCI(2),stdCropRowAvg,stdCropRowAvgCI(1),stdCropRowAvgCI(2),varCropRowAvg)
fclose(fileID);


%% Plot cropped rowAvg
figureRowAvgCrop = figure;

% Create axes
axes1 = axes('Parent',figureRowAvgCrop,...
    'Position',[0.13 0.11 0.611666666666667 0.815]);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(cropTime,cropSPR,'Parent',axes1,...
    'Color',[14/255 65/255 160/255],...
    'DisplayName','Row Average');

xlabel('Time (ms)','FontName','Times New Roman');
ylabel({'Average Intensity (a.u.)'},'FontName','Times New Roman');

% Create axes
axes2 = axes('Parent',figureRowAvgCrop,...
    'Position',[0.784895833333333 0.11 0.120104166666666 0.815],...
    'CLim',[1 2]);
view(axes2,[270 90]);
box(axes2,'on');
hold(axes2,'all');

% Create plot
histfit(cropSPR,nbins);
ylabel({'Count'},'FontName','Times New Roman');
grid on


% Get distribution curve from histfit
pdCropSPR = fitdist(cropSPR,'normal');
[bincounts,bincenters]=hist(cropSPR,nbins);
hhS = bar(bincenters,bincounts,[min(cropSPR),max(cropSPR)],'hist');
xdS = get(hhS,'Xdata');             % Gets the x-data of the bins.
rangexS = max(xdS(:)) - min(xdS(:)); % Finds the range of this data.
binwidthS = rangexS/nbins;          % Finds the width of each bin.
nS = numel(cropSPR);
areaCropSPR = nS * binwidthS;
qS = icdf(pdCropSPR,[0.0013499 0.99865]); % three-sigma range for normal distribution
xS = linspace(qS(1),qS(2));
yS = areaCropSPR * pdf(pdCropSPR,xS);


% compare to dist of just noise
% noise same variance as rowAvg
ydb = pow2db(varCropRowAvg);
 y1 = wgn(length(cropSPR),1,ydb);
 
varNoise = var(y1);

figureNoise = figure;

% Create axes
axes1 = axes('Parent',figureNoise,...
    'Position',[0.13 0.11 0.611666666666667 0.815]);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(cropTime,y1,'Parent',axes1,...
    'Color',[14/255 65/255 160/255],...
    'DisplayName','Row Average');

xlabel('Time (ms)','FontName','Times New Roman');
ylabel({'Average Intensity (a.u.)'},'FontName','Times New Roman');

% Create axes
axes2 = axes('Parent',figureNoise,...
    'Position',[0.784895833333333 0.11 0.120104166666666 0.815],...
    'CLim',[1 2]);
view(axes2,[270 90]);
box(axes2,'on');
hold(axes2,'all');

% Create plot 
histfit(y1,nbins);
ylabel({'Count'},'FontName','Times New Roman');
grid on


% Get distribution curve from histfit for noise
pdNoise = fitdist(y1,'normal');
[bincounts,bincenters]=hist(y1,nbins);
hhN = bar(bincenters,bincounts,[min(y1),max(y1)],'hist');
xdN = get(hhN,'Xdata');             % Gets the x-data of the bins.
rangexN = max(xdN(:)) - min(xdN(:)); % Finds the range of this data.
binwidthN = rangexN/nbins;          % Finds the width of each bin.
nN = numel(y1);
areaN = nN * binwidthN;
qN = icdf(pdNoise,[0.0013499 0.99865]); % three-sigma range for normal distribution
xN = linspace(qN(1),qN(2));
yN = areaN * pdf(pdNoise,xN);


% Plot to compare
figure;
plot(xS,yS,'Color',[14/255 65/255 160/255],...
    'DisplayName','SPR Row Average','LineWidth',2)
hold on
plot(xN,yN,'Color',[59/255 113/255 86/255],...
    'DisplayName','White Gaussian Noise','LineWidth',2)
ylabel({'Count'},'FontName','Times New Roman');
xlabel({'Average Intensity (a.u.)'},'FontName','Times New Roman');
legend show
grid on

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xCor
[C1,lag1] = xcorr(cropVm,cropSPR);        

figure('units','normalized','outerposition',[pos]);
bar((lag1/fs_spr)*1000,C1,'EdgeColor',[14/255 160/255 109/255],'FaceColor',[14/255 160/255 109/255])
box on
grid on
ylabel('Amplitude') 
xlabel('Time (ms)') 


%%%%%%%%%%%%%%%%%
%% FFT & PSD Analysis
dFs = fs_spr/length(cropTime);      % freq increment 
f_spr = 0:dFs:fs_spr/2-dFs;       % frequency domain

dFe = fs_ephys/length(cropTime);      % freq increment 
f_ephys = 0:dFs:fs_ephys/2-dFe;       % frequency domain

fftRowAvg = abs(fft(cropSPR));
fftVmAvg = abs(fft(cropVm));
fftVmAvg32 = abs(fft(Vm_avg));

figure('units','normalized','outerposition',[pos]);
subplot(3,1,1)
plot(f_ephys,fftVmAvg32(1:length(f_ephys)),'LineWidth',2,...
    'Color',[0.627450980392157 0.0549019607843137 0.729411764705882],...
    'DisplayName','Original Vm Average');
title('FFT of Electrophysiology Average (32 kHz)')

subplot(3,1,2)
plot(f_spr,fftVmAvg(1:length(f_spr)),'LineWidth',2,...
    'Color',[0.0509803921568627 0.647058823529412 0.541176470588235],...
    'DisplayName','Resampled');
title('FFT of Resampled (8 kHz) Electrophysiology Average')

subplot(3,1,3)
plot(f_spr,fftRowAvg(1:length(f_spr)),'LineWidth',2,...
    'Color',[14/255 65/255 160/255],...
    'DisplayName','SPR Row Average');
title('FFT of SPR Row Average')
xlabel('Frequency (Hz)','FontName','Times New Roman');

% same plot
figure('units','normalized','outerposition',[pos]);
plot(f_ephys,fftVmAvg32(1:length(f_ephys)),'LineWidth',2,...
    'Color',[0.627450980392157 0.0549019607843137 0.729411764705882],...
    'DisplayName','Original Vm Average');
hold on
plot(f_spr,fftVmAvg(1:length(f_spr)),'LineWidth',2,...
    'Color',[0.0509803921568627 0.647058823529412 0.541176470588235],...
    'DisplayName','Resampled');
xlabel('Frequency (Hz)','FontName','Times New Roman');
grid on
legend show

% Normalising FFTs
normData = [fftVmAvg',fftRowAvg];
normcData = normc(normData);
figure('units','normalized','outerposition',[pos]);
plot1=plot(f_spr,normcData(1:length(f_spr),:),'LineWidth',2);
set(plot1(1),...
    'Color',[0.0509803921568627 0.647058823529412 0.541176470588235],...
    'DisplayName','Vm Average');
set(plot1(2),...
    'Color',[0.0549019607843137 0.254901960784314 0.627450980392157],...
    'DisplayName','SPR Row Average');
xlabel('Frequency (Hz)','FontName','Times New Roman');
grid on
legend show

differenceFFT = normcData(:,1) - normcData(:,2);

figure('units','normalized','outerposition',[pos]);
plot(f_spr,differenceFFT(1:length(f_spr)),'LineWidth',2,...
    'Color',[255/255 100/255 100/255])
xlabel('Frequency (Hz)','FontName','Times New Roman');
ylabel({'Amplitude (norm)'},'FontName','Times New Roman');
grid on


%% Spetral Coherence
%PSD
[P1,f1] = periodogram(cropVm,[],[],fs_spr,'power');
[P2,f2] = periodogram(cropSPR,[],[],fs_spr,'power');

figure('units','normalized','outerposition',[pos]);
subplot(2,2,1)
plot(cropTime,cropVm,'Color',[0.0509803921568627 0.647058823529412 0.541176470588235],...
    'LineWidth',2)
ylabel('Voltage (mV)','FontName','Times New Roman')
grid on
title('Time Series')
subplot(2,2,3)
plot(cropTime,cropSPR,'Color',[0.0549019607843137 0.254901960784314 0.627450980392157],...
    'LineWidth',1)
ylabel('Average Light Intensity (a.u.)')
grid on
xlabel('Time (ms)','FontName','Times New Roman')
subplot(2,2,2)
semilogy(f1,P1,'Color',[0.0509803921568627 0.647058823529412 0.541176470588235],...
    'LineWidth',1)
ylabel('Amplitude','FontName','Times New Roman')
grid on
axis tight
title('Power Spectrum')
subplot(2,2,4)
semilogy(f2,P2,'Color',[0.0549019607843137 0.254901960784314 0.627450980392157],...
    'LineWidth',1)
ylabel('Amplitude','FontName','Times New Roman')
grid on
axis tight
xlabel('Frequency (Hz)','FontName','Times New Roman');


[Cxy,f] = mscohere(cropVm,cropSPR,[],[],[],fs_spr);
Pxy     = cpsd(cropVm,cropSPR,[],[],[],fs_spr);
phase   = -angle(Pxy)/pi*180;
[pks,locs] = findpeaks(Cxy,'MinPeakHeight',0.8);

figure('units','normalized','outerposition',[pos]);
subplot(2,1,1)
plot(f,Cxy,'Color',[89/255 51/255 84/255],'LineWidth',2)
title('Coherence Estimate')
grid on
hold on

subplot(2,1,2)
plot(f,phase,'Color',[89/255 51/255 84/255],'LineWidth',2)
title('Cross-spectrum Phase (deg)')
grid on
hold on
xlabel('Frequency (Hz)')

