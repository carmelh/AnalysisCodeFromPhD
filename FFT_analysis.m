function[rowAvg,newTimeSPR] = FFT_analysis(APs, fs_spr, num_images_pt, Npx, numTriggers, dt_offset, timeSPR, timeAP,nbins)
%FFT_ANALYSIS - 1. applys hanning window to the image data. 
%2. aligns the images using an FFT. 
%3. average over all the repeats

%APs             3d matrix   pixels x time x separate AP
%fs_spr          scalar      sampling frequency of imaging
%num_images_pt   scalar      number of images per separate trigger
%numTriggers     scalar      number of triggers
%Npx         scalar      number of pixels
%dt_offset   vector      number of pixels
%timeSPR     vector      time points
%timeAP      scalar      time at which the AP occured
%nbins       scalar      number of bins for the histograms


% hanning window 
winlen = round(num_images_pt); % length of the analysis window
temp_hann = hanning(fs_spr*0.01); % for the edges - arbitary
winfun = [temp_hann(1:length(temp_hann)/2); ones(winlen-length(temp_hann),1); ...
    temp_hann(length(temp_hann)/2+1:end)];

APs_han = bsxfun(@times,winfun,APs); % apply hanning window along t-axis of each AP

% harmonics for FFT
h_ord=cat(2,linspace(0,num_images_pt/2-1,num_images_pt/2),num_images_pt/2,-fliplr(linspace(0,num_images_pt/2-1,num_images_pt/2)));
h_ord=h_ord(1,1:end-1);

Images_aligned = zeros(num_images_pt, Npx, numTriggers); % assigning space
Images_aligned(:,:,1) = APs_han(:,:,1);
 
% FFT
for ev=2:numTriggers
    dt=dt_offset(1,ev-1);
    Fg=fft2(squeeze(APs_han(:,:,ev)));

    Images_aligned(:,:,ev)=real(ifft2(Fg.*repmat(exp(-1j*2*pi*fs_spr*h_ord*dt/num_images_pt)',[1 Npx])));
end;

% moving time points so AP is at time 0
newTimeSPR = ([0-timeAP :max(timeSPR)/(num_images_pt-1): (max(timeSPR) - timeAP)])*1000;

%% Single AP
imagesAligned1 = Images_aligned(:,:,1); % all pixels for 1st AP
AP1 = mean(imagesAligned1,2);
figureAP1 = figure;

% Create axes
axes1 = axes('Parent',figureAP1,...
    'Position',[0.13 0.11 0.611666666666667 0.815]);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(newTimeSPR,AP1,'Parent',axes1,...
    'Color',[0.0549019607843137 0.552941176470588 0.643137254901961],...
    'DisplayName','Row Average');

xlabel('Time (ms)','FontName','Times New Roman');
ylabel({'Intensity (a.u.)'},'FontName','Times New Roman');

% Create axes
axes2 = axes('Parent',figureAP1,...
    'Position',[0.784895833333333 0.11 0.120104166666666 0.815],...
    'CLim',[1 2]);
view(axes2,[270 90]);
box(axes2,'on');
hold(axes2,'all');

% Create plot
histfit(AP1,nbins);
ylabel({'Count'},'FontName','Times New Roman');
grid on

[mean1AP,std1AP,mean1APCI,std1APCI] = normfit(AP1);
var1AP= var(AP1);


%% Average APs
Images_avg = mean(Images_aligned,3);  % Average across 3rd dimension
pixel_num = [1:Npx];
figure; pcolor(pixel_num,(newTimeSPR),Images_avg); shading flat; xlabel('Pixel Number'); ylabel('Time (ms)'); %pcolor plot of all averaged pixels 

% average of each row
rowAvg = mean(Images_avg,2);

% Stats
[meanRowAvg,stdRowAvg,meanRowAvgCI,stdRowAvgCI] = normfit(rowAvg);
varRowAvg = var(rowAvg);

%save mean, stddev and variance
fileID = fopen('statistics.txt','w');
formatSpec = '1 AP: mean = %8.4e (confidence intervals %8.4e and %8.4e), standard deviation = %8.4e (confidence intervals %8.4e and %8.4e) and variance = %8.4e \n';
fprintf(fileID,formatSpec,mean1AP,mean1APCI(1),mean1APCI(2),std1AP,std1APCI(1),std1APCI(2),var1AP)
formatSpec = 'rowAvg: mean = %8.4e (confidence intervals %8.4e and %8.4e), standard deviation = %8.4e (confidence intervals %8.4e and %8.4e) and variance = %8.4e \n';
fprintf(fileID,formatSpec,meanRowAvg,meanRowAvgCI(1),meanRowAvgCI(2),stdRowAvg,stdRowAvgCI(1),stdRowAvgCI(2),varRowAvg)
fclose(fileID);


%% Plot rowAvg
figureRowAvg = figure;

% Create axes
axes1 = axes('Parent',figureRowAvg,...
    'Position',[0.13 0.11 0.611666666666667 0.815]);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(newTimeSPR,rowAvg,'Parent',axes1,...
    'Color',[14/255 65/255 160/255],...
    'DisplayName','Row Average');

xlabel('Time (ms)','FontName','Times New Roman');
ylabel({'Average Intensity (a.u.)'},'FontName','Times New Roman');

% Create axes
axes2 = axes('Parent',figureRowAvg,...
    'Position',[0.784895833333333 0.11 0.120104166666666 0.815],...
    'CLim',[1 2]);
view(axes2,[270 90]);
box(axes2,'on');
hold(axes2,'all');

% Create plot
histfit(rowAvg,nbins);
ylabel({'Count'},'FontName','Times New Roman');
grid on

