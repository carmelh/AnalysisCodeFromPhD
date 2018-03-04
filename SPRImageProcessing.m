function [ APs , Npx] = SPRImageProcessing( imagesDbl, num_images, num_images_pt, fs_spr, timeSPR, Npx )
%SPRIMAGEPROCESSING processes the optical images from the SPR data
%   options to crop the images
%   options to bandpass filter the images or to remove control area from cell ROI

%imagesDbl       3d matrix   image data 3D matrix - pixel x pixel x time
%num_images      scalar      total number of images
%num_images_pt   scalar      number of images per separate trigger
%fs_spr      scalar      frequency of the imaging
%timeSPR     vector      time points
%Npx         scalar      number of pixels

figure; imshow(imagesDbl(:,:,1));

reply = input('Do you want to crop the Images? (Y/N): ','s');
if strcmp(reply,'Y')
  disp('Please Select ROI')
  imshow(imagesDbl(:,:,1))
  cellROIh = imrect;
  positionCell = wait(cellROIh);
  close 
  
  disp('Please Select control ROI') 
  imshow(imagesDbl(:,:,1))
  contROIh = imrect;
  positionCont = wait(contROIh);
  close
      for k = 1:num_images
        cellROI(:,:,k) = imcrop(imagesDbl(:,:,k), positionCell);
        contROI(:,:,k) = imcrop(imagesDbl(:,:,k), positionCont);
      end
      
    images2D = (reshape(cellROI,[],num_images))'; % Reshape imagestack to 1D vector x_pix * y_pix by time
    Npx = size(images2D,2);
else
    images2D = (reshape(imagesDbl,[],num_images))';
end 
  

reply = input('How do you want to process the Images? (F)ilter, (A)verage, (N)othing: ','s');
if strcmp(reply,'F')
% Apply filter to columns
    imagesFiltered = bandpass(images2D, fs_spr, timeSPR);

% reshape to 3d matrix - pixel x time x separate AP
    APs = permute(reshape(imagesFiltered.', Npx,num_images_pt,[]),[2 1 3]);

elseif strcmp(reply,'A')
    % Need to scale the averages to subtract control area from cell ROI
    % Taking average of all pixels over time across both the control and cell area 
    meanContROI = mean(mean(contROI));                % get mean over entire control ROI over time
    mean2ContROI = mean2(contROI);
    meanCellROI = mean(mean(cellROI));                % get mean over entire cell ROI over time
    mean2CellROI = mean2(cellROI);

    % Scaling the control ROI
    scaledContROI = meanContROI + (mean2CellROI - mean2ContROI);
    
    % need to make scaled control ROI same size as cell region of interest to
    % subtract it 
    [xCell, yCell, zCell] = size(cellROI);
    scaledContROIResize=zeros([xCell yCell zCell]);
    
    % for loop to fill matrix with values from scaled control ROI  
    for k = 1:num_images
        scaledContROIResize(:,:,k) = scaledContROI(:,:,k); 
    end

    imagesSubtracted = cellROI - scaledContROIResize;  % subtracting control area from cell area    

    imagesAveraged = (reshape(imagesSubtracted,[],num_images))';
    APs = permute(reshape(imagesAveraged.', Npx,num_images_pt,[]),[2 1 3]);
else
% No filter
    APs = permute(reshape(images2D.', Npx,num_images_pt,[]),[2 1 3]);
end

end

