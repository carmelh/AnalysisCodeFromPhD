function[num_images,Npx, imagesDbl] = ImportTiff_v4()
% Import Tiff Stack returns an array with xy in first two dimensions, stack
% in third dimension (time)

[filename_spr, pathname_spr] = uigetfile({'*.tif', 'Image Files (*.tif)'}, 'Choose the image to open.'); % Open Image File
name_spr=fullfile(pathname_spr,filename_spr);


info = imfinfo(name_spr);
[img_x, img_y]=size(imread(name_spr, 1, 'Info', info));
img_type=['uint' num2str(info(1).BitDepth)];
num_images = numel(info);
Npx = img_x * img_y;

%preallocate array of correct type and size to improve speed
xyt_stack=zeros([img_x img_y num_images], img_type);
    
for k = 1:num_images
    A = imread(name_spr, k, 'Info', info);
    xyt_stack(:,:,k)=A(:,:);
end

imagesDbl = im2double(xyt_stack);         % Converting current image to double


