%% *IMAGE FILTERING BPN LAB*
% 1.0.0
%
% By José Ossandón (<jose.ossandon@uni-hamburg.de>)
%% Image preparation
%
% Let's first open an example image and show it. The original image is 400x400
% pixels but we are going to reshape it to a non-square size for the
% examples below.

imagePath          = '/Users/jossando/trabajo/CSF/KDEF_400_final_cut/400px_600_AF31NEHL.bmp';
thisImage          = imread(imagePath);
[m,n,p]            = size(thisImage);
imshow(thisImage)
title(sprintf('%dx%d image',m,n))
%%
% Rehsaping of the images to 600x650. This makes it bigger but also changes
% a bit the aspect ratio. In generla, to change the aspect ration will not
% be recommendable, but here we are doing to learn how to deal with
% non-square images
targetSizeInPix    = [600 650];      % be careful with the ordering of numbers, MATLAB works with (rows,columns) but screen resolution are usually expressed in width,height    
thisImage          = imresize(thisImage,targetSizeInPix);
[m,n,p]            = size(thisImage);
imshow(thisImage)
title(sprintf('%dx%d image',m,n))
%%
% There is two thing to have in mind when preparing images to be filtered:
%
% # Image background: if the images are not going to be presented to the
% full size of the screen, it is important to take care which will be the
% background color, usually set to black or the midlle gray level (level of 0.5 in a scale of 0 to 1 
% or 127/128 in a scale between 0 and 5). It can also be that the original
% image has already a background color that need to be taken in account.
% Our example image has a gray background (138). 
% # Image aspect ratio: to use the code presented here images need to be a
% square. Otherwise the butterworth and CSF filters with not work
% correctly. This can be achieved in two main ways, by preparing the image
% before filtering or by set the filtering with padding.
% 
% Here, we are going to first put our 600x650 image in a larger square (800x800) backgound
% of the same color
bkg_size            = [800 800];
bkg_color           = 138;
square_image        = uint8(repmat(bkg_color.*ones(bkg_size(1),bkg_size(2)),1,1,3));
colIndxs            = bkg_size(1)./2+1-n/2:bkg_size(1)./2+n/2;
rowIndxs            = bkg_size(2)./2+1-m/2:bkg_size(2)./2+m/2;

for rbg = 1:p
    square_image(rowIndxs,colIndxs,rbg) = thisImage(:,:,rbg);
end
imshow(square_image)
title(sprintf('%dx%d image centered in a %dx%d gray background',m,n,bkg_size(1),bkg_size(2)))
%% 
% Alternatively, this is how it would look with a black background
bkg_color            = 0;
square_image_black  = uint8(repmat(bkg_color.*ones(bkg_size(1),bkg_size(2)),1,1,3));

for rbg = 1:p
    square_image_black(rowIndxs,colIndxs,rbg) = thisImage(:,:,rbg);
end
imshow(square_image_black)
title(sprintf('%dx%d image centered in a %dx%d black background',m,n,bkg_size(1),bkg_size(2)))
