%%
%  imagePath       = '/Users/jossando/trabajo/FaceEEG/2 experiments scripts/1_FO_loc/stimuli/Grey/face_41.jpg';
% screenSize          = [1080 1080];
% options.pixxgrade   = 50;
% 
%  image           = .5.*ones(screenSize(2),screenSize(1));
%   image(screenSize(2)/2-200+1:screenSize(2)/2+200,screenSize(1)/2-200+1:screenSize(1)/2+200) = imresize(double(imread(imagePath))./255,2);
 %%
 imagePath          = '/Users/jossando/trabajo/CSF/KDEF_400_final_cut/400px_600_AF31NEHL.bmp';
 % imagePath          = '/Users/jossando/Desktop/Screenshot 2024-09-05 at 10.15.54.png';
finalSize          = [1100 1100];
 targetSizeInPix     = [1050 1050];
bkg_color           = 138;
 thisImage          = imresize(imread(imagePath),targetSizeInPix );
image               = uint8(repmat(bkg_color.*ones(finalSize(2),finalSize(1)),1,1,3));
colIndxs            = finalSize(2)./2+1-targetSizeInPix(1)/2:finalSize(2)./2+targetSizeInPix(1)/2;
rowIndxs            = finalSize(1)./2+1-targetSizeInPix(2)/2:finalSize(1)./2+targetSizeInPix(2)/2;;

for rbg = 1:3
    image(colIndxs,rowIndxs,rbg) = thisImage(:,:,rbg);
   
end
options.pixxgrade   = 70;
fh = figure;
imshow(image)
set(fh,'Position',[0 0 finalSize])
set(gca,'Position',[0 0 1 1])
 %%
options.pixxgrade   = 70*2;
options.padding    = [];
options.cutoff      = 3;
options.filterType  = 'gaussian';% 'butterworth';%
options.toplot      = 1; %   - 1 to toplot filtered images and the filters
options.computeRadialSpectra      = 1;
[hpimage, lpimage, filter, radialAvgFrq] =image_filter(image,options);
fh = figure,imshow(lpimage)
set(fh,'Position',[0 0 finalSize])
toc
%%
 options.pixxgrade   = 70*2;
 options.padding    = [];

options.cutoff      = [3 .1];
options.filterType  = 'gaussian_custom_cutoff';% 'butterworth';%
[hpimage, lpimage, filter, radialAvgFrq] =image_filter(image,options);

fh=figure,imshow(lpimage)
set(fh,'Position',[0 0 finalSize])
%
set(gca,'Position',[0 0 1 1])
%%
options.order       = 6;
options.filterType  = 'butterworth_custom_cutoff';% 'butterworth';%
[hpimage, lpimage, filter, radialAvgFrq] =image_filter(image,options);
fh=figure,imshow(lpimage)
set(fh,'Position',[0 0 finalSize])
%
set(gca,'Position',[0 0 1 1])
%%
options.cutoff      = [6 .1];
options.filterType  = 'gaussian_custom_cutoff';% 'butterworth';%
[hpimage, lpimage, filter, radialAvgFrq] =image_filter(image,options);
figure,imshow(lpimage)
%%
options.cutoff      = 3;
options.order       = 3;
options.filterType  = 'butterworth';% 'butterworth';%
[hpimage, lpimage, filter, radialAvgFrq] =image_filterbk(uint8(image.*255),options);

options.cutoff      = 3;
options.order       = 5;
options.filterType  = 'butterworth';% 'butterworth';%
[hpimage, lpimage, filter, radialAvgFrq] =image_filterbk(uint8(image.*255),options);

%%
options.cutoff      = [3 .1];

options.filterType  = 'gaussian_custom_cutoff';% 'butterworth';%
[hpimage, lpimage, filter, radialAvgFrq] =image_filterbk(uint8(image.*255),options);

options.order       = 4;
options.filterType  = 'butterworth_custom_cutoff';% 'butterworth';%
[hpimage, lpimage, filter, radialAvgFrq] =image_filter(uint8(image.*255),options);

options.order       = 6;
options.filterType  = 'butterworth_custom_cutoff';% 'butterworth';%
[hpimage, lpimage, filter, radialAvgFrq] =image_filter(uint8(image.*255),options);







%%
options.pixxgrade   = 40;
options.cutoff      = [6 .1];
options.filterType  = 'guassian_custom_cutoff';% 'butterworth';%
options.toplot      = 1; %   - 1 to toplot filtered images and the filters

n = 1920;
m = 1080;
sinf = 6;
% degperim  = 40;
x = 0:1./options.pixxgrade:n./options.pixxgrade-1./options.pixxgrade;
y = 0:1./options.pixxgrade:m./options.pixxgrade-1./options.pixxgrade;
ima = repmat(0.25*(1+sin(2*pi*sinf*x)),m,1)+repmat(0.25*(1+sin(2*pi*.5*sinf*y))',1,n);
figure,
imshow(ima,[])
set(gca, 'XTick',0:options.pixxgrade:n,'YTick',0:options.pixxgrade:m)
axis on, grid on

[hpimage, lpimage, filter, radialAvgFrq] =image_filter(uint8(ima.*255),options);

%%
