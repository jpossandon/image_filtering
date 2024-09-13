%% *IMAGE FILTERING BPN LAB*
% 1.0.0
%
% By José Ossandón (<jose.ossandon@uni-hamburg.de>)
%% Image preparation
%
% Let's first open an example image and show it. The original image is 400x400
% pixels but we are going to reshape it to a non-square size for the
% examples below.
clear
imagePath          = '/Users/jossando/trabajo/CSF/KDEF_400_final_cut/400px_600_AF31NEHL.bmp';
thisImage          = imread(imagePath);
[m,n,p]            = size(thisImage);
imshow(thisImage)
title(sprintf('%dx%d image',m,n))
%%
% Rehsaping of the images to 600x650. This makes it bigger but also changes
% a bit the aspect ratio. In generla, to change the aspect ratio is not
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
% image has already a background color that needs to be taken in account.
% Our example image has a gray background (138). 
% # Image aspect ratio: to use the code presented here images need to be a
% square. Otherwise the butterworth and CSF filters with not work
% correctly. This can be achieved in two main ways, by preparing the image
% before filtering or by setting the filtering with padding.
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

%% Image filtering: guassian filter with image_filter.m
%
% The function image_filter.m requires two inputs and image to fitler  and a 'options' structure 
% which defines the different parameters of the function (see documentation).
% To filter an image is first necessary to define the image resolution in
% pixels per degree. This need to be calculated from the image (screen)
% size, its size in pixels, and the distance of observation. For
% exammple:

screen_size_horizontal  = 52.7; % in some unit, here centimeters
screen_distance         = 80;   % in the same units than above
horizontal_resolution   = 1920; % in pixels
screen_vis_degree       = atan(screen_size_horizontal/screen_distance); % atan is the inverse tangent given the degrees in radian
screen_vis_degree       = screen_vis_degree*180/pi; % which can be transformed to degrees by multypling with (180/pi)

pixels_per_degree       = round(horizontal_resolution./screen_vis_degree); % and the pixells per degree can be caluclated by dividing the screen/image resolution by how many visual degrees it spans (which we round to the next integer) 
display(sprintf('At a distance of %1.1f cms, %1.2f cms cover %1.2f visual degrees.\n',screen_distance,screen_size_horizontal,screen_vis_degree))
display(sprintf('At that distance and screen sixe, %d pixels result in a resolution of %d pixels per degreee (ppd).',horizontal_resolution,pixels_per_degree ))

%%
% The filter gaussian can be set as options.filterType 'gaussian' or
% 'gaussian_custom_cutoff'. Let's look first at 'gaussian', for which we
% only need to define the cutoff frequency:

optionsG3.pixxgrade       = pixels_per_degree;
optionsG3.cutoff          = 3;
optionsG3.filterType      = 'gaussian';

[~,lpimageG3] =image_filter(square_image,optionsG3);
imshow(lpimageG3)
title(sprintf('%dx%d gray background image, %s filter cutoff %d cyc/deg cutoff',bkg_size(1),bkg_size(2),optionsG3.filterType,optionsG3.cutoff))

%%
% The cutoff for a guassian filter withoout additional specification meens
% a reduction of ~0.607 in amplitude at that frequency and with a very slow
% 'roll-off'. To see this we can add the option of computing the image and
% filtered image (radial) spectra and plotting it

optionsG3.computeRadialSpectra      = 1;
optionsG3.toplot                    = 1; %  
[~,lpimageG3,filtersG3,radialAvgFrqG3] =image_filter(square_image,optionsG3);

%%
% A reduction of only 0.607 amplitude at a given frequency might be too restricting or little. 
% Therefore we can only specifiy a desired amplitude reduction at agiven frequency. For this
% we use the filterType 'gaussian_custom_cutoff'. For example, let's filter with a reduction of 10% amplitude
% at 3 cyc/degree, and compare the images and spectra to the standard gaussian filter described above:

optionsG3_1.pixxgrade               = pixels_per_degree;
optionsG3_1.filterType              = 'gaussian_custom_cutoff';
optionsG3_1.cutoff                  = [3 .1];        % for custom_cutoff the first value is the cutoff requency and the second one the amplitude reduction           
optionsG3_1.computeRadialSpectra    = 1;
optionsG3_1.toplot                  = 0; %  
[~,lpimageG3_1,filtersG3_1,radialAvgFrqG3_1] =image_filter(square_image,optionsG3_1);

fh = figure;
fh.Position(3:4) = fh.Position(3:4)*2; % just making the figure bigger 
subplot('Position',[.05 .55 .4 .4])
imshow(lpimageG3_1)
title(sprintf('%s filter cutoff %d cyc/deg %1.1f cutoff',optionsG3_1.filterType,optionsG3_1.cutoff(1),optionsG3_1.cutoff(2)))
subplot('Position',[.55 .55 .4 .4])
imshow(lpimageG3)
title(sprintf('%s filter cutoff %d cyc/deg std cutoff',optionsG3.filterType,optionsG3.cutoff(1)))

% this plots the spectra of the two filters
subplot('Position',[.05 .05 .4 .4])
freq_spacing=.1;
xx = 0:1:length(radialAvgFrqG3{1})-1;
h(1) = plot(xx,10*log10(radialAvgFrqG3{1}.^2),'k','LineWidth',2);   % the original image spectra
hold on
h(2) = plot(xx,10*log10(radialAvgFrqG3{1}.^2.*.61.^2),'b:');        % the original image spectra with an amplitude reudction of 0.61
h(3) = plot(xx,10*log10(radialAvgFrqG3{1}.^2.*.1.^2),'r:');         % the original image spectra with an amplitude reudction of 0.1
% plot(0:1:length(radialProfile)-1,10*log10(radialProfile.^2.*ar.^2));
h(4) = plot(xx,10*log10(radialAvgFrqG3{2}.^2));                     % the image filtered with 'gauss' filter
h(5) = plot(xx,10*log10(radialAvgFrqG3_1{2}.^2));                   % the image filtered with 'gauss)custom_cutoff' filter
axis([0 8/freq_spacing -75 10])
vline([1.5 3]./freq_spacing)
legend(h,{'Original','0.61*original','0.1*original','Gauss60%','Gauss10%'},'box','off')
set(gca,'XTick',0:1./freq_spacing:length(radialAvgFrqG3{1})-1,'XTickLabel',0:1:(length(radialAvgFrqG3{1})-1)*freq_spacing)
xlabel('Spatial Frequency (cyc/deg)')
ylabel('Power (dB)')
title('Example filtering cut-off 3 cyc')

% a cross section of the filters
subplot('Position',[.55 .05 .4 .4])
hh(1) = plot(filtersG3(size(filtersG3,1)/2+1,:));, hold on
hh(2) = plot(filtersG3_1(size(filtersG3_1,1)/2+1,:));
legend(hh,{'gauss','gauss_custom_cutoff'})
title('filters cross-section')

%% Image filtering: butterwoA different type of filter
% The second type of filter available in image_filter.m is a butterworth
% filter. As the guassian filter it can be specified as a standard 'butterworth', in which
% the cutoff freuqency is where the a,plitude is reduced in half, o as a'butterwoth_custom_cutoff'
% in which the cutoff frequency amplitude reduction can be specified.
% In addition,these filters need their 'order' to be specified. Filters with a lower order (e.g.,1) have
% a very slow roll-off similar to the gaussian filter. Higer order filters
% have a faster roll-off at the expense of having 'ripples' in the passband
% and stopband, this translates to 'riginng artifacts' in the spatial
% domain. 

% Let's see first a standard 'butterworth', with thre different orders
% 1,3,6:
optionsG3_butt.pixxgrade               = pixels_per_degree;
optionsG3_butt.filterType              = 'butterworth';
optionsG3_butt.cutoff                  = [3];        % for custom_cutoff the first value is the cutoff requency and the second one the amplitude reduction  
optionsG3_butt.order                   = 1;
optionsG3_butt.computeRadialSpectra    = 1;
optionsG3_butt.toplot                  = 0; %  
[~,lpimageG3_butt_ord1,filtersG3_butt_ord1,radialAvgFrqG3_butt_ord1] =image_filter(square_image,optionsG3_butt);
figure,imshow(lpimageG3_butt_ord1)
optionsG3_butt.order                   = 5;
[~,lpimageG3_butt_ord5,filtersG3_butt_ord5,radialAvgFrqG3_butt_ord5] =image_filter(square_image,optionsG3_butt);
figure,imshow(lpimageG3_butt_ord5)
optionsG3_butt.order                   = 10;
[~,lpimageG3_butt_ord10,filtersG3_butt_ord10,radialAvgFrqG3_butt_ord10] =image_filter(square_image,optionsG3_butt);
figure,imshow(lpimageG3_butt_ord10)


% The cutoff for a guassian filter withoout additional specification meens
% a reduction of ~0.607 in amplitude at that frequency and with a very slow
% 'roll-off'. To see this we can add the option of computing the image and
% filtered image (radial) spectra and plotting it

% options.padding     = [];
% options.cutoff      = 3;
% options.filterType  = 'gaussian';% 'butterworth';%
% options.toplot      = 1; %   - 1 to toplot filtered images and the filters
% options.computeRadialSpectra      = 1;
% [hpimage, lpimageG3, filterG3, radialAvgFrqG3] =image_filter(image,options);

%%
% close all