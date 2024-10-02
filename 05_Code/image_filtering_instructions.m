%% *IMAGE FILTERING BPN LAB*
% 1.0.0
%
% By José Ossandón (<jose.ossandon@uni-hamburg.de>)
%
% To reproduce this script it is necessary to have the filtering function
% contianed in the following github repository https://github.com/jpossandon/image_filtering.git
% and https://github.com/jpossandon/csf_psi.git
%% Image preparation
%
% Let's first open an example image and show it. The original image is 400x400
% pixels but we are going to reshape it to a non-square size for the
% examples below.
clear
thisScriptDir      = fileparts(mfilename('fullpath'))
imagePath          = fullfile(thisScriptDir,'400px_600_AF31NEHL.bmp');
thisImage          = imread(imagePath);
[m,n,p]            = size(thisImage);
imshow(thisImage)
title(sprintf('%dx%d image',m,n))
%%
% Reshaping the image to 600x650. This makes it bigger but also changes
% a bit the aspect ratio. In general, to change the aspect ratio is not
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
% or 127/128 in a scale between 0 and 255). 
% 
% It can also be that the original
% image has already a background color that needs to be taken in account.
% Our example image has a gray background (138/255). 
%
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
% Alternatively, this is how it would look like with a black background
bkg_color            = 0;
square_image_black  = uint8(repmat(bkg_color.*ones(bkg_size(1),bkg_size(2)),1,1,3));

for rbg = 1:p
    square_image_black(rowIndxs,colIndxs,rbg) = thisImage(:,:,rbg);
end
imshow(square_image_black)
title(sprintf('%dx%d image centered in a %dx%d black background',m,n,bkg_size(1),bkg_size(2)))

%% Image filtering: guassian filter with image_filter.m
%
% The function image_filter.m requires two inputs, an image to filter and an 'options' structure 
% which defines the different parameters of the function (see  image_filter.m documentation for details).
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
display(sprintf('At that distance and screen size, %d pixels result in a resolution of %d pixels per degreee (ppd).',horizontal_resolution,pixels_per_degree ))

%%
% The filter gaussian can be set  options.filterType 'gaussian' or
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
% filtered image (radial) spatial spectra and plotting it

optionsG3.computeRadialSpectra      = 1;
optionsG3.toplot                    = 1; %  
[~,lpimageG3,filtersG3,radialAvgFrqG3] =image_filter(square_image,optionsG3);

%%
% A reduction of only 0.607 amplitude at a given frequency might be too restricting or little. 
% Therefore we can only specifiy a desired amplitude reduction at a given frequency. For this
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
title(sprintf('%s filter cutoff %d cyc/deg %1.1f amp. red.',optionsG3_1.filterType,optionsG3_1.cutoff(1),optionsG3_1.cutoff(2)))
subplot('Position',[.55 .55 .4 .4])
imshow(lpimageG3)
title(sprintf('%s filter cutoff %d cyc/deg std 0.607 amp. red.',optionsG3.filterType,optionsG3.cutoff(1)))

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

%%
% Note that the filter with a reduction of 0.1 ampitude, not only reduce
% much more the amplitude at the cut-off frequency but also at lower
% frequencies, resulting in a much more blurred image.

%% Image filtering: butterworth filter with image_filter.m
% The second type of filter available in image_filter.m is a butterworth
% filter. As with the guassian filter, it can be specified as a standard 'butterworth', in which
% the cutoff frequency is where the amplitude is reduced 0.71, or as a 'butterwoth_custom_cutoff'
% in which the cutoff frequency amplitude reduction can be specified.
% In addition,these filters need that their 'order' to be specified. Filters with a lower order (e.g.,1) have
% a very slow roll-off, similar to the gaussian filter. Higer order filters
% have a faster roll-off at the expense of having 'ripples' in the passband
% and stopband, this translates to 'ringinng artifacts' in the spatial
% domain. 

% Let's see first a standard 'butterworth', with thre different orders
% 1,3,6:
optionsButt3.pixxgrade               = pixels_per_degree;
optionsButt3.filterType              = 'butterworth';
optionsButt3.cutoff                  = [3];       
optionsButt3.order                   = 1;
optionsButt3.computeRadialSpectra    = 1;
optionsButt3.toplot                  = 0; %  
[~,lpimageButt3_ord1,filtersButt3_ord1,radialAvgFrqButt3_ord1] =image_filter(square_image,optionsButt3);
figure,imshow(lpimageButt3_ord1)
title(sprintf('%s filter cutoff %d cyc/deg %d order std 0.5 amp. red.',optionsButt3.filterType,optionsButt3.cutoff(1),optionsButt3.order))
pause(1) % for some reason the cript need to pause a bit otherwise the titles go to the wrong figures

optionsButt3.order                   = 5;
[~,lpimageButt3_ord5,filtersButt3_ord5,radialAvgFrqButt3_ord5] =image_filter(square_image,optionsButt3);
figure,imshow(lpimageButt3_ord5)
title(sprintf('%s filter cutoff %d cyc/deg %d order std 0.5 amp. red.',optionsButt3.filterType,optionsButt3.cutoff(1),optionsButt3.order))
pause(1)

optionsButt3.order                   = 10;
[~,lpimageButt3_ord10,filtersButt3_ord10,radialAvgFrqButt3_ord10] =image_filter(square_image,optionsButt3);
figure,imshow(lpimageButt3_ord10)
title(sprintf('%s filter cutoff %d cyc/deg %d order std 0.5 amp. red.',optionsButt3.filterType,optionsButt3.cutoff(1),optionsButt3.order))

%%
% As with the gaussian filter, the butterwoth filter can also be specified
% to have a desired amplitude reduction at a given frequency. For this
% we use the filterType 'butterworth_custom_cutoff'. Let's compare the
% filter above order 5 withone with the same order but a reduction of
% amplitude of 0.1 at the cutoff frequency:

optionsButt3_1.pixxgrade               = pixels_per_degree;
optionsButt3_1.filterType              = 'butterworth_custom_cutoff';
optionsButt3_1.cutoff                  = [3 .1];        % for custom_cutoff the first value is the cutoff requency and the second one the amplitude reduction  
optionsButt3_1.order                   = 5;
optionsButt3_1.computeRadialSpectra    = 1;
optionsButt3_1.toplot                  = 0; %  

% plotting the already fitlered image
fh = figure;
fh.Position(3:4) = fh.Position(3:4)*2; % just making the figure bigger 
subplot('Position',[.05 .55 .4 .4])
imshow(lpimageButt3_ord5)
title(sprintf('%s filter cutoff %d cyc/deg %d order std 0.5 amp. red.',optionsButt3.filterType,optionsButt3.cutoff(1),optionsButt3.order))
pause(1)

[~,lpimageButt3_1_ord5,filtersButt3_1_ord5,radialAvgFrqButt3_1_ord5] =image_filter(square_image,optionsButt3_1);
subplot('Position',[.55 .55 .4 .4])
imshow(lpimageButt3_1_ord5)
title(sprintf('%s filter cutoff %d cyc/deg %d order %1.1f amp. red.',optionsButt3_1.filterType,optionsButt3_1.cutoff(1),optionsButt3_1.order,optionsButt3_1.cutoff(2)))
pause(1)

% this plots the spectra of the two filters
subplot('Position',[.05 .05 .4 .4])
freq_spacing=.1;
xx = 0:1:length(radialAvgFrqButt3_ord5{1})-1;
h(1) = plot(xx,10*log10(radialAvgFrqButt3_ord5{1}.^2),'k','LineWidth',2);   % the original image spectra
hold on
h(2) = plot(xx,10*log10(radialAvgFrqButt3_ord5{1}.^2.*.71.^2),'b:');        % the original image spectra with an amplitude reudction of 0.71
h(3) = plot(xx,10*log10(radialAvgFrqButt3_ord5{1}.^2.*.1.^2),'r:');         % the original image spectra with an amplitude reudction of 0.1
% plot(0:1:length(radialProfile)-1,10*log10(radialProfile.^2.*ar.^2));
h(4) = plot(xx,10*log10(radialAvgFrqButt3_ord5{2}.^2));                     % the image filtered with 'butterworth' filter order 5
h(5) = plot(xx,10*log10(radialAvgFrqButt3_1_ord5{2}.^2));                   % the image filtered with 'butterworth_custom_cutoff' filter order 5
axis([0 8/freq_spacing -75 10])
vline([1.5 3]./freq_spacing)
legend(h,{'Original','0.7*original','0.1*original','Butt0.7%Ord5','Butt0.1%Ord5'},'box','off')
set(gca,'XTick',0:1./freq_spacing:length(radialAvgFrqButt3_ord5{1})-1,'XTickLabel',0:1:(length(radialAvgFrqButt3_ord5{1})-1)*freq_spacing)
xlabel('Spatial Frequency (cyc/deg)')
ylabel('Power (dB)')
title('Example filtering Butterwoth cut-off 3 cyc')
pause(1)

% a cross section of the filters
subplot('Position',[.55 .05 .4 .4])
hh(1) = plot(filtersButt3_ord5(size(filtersButt3_ord5,1)/2+1,:));, hold on
hh(2) = plot(filtersButt3_1_ord5(size(filtersButt3_1_ord5,1)/2+1,:));
legend(hh,{'Butterworth','Butterworth_custom_cutoff'})
title('filters cross-section')
pause(1)

%% CSF filter with im_CSF_filter.m
% The last filtering option is the experimental CSF filter which
% 'filters' an image following a Contrast Sensitive Function
% specification according to the truncated log-parabola model 
% (for details of this CSF model see Lesmes, L. A. (2010). Journal of
% Vision, 10(3), 1–21).
%
% Let's see how the parametrization works, t requires you to install
% our CSF implementation you can find in https://github.com/jpossandon/csf_psi.git)
% The truncated log-parabola model requires 4 parameters: peak frequency
% (p_f), peak sensitivity (gamma), bandwidth(bw) and truncation (delta)

% 'normal human vision' csf parameters
CSFparams.p_f    = 4;
CSFparams.gamma  = 200;
CSFparams.delta  = 1;
CSFparams.bw     = 2;

% To calcualte the contrast sensitive function we require additionally a
% set of spatial frequency and csf.m
spfreqs                         = 2.^[-2:.25:5];
[S]             = csf(CSFparams.p_f ,CSFparams.gamma,CSFparams.delta,CSFparams.bw ,spfreqs);

fh = figure;
fh.Position(3) = fh.Position(3)*2; 
subplot('Position',[.05 .1 .4 .8])
plotCSF(spfreqs,S);
title('CSF ''normal'' vision')
% and 'lower human vision' csf parameters
CSFparams.p_f    = 1;
CSFparams.gamma  = 100;
CSFparams.delta  = .5;
CSFparams.bw     = 1;

% To calcualte the contrast sensitive function we require additionally a
% set of spatial frequency and csf.m

[Slow]             = csf(CSFparams.p_f ,CSFparams.gamma,CSFparams.delta,CSFparams.bw ,spfreqs);
subplot('Position',[.55 .1 .4 .8])
plotCSF(spfreqs,Slow);
title('CSF ''low'' vision')
pause(2)

%%
% Now, let's filter the image according to the 'low vision' CSF. For this
% we need to use the function im_CSF_filter.m and specify the CSF parameters (alrady done above) and an option
% structure (see im_CSF_filter.m documentation for details)
% (note that the function only takes in account contrast sensitivties for
% spatial frequencies above 0.5 cycles per degree)
optionsCSF.pixdegree        = pixels_per_degree;
optionsCSF.padding          = [];
optionsCSF.toplot           = {};
optionsCSF.BPfilter_spacing = 1;
optionsCSF.filtertype       = 'cosine';
[mod_Image] = im_CSF_filter(square_image,CSFparams,optionsCSF);
title(sprintf('Image filtered according to CSF',optionsG3.filterType,optionsG3.cutoff(1)))

fh = figure;
fh.Position(3) = fh.Position(3)*2; 
subplot('Position',[.05 .1 .4 .8])
imshow(mod_Image)
% to check the resulting spectra we use the image_filter.m (this will also
% filter the modified image further with the specifcations in the option
% structure, but we just do not look at this)
[~, ~, ~, radialSpectraCSFlowvis]      = image_filter(mod_Image,optionsButt3_1);

subplot('Position',[.55 .1 .4 .8]), hold on
freq_spacing=.1;
xx = 0:1:length(radialAvgFrqButt3_ord5{1})-1;
hh(1) = plot(xx,10*log10(radialAvgFrqButt3_ord5{1}.^2),'k','LineWidth',2);   % the original image spectra
hh(2) = plot(xx,10*log10(radialSpectraCSFlowvis{1}.^2),'r','LineWidth',2);   % CSF filtereed image spectra
legend(hh,'Original','CSF modified')
axis([0 8/freq_spacing -75 10])
set(gca,'XTick',0:1./freq_spacing:length(radialAvgFrqButt3_ord5{1})-1,'XTickLabel',0:1:(length(radialAvgFrqButt3_ord5{1})-1)*freq_spacing)


%%
close all