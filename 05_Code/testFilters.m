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
finalSize          = [450 450];
 targetSizeInPix     = [400 400];
bkg_color           = 138;
 thisImage          = imresize(imread(imagePath),targetSizeInPix );
image               = uint8(repmat(bkg_color.*ones(finalSize(2),finalSize(1)),1,1,3));
colIndxs            = finalSize(2)./2+1-targetSizeInPix(1)/2:finalSize(2)./2+targetSizeInPix(1)/2;
rowIndxs            = finalSize(1)./2+1-targetSizeInPix(2)/2:finalSize(1)./2+targetSizeInPix(2)/2;;

for rbg = 1:3
    image(colIndxs,rowIndxs,rbg) = thisImage(:,:,rbg);
   
end
options.pixxgrade   = 70;
% fh = figure;
% imshow(image)
% set(fh,'Position',[0 0 finalSize])
% set(gca,'Position',[0 0 1 1])
 %%
options.pixxgrade   = 70;
options.padding    = [];
options.cutoff      = 3;
options.filterType  = 'gaussian';% 'butterworth';%
options.toplot      = 1; %   - 1 to toplot filtered images and the filters
options.computeRadialSpectra      = 1;
[hpimage, lpimageG3, filterG3, radialAvgFrqG3] =image_filter(image,options);
% fh = figure,imshow(lpimageG3)
% set(fh,'Position',[0 0 finalSize])

%%
 options.pixxgrade   = 70;
 options.padding    = [];
options.toplot      = 0; %   - 1 to toplot filtered images and the filters

options.cutoff      = [3 .1];
options.filterType  = 'gaussian_custom_cutoff';% 'butterworth';%
[hpimage, lpimageG3_1, filterG3_1, radialAvgFrqG3_1] =image_filter(image,options);

options.cutoff      = [1.5 .1];
[hpimage, lpimageG1_5_1, filterG1_5_1, radialAvgFrqG1_5_1] =image_filter(image,options);

% 
% fh=figure,imshow(lpimageG3_1)
% set(fh,'Position',[0 0 finalSize])
% %
% set(gca,'Position',[0 0 1 1])
%%
options.cutoff      = [3 .1];
options.order       = 3;
options.filterType  = 'butterworth_custom_cutoff';% 'butterworth';%
[hpimage, lpimageBW3_3, filterBW3_3, radialAvgFrqBW3_3] =image_filter(image,options);

options.order       = 6;
options.filterType  = 'butterworth_custom_cutoff';% 'butterworth';%
[hpimage, lpimageBW3_6, filterBW3_6, radialAvgFrqBW3_6] =image_filter(image,options);



%%
freq_spacing=.1;
xx = 0:1:length(radialAvgFrqG3{1})-1;
figjp,hold on
plot(xx,10*log10(radialAvgFrqG3{1}.^2),'k','LineWidth',2)
plot(xx,10*log10(radialAvgFrqG3{1}.^2.*.61.^2),'b:')
plot(xx,10*log10(radialAvgFrqG3{1}.^2.*.1.^2),'r:')
% plot(0:1:length(radialProfile)-1,10*log10(radialProfile.^2.*ar.^2));
plot(xx,10*log10(radialAvgFrqG3{2}.^2))
plot(xx,10*log10(radialAvgFrqG3_1{2}.^2))
plot(xx,10*log10(radialAvgFrqBW3_3{2}.^2))
plot(xx,10*log10(radialAvgFrqBW3_6{2}.^2))
plot(xx,10*log10(radialAvgFrqG1_5_1{2}.^2))
    axis([0 8/freq_spacing -75 10])
    vline([1.5 3]./freq_spacing)
legend({'Original','0.61*original','0.1*original','Gauss60%','Gauss10%','Butter10%Ord3','Butter10%Ord6','GuassHalfCO10%'},'box','off')
set(gca,'XTick',0:1./freq_spacing:length(radialAvgFrqG3{1})-1,'XTickLabel',0:1:(length(radialAvgFrqG3{1})-1)*freq_spacing)
xlabel('Spatial Frequency (cyc/deg)')
ylabel('Power (dB)')
title('Example filtering cut-off 3 cyc')

%%
fh = figure;
fh.Units = 'centimeter';
fh.Position = [0 0 17.6 17.6*2/3]*3;
nrows = 2;
ncols = 3;
hmargin =.1;
vmargin =.05;
adjust = 0;
subplot('Position',subplotFull(1,1,nrows,ncols,hmargin,vmargin,adjust))
imshow(image)
title('Original')
subplot('Position',subplotFull(1,2,nrows,ncols,hmargin,vmargin,adjust))
imshow(lpimageG3)
title('Gauss. cutoff 3 .6Amp')
subplot('Position',subplotFull(1,3,nrows,ncols,hmargin,vmargin,adjust))
imshow(lpimageG3_1)
title('Gauss. cutoff 3 .1Amp')
subplot('Position',subplotFull(2,1,nrows,ncols,hmargin,vmargin,adjust))
imshow(lpimageBW3_3)
title('BW. cutoff 3 .1Amp Ord. 3')
subplot('Position',subplotFull(2,2,nrows,ncols,hmargin,vmargin,adjust))
imshow(lpimageBW3_6)
title('BW. cutoff 3 .1Amp Ord. 6')
subplot('Position',subplotFull(2,3,nrows,ncols,hmargin,vmargin,adjust))
imshow(lpimageG1_5_1)
title('Gauss. cutoff 1.5 .1Amp')

%%
% image change checks
close all
options.pixxgrade   = 70;
options.padding    = [];
options.cutoff      = [4 .1];
options.filterType  = 'butterworth_custom_cutoff';% 'butterworth';%
options.toplot      = 1; %   - 1 to toplot filtered images and the filters
options.computeRadialSpectra      = 1;
[~, lpimageG3, ~, radialAvgFrqG3] =image_filter(image,options);

imageBigunfilt = imresize(image,[900 900]);
options.toplot      = 0;
[~, lpimageG3big, ~, radialAvgFrqG3big] =image_filter(imageBigunfilt,options);
hold on, plot(10*log10(radialAvgFrqG3big{1}.^2))
hold on, plot(10*log10(radialAvgFrqG3big{2}.^2))

imageBig = imresize(lpimageG3,[900 900]);
options.toplot      = 0;
[~, lpimagelpG3big, ~, radialAvgFrqlpG3big] =image_filter(imageBig,options);
hold on, plot(10*log10(radialAvgFrqlpG3big{1}.^2))

figure,imshow(image), title('original')
figure,imshow(lpimageG3), title('original filt')
figure,imshow(imageBigunfilt),title('original bigger')
figure,imshow(lpimageG3big),title('original bigger filt')
figure,imshow(imageBig),title('original filt bigger')

%%
% image change checks, other way around
close all
options.pixxgrade   = 70;
options.padding    = [];
options.cutoff      = [4 .1];
options.filterType  = 'butterworth_custom_cutoff';% 'butterworth';%
options.toplot      = 1; %   - 1 to toplot filtered images and the filters
options.computeRadialSpectra      = 1;
[~, lpimageG3, ~, radialAvgFrqG3] =image_filter(image,options);

imageBigunfilt = imresize(image,[900 900]);
options.toplot      = 0;
[~, lpimageG3big, ~, radialAvgFrqG3big] =image_filter(imageBigunfilt,options);
hold on, plot(10*log10(radialAvgFrqG3big{1}.^2))
hold on, plot(10*log10(radialAvgFrqG3big{2}.^2))

imageSmall = imresize(lpimageG3big,[450 450]);
options.toplot      = 0;
[~, lpimagelpG3small, ~, radialAvgFrqlpG3small] =image_filter(imageSmall,options);
hold on, plot(10*log10(radialAvgFrqlpG3small{1}.^2))

figure,imshow(image), title('original')
figure,imshow(lpimageG3), title('original filt')
figure,imshow(imageBigunfilt),title('original bigger')
figure,imshow(lpimageG3big),title('original bigger filt')
figure,imshow(imageSmall),title('original bigger filt smaller')
