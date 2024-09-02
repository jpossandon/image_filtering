% imagePath       = '/Users/jossando/Desktop/hamburg.png';
 % image           = double(imread(imagePath))./255;

%  imagePath       = '/Users/jossando/trabajo/FaceEEG/2 experiments scripts/1_FO_loc/stimuli/Grey/face_41.jpg';
% 
%  screenSize      = [1920 1080];
%  image           = .5.*ones(screenSize(2),screenSize(1));
%  image(screenSize(2)/2-200+1:screenSize(2)/2+200,screenSize(1)/2-200+1:screenSize(1)/2+200) = imresize(double(imread(imagePath))./255,2);
%  % image = imresize(double(imread(imagePath))./255,4);
% pixdegree       = 40;

 imagePath          = '/Users/jossando/trabajo/CSF/KDEF_400_final_cut/400px_600_AF31NEHL.bmp';
 screenSize          = [1000 1000];

% imagePath           = '/Users/jossando/Desktop/hamburg.png'
% screenSize          = [1920 1080];
% targetSizeInPix     = [];
targetSizeInPix     = [784 784];

 thisImage          = imread(imagePath);
if ~isempty(targetSizeInPix)
    thisImage          = imresize(thisImage,targetSizeInPix );
end
[m,n,p]             = size(thisImage);
imagen               = uint8(repmat(138.*ones(screenSize(2),screenSize(1)),1,1,3));
colIndxs            = screenSize(2)./2+1-m/2:screenSize(2)./2+m/2;
rowIndxs            = screenSize(1)./2+1-n/2:screenSize(1)./2+n/2;
pixdegree           = 50;
for rbg = 1:3
    imagen(colIndxs,rowIndxs,rbg) = thisImage(:,:,rbg);
   
end
clearvars -except CSDparams pixdegree screenSize imagen
% 'normal' params
CSFparams.p_f    = 5;
CSFparams.gamma  = 200;
CSFparams.delta  = 1;
CSFparams.bw     = 2;

% 'catarct' params
% CSFparams.p_f    = 1.8;
% CSFparams.gamma  = 150;
% CSFparams.delta  = 1;
% CSFparams.bw     = 1.5;

% padding          = 138./255;

%%
optionsCSF.pixdegree        = pixdegree;
optionsCSF.padding          = [];
optionsCSF.toplot           = {'mod_image','filters'};
optionsCSF.BPfilter_spacing = 1;
optionsCSF.filtertype       = 'cosine';
[mod_Image] = im_CSF_filter(imagen,CSFparams,optionsCSF);

% laternative with gaussian filter
%%
optionsFilter.pixxgrade             = pixdegree;
optionsFilter.toplot                = 1; %   - 1 to toplot filtered images and the filters
optionsFilter.computeRadialSpectra  = 1;
optionsFilter.padding               = 138;
optionsFilter.cutoff                = [4 .1];
optionsFilter.order                 = [];
optionsFilter.filterType            = 'gaussian_custom_cutoff';
[~, lpimage, ~, radialSpectra]      = image_filter(uint8(imagen.*255),optionsFilter);
figure, imshow(lpimage)
%
optionsFilter.toplot                = 0; 
[~, ~, ~, radialSpectra]      = image_filter(uint8(mod_Image.*255),optionsFilter);

hold on, plot(10.*log10(radialSpectra{1}.^2))
