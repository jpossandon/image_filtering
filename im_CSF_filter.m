function [mod_image,bpfreqs_cycperdeg] = im_CSF_filter(imagen,CSFparams,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [mod_image,bpfreqs_cycperdeg] = im_CSF_filter(imagen,CSFparams,options)
% Implements the filtering of images according to a contrast sensitive
% function, as proposes by Peli (J. Opt. Soc. Am. A,(1990) 7(10):2032-2040)

[m,n,p]         = size(imagen);

if ~isempty(options.padding)
    PQ          = paddedsize([m n],'PWR2');
    image       = uint8(options.padding.*ones(PQ(1),PQ(2),p));
    padrows     = PQ(1)./2-m/2+1:PQ(1)./2+m/2;
    padcols     = PQ(2)./2-n/2+1:PQ(2)./2+n/2;
    for im_chan = 1:p
        image(padrows,padcols,im_chan) = imagen(:,:,im_chan);
    end
else
    image       = imagen;
    padrows     = 1:m;
    padcols     = 1:n;
end
for im_chan = 1:p
    F{im_chan}      = fft2(image(:,:,im_chan));
    imDC(im_chan)   = mean(mean(image(:,:,im_chan)));
end
[mm,nn,p]           = size(image);
degperim            = mm/options.pixdegree;

%%
% %%%%%%%%%%%%%%%%%%%%%%%%
% % log2 filter bank
% %the filters (filters are 2D matrices to be multiplied with the frequency
% % domain representation of the image (each channel)
% %%%%%%%%%%%%%%%%%%%%%%%%
exp_space           = options.BPfilter_spacing;
max_central_exp     = floor(log2(min(mm/2,nn/2)));          % the max frequency possible (in cycles per image)
exp_i               = 0:exp_space:max_central_exp;          % centre frequency exponents for center frequency octaves (2.^exp_i)
bpfreqs_cycperdeg   = 2.^exp_i./degperim;

% this is needed to create the filter according to the position in the frequency domain in polar coordinates
[X,Y]           = meshgrid(-nn/2:1:nn/2-1,-mm/2:1:mm/2-1);
[TH,R]          = cart2pol(X,Y);
clear Y TH

filters         = [];
if strcmp(options.filtertype,'cosine')
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Cosine log filter bank
    %%%%%%%%%%%%%%%%%%%%%%%%
    for it = 1:length(exp_i)
        filter      = .5*(1+cos(pi*log2(R)/exp_space-pi*exp_i(it)./exp_space));              % the cosine log filter
        filter(R<2.^(exp_i(it)-exp_space) | R>2.^(exp_i(it)+exp_space)) = 0;                % wee keepn only one octave around the center frequency
        filters     = cat(3,filters,filter);
    end
elseif strcmp(options.filtertype,'gaussian')
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Gaussian log filter bank
    %%%%%%%%%%%%%%%%%%%%%%%%
    for it = 1:length(exp_i)
        sd          = exp_space./ (2.*sqrt(2*(log(2)))); % FWHM is set to exp_space
        mean1       = exp_i(it);
        filter      = exp(-((log2(R)-mean1).^2)./2./sd.^2);
        filters     = cat(3,filters,filter);
    end
 end
clear filter R
if ismember('filters',options.toplot)
    figure, hold on
    plot([log2(X(mm/2+1,nn/2+2:end))],squeeze(filters(mm/2+1,nn/2+2:end,1:end)),'.-')
    plot([log2(X(mm/2+1,nn/2+2:end))],sum(squeeze(filters(mm/2+1,nn/2+2:end,1:end)),2),'.-')
    set(gca,'XTick',[0:exp_space:max_central_exp],'XTickLabels',2.^[0:exp_space:max_central_exp])
    xlabel('Center Frequency (cyc/image)')
    ylabel('Magnitude')

    figure, hold on
    plot([X(mm/2+1,nn/2+2:end)],squeeze(filters(mm/2+1,nn/2+2:end,1:end)),'.-')
    plot([X(mm/2+1,nn/2+2:end)],sum(squeeze(filters(mm/2+1,nn/2+2:end,1:end)),2),'.-')
    % set(gca,'XTick',[0:exp_space:max_central_exp],'XTickLabels',[0:exp_space:max_central_exp])
    xlabel('Center Frequency (cyc/image)')
    ylabel('Magnitude')
end
clear X
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creation of bandpass and contrast images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contrast_images  = {};
bandpass_images  = {};
for im_chan = 1:p
    bp_image  = [];
    for it = 1:length(exp_i)
        % here we first center the image spectrum, then multiply it with the cosine torus shape filter,
        % shift it back and inverse forueir transform it
        bp_image = cat(3,bp_image,ifft2(ifftshift(fftshift(F{im_chan}).*filters(:,:,it))));
        % contrast images are created by diving the bandpass image by the
        % low pass fitlered version at the same scale ( which is equalt to
        % the sum of all bandpass images below this level plus the DC
        % component
        if it>1
            thisContIm = bp_image(:,:,it)./(imDC(im_chan)+sum(bp_image(:,:,1:it-1),3)); % dont forget DC checl this and extreme values
            contrast_images{im_chan}(:,:,it-1) = thisContIm(padrows,padcols);
        end
    end
    bandpass_images{im_chan} = bp_image(padrows,padcols,:);
end
clear bp_image F filters  thisContIm

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modifying image accoriding to a CSF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
low_limit       = .5;                               % in cycle per degree, images are not odified according to frequencies below this number, set above 0.5 c/degree to avoid ringing artifact of the lwoer band pass filters
[S]             = csf(CSFparams.p_f ,CSFparams.gamma,CSFparams.delta,CSFparams.bw ,bpfreqs_cycperdeg);
S(bpfreqs_cycperdeg<low_limit) = log10(10000);      % not realy necesay since these frequencioes are not evaluated
Sens            = 10.^S;
Sens(Sens<.1)   = .1;
contS           = 1./Sens;

mod_image       = uint8([]);
thresholdedBP   = {};
for im_chan = 1:p
    chan_mod_image  = [];
    for it = 1:length(exp_i)
        if it ==1
            chan_mod_image = bandpass_images{im_chan}(:,:,1)+imDC(im_chan);
        else
            thisBP = bandpass_images{im_chan}(:,:,it);
            if 2.^exp_i(it)/degperim>low_limit % TODO only do above some cyclesperdegree
                thisBP(abs(contrast_images{im_chan}(:,:,it-1))<contS(it)) = 0;
            end
            chan_mod_image = chan_mod_image+thisBP;
            thresholdedBP{im_chan}(:,:,it-1) = thisBP;
        end
    end
    mod_image(:,:,im_chan) = uint8(chan_mod_image);
end
if ismember('mod_image',options.toplot)
    figure,imshow(mod_image)
end


%%
if ismember('filtering_process',options.toplot)
    nrows = 5;
    ncols = ceil(length(exp_i)/2);
    hmargin = .05;
    vmargin = 0.05;
    adj     = 0;
    fh = figure;
    % fh.Units = 'centimeters'
    fh.Position = [1 1 1700 1700*m/n/ncols*(nrows+1)];

    subplot('Position',subplotFull(1,1,nrows,ncols,hmargin,vmargin,adj))
    imshow(imagen)
    ylabel('Original')
    itt = 1;
    for it=1:2:length(exp_i)
        subplot('Position',subplotFull(2,itt,nrows,ncols,hmargin,vmargin,adj))
        if p>1
            imshow(cat(3,imDC(1)+bandpass_images{1}(:,:,it),imDC(2)+bandpass_images{2}(:,:,it),imDC(3)+bandpass_images{3}(:,:,it)))
        else
            imshow(imDC(1)+bandpass_images{1}(:,:,it))
        end
        title(sprintf('%1.2f / %2.2f',2.^exp_i(it),2.^exp_i(it)./degperim))
        if itt==1, ylabel('Bandpass'),end
        if itt>1
            subplot('Position',subplotFull(3,itt,nrows,ncols,hmargin,vmargin,adj))
            if p>1
                imshow(.5+cat(3,contrast_images{1}(:,:,it-1),contrast_images{2}(:,:,it-1),contrast_images{3}(:,:,it-1)))
            else
                imshow(.5+contrast_images{1}(:,:,it-1))
            end
            if itt==2, ylabel('Contrast Im'),end

            subplot('Position',subplotFull(4,itt,nrows,ncols,hmargin,vmargin,adj))
            if p>1
                imshow(.5+cat(3,thresholdedBP{1}(:,:,it-1),thresholdedBP{2}(:,:,it-1),thresholdedBP{3}(:,:,it-1)))
            else
                imshow(.5+thresholdedBP{1}(:,:,it-1))
            end
            if itt==2, ylabel('Bandpass >  CSF'),end
        end
        itt = itt+1;
    end
    subplot('Position',subplotFull(5,1,nrows,ncols,hmargin,vmargin,adj))
    imshow(mod_image)
    ylabel('Modified')
end

   