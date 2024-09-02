function [hpimage, lpimage, filter, radialAvgFrq] = image_filter(imagen,options);%pixxgrade,cutoff,filterType,toplot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% jpimfourier, low and high image filtering with a gaussian filter
% characterized by its standard deviation stdgaus (in cycles per visual
% degree) that correspond to the filter frequency cutoff (0.607 reduction)
%   Inputs
%               imagen                  - either the path and name of animage archive of the name of an unit8 image variable
%                                           in the workspace '
%               options.pixxgrade       - pixels per visual degree, accordint to screen resolution and observer distance
%               options.cutoff          - for filterType 'gaussian' the cutoff is one number in cycles per degree resulting in a reduction of ampitude of ~0.607 and
%                                           equalt to the standard deviation of the gaussain filter in cycles/degree
%                                         for filterType 'gaussian_custom_cutoff' two numbers are needed [fc ar], the frequency
%                                         cutoff(fc) and the desired amplitude reduction(ar)
%               options.order           - for butterworth filter, lower order are
%                                       similar to gaussian, higher order produce
%                                       sharper filter at the expense of more
%                                       ringing
%               options.filterType      - 'gaussian', 'gaussian_custom_cutoff' or 'butterworth'
%               options.padding         - if set to empty, images are not padded, this only work well 
%                                       for guassian filters or for butterwoth filter with square images. 
%                                       - if set to a number (0-255) the
%                                       images are symmetrically padded
%                                       with the value to a square size equalt to the next power of
%                                       2 of the larger image size
%               options.toplot          - 1 to toplot filtered images and the filters
%
%  Outputs
%               hpimage, lpimage: corresponding images filtered

% 23/07/07, JPO, from Gonzales RC. et al. -Digital Image Processing
% Using MATLAB (2004). Usar bajo su exclusiva responabilidad
% 02/2024, english documantation
% 08/2024, changed to a more gneral function using two different ttypes of
% filters and cutoff
% 8/2024, added padding so butterworth filter can be used in non-suqare
% images, fixed the radially averaged spectra problems.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD padding
if ~isnumeric(imagen)
    imagen=imread(imagen);            % abre imagen si es el nombre de una imagen
end
if strcmp(options.filterType,'gaussian')
    cutoff      = options.cutoff;
    ar          = exp(-.5);
elseif strcmp(options.filterType,'gaussian_custom_cutoff')
    cutoff      = options.cutoff(1)./(sqrt(2.*log(1./options.cutoff(2))));
    ar          = options.cutoff(2);
end

pixxgrade   = options.pixxgrade;
if strcmp(options.filterType,'butterworth') | strcmp(options.filterType,'butterworth_custom_cutoff')
    if isfield(options,'order')
        order_n = options.order;
    else
        order_n = 3;
        warning('Butterwoth filter order was not specified so it is set to 3')
    end
end
if strcmp(options.filterType,'butterworth')
    cutoff      = options.cutoff;
    ar          = 1/sqrt(2); % butterworth cutoff is ~3dB
elseif strcmp(options.filterType,'butterworth_custom_cutoff')
    cutoff      = options.cutoff(1);
    ar          = options.cutoff(2); 
end

[m,n,p]=size(imagen);

for depth=1:p
    % Transformacion de imagen al dominio de la frecuencia
    if options.padding 
        PQ = paddedsize([m n],'PWR2'); % f is floating point.
        thisIm = options.padding.*ones(PQ(1),PQ(2));
        thisIm(PQ(1)./2-m/2+1:PQ(1)./2+m/2,PQ(2)./2-n/2+1:PQ(2)./2+n/2) = imagen(:,:,depth);
        F  = fft2(single(thisIm)); 
    else
         F  = fft2(single(imagen(:,:,depth))); 
    end
    [mm,nn] = size(F);
    Fc = fftshift(F);                       % pone le frecuencia 0 al centro
    S2 = log(1+abs(Fc));                    % poder espectral en escala logaritmica
    
    % Creacion de los filtros
    DOU = cutoff*(mm/pixxgrade);            % frequencies to samples (multiplity for N/Fs)
    DOV = cutoff*(nn/pixxgrade);
    [U, V] = dftuv(mm,nn);
     D = hypot(U,V);
    if strcmp(options.filterType,'gaussian') | strcmp(options.filterType,'gaussian_custom_cutoff')
        % gaussian
        Hlp = exp(-0.5*((U.^2)./(DOU^2)+(V.^2)./(DOV^2)));
           % Hlp = exp(-(D.^2)./(2*(D0^2)));

    elseif strcmp(options.filterType,'butterworth')
        % butterworth
        % this works correctly only for an square images so use padding
        % otherwise
         % % determining the filtering mask
          Hlp = sqrt(1./(1 + (D./DOU).^(2*order_n)));
 
    elseif strcmp(options.filterType,'butterworth_custom_cutoff')
         DOU = DOU./power(-1+options.cutoff(2).^(-2),1./2./order_n); % adjustment of the cutoff (where power is reduced 0.5 and amplitude 0.708) so the desired reduction ar occur at the desired cutoff
         Hlp = sqrt(1./(1 + (D./DOU).^(2*order_n)));
    end
    Hhp = 1-Hlp;
    
    lpFiltered         =(real(ifft2(Hlp.*F)));
   hpFiltered         =(real(ifft2(Hhp.*F)));
     if options.padding
        lpFiltered         = lpFiltered(PQ(1)./2-m/2+1:PQ(1)./2+m/2,PQ(2)./2-n/2+1:PQ(2)./2+n/2);
         hpFiltered         = hpFiltered(PQ(1)./2-m/2+1:PQ(1)./2+m/2,PQ(2)./2-n/2+1:PQ(2)./2+n/2);
       end
    
    hpimage(:,:,depth) = hpFiltered;
    lpimage(:,:,depth) = lpFiltered;
  
    hpimage(:,:,depth) = (hpimage(:,:,depth)+mean2(imagen (:,:,depth)));
end
hpimage = uint8(hpimage);
lpimage = uint8(lpimage);
filter  = fftshift(Hlp);
%%
% generate redial averaged spectra
freq_spacing = .1; % this is arbitrary, small enough to have good resoltuion, bigh enoug to have good smoothing and to be over the minimal frequency resolution
if options.computeRadialSpectra == 1
    for ee = 1:3
        if ee ==1
            if p==3
                amplOrig = abs(fftshift(fft2(double(rgb2gray(imagen)))));
            else
                amplOrig = abs(fftshift(fft2(double(imagen))));
            end
        elseif ee==2
            if p==3
                amplOrig = abs(fftshift(fft2(double(rgb2gray(lpimage)))));
            else
                amplOrig = abs(fftshift(fft2(double(lpimage))));
            end
        elseif ee==3
            if p==3
                amplOrig = abs(fftshift(fft2(double(rgb2gray(hpimage)))));
            else
                amplOrig = abs(fftshift(fft2(double(hpimage))));
            end
        end
        % amplOrig = amplOrig.^2;
        % Get the average radial profile
        midRow = m/2+1;
        midCol = n/2+1;
        shorterSiz = min(m/2,n/2);
        % maxRadius = ceil(sqrt((shorterSiz./pixxgrade)^2 + (shorterSiz./pixxgrade)^2));
        maxRadius = ceil(sqrt((midRow./m*pixxgrade)^2 + (midCol./n.*pixxgrade)^2));   % in degrees to take care of non square images
        freqSize          = maxRadius./freq_spacing;%.*pixxgrade;
        radialProfile = zeros(freqSize+1, 1);  % in cycles per image
        count = zeros(freqSize+1, 1);
        for col = midCol-n/2 : midCol+n/2-1
         for row = midRow-m/2 : midRow+m/2-1
                radius = sqrt(((row - midRow)./m.*pixxgrade) ^ 2 + ((col - midCol)./n.*pixxgrade) ^ 2); % radius in visual degrees
                thisIndex = ceil(radius./freq_spacing) + 1;
                radialProfile(thisIndex) = radialProfile(thisIndex) + amplOrig(row, col);
                count(thisIndex) = count(thisIndex) + 1;
            end
        end
        % Get average
        ixr = radialProfile==0;
        radialProfile = radialProfile ./ count ./ m./n;  % note that by dividing by count this means that the resultion ampokitude is the averages across orientaitons so for a test pattern in one direction the resulting amplitude will be expected./count
        radialProfile(ixr) = 0.000000001;
        radialAvgFrq{ee}   = radialProfile;
    end
end


% Plot
if options.toplot==1

    % figure
    % subplot(2,2,1)
    % imshow(imagen)
    % title('Original','FontWeight','Bold')
    % subplot(2,2,2)
    % imshow(S2,[])
    % title('Power Spectrum(log)','FontWeight','Bold')
    % subplot(2,2,3)
    % imshow(fftshift(Hlp),[])
    % title(sprintf('LP Gaussian Filter std=%2.1f cycles/degree',cutoff),'FontWeight','Bold')
    % subplot(2,2,4)
    % imshow(fftshift(Hhp),[])
    % title(sprintf('HP Gaussian Filter std=%2.1f cycles/degree',cutoff),'FontWeight','Bold')
    % figure
    % imshow(hpimage)
    % title('High Pass Filtered Image','FontWeight','Bold')
    % figure
    % imshow(lpimage)
    % title('Low Pass Filtered Image','FontWeight','Bold')

    figure, hold on
    for ee = 1:2
        radialProfile = radialAvgFrq{ee};
        plot(0:1:length(radialProfile)-1,10*log10(radialProfile.^2)); % dB as a ratio of powers
        if ee==1
                % for 'gaussian' change of -4.3 dB in power equivalent to a reduction
                % of amplitude of 0.607 for standard cutoff equal to the sd of
                % the gausian.
                % for 'gaussian_custom_cutoff' change of x dB in power equivalent to a reduction
                % of amplitude of ar for a cutoff equal to
                % sqrt(2.*log(1./ar))*gaussianSD.

                % for 'butterworth' change of ~-3.1 dB in power equivalent to a reduction
                % of amplitude of 0.708 for teh filter standard cutoff
                % for 'butterworth_custom_cutoff' change of x dB in power equivalent to a reduction
                % of amplitude of ar for a cutoff that is adjusted by DOU./power(-1+options.cutoff(2).^(-2),1./2./order_n);
                plot(0:1:length(radialProfile)-1,10*log10(radialProfile.^2.*ar.^2));
        end
         vline(options.cutoff(1)./freq_spacing)
        text(options.cutoff(1)./freq_spacing+.25./freq_spacing,0,'Cut-off','Color',[1 0 0])
    end
    if ee==2
        legend({'Original',sprintf('%1.2f*Original/%1.1fdB Power',ar,10*log10(ar.^2)),'Low-pass'})
    elseif ee==3
        legend({'Original','0.607 Original','Low-pass','high-pass'})
    end
    axis([0 10/freq_spacing -100 20])

    xlabel('Spatial Frequency (cycles/degree)')
    ylabel('Power (dB)')
    if ismember(options.filterType,{'butterworth','butterworth_custom_cutoff'})
        title(sprintf('Radially avg. spectrum %s order %d',options.filterType,options.order))
    else
        title(sprintf('Radially avg. spectrum %s',options.filterType))
    end
    % set(gca,'XTick',0:m/pixxgrade:length(radialProfile)-1,'XTickLabel',0:1:(length(radialProfile)-1)*pixxgrade/m) % in the plot integer frequencies are at mutltiples of N/Fx the number of samples of frequencies 1,2,3,
    set(gca,'XTick',0:1./freq_spacing:length(radialProfile)-1,'XTickLabel',0:1:(length(radialProfile)-1)*freq_spacing)
end


%%
% checks radial averagedspectrum of the grayscale image
% figure,
% hold on


%%
% CHECK, this works
% pixdegree = 40;
% degperim  = 48;
% x = 0:1./pixdegree:degperim-1./pixdegree;
% ima = repmat(0.5*(1+sin(2*pi*3*x)),m,1);
% figure,
% imshow(ima)
% set(gca, 'XTick',0:40:1920,'YTick',0:pixdegree:1080)
% axis on, grid on
% 
% [hpimage, lpimage]=jpimfourier(uint8(ima.*255),pixdegree,3,0);
% figure, imshow(double(lpimage)./255)
% set(gca, 'XTick',0:40:1920,'YTick',0:pixdegree:1080)
% axis on, grid on
% (max(lpimage(:))-min(lpimage(:)))/255 % this should be  ~ 0.607