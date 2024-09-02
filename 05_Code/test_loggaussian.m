x = 0:.01:30;
sd = 2;
fwhm =  2.*sqrt(2*(log(2)))*sd;
figjp, hold on
mean1      = 15;
gaussian1  = exp(-((x-mean1).^2)./2/sd.^2);
plot(x,gaussian1,'.-r')
vline([mean1 mean1-fwhm/2 mean1+fwhm/2])
hline(.5)
mean2 = mean1+fwhm;
gaussian2  = exp(-((x-mean2).^2)./2/sd.^2);
plot(x,gaussian2,'.-b')
mean3 = mean2+fwhm;
gaussian3  = exp(-((x-mean3).^2)./2/sd.^2);
plot(x,gaussian3,'.-g')
mean4 = mean3+fwhm;
gaussian4  = exp(-((x-mean4).^2)./2/sd.^2);
plot(x,gaussian4,'.-g')
plot(x,gaussian1+gaussian2+gaussian3+gaussian4)
%%
x           = 0:.01:50;
sd          = 2%1.5;
fwhm        =  2.*sqrt(2*(log(2)))*sd;
figjp, hold on
mean1       = 8;
gaussian1   = exp(-((log2(x)-log2(mean1)).^2)./2/log2(sd).^2);
plot(log2(x),gaussian1,'.-r')

%%
    % sd          = 1.5;
    fwhm        = 1;% 2.*sqrt(2*(log(2)))*sd;
    sd          = 2.^(fwhm./ (2.*sqrt(2*(log(2)))));
    mean1       = 3%exp_i(it);
    filter   = exp(-((log2(x)-mean1).^2)./2/log2(sd).^2);

    figjp,plot(log2(x),filter,'.-r')
