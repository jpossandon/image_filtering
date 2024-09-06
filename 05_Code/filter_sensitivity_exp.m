%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% image sensitivity experiment

beep OFF
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paths, experiment parameter and inputs

if ismac
     paths.parent    = fullfile(filesep,'Users','jossando','trabajo','image_filtering');
else
    paths.parent    = fullfile('D:','Jose Ossandon','image_filtering');
    addpath('D:\Jose Ossandon\csf_psi')
    addpath('D:\Jose Ossandon\image_filtering\05_Code\external')
end
paths.imDIR     = fullfile(paths.parent,'KDEF_400_final_cut');
paths.result    = fullfile(paths.parent, '06_RawData'); 

answer         = inputdlg({'Participant ID:','Participant decimal acuity:','CSF peak frequency:','CSF gamma:','CSF delta:','CSF bw:','Screen distance:'},'Filter Sens.');
fileMAT        = sprintf('FS%03d.mat',str2num(answer{1}));

if isempty(answer{3})

    PARAMS.conditions      = struct( 'name',{'same','gaussian','gaussian_01cutoffhalf','gaussian_01cutoff','butterworth_custom_01cutoff_order3','butterworth_custom_01cutoff_order6'},...
                                'type',{'same','gaussian','gaussian_custom_cutoff','gaussian_custom_cutoff','butterworth_custom_cutoff','butterworth_custom_cutoff'},...
                                     'cutoff',{NaN,NaN,0.1,0.1,0.1,0.1},...
                                     'order',{NaN,NaN,NaN,NaN,3,6});
else
PARAMS.conditions      = struct( 'name',{'same','gaussian','gaussian_01cutoffhalf','gaussian_01cutoff','butterworth_custom_01cutoff_order3','butterworth_custom_01cutoff_order6','CSF_filter'},...
                                'type',{'same','gaussian','gaussian_custom_cutoff','gaussian_custom_cutoff','butterworth_custom_cutoff','butterworth_custom_cutoff','CSF_filter'},...
                                     'cutoff',{NaN,NaN,0.1,0.1,0.1,0.1,NaN},...
                                     'order',{NaN,NaN,NaN,NaN,3,6,NaN});
end
% my screen
% PARAMS.setup.scr_wdth               = 51.7; % centimeters
% PARAMS.setup.scr_hght               = 51.7./1.6;
% screen RAUM 305
PARAMS.setup.scr_wdth               = 52.7; % centimeters
PARAMS.setup.scr_hght               = 29.5;


PARAMS.setup.Dist_from_screen       = str2double(answer(7));
PARAMS.setup.centerdotsize          = 15;

PARAMS.images.targetStimSize        = [21 21]; % width,height in cms
PARAMS.images.backgroundColor       = 138;



PARAMS.subject.CSFparams.p_f        = str2double(answer(3));
PARAMS.subject.CSFparams.gamma      = str2double(answer(4));
PARAMS.subject.CSFparams.delta      = str2double(answer(5));
PARAMS.subject.CSFparams.bw         = str2double(answer(6));
if ~isempty(answer{2})
    PARAMS.subject.visualAcuity     = str2double(answer(2));
end
    PARAMS.subject.spfreqs                         = .1:.05:50;
    [PARAMS.subject.S]                             = csf(PARAMS.subject.CSFparams.p_f ,PARAMS.subject.CSFparams.gamma,PARAMS.subject.CSFparams.delta,PARAMS.subject.CSFparams.bw ,PARAMS.subject.spfreqs);
if isempty(answer{2})
    PARAMS.subject.visualAcuity     = PARAMS.subject.spfreqs(find(10.^PARAMS.subject.S<1,1,'first'))/30;
end

PARAMS.subject.gratingCutoff        = PARAMS.subject.visualAcuity*30;
 
 % answer         = inputdlg({'Participant ID:'},'Filt Sens.');
 % fileMAT        = sprintf('FS%03d.mat',str2num(answer{1}));

%%

KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE');
Screen('Preference', 'SkipSyncTests', 1);
Screen('Preference', 'ConserveVRAM',  4096)
% ListenChar(2)
AssertOpenGL;

screens         = Screen('Screens');
screenNumber    = max(screens);
bkg_gray        = PARAMS.images.backgroundColor;

% open full window
[win, windowRect] = PsychImaging('OpenWindow', screenNumber,bkg_gray);
[PARAMS.setup.scr_pixwdth, PARAMS.setup.scr_pixhgt] = Screen('WindowSize', win);
PARAMS.setup.screen_ratio = PARAMS.setup.scr_pixwdth/PARAMS.setup.scr_pixhgt;
Screen('TextSize', win, 30);
rectOval = [PARAMS.setup.scr_pixwdth/2-PARAMS.setup.centerdotsize/2,...
            PARAMS.setup.scr_pixhgt/2-PARAMS.setup.centerdotsize/2,...
            PARAMS.setup.scr_pixwdth/2+PARAMS.setup.centerdotsize/2,...
            PARAMS.setup.scr_pixhgt/2+PARAMS.setup.centerdotsize/2];

% HideCursor(win);
clear TRIALS

%%
% trial randomization
PARAMS.TRIALS.trialsPerCondition    = 30;
PARAMS.TRIALS.trial_per_block       = 15;
PARAMS.TRIALS.totalTrials           = PARAMS.TRIALS.trialsPerCondition.*length(PARAMS.conditions);
PARAMS.TRIALS.imageDuration         = .25; %in seconds
PARAMS.TRIALS.ISI                   = .1; %period between the two images, in seconds
PARAMS.TRIALS.centerJitter          = 20; %jitter around the center for image presentation
TRIALS.trialCondition               = repmat({PARAMS.conditions.name},1,PARAMS.TRIALS.trialsPerCondition);
TRIALS.trialConditionFilterType     = repmat({PARAMS.conditions.type},1,PARAMS.TRIALS.trialsPerCondition);
permOrder                           = randperm(length(TRIALS.trialCondition));
TRIALS.trialCondition               = TRIALS.trialCondition(permOrder); %random order
TRIALS.trialConditionFilterType     = TRIALS.trialConditionFilterType(permOrder); %random order

%%
% calcualte image target dimension
PARAMS.setup.pixels_per_degree  = round(PARAMS.setup.scr_pixwdth/(2*atan(PARAMS.setup.scr_wdth/2/PARAMS.setup.Dist_from_screen)*180/pi)); 
% PARAMS.images.targetSizeInPix   = round(PARAMS.images.targetStimSize.* PARAMS.setup.pixels_per_degree);
PARAMS.images.targetSizeInPix   = round(PARAMS.images.targetStimSize.* PARAMS.setup.scr_pixwdth./PARAMS.setup.scr_wdth);
PARAMS.images.targetSizeInPix   = PARAMS.images.targetSizeInPix+mod(PARAMS.images.targetSizeInPix,2); % get evensize images
PARAMS.images.filtSize           = [1000 1000];
% images
imagesNames = dir(paths.imDIR);
imagesNames([imagesNames.isdir] | ismember({imagesNames.name},{'Thumbs.db'}) |  startsWith({imagesNames.name},'._')) = [];

TRIALS.imagesList           = {imagesNames(randi(length(imagesNames),1,PARAMS.TRIALS.totalTrials)).name};
TRIALS.imageOrder           = [TRIALS.trialCondition; repmat({'same'},1,length(TRIALS.trialCondition))];
for ee = 1:size(TRIALS.imageOrder,2),if rand(1)>.5, TRIALS.imageOrder(:,ee) = flipud(TRIALS.imageOrder(:,ee));,end,end
TRIALS.xjitter     = [round(PARAMS.TRIALS.centerJitter.*(rand(1,PARAMS.TRIALS.totalTrials)-.5));round(PARAMS.TRIALS.centerJitter.*(rand(1,PARAMS.TRIALS.totalTrials)-.5))];
TRIALS.yjitter     = [round(PARAMS.TRIALS.centerJitter.*(rand(1,PARAMS.TRIALS.totalTrials)-.5));round(PARAMS.TRIALS.centerJitter.*(rand(1,PARAMS.TRIALS.totalTrials)-.5))];
colIndxs1         = PARAMS.images.filtSize(1)./2+1-PARAMS.images.targetSizeInPix(1)/2:PARAMS.images.filtSize(1)./2+PARAMS.images.targetSizeInPix(1)/2;
rowIndxs1         = PARAMS.images.filtSize(2)./2+1-PARAMS.images.targetSizeInPix(2)/2:PARAMS.images.filtSize(2)./2+PARAMS.images.targetSizeInPix(2)/2;
colIndxs         = PARAMS.setup.scr_pixhgt./2+1-PARAMS.images.filtSize(1)/2:PARAMS.setup.scr_pixhgt./2+PARAMS.images.filtSize(1)/2;
rowIndxs         = PARAMS.setup.scr_pixwdth./2+1-PARAMS.images.filtSize(2)/2:PARAMS.setup.scr_pixwdth./2+PARAMS.images.filtSize(2)/2;

%colIndxs         = PARAMS.setup.scr_pixhgt./2+1-PARAMS.images.targetSizeInPix(1)/2:PARAMS.setup.scr_pixhgt./2+PARAMS.images.targetSizeInPix(1)/2;
%rowIndxs         = PARAMS.setup.scr_pixwdth./2+1-PARAMS.images.targetSizeInPix(2)/2:PARAMS.setup.scr_pixwdth./2+PARAMS.images.targetSizeInPix(2)/2;
 
%%
% setup common filter options
% gaussian/buttherworth filters
optionsFilter.pixxgrade             = PARAMS.setup.pixels_per_degree;
optionsFilter.toplot                = 0; %   - 1 to toplot filtered images and the filters
optionsFilter.computeRadialSpectra  = 0;
optionsFilter.padding               = [];%138;%[];
% CSF filters
optionsCSF.pixdegree                = PARAMS.setup.pixels_per_degree;
optionsCSF.padding                  = [];%138;
optionsCSF.toplot                   = {''};
optionsCSF.BPfilter_spacing         = 1;%.5;
optionsCSF.filtertype               = 'cosine';
%%

Screen('FillRect',win,bkg_gray); 
DrawFormattedText(win, 'Image filtering sensitivity experiment\n\n (Press any key to start)','center','center');
Screen('Flip', win);
KbWait([], 2);

flagscape =0;

for trl = 1:length(TRIALS.trialCondition)
    % preload images of this block
    if rem(trl,PARAMS.TRIALS.trial_per_block)==1
         if trl>1
            Screen('FillRect',win,bkg_gray); 
            DrawFormattedText(win, sprintf('Block %d of %d finished\n\n (Press any key to continue)',floor(trl/PARAMS.TRIALS.trial_per_block),ceil(length(TRIALS.trialCondition)/PARAMS.TRIALS.trial_per_block)),'center','center');
            Screen('Flip', win);
            WaitSecs(.1)
            KbWait([], 2);
        end
        DrawFormattedText(win, 'Wait for images upload ...\n\n ','center','center');
        Screen('Flip', win);
        theseImages = {};
        for ims = 0:PARAMS.TRIALS.trial_per_block-1
            thisImage        = imresize(imread(fullfile(paths.imDIR,TRIALS.imagesList{trl+ims})),PARAMS.images.targetSizeInPix);
            thisImage1       = uint8(repmat(PARAMS.images.backgroundColor.*ones(PARAMS.setup.scr_pixhgt,PARAMS.setup.scr_pixwdth),1,1,3));
            thisImage2       = thisImage1;
            thisImage11       = uint8(repmat(PARAMS.images.backgroundColor.*ones(PARAMS.images.filtSize),1,1,3));
            thisImage22       = thisImage11;
            for rbg = 1:3
                thisImage11(colIndxs1,rowIndxs1,rbg) = thisImage(:,:,rbg);
                thisImage22(colIndxs1,rowIndxs1,rbg) = thisImage(:,:,rbg);
            end
            
            if ~strcmp(TRIALS.trialCondition{trl+ims},'same')
                
                 optionsFilter.filterType      = TRIALS.trialConditionFilterType{trl+ims};%
                 iXfilt                        = find(ismember({PARAMS.conditions.name},TRIALS.trialCondition{trl+ims}));
                if ismember(optionsFilter.filterType,{'gaussian','gaussian_custom_cutoff','butterworth','butterworth_custom_cutoff'})
                    if strcmp(optionsFilter.filterType,'gaussian_custom_cutoff')
                    % setup this trial filter options
                        if strcmp(TRIALS.trialCondition{trl+ims},'gaussian_01cutoffhalf')
                             optionsFilter.cutoff    = [PARAMS.subject.gratingCutoff./2 PARAMS.conditions(iXfilt).cutoff];
                        else
                            optionsFilter.cutoff    = [PARAMS.subject.gratingCutoff PARAMS.conditions(iXfilt).cutoff];
                        end
                        optionsFilter.order     = [];
                    elseif strcmp(optionsFilter.filterType,'gaussian')
                    % setup this trial filter options
                        optionsFilter.cutoff    = [PARAMS.subject.gratingCutoff];
                        optionsFilter.order     = [];
                    elseif strcmp(optionsFilter.filterType,'butterworth')
                        optionsFilter.cutoff    = PARAMS.subject.gratingCutoff;
                        optionsFilter.order     = PARAMS.conditions(iXfilt).order;
                    elseif strcmp(optionsFilter.filterType,'butterworth_custom_cutoff')
                        optionsFilter.cutoff    = [PARAMS.subject.gratingCutoff PARAMS.conditions(iXfilt).cutoff];
                        optionsFilter.order     = PARAMS.conditions(iXfilt).order;
                    end
                    
                    if strcmp(TRIALS.imageOrder{1,trl+ims},'same')
                        [~, thisImage22] = image_filter(thisImage22,optionsFilter);
                    else
                        [~, thisImage11] = image_filter(thisImage11,optionsFilter);
                    end
                elseif ismember(optionsFilter.filterType,{'CSF_filter'})
                    if strcmp(TRIALS.imageOrder{1,trl+ims},'same')
                        [thisImage22] = im_CSF_filter(thisImage22,PARAMS.subject.CSFparams,optionsCSF);
                    else
                        [thisImage11] = im_CSF_filter(thisImage11,PARAMS.subject.CSFparams,optionsCSF);
                    end
                end
                    
            end
            for rbg = 1:3
                thisImage1(colIndxs-TRIALS.xjitter(1,trl+ims),rowIndxs-TRIALS.yjitter(1,trl+ims),rbg) = thisImage11(:,:,rbg);
                thisImage2(colIndxs-TRIALS.xjitter(2,trl+ims),rowIndxs-TRIALS.yjitter(2,trl+ims),rbg) = thisImage22(:,:,rbg);
            end
            theseImages{1,trl+ims} = thisImage1;
            theseImages{2,trl+ims} = thisImage2;
        end
        DrawFormattedText(win, 'Get Ready!\n\n(Press any key to start) ','center','center');
        Screen('Flip', win);
        WaitSecs(.1)
         KbWait([], 2);
         WaitSecs(2)  
    end
   
    Screen('FillRect',win,bkg_gray);
    Screen('FillOval', win,[255 0 0] , rectOval);
    Screen('Flip', win);
    WaitSecs(.5)
    textureIndex(1) = Screen('MakeTexture', win,theseImages{1,trl});
    textureIndex(2) = Screen('MakeTexture', win,theseImages{2,trl});
    Screen('DrawTextures', win, textureIndex(1)); % background
    Screen('FillOval', win,[255 0 0] , rectOval);
    Screen('Flip', win); 
    WaitSecs(PARAMS.TRIALS.imageDuration);
    Screen('FillRect',win,bkg_gray);
    Screen('FillOval', win,[255 0 0] , rectOval);
    Screen('Flip', win);
    WaitSecs(PARAMS.TRIALS.ISI);
    Screen('DrawTextures', win, textureIndex(2)); % background
    Screen('FillOval', win,[255 0 0] , rectOval);
    Screen('Flip', win); 
    WaitSecs(PARAMS.TRIALS.imageDuration);

    Screen('FillRect',win,bkg_gray);
    DrawFormattedText(win, 'Same(s) or different(d)?\n(Press Esc to abort the experiment)','center',round(PARAMS.setup.scr_pixhgt*.25),[0 0 0]);
    Screen('Flip', win);

    while KbCheck; end % Wait until all keys are released.

    while 1
        % Check the state of the keyboard.
        [ keyIsDown, seconds, keyCode ] = KbCheck;

        % If the user is pressing a key, then display its code number and name.
        if keyIsDown
            if strcmp(KbName(keyCode),'s')
               TRIALS.response{trl} = 's';
               if  strcmp(TRIALS.trialCondition(trl),'same')
                    TRIALS.result{trl} ='H';
               else
                    TRIALS.result{trl} ='FA';
               end
                break
            elseif strcmp(KbName(keyCode),'d')   
               TRIALS.response{trl} = 'd';
               if  strcmp(TRIALS.trialCondition(trl),'same')
                    TRIALS.result{trl} ='M';
               else
                    TRIALS.result{trl} ='CR';
               end
               
                break
            elseif keyCode(escapeKey)
%                 sca
                flagscape = 1;
                break
                %error(sprintf('Experiment aborted by experimenter on trial %d',trl))
            end
        end
    end
    if flagscape
        break
    else
        
      save(fullfile(paths.result,fileMAT),'PARAMS','TRIALS')
        
    end
    
    Screen('FillRect',win,127);
    Screen('FillOval', win,[255 0 0] , rectOval);
    Screen('Flip', win);
  
end
sca

%%
% result summary


 iX       = ismember(TRIALS.trialCondition,'same');
 hits     = sum(ismember(TRIALS.result(iX),'H'));
 misses   = sum(ismember(TRIALS.result(iX),'M'));
 Hp       = hits./(hits+misses);
 if Hp==1
     Hp=1-1/2/sum(iX);
 elseif Hp==0
     Hp= 1/2/sum(iX);
 end
zHp       =  norminv(Hp,0,1);
for ee = find(~ismember({PARAMS.conditions.name},{'same'}))
    thisCond = PARAMS.conditions(ee).name;
    iX       = ismember(TRIALS.trialCondition,thisCond);
    FA       = sum(ismember(TRIALS.result(iX),'FA'));
    CR       = sum(ismember(TRIALS.result(iX),{'CR'}));
    FAp      = FA./(FA+CR);
     if FAp==1
        FAp=1-1/2/sum(iX);
    elseif FAp==0
        FAp= 1/2/sum(iX);
    end
    zFAp     =  norminv(FAp,0,1);
    dprime   = zHp-zFAp;
    % fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    fprintf('\n\n%s\n',thisCond)
    fprintf('%%%%%%%%%%%%%%%%%%Response%%%%%%%%%%%%%%%%\n')
    fprintf('        ''S'' | ''D'' | Total\n')
    fprintf('--------------------------\n')
    fprintf('Orig.    %02d | %02d  | %d \n',hits,misses,hits+misses)
    fprintf('Diff.    %02d | %02d  | %d \n',FA,CR,FA+CR)
    fprintf('\nd'' = %1.2f\n',dprime)
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    TRIALS.([thisCond '_dprime']) = dprime;

end
figure,plotCSF(PARAMS.subject.spfreqs,PARAMS.subject.S)
save(fullfile(paths.result,fileMAT),'PARAMS','TRIALS')

