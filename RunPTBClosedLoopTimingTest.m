clear; clc; close all;

% Written by Matt Anderson on 2/12/2025. Submit bug reports via the dedicated
% github page: https://github.com/mattanderson94/Method-for-evaluating-closed-loop-latency-of-gaze-contingent-rendering-

% the following test works on this system:

% -----------------------------------------------------------------------------------------------------
% MATLAB Version: 23.2.0.2668659 (R2023b) Update 9
% MATLAB License Number: xxxxxxx
% Operating System: Microsoft Windows 10 Pro Version 10.0 (Build 19045)
% Java Version: Java 1.8.0_202-b08 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
% -----------------------------------------------------------------------------------------------------
% MATLAB                                                Version 23.2        (R2023b)   
% Computer Vision Toolbox                               Version 23.2        (R2023b)   
% Curve Fitting Toolbox                                 Version 23.2        (R2023b)   
% DatapixxToolbox                                       Version 0.9,        Aug        
% Image Processing Toolbox                              Version 23.2        (R2023b)   
% Optimization Toolbox                                  Version 23.2        (R2023b)   
% Psychtoolbox                                          Version 3.0.19      17 February
% Signal Processing Toolbox                             Version 23.2        (R2023b)   
% Statistics and Machine Learning Toolbox               Version 23.2        (R2023b)   

% We are using a desktop-mounted Eyelink 1000+, with the latest Eyelink
% Devkit API downloaded (as of 2/12/2025). See https://www.sr-research.com/support/forum-9.html

%% -------    PRELIMINARY PUPIL PROPERTIES    ------- %%

% edf/mat file name?
fname = 'TestMeasurementCR';

% render corneal reflection or no?
renderCR = true;

% define pupil-motion parameters
pupilRadiusDeg = .4;
slowPhaseAmplitudeDeg = 1;
slowPhaseSpeedDegPerSec = 3;
addedLatencySecs = 0; % if not a multiple of the frame rate, this is rounded up to the nearest frame
trialDurationSecs = 60;

% screen dimensions
screenWidth = 170;
screenHeight = 96;
screenDistance = 100;

% pupil and bg cols
pupilColor = [0,0,0];
cornealReflectionColor = [255,255,255];
foregroundColor = [255,255,255];
backgroundColor = [100,100,100];

%% -------    PSYCHTOOLBOX SETUP    ------- %%

% Ideally, we should have no sync issues (timing really matters here)
Screen('Preference', 'SkipSyncTests', 1);

% boilerplate PTB stuff...
PsychDefaultSetup(2);
screens = Screen('Screens');
screenNumber = max(screens);

[window, windowRect] = Screen('OpenWindow', screenNumber, backgroundColor);
[pxWidth, pxHeight] = Screen('WindowSize', window);
refreshRate = Screen('NominalFrameRate', window);

%% -------    DEFINE OTHER PUPIL PROPERTIES BASED ON SCREEN RESOLUTION AND REFRESH RATE    ------- %%

% convert degrees to pixels
wholescreendeg = rad2deg(atan((screenWidth/2)/screenDistance)*2);
pxperdeg = pxWidth/wholescreendeg;
pupilradpx = round(pupilRadiusDeg*pxperdeg);
amplitudepx = round(slowPhaseAmplitudeDeg*pxperdeg);
speedpx = slowPhaseSpeedDegPerSec*pxperdeg;

% if we draw the corneal reflection, we make sure it's smaller than the
% pupil by same fraction
CRradpx = pupilradpx/9; % radius
CRxd = pupilradpx/2; % x-displacement of CR, relative to CENTER of pupil

% place pupils at the center of the screen.
Lpupilx = round(pxWidth*.450);
Rpupilx = round(pxWidth*.550);
startpupily = round(pxHeight*.3); % y-location a bit higher than the center. Works better for our setup

% based on the nominal refresh rate, calculate how much we want to move the
% pupil on each frame
dy = speedpx/refreshRate; % y-displacement per frame
dyrect = -[[0;dy;0;dy],[0;dy;0;dy]]; % reformat to a rect we can add to the pupil position rects

% pupil motion is a triangular waveform, so it moves up and down
% periodically. The amplitude is fixed, and since we define dy per frame,
% we can just count the frames, and then reverse the pupil direction based
% on the frame count
nframesAmp = round(amplitudepx/speedpx*refreshRate);
framectr = 1; % start ctr at 1

% build initial pupil and CR rects. The two columns are the left and right
% pupils respectively
startpupilrects = [[Lpupilx-pupilradpx; startpupily-pupilradpx; Lpupilx+pupilradpx; startpupily+pupilradpx],...
    [Rpupilx-pupilradpx; startpupily-pupilradpx; Rpupilx+pupilradpx; startpupily+pupilradpx]];

startCRrects = [[Lpupilx-CRradpx-CRxd; startpupily-CRradpx; Lpupilx+CRradpx-CRxd; startpupily+CRradpx],...
    [Rpupilx-CRradpx-CRxd; startpupily-CRradpx; Rpupilx+CRradpx-CRxd; startpupily+CRradpx]];

% put instructions on screen at the bottom
% txtRectX = pxWidth*.5;
txtRectY = pxHeight*.9;

%% -------    PUPIL-ONLY SETUP FOR EYELINK    ------- %%

% the main thing we need to do here is just enable pupil-only tracking if
% we don't want to draw the CR. To do this, we need to connect the tracker,
% reconfigure it, and then disconnect. Then the Pupil tracking mode should
% be visible on the monitor connected to the eyetracker...
el = EyelinkInitDefaults(window);

if ~renderCR
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    Eyelink('Shutdown');
    if Eyelink('Initialize',  'PsychEyelinkDispatchCallback') ~= 0
        fprintf('Eyelink Init aborted.\n');
        Eyelink('Shutdown');  % cleanup function
        error('EYELINK could not be initialized');
    end
    
    Eyelink('command', 'force_corneal_reflection = FALSE')
    Eyelink('command', 'corneal_mode = FALSE')
    Eyelink('command', 'allow_pupil_without_cr = TRUE')
    Eyelink('command','inputword_is_window = ON');

    Eyelink('Shutdown');

end

%% -------    PUT PUPIL ON-SCREEN FOR EXPERIMENTER TO SET UP EYETRACKER    ------- %%

% Now let's put our pupil on the screen, and allow the experimenter to
% adjust the camera to make sure it is captured correctly

Screen('FillOval', window, pupilColor, startpupilrects);

if renderCR
    Screen('FillOval', window, cornealReflectionColor, startCRrects);
end
DrawFormattedText(window,'Set up your camera now so that it can capture the two pupils simultaneously\nPress any key to continue when you are ready','Center',txtRectY,foregroundColor);
Screen('Flip', window);

KbStrokeWait

%% -------    NOW SIMULATE THE PUPIL MOTION AND MAKE SURE THE CAMERA IS PICKING UP THE CHANGES IN GAZE POSITION    ------- %%

checkingSetup = true;

% Define key codes
escKey = KbName('Escape');  % Escape key
cKey = KbName('c');         % 'C' key

while checkingSetup

    % shift pupils up or down
    newpupilrects = startpupilrects-dyrect;
    startpupilrects = newpupilrects;

    % shift CR up or down
    newCRrects = startCRrects-dyrect;
    startCRrects = newCRrects;

    % measure if max amplitude is met. Switch direction if so.
    framectr = framectr+1;
    if framectr>nframesAmp
        framectr = 1;
        dyrect = -dyrect;
    end

    Screen('FillOval', window, pupilColor, newpupilrects);

    if renderCR
        Screen('FillOval', window, cornealReflectionColor, newCRrects);
    end

    DrawFormattedText(window,'Check that gaze position is smoothly tracked on the eyelink monitor\nIf it is not, press "Esc" and do setup again\nIf tracking looks good, press "c" to continue','Center',txtRectY,foregroundColor);
    Screen('Flip', window);

    [keyIsDown, ~, keyCode] = KbCheck();
    
    if keyIsDown && keyCode(escKey)
        sca;
        return;
    end

    if keyIsDown && keyCode(cKey)
        break;  
    end

end


%% -------    NOW WE BEGIN PHASE 1 MEASUREMENT    ------- %%

% continue to move the pupil along the same trajectory as before.
% Measurement here will take place in two phases:

% PHASE 1: We measure pupil position for calibration purposes only. In
% gaze-contingent rendering, we need to know the mapping between the sensor
% (camera) coordinates of the pupil (-CR), and the screen coordinates. A
% calibration model is typically built for this purpose. For our current
% setup, it's not easy to rely on the standard 9-point calibration grid. But, 
% a simple solution here is to collect a bunch of measurements of our
% pupils moving up and down, measure the raw coordinates reported by the
% eye tracker, and then build a polynomial calibration model. Since X is fixed in all
% cases, we only need to do this for the y dimension. Note that the
% accuracy of this model is not super critical. It just makes the analysis
% a bit easier further down the line. 

% first some more eyelink setup
if Eyelink('Initialize',  'PsychEyelinkDispatchCallback') ~= 0
    fprintf('Eyelink Init aborted.\n');
    Eyelink('Shutdown');  % cleanup function
    error('EYELINK could not be initialized');
end

Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA,GAZERES,HREF,PUPIL,STATUS,INPUT,HMARKER');
Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,AREA,GAZERES,HREF,PUPIL,STATUS,INPUT,HMARKER');
el.edfFile = sprintf('ELTmp.edf'); 
Eyelink('openfile', el.edfFile);
Eyelink('StartRecording');

% now we define our variables for the calibration loop
calibrateDurationSecs = 10;

nrows = ceil(refreshRate*calibrateDurationSecs*1.5);
emptnum = nan(nrows,1);
calibgazedata = table(emptnum,emptnum,emptnum,emptnum,emptnum,emptnum,emptnum,emptnum,emptnum,emptnum,...
    'variablenames', {'ELtime','PTBtime','LGazeX','LGazeY','RGazeX','RGazeY','PTBLGazeX','PTBLGazeY','PTBRGazeX','PTBRGazeY'});
tblctr = 1;

currtime = GetSecs;
while (GetSecs-currtime) < calibrateDurationSecs

    % get eyelink gaze position data in px if available
    evt = GetCurrentDataRaw(el);
    gazeX = evt.px;
    gazeY = evt.py;

    % shift pupils up or down
    newpupilrects = startpupilrects-dyrect;
    startpupilrects = newpupilrects;

    % shift CR up or down
    newCRrects = startCRrects-dyrect;
    startCRrects = newCRrects;

    % measure if max amplitude is met. Switch direction if so.
    framectr = framectr+1;
    if framectr>nframesAmp
        framectr = 1;
        dyrect = -dyrect;
    end

    % Draw pupils
    Screen('FillOval', window, pupilColor, newpupilrects);

    if renderCR
        Screen('FillOval', window, cornealReflectionColor, newCRrects);
    end

    DrawFormattedText(window,'Calibrating...','Center',txtRectY,foregroundColor);
    Screen('Flip', window);

    % also, save all the gaze-contingent pos and flip-time data
    calibgazedata.ELtime(tblctr) = evt.time;

    PTBgx = mean(newpupilrects([1,3],:));
    PTBgy = mean(newpupilrects([2,4],:));

    calibgazedata.PTBLGazeX(tblctr) = PTBgx(1);
    calibgazedata.PTBRGazeX(tblctr) = PTBgx(2);
    calibgazedata.PTBLGazeY(tblctr) = PTBgy(1);
    calibgazedata.PTBRGazeY(tblctr) = PTBgy(2);

    calibgazedata.LGazeX(tblctr) = evt.px(1);
    calibgazedata.RGazeX(tblctr) = evt.px(2);
    calibgazedata.LGazeY(tblctr) = evt.py(1);
    calibgazedata.RGazeY(tblctr) = evt.py(2);
    tblctr = tblctr+1;

end

Eyelink('StopRecording');

% The relationship between on-screen pixel position and sensor-based gaze y
% position is highly linear. We can capture it with just a two-term model I think
idxs = 1:(tblctr-1); 
p = polyfit(calibgazedata.LGazeY(idxs),calibgazedata.PTBLGazeY(idxs),1);  

% confirm yourself with some plots (note that the resulting data may be
% rectangular, and this reflects the fact that there is of course a latency
% (which we are trying to measure!!!) between reported on-screen position
% from PTB, and reported on-screen position from the Eyelink. 
% figure; hold on;
% scatter(calibgazedata.PTBLGazeY,calibgazedata.LGazeY)
% scatter(calibgazedata.PTBRGazeY,calibgazedata.RGazeY)

%% -------    NOW WE BEGIN PHASE 2 MEASUREMENT    ------- %%

% now we collect the gaze data for real, to extract our latency estimates. 

% Phase 2: Here, we use our calibration model built above, but this time we
% continuously draw the second pupils' location based on the measured
% location of the first pupil. 

% define fresh table for storing data if needed
nrows = ceil(refreshRate*trialDurationSecs*1.5);

emptnum = nan(nrows,1);
emptcell = cell(nrows,1);
gazedata = table(emptnum,emptnum,emptnum,emptnum,emptnum,emptnum,emptnum,emptnum,emptnum,emptnum,emptnum,emptcell,emptcell,...
    'variablenames', {'ELtime','PTBtime','LGazeX','LGazeY','RGazeX','RGazeY','PTBLGazeX','PTBLGazeY','PTBRGazeX','PTBRGazeY','RGazeYPred','allPupilRects','allCRRects'});

tblctr = 1;

nframesAddedLatency = ceil(addedLatencySecs*refreshRate);

% overwrite the calibration file
el.edfFile = sprintf('ELTmp.edf'); 
Eyelink('openfile', el.edfFile);
Eyelink('StartRecording');

currtime = GetSecs;

while (GetSecs-currtime) < trialDurationSecs

    % get eyelink gaze position data in px if available
    evt = GetCurrentDataRaw(el);
    gazeX = evt.px;
    gazeY = evt.py;

    % use calibration model to get right-eye Y position
    evt.predgy = gazeY(1)*p(1)+p(2); 

    % shift up or down as usual
    newpupilrects = startpupilrects-dyrect;

    % but replace right-eye y position with the measured position of the
    % left eye (after applying calibration model)
    newpupilrects(2,2) = evt.predgy-pupilradpx;
    newpupilrects(4,2) = evt.predgy+pupilradpx;
    startpupilrects = newpupilrects;

    % do same for corneal reflection rect
    newCRrects = startCRrects-dyrect;
    newCRrects(2,2) = evt.predgy-CRradpx;
    newCRrects(4,2) = evt.predgy+CRradpx;
    startCRrects = newCRrects;

    % measure if max amplitude is met. Switch direction if so.
    framectr = framectr+1;
    if framectr>nframesAmp
        framectr = 1;
        dyrect = -dyrect;
    end

    drawpupilrects = newpupilrects;
    drawCRrects = newCRrects;

    % we can simulate an 'added' latency by taking a right pupil
    % location that corresponds to n frames back in time...
    pupilrectidx = tblctr-nframesAddedLatency;

    % the first few frames will have zero added latency, so should be
    % removed based on the latency setting
    if pupilrectidx>0 && addedLatencySecs>0

        % we only replace the right pupil because the right pupil
        % location is drawn based on the recorded left pupil
        pastpupilrects = gazedata.allPupilRects{pupilrectidx};
        drawpupilrects(:,2) = pastpupilrects(:,2);

        pastCRrects = gazedata.allCRRects{pupilrectidx};
        drawCRrects(:,2) = pastCRrects(:,2);

    end

    % Draw pupils
    Screen('FillOval', window, pupilColor, drawpupilrects);

    if renderCR
        Screen('FillOval', window, cornealReflectionColor, drawCRrects);
    end

    DrawFormattedText(window,'Measuring Latency...','Center',txtRectY,foregroundColor);
    Screen('Flip', window);

    % also, save all the gaze-contingent pos and flip-time data
    gazedata.ELtime(tblctr) = evt.time;

    PTBgx = mean(newpupilrects([1,3],:));
    PTBgy = mean(newpupilrects([2,4],:));

    gazedata.PTBLGazeX(tblctr) = PTBgx(1);
    gazedata.PTBRGazeX(tblctr) = PTBgx(2);
    gazedata.PTBLGazeY(tblctr) = PTBgy(1);
    gazedata.PTBRGazeY(tblctr) = PTBgy(2);

    gazedata.LGazeX(tblctr) = evt.px(1);
    gazedata.RGazeX(tblctr) = evt.px(2);
    gazedata.LGazeY(tblctr) = evt.py(1);
    gazedata.RGazeY(tblctr) = evt.py(2);

    gazedata.RGazeYPred(tblctr) = evt.predgy;

    gazedata.allPupilRects{tblctr} = newpupilrects;
    gazedata.allCRRects{tblctr} = newCRrects;

    tblctr = tblctr+1;

end

Eyelink('StopRecording');

Screen('CloseAll')

% retreive Eyelink file and store in current 'Data' folder
destpath = [fileparts(mfilename('fullpath')),filesep,'EyelinkData',filesep];
DownloadFile(el,destpath,[fname,'.edf'])

% truncate table to correct sz
gazedata = gazedata(1:tblctr-1,:);

% and save if needed
save([destpath,fname,'.mat'],"gazedata")



%% -------    USEFUL FUNCTIONS    ------- %%

function evt = GetCurrentDataRaw(el)

    eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
    if eye_used == el.BINOCULAR % if both eyes are tracked
        eye_used = el.LEFT_EYE; % use left eye
    end

    % get the sample in the form of an event structure
    evt = Eyelink('NewestFloatSampleRaw',eye_used);
    
    if eye_used ~= -1 % do we know which eye to use yet?

        % if we do, get current gaze position from sample
        x = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
        y = evt.gy(eye_used+1);

        % do we have valid data and is the pupil visible?
        if x~=el.MISSING_DATA && y~=el.MISSING_DATA && evt.pa(eye_used+1) > 0
            evt.mx=x;
            evt.my=y;
        end
    end

end

function [files]= DownloadFile(el, path, newFileName)

    % download data file
    try

        % close EL file and copy file over to local directory
        Eyelink('CloseFile');

        fprintf('Receiving data file ''%s''\n', newFileName);

        % do not overwrite existing file
        [~,name,ext] = fileparts(newFileName);
        ctr = 1;
        while exist(fullfile(path, newFileName),'file')
            ctr = ctr+1;
            newFileName = [name, sprintf('%02d',ctr), ext];
        end

        % send that data over boi!
        status=Eyelink('ReceiveFile',el.edfFile, ...
            convertStringsToChars(fullfile(path, newFileName)));
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        
        files = {newFileName};

    catch ex
        getReport(ex)
        cprintf('red', sprintf('++ EYELINK :: Problem receiving data file ''%s''\n', this.el.edfFile));
        files = {};
    end
end