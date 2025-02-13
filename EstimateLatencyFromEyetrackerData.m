clear; clc; close all;

%% Specify system settings and analysis method

addpath(genpath('AnalysisTools/'));

fname = 'TEST1000HZ_CR_0MsecsLatency';

% manually specify our setup 
refreshRate = 144;
sampleRate = 1000; % per second
emAmplitude = 1; % in degrees


% latency analysis method #1: 
% find matching pairs of left-pupil and right-pupil positions. If the
% sample rate is higher than the screen refresh rate (e.g., 1000 Hz versus
% 144 Hz), this is the preferred method.

% latency analysis method #2: 
% estimate velocity of the first pupil and spatial offset of the left and 
% right pupils at each sample. Then predict time-offset from these values.
% Mainly useful for lower eye tracker sample rates. (e.g., 250 Hz)
latencymeasurementmethod = 1; 

if latencymeasurementmethod == 1

    % super important: 
    % this parameter is used to locate the eye tracker samples corresponding to the time when
    % the on-screen pupil changes position (versus when it was static). You may need to mess around
    % with this value. The way to do this is to add a break-point the line
    % before this value is used, and just look at the velocity values when
    % (i) the pupil changes position -- i.e., signal, versus (ii) when the pupil
    % is static -- noise. I did this by eye, but use whatever method you like
    velocityChangeThreshold = .0025;

else
    msecsFilter = 40; 
    nsampsFilter = msecsFilter/(1000/sampleRate);
end

%% Load in files
    
edffname = sprintf('EyelinkData%s%s.edf',filesep,fname);

% if edf2mat does not work, type in 'open Edf2Mat' and follow the
% instructions on the top of the script, or visit https://github.com/uzh/edf-converter
emdata = Edf2Mat(edffname);
    
% load in corresponding matlab table
matfname = sprintf('EyelinkData%s%s.mat',filesep,fname);
data = load(matfname);

%% make sure we only take the data during on-screen pupil motion

% remove samples outside trial
t = emdata.Samples.time;
idxs = t>data.gazedata.ELtime(1) & t<data.gazedata.ELtime(end);
t = t(idxs)/1000; % to seconds

% grab raw pupil data
gx = emdata.Samples.px(idxs,:);
gy = emdata.Samples.py(idxs,:);

%% Now we rescale the data so we can match up the y positions accurately

% normalize the left and right eye data to mean=0, and std=1
gy(:,1) = (gy(:,1)-mean(gy(:,1)))./std(gy(:,1));
gy(:,2) = (gy(:,2)-mean(gy(:,2)))./std(gy(:,2));

% although we have normalized the left and right-eye pupil data,
% so that the kinks should match up (y-position should be similar for a
% given frame), there may be a slow drift in measured y-position due to
% lighting, small mirror position changes (I used a hinged mirror,
% which may unfold/fold on its own), etc. Basically, I'm trying to
% account for the lack of control here. To correct for this, I simply
% fit a polynomial to the data, and undo the fitted effects. This way,
% everything is nice and flat over time
p = polyfit(t,gy(:,1),3);
gy(:,1) = gy(:,1) + polyval(-p,t);

p = polyfit(t,gy(:,2),3);
gy(:,2) = gy(:,2) + polyval(-p,t);

% and we specified the eye-movement amplitude in the script, so we can
% enforce that here
gy = gy./range(gy(:)).*emAmplitude;

%% Estimate frame-by-frame latency!

if latencymeasurementmethod == 2

    gy = movmean(gy,nsampsFilter,1);

    % note that, for very long latencies, more data need to be removed
    lthresh = prctile(gy(:,1),5);
    uthresh = prctile(gy(:,1),95);
    
    % this is the analysis method that only works for periods inbetween the
    % reversals. Here, we find the locations where velocities tend to flatten out
    vidxs = gy(:,1)>lthresh & gy(:,2)>lthresh & gy(:,1)<uthresh & gy(:,2)<uthresh;
    
    % how much does the left eye move each frame? 
    dyL = median(abs(diff(gy(vidxs,1))))/median(diff(t(vidxs))); % in pixels per second
    
    % now subtract the L and R signals to get the estimated offset per-frame
    Rdelay = abs((gy(vidxs,1)-gy(vidxs,2))./dyL)*1000; % convert back to secs

else

    % This method only really works when the sample rate far-exceeds
    % the refresh rate of the display.

    % First, we find matches in y-gaze positions. To do this, we 
    % find the kinks in the triangular waveform (when same y pos is
    % measured repeatedly). 
    dyL = abs([0; diff(gy(:,1))]); 
    locsL = find(dyL<velocityChangeThreshold); % arbitrary threshold at 1000Hz. works well. May not at other sampling rates or noise levels
    locsLflat = locsL(diff([locsL(1);locsL])>1)-1; % take first in flat period

    % now same for Right-pupil
    dyR = abs([0;diff(gy(:,2))]); 
    locsR = find(dyR<velocityChangeThreshold);
    locsRflat = locsR(diff([locsR(1);locsR])>1)-1; 

    % now we can see what function we can apply to best align the left and
    % right-eye pupil data. First let's find the right-eye position closest
    % to the left-eye position within a certain time-frame
    nsampwindow = 15; % time-frame to search for matches. This is in eyelink samples, so may be sensitive to sample rate and latency

    n = length(locsLflat);
    locsLflatMatch = zeros(n,1);

    nnomatch = 0;
    nmatch = 0;

    for j = 1:n
        
        gyL = gy(locsLflat(j),1);

        windowidxs = locsRflat(j:min(length(locsRflat),(j+nsampwindow)));

        gyRall = gy(windowidxs,2);

        ydiff = abs(gyL-gyRall);

        % if there are multiple peaks (at the reversals), take the first
        % peak that we observe
        if length(ydiff)>3

            [pks,locs] = findpeaks(-abs(ydiff));

            % never re-use the same match:
            alreadyusedidxs = ismember(windowidxs(locs),locsLflatMatch);
            locs(alreadyusedidxs) = [];
        
            if isempty(locs)
                fprintf('No match/peak found. Replacing with NaN...\n')
                nnomatch = nnomatch+1;
                locsLflatMatch(j) = NaN;
            else
                matchidx = locs(1);
                nmatch = nmatch+1;
                locsLflatMatch(j) = windowidxs(matchidx);
            end
        else
            locsLflatMatch(j) = NaN;
        end

    end

    % remove the first n samples because we need a 'burn-in' period
    locsLflatMatch(1:nsampwindow) = NaN;

    % an estimate of the delay can be calculated by simply the
    % time-difference between the matches
    tdelay = locsLflatMatch-locsLflat;
    Rdelay = tdelay;

    % for plotting, let's show one second starting at the 500th frame
    locsLflattidxs = locsLflat(locsLflat>=500 & locsLflat<=(sampleRate+500));
    locsRflattidxs = locsRflat(locsRflat>=500 & locsRflat<=(sampleRate+500));

    if nnomatch>1
        msg = sprintf('\n\nTotal of %i non-matches found. This represents %.2f%% of total samples.\nIf this percentage exceeds 5.00%%, then thoroughly check your data, and try some different sample-matching parameters.\n',nnomatch,nnomatch/(nnomatch+nmatch)*100);
        warning(msg);
    end

end



%% And perform cross-correlation analysis too

% if we have isolated single points per on-screen pupil stimulus, use those
% points. If not, use all available points
if latencymeasurementmethod == 2
    tdelay = median(diff(t))*1000; % median delay between samples
else
    tdelay = median(diff(locsLflat)); % median delay between matched pts
end

msecdelays = 150; % range of delays to test
ndelays = floor(msecdelays/tdelay);
poscorrs = zeros(ndelays+1,1);
velcorrs = zeros(ndelays+1,1);
poserr = zeros(ndelays+1,1);

for j = 0:ndelays

    if latencymeasurementmethod == 2
        tdelayidxs = (j+1):length(t);
        x = gy(tdelayidxs,2);
        y = gy(:,1);
        y = y(1:length(x)); % ensure same length
    else
        tdelayidxs = locsLflat((j+1):end);
        x = gy(tdelayidxs,2);
        y = gy(locsLflat,1);
        y = y(1:length(x)); % ensure same length
    end

    % spatial correlation
    idx = ~isnan(y) & ~isnan(x);
    poscorrs(j+1) = corr(x(idx),y(idx));

    % difference
    poserr(j+1) = sum((x-y).^2);

    % velocity
    x = diff(x);
    y = diff(y);

    % velocity correlation
    idx = ~isnan(y) & ~isnan(x);
    velcorrs(j+1) = corr(x(idx),y(idx));
    
end

tdelays = 0:tdelay:msecdelays;

% peak of crosscorrelation function
[cmax,maxidx] = max(velcorrs);
tmax = tdelays(maxidx);
        

%% Plot eye-trace

tidxs = 500:(sampleRate+500); 


figure; 
hold on; 

plt1 = plot(tidxs,gy(tidxs,1),'LineWidth',2,'Color',[1.0000    0.4118    0.1608]);

% if we chose second method, we can indicate the matching points we
% found with markers. First for the left eye
if latencymeasurementmethod == 1
    scatter(locsLflattidxs,gy(locsLflattidxs,1),15,'MarkerFaceColor',[1.0000    0.4118    0.1608]);
end

plt2 = plot(tidxs,gy(tidxs,2),'LineWidth',2,'Color',[0.4667    0.6745    0.1882]);

% and again for the right eye
if latencymeasurementmethod == 1
    scatter(locsRflattidxs,gy(locsRflattidxs,2),15,'MarkerFaceColor',[0.4667    0.6745    0.1882]);
end

% and if we chose the first method, we can show those percentiles we
% used to threshold the data
if latencymeasurementmethod == 2
    plot(tidxs,repmat(uthresh,size(tidxs)),'LineWidth',2,'Color',[.5,.5,.5],'LineStyle','--')
    plot(tidxs,repmat(lthresh,size(tidxs)),'LineWidth',2,'Color',[.5,.5,.5],'LineStyle','--')
end

legend([plt1,plt2],'Left (drawn first)','Right')
legend boxoff
set(gca,'LineWidth',2,'FontSize',18)
box off;
ylabel('Y Gaze Position (Deg)')
xlabel('Time (msecs)')

% show from .5 secs in (exclude burn in period), and show 1 sec
xlim([500,500+sampleRate])
xticks(500:sampleRate/4:(sampleRate+500))

xticklabels(0:sampleRate/4:sampleRate)
ylim([-emAmplitude.*.55,emAmplitude.*.55])


%% And show the matched-up points with some arrows

if latencymeasurementmethod == 1
    figure
    hold on;

    % we apply a scaling factor to the y-axis here so that X and Y are on
    % a similar scale (quiver is annoyingly sensitive to scale differences)
    sfac = 400;

    plt1 = plot(tidxs,gy(tidxs,1).*sfac,'LineWidth',3,'Color',[1.0000    0.4118    0.1608]);
        scatter(locsLflattidxs,gy(locsLflattidxs,1).*sfac,75,'MarkerFaceColor',[1.0000    0.4118    0.1608]);

    plt2 = plot(tidxs,gy(tidxs,2).*sfac,'LineWidth',3,'Color',[0.4667    0.6745    0.1882]);
    scatter(locsRflattidxs,gy(locsRflattidxs,2).*sfac,75,'MarkerFaceColor',[0.4667    0.6745    0.1882]);

    valididxs = ~isnan(locsLflatMatch) & locsLflatMatch>=500 & locsLflatMatch<=(sampleRate+500);
    U = locsLflatMatch(valididxs)-locsLflat(valididxs);
    V = (gy(locsLflatMatch(valididxs),2)-gy(locsLflat(valididxs),1)).*sfac;
    X = locsLflat(valididxs);
    Y = gy(locsLflat(valididxs),1).*sfac;
    quiver(X,Y,U,V, 'off','LineWidth',1.5,'Color',[.7,.7,.7]);

    box off;
    % xlim([920,1120])
    %ylim([250,520])
    axis off;
end

%% Show distribution of latencies, based on applied analysis

edges = 0:3:100;

figure;
histogram(Rdelay,edges,'LineWidth',2,'FaceColor',[.6,.6,.6])
set(gca,'LineWidth',2,'FontSize',18)
box off;
axis square
ylabel('Frequency')
xlabel('Est. Time-Lag (msecs)')
xlim([0,100])
% ylim([0,10000])


%% and show crosscorrelation results, for position and velocity

figure;
hold on;
plot(tdelays,poscorrs,'Color','k','LineWidth',4)
plot(tdelays,velcorrs,'Color',[0.8196    0.4627    0.8902],'LineWidth',4)
plot([tmax,tmax],[0,cmax],'Color',[0.8196    0.4627    0.8902],'LineWidth',2,'LineStyle','-.')
set(gca,'LineWidth',2,'FontSize',18)
box off;
axis square
ylabel('Cross-correlation')
xlabel('Time-Lag (msecs)')
% ylim([.7,1])
xlim([0,max(tdelays)])

%% And print out basic descriptives if needed

mu = mean(Rdelay,'omitmissing');
mdn = median(Rdelay,'omitmissing');
sigma = std(Rdelay,[],'omitmissing');

fprintf('File %s, Mean = %.2f, Median = %.2f, Std = %.2f, Max of CrossCorr = %.2f \n\n',fname,mu,mdn,sigma,tmax)
