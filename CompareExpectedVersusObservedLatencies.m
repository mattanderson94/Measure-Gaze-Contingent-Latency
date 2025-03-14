clear; clc; close all;

% This script is useful for validating our method. What we do here is use
% the latency estimate extracted by running
% EstimateLatencyFromEyetrackerData(), to find our 'baseline'. Then, we
% hard-code the latency value (in msecs) here, and check whether adding artificial delays
% produce the baseline plus the added delay. If so, we can show that our
% method produces measured latencies that match what we would expect given
% the frame rate of our monitor/setup. 

%% Specify system settings and analysis method

% add path to directory where @Edf2mat is located. For me, I created a
% directory called ./Analysistools/ and copied a fresh version into there
addpath(genpath('AnalysisTools/'));

% list all the files with artificially added latencies 
fnames = {'TEST1000HZ_7MsecsLatency','TEST1000HZ_14MsecsLatency',...
    'TEST1000HZ_21MsecsLatency','TEST1000HZ_56MsecsLatency'};
nfiles = length(fnames);

% hard-code baseline
baselinelatency = 21; % msecs

% hard-code the expected delay in each file. 
refreshRate = 144;
msecsdelayperfile = [7,14,21,56];
nframesdelay = round(msecsdelayperfile ./ ((1/refreshRate)*1000));

% manually specify our setup 
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


%% Set up our plots. There will be 4 in total: 
% (1) histogram of latencies per file
% (2) cross-correlation results per file
% (3) expected versus observed results from latency distribution analyses
% (4) expected versus observed results from crosscorrelation analyses

figure(1); hold on;
figure(2); hold on;
figure(3); hold on;
plot([0,100],[0,100],'Color',[.7,.7,.7],'LineStyle','-.','LineWidth',2)
figure(4); hold on;
plot([0,100],[0,100],'Color',[.7,.7,.7],'LineStyle','-.','LineWidth',2)

cols = [0.0745    0.6235    1.0000;
    0.8196    0.4627    0.8902;
    0.4667    0.6745    0.1882;
    0.9294    0.6941    0.1255];

for i = 1:nfiles
    
    %% Load in files
        
    edffname = sprintf('EyelinkData%s%s.edf',filesep,fnames{i});
    
    % if edf2mat does not work, type in 'open Edf2Mat' and follow the
    % instructions on the top of the script, or visit https://github.com/uzh/edf-converter
    emdata = Edf2Mat(edffname);
        
    % load in corresponding matlab table
    matfname = sprintf('EyelinkData%s%s.mat',filesep,fnames{i});
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
            
        
    %% distribution of latencies, based on applied analysis

    edges = 0:3:100;
    
    figure(1);
    histogram(Rdelay,edges,'LineWidth',2,'FaceColor',cols(i,:))
    
    %% crosscorrelation results, for velocity only
    
    figure(2);
    hold on;
    % plot(tdelays,poscorrs,'Color','k','LineWidth',4)
    plot(tdelays,velcorrs,'Color',cols(i,:),'LineWidth',4)
    plot([tmax,tmax],[0,cmax],'Color',cols(i,:),'LineWidth',2,'LineStyle','-.')
    
    %% And compare expected versus observed for both metrics

    expectedLatency = baselinelatency+(nframesdelay(i)*(1/refreshRate*1000));
    
    mdn = median(Rdelay,'omitmissing');
    % mu = mean(Rdelay,'omitmissing');
    % sigma = std(Rdelay,[],'omitmissing');

    figure(3);
    scatter(expectedLatency,mdn,175,'filled','MarkerFaceColor',cols(i,:),'MarkerEdgeColor','k')

    figure(4);
    scatter(expectedLatency,tmax,175,'filled','MarkerFaceColor',cols(i,:),'MarkerEdgeColor','k')
    
end

%% now apply global formatting
figure(1);
plot([baselinelatency,baselinelatency],[0,10000],'LineStyle',':','LineWidth',1.5,'Color','k')
set(gca,'LineWidth',2,'FontSize',18)
box off;
axis square
ylabel('Frequency')
xlabel('Est. Time-Lag (msecs)')
xlim([0,100])

figure(2);
set(gca,'LineWidth',2,'FontSize',18)
box off;
axis square
ylabel('Cross-correlation')
xlabel('Time-Lag (msecs)')
ylim([.7,1])
xlim([0,max(tdelays)])

figure(3)
axis square;
xlim([0,100])
ylim([0,100])
xlabel('Expected Latency (msecs)')
ylabel('Measured Latency (msecs)')
set(gca,'LineWidth',2,'FontSize',18)
box off;

figure(4)
axis square;
xlim([0,100])
ylim([0,100])
xlabel('Expected Latency (msecs)')
ylabel('Measured Latency (msecs)')
set(gca,'LineWidth',2,'FontSize',18)
box off;
