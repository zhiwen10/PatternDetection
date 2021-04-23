% add path for dependencies 
githubDir1 = 'C:\Users\Steinmetz lab\Documents\MATLAB';
addpath(genpath(fullfile(githubDir1, 'HS')))
%%
githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'NeuroPattToolbox')))
%%
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\flowfield'))
%%
mn = 'ZYE_0012';
td = '2020-10-16';
en = 5;
serverRoot = expPath(mn, td, en);
%%
wf_svd
%%
newTrialTimes1 = 1681;
foldername = [mn '_' td '_' num2str(newTrialTimes1) 's_frame24']
% mkdir(foldername);
duration = [0 2];
ntime = 2;
ntrials = length(newTrialTimes1);
fs = 1/mean(diff(t));
eventID = 1:size(newTrialTimes1,1);
eventID = eventID';
[avgPeriEventV, winSamps] = eventLockedAvgSVD(U, dV, t, newTrialTimes1, eventID, duration);
Ur = reshape(U, size(U,1)*size(U,2), size(U,3));

%% mean trace
avgPeriEventV1 = squeeze(mean(avgPeriEventV,1));
% avgPeriEventV1 = squeeze(avgPeriEventV(183,:,:));
meanTrace1 = Ur*avgPeriEventV1;
meanTrace1 = double(meanTrace1);
meanTrace = reshape(meanTrace1, size(U,1), size(U,2), size(meanTrace1,2));
%%
for k = 1:size(meanTrace,3)
    meanTraceDown(:,:,k) = imresize(meanTrace(:,:,k),0.12,'bicubic');
end
meanTraceDown2 = reshape(meanTraceDown, size(meanTraceDown,1)*size(meanTraceDown,2), size(meanTraceDown,3));
clear meanTrace1
meanTrace1 = meanTraceDown2;
dsize = size(meanTraceDown,1);
%% filter 2-8Hz
npoints = 35*ntime+1;
traceTemp = meanTrace1-repmat(mean(meanTrace1,2),1,npoints);
% filter and hilbert transform work on each column
traceTemp = traceTemp';
[f1,f2] = butter(2, [2 8]/(fs/2), 'bandpass');
traceTemp = filtfilt(f1,f2,traceTemp);
traceHilbert =hilbert(traceTemp);
tracePhase = angle(traceHilbert);
traceAmp = abs(traceHilbert);
traceTemp1 = reshape(traceTemp, size(traceTemp,1),dsize, dsize);
tracePhase1 = reshape(tracePhase,size(tracePhase,1),dsize, dsize);
traceHilbert1 = reshape(traceHilbert,size(traceHilbert,1),dsize, dsize);
traceAmp1 = abs(traceHilbert1);
AmpMax = max(max(max(traceAmp1)));
traceAmp1 = traceAmp1/AmpMax;
%%
params = setNeuroPattParams(35);
%%
for i = 24:30
    A1 = squeeze(tracePhase1(i,:,:));
    A2 = squeeze(tracePhase1(i+1,:,:));
    [u, v] = HS_mod(A1, A2);
    %%
    uvsum = sqrt(u.^2+v.^2);
    u1 = u./uvsum;
    v1 = v./uvsum;
    %%
    [critpointLocs] = findAllPatterns2(u, v, params);
    %%
    scale1 = 2;
    figure; quiver(u1(2:2:end,2:2:end),v1(2:2:end,2:2:end),1)
    set(gca,'YDir','reverse');
    axis image
    hold on; scatter(critpointLocs{1,1}(:,2)/scale1,critpointLocs{1,1}(:,1)/scale1,'m','lineWidth',2)
    hold on; scatter(critpointLocs{1,2}(:,2)/scale1,critpointLocs{1,2}(:,1)/scale1,'r','lineWidth',2)
    hold on; scatter(critpointLocs{1,3}(:,2)/scale1,critpointLocs{1,3}(:,1)/scale1,'g','lineWidth',2)
    hold on; scatter(critpointLocs{1,4}(:,2)/scale1,critpointLocs{1,4}(:,1)/scale1,'c','lineWidth',2)
    hold on; scatter(critpointLocs{1,5}(:,2)/scale1,critpointLocs{1,5}(:,1)/scale1,'k','lineWidth',2)
    % legend('flow field','Sink','Source','Spiral-in','Spiral-out','saddle')
    title(['frame' num2str(i)])
    saveas(gcf,['frame' num2str(i) '.png'])
end