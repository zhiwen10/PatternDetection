% add path for dependencies 
githubDir1 = 'C:\Users\Steinmetz lab\Documents\MATLAB';
addpath(genpath(fullfile(githubDir1, 'HS')))
% https://www.mathworks.com/matlabcentral/fileexchange/22756-horn-schunck-optical-flow-method
%%
githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines')))
addpath(genpath(fullfile(githubdir2, 'NeuroPattToolbox')))
% https://github.com/BrainDynamicsUSYD/NeuroPattToolbox
%%
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\flowfield'))
%%
mn = 'ZYE_0012';
td = '2020-10-16';
en = 5;
serverRoot = expPath(mn, td, en);
%%
wf_svd
%% get 2 seconds of trace for all pixels from TrialTime
newTrialTimes1 = 1681; duration = [0 2];ntime = 2;
Fs = 1/mean(diff(t));
eventID = 1;
[PeriEventV, winSamps] = eventLockedAvgSVD(U, dV, t, newTrialTimes1, eventID, duration);
Ur = reshape(U, size(U,1)*size(U,2), size(U,3));
PeriEventV1 = squeeze(PeriEventV);
Trace1d = Ur*PeriEventV1;
Trace1d = double(Trace1d);
Trace2d = reshape(Trace1d, size(U,1), size(U,2), size(Trace1d,2));
%% downsize trace by 2D smoothing
for k = 1:size(Trace2d,3)
    Trace2dDown(:,:,k) = imresize(Trace2d(:,:,k),0.25,'bicubic');
end
Trace2dDown2 = reshape(Trace2dDown, size(Trace2dDown,1)*size(Trace2dDown,2), size(Trace2dDown,3));
dsize = size(Trace2dDown,1);
%% filter 2-8Hz
npoints = 35*ntime+1;
traceTemp = Trace2dDown2-repmat(mean(Trace2dDown2,2),1,npoints);
% filter and hilbert transform work on each column
traceTemp = traceTemp';
[f1,f2] = butter(2, [2 8]/(Fs/2), 'bandpass');
traceTemp = filtfilt(f1,f2,traceTemp);
traceHilbert =hilbert(traceTemp);
tracePhase = angle(traceHilbert);
traceAmp = abs(traceHilbert);
traceTemp1 = reshape(traceTemp, size(traceTemp,1),dsize, dsize);
tracePhase1 = reshape(tracePhase,size(tracePhase,1),dsize, dsize);
traceHilbert1 = reshape(traceHilbert,size(traceHilbert,1),dsize, dsize);
traceAmp1 = abs(traceHilbert1);
%%
frame = 24;
A1 = squeeze(tracePhase1(frame ,:,:));
A2 = squeeze(tracePhase1(frame +1,:,:));
%% Horn and Shunck method
[fx, fy, ft] = computeDerivatives(A1,A2);
[u, v] = HS_mod(A1, A2);
% normalize flow field vector to length of 1
uvsum = sqrt(u.^2+v.^2);
u1 = u./uvsum;
v1 = v./uvsum;
%% use our filtered data, but GONG's flow field
Fs = 35;
params = setNeuroPattParams_mod(Fs);
wvcfs1 = permute(traceHilbert1,[2,3,1]);
[vx, vy, csteps] = opticalFlow(wvcfs1(:,:,24:25),[],params.opAlpha, params.opBeta, ~params.useAmplitude);
% normalize flow field vector to length of 1
uvsum1 = sqrt(vx.^2+vy.^2);
vx1 = vx./uvsum1;
vy1 = vy./uvsum1;
%% full piepeline from GONG et al., including filtering and hilbert transform
tic
data = Trace2dDown;
data = bsxfun(@minus, data, mean(data,3));
timeDim = 3;
ntimesteps = size(timeDim,3) - 1;
wvcfs = filterSignal(data,params.hilbFreqLow, params.hilbFreqHigh,Fs, timeDim, true);
[vxGong, vyGong, csteps] = opticalFlow(wvcfs(:,:,24:25),[],params.opAlpha, params.opBeta, ~params.useAmplitude);
% normalize flow field vector to length of 1
uvsumGong = sqrt(vxGong.^2+vyGong.^2);
vxGong1 = vxGong./uvsumGong;
vyGong1 = vyGong./uvsumGong;
T = toc;
%% plot normalized flow field
figure; 
h1 = subplot(3,3,1)
imagesc(A1(2:2:end,2:2:end)+pi);
colormap(hsv);axis image
originalSize1 = get(gca, 'Position');
colorbar
set(h1, 'Position', originalSize1);
title('PhaseMap1')
h2 = subplot(3,3,2)
imagesc(A2(2:2:end,2:2:end)+pi);
colormap(hsv);axis image
originalSize1 = get(gca, 'Position');
colorbar
set(h2, 'Position', originalSize1);
title('Phasemap2')

h4 = subplot(3,3,4)
imagesc(fx(2:2:end,2:2:end));
colormap(hsv);axis image
originalSize1 = get(gca, 'Position');
colorbar
set(h4, 'Position', originalSize1);
title('X derivative')
h5 = subplot(3,3,5)
imagesc(fy(2:2:end,2:2:end));
colormap(hsv);axis image
originalSize1 = get(gca, 'Position');
colorbar
set(h5, 'Position', originalSize1);
title('Y derivative')
h6 = subplot(3,3,6)
imagesc(ft(2:2:end,2:2:end));
colormap(hsv);axis image
originalSize1 = get(gca, 'Position');
colorbar
set(h6, 'Position', originalSize1);
title('T derivative')

% Horn and Shunck method
h7 = subplot(3,3,7)
quiver(u1(2:2:end,2:2:end),v1(2:2:end,2:2:end),1)
set(gca,'YDir','reverse'); axis image
originalSize2 = get(gca, 'Position');
title({'Order2 filter','Horn & Shunck'})
% Charbonnier penalty
h8 = subplot(3,3,8)
quiver(vx1(2:2:end,2:2:end),vy1(2:2:end,2:2:end),1)
set(gca,'YDir','reverse'); axis image
originalSize2 = get(gca, 'Position');
title({'Order2 filter','Charbonnier penalty'})
% Charbonnier penalty
h9 = subplot(3,3,9)
quiver(vxGong1(2:2:end,2:2:end),vyGong1(2:2:end,2:2:end),1)
set(gca,'YDir','reverse'); axis image
originalSize2 = get(gca, 'Position');
title({'Order8 filter','Charbonnier penalty'})