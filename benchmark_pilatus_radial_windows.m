clearvars
close all
clc

%% PILATUS RADIAL-ONLY WINDOW BENCHMARK
%
% Purpose:
%   Benchmark the pipeline you actually want:
%
%       TIFF
%         -> bad-pixel mask
%         -> per-frame normalization from hard-coded ROI
%         -> normalized even-odd differences
%         -> average W pair-differences per output
%         -> radial DeltaI(q)
%
% Outputs kept during the benchmark:
%   1) one final 2-D average difference over the run
%   2) DeltaI(q,t/window-index)
%   3) per-frame normalization trace
%
% NO q/phi maps are calculated.
% NO per-window 2-D detector differences are saved.
%
% The benchmark compares different window sizes while keeping everything
% else fixed.  Every TIFF is still read once, so larger W mainly reduces
% the number of radial integrations and the number of DeltaI(q) outputs.
%
% MATLAB R2026a + Parallel Computing Toolbox
%
% RUN THIS FILE DIRECTLY.

%% ========================================================================
% USER SETTINGS
% ========================================================================

inDir = "X:\2026-3\RuettSept26\data\fresh_ag_1\run3_UV\";

benchRoot = ...
    "C:\Users\11idbuser\Downloads\Natan\pilatus_fast_diff_matlab\analysis\benchmark_radial_windows";

% Compare these window sizes.
% W=1  -> one pair per DeltaI(q)
% W=3  -> mean of 3 pair-differences per DeltaI(q)
% W=10 -> mean of 10 pair-differences per DeltaI(q)
WINDOWS_TO_TEST = [1 3 10];

% Use the CPU worker count you currently use.
CPU_WORKERS = 2;

% 180 TIFFs gives:
%   W=1  -> 90 radial outputs
%   W=3  -> 30 radial outputs
%   W=10 -> 9 radial outputs
BENCH_FILES_REQUESTED = 180;

% Radial resolution.
% There is little reason to reduce this aggressively because q-only
% integration is much smaller than q/phi. Keep 1000 unless you know you
% need less radial resolution.
NQ = 1000;

% Geometry/polarization convention.
% "divide"   : Icorr = I ./ (geometry .* polarization)
% "multiply" : Icorr = I .* (geometry .* polarization)
CorrectionMode = "divide";

DELETE_TEMP_OUTPUTS = true;

%% ========================================================================
% FILE LIST
% ========================================================================

D = dir(fullfile(inDir,"*.tif"));

if isempty(D)
    error('No TIFF files found in %s',inDir);
end

names     = string({D.name})';
paths     = fullfile(string({D.folder})',names);
fileBytes = double([D.bytes])';

frameIDs = nan(numel(names),1);

for k = 1:numel(names)
    tok = regexp(names(k),'\d+','match');
    if isempty(tok)
        error('No numeric frame ID found in filename: %s',names(k));
    end
    frameIDs(k) = str2double(tok{end});
end

[frameIDs,ord] = sort(frameIDs);
paths     = paths(ord);
fileBytes = fileBytes(ord);

nFiles = numel(paths);

fprintf('\n============================================================\n');
fprintf('PILATUS RADIAL-ONLY WINDOW BENCHMARK\n');
fprintf('============================================================\n');
fprintf('Files found : %d\n',nFiles);
fprintf('Input size  : %.3f TB\n',sum(fileBytes)/1e12);
fprintf('Frame IDs   : %d ... %d\n',frameIDs(1),frameIDs(end));
fprintf('Workers     : %d\n',CPU_WORKERS);
fprintf('q bins      : %d\n',NQ);
fprintf('Windows     : %s\n',mat2str(WINDOWS_TO_TEST));

%% ========================================================================
% CALIBRATION / MASK / NORMALIZATION
% ========================================================================

load("qgp.mat");
load("lab6avg_tz400.mat");

qMap         = QGP.qMap;
geometry     = QGP.geometryCorr;
polarization = QGP.polarizationCorr;

badMask = lab6_avg <= 0;
detectorValid = ~badMask;

% Hard-coded normalization ROI.
normMask = false(size(qMap));

NORM_ROWS = 1:408;
NORM_COLS = 990:1475;

normMask(NORM_ROWS,NORM_COLS) = true;
normMask(badMask) = false;

normIdx = find(normMask);

% Fixed reciprocal-space correction.
corrMap = geometry .* polarization;

if CorrectionMode == "divide"
    qStaticValid = detectorValid & isfinite(qMap) & ...
                   isfinite(corrMap) & corrMap ~= 0;
elseif CorrectionMode == "multiply"
    qStaticValid = detectorValid & isfinite(qMap) & isfinite(corrMap);
else
    error('CorrectionMode must be "divide" or "multiply".');
end

% Detector dimensions.
A0 = read_tiff_fast(paths(1));
[ny,nx] = size(A0);

if ~isequal(size(A0),size(qMap))
    error('TIFF is %dx%d but qMap is %dx%d.', ...
        ny,nx,size(qMap,1),size(qMap,2));
end

clear A0

fprintf('Detector     : %d x %d\n',ny,nx);
fprintf('Norm ROI     : %d good pixels\n',numel(normIdx));

%% ========================================================================
% PRECOMPUTE q-BIN MEMBERSHIP ONCE
% ========================================================================

qGood = qMap(qStaticValid);

qEdges = linspace(min(qGood),max(qGood),NQ+1);
qCenter = (qEdges(1:end-1)+qEdges(2:end))/2;

qBinMap = discretize(qMap,qEdges);

qMask = qStaticValid & ~isnan(qBinMap);
qBins = double(qBinMap(qMask));

if CorrectionMode == "divide"
    qFactor = single(1 ./ corrMap(qMask));
else
    qFactor = single(corrMap(qMask));
end

fprintf('q pixels     : %d\n',nnz(qMask));

%% ========================================================================
% PROCESS POOL SETUP - NOT TIMED
% ========================================================================

p = gcp('nocreate');

if CPU_WORKERS <= 1

    if ~isempty(p)
        delete(p);
    end

else

    if isempty(p) || p.NumWorkers ~= CPU_WORKERS

        if ~isempty(p)
            delete(p);
        end

        fprintf('Starting %d-process pool (not timed)...\n',CPU_WORKERS);
        parpool("Processes",CPU_WORKERS);
    end
end

%% ========================================================================
% BENCHMARK BLOCKS
% ========================================================================

nTests = numel(WINDOWS_TO_TEST);

% Make the timed block a multiple of ALL requested 2*window group sizes,
% so every test sees the same number of TIFFs with no discarded tail.
groupSizes = 2*WINDOWS_TO_TEST;
commonMultiple = groupSizes(1);

for k = 2:numel(groupSizes)
    commonMultiple = lcm(commonMultiple,groupSizes(k));
end

benchN = floor(BENCH_FILES_REQUESTED/commonMultiple)*commonMultiple;

if benchN < commonMultiple
    benchN = commonMultiple;
end

% One independent block per test to reduce cache reuse.
if nFiles < nTests*benchN
    benchN = floor((nFiles/nTests)/commonMultiple)*commonMultiple;
end

if benchN < commonMultiple
    error('Not enough files for the requested independent benchmark blocks.');
end

starts = round(linspace(1,nFiles-benchN+1,nTests));
starts = 1 + floor((starts-1)/commonMultiple)*commonMultiple;

for k = 2:nTests
    if starts(k) < starts(k-1)+benchN
        starts(k) = starts(k-1)+benchN;
    end
end

if starts(end)+benchN-1 > nFiles
    starts(end) = nFiles-benchN+1;
    starts(end) = 1 + floor((starts(end)-1)/commonMultiple)*commonMultiple;
end

fprintf('Timed/test   : %d TIFFs (~%.2f GB)\n', ...
    benchN,median(fileBytes)*benchN/1e9);

if isfolder(benchRoot)
    try
        rmdir(benchRoot,'s');
    catch
    end
end

mkdir(benchRoot);

%% ========================================================================
% BENCHMARK EACH WINDOW
% ========================================================================

elapsed      = nan(nTests,1);
tiffPerSec   = nan(nTests,1);
inputMBps    = nan(nTests,1);
projectedMin = nan(nTests,1);
nOutputsFull = nan(nTests,1);
outputMB     = nan(nTests,1);

for t = 1:nTests

    W = WINDOWS_TO_TEST(t);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('WINDOW = %d\n',W);
    fprintf('------------------------------------------------------------\n');

    idx = starts(t):starts(t)+benchN-1;

    timedPaths = paths(idx);
    timedBytes = sum(fileBytes(idx));

    nOutBench = benchN/(2*W);
    nOutFull  = floor(nFiles/(2*W));

    fprintf('Benchmark outputs : %d DeltaI(q) curves\n',nOutBench);
    fprintf('Full-run outputs  : %d DeltaI(q) curves\n',nOutFull);

    drawnow
    t0 = tic;

    R = run_radial_only_pipeline( ...
        timedPaths,W,ny,nx,NQ,CPU_WORKERS, ...
        detectorValid,normIdx,qMask,qBins,qFactor);

    % Include representative uncompressed output write in timing.
    R.q = qCenter;
    R.window = W;

    outFile = fullfile(benchRoot,sprintf('window_%d.mat',W));

    save(outFile,'R','-v7','-nocompression');

    elapsed(t) = toc(t0);

    info = dir(outFile);

    outputMB(t)     = info.bytes/1e6;
    tiffPerSec(t)   = benchN/elapsed(t);
    inputMBps(t)    = timedBytes/1e6/elapsed(t);
    projectedMin(t) = nFiles/tiffPerSec(t)/60;
    nOutputsFull(t) = nOutFull;

    fprintf('Elapsed           : %.2f s\n',elapsed(t));
    fprintf('Throughput        : %.2f TIFF/s\n',tiffPerSec(t));
    fprintf('Input throughput  : %.1f MB/s\n',inputMBps(t));
    fprintf('Output file       : %.1f MB\n',outputMB(t));
    fprintf('Projected full run: %.1f min\n',projectedMin(t));

    if DELETE_TEMP_OUTPUTS && isfile(outFile)
        delete(outFile);
    end
end

%% ========================================================================
% RESULTS
% ========================================================================

T = table( ...
    WINDOWS_TO_TEST(:), ...
    nOutputsFull, ...
    elapsed, ...
    tiffPerSec, ...
    inputMBps, ...
    outputMB, ...
    projectedMin, ...
    'VariableNames', ...
    {'Window','FullRun_DeltaIq_Curves','Seconds','TIFF_per_s', ...
     'Input_MB_per_s','Benchmark_Output_MB','Projected_full_min'});

T = sortrows(T,'Seconds','ascend');

fprintf('\n\n============================================================\n');
fprintf('FINAL WINDOW BENCHMARK\n');
fprintf('============================================================\n');
disp(T);

fprintf('\nInterpretation:\n');

fastestW = T.Window(1);

fprintf('Fastest measured window: W=%d\n',fastestW);

if height(T) >= 2

    speedGain = 100*(T.Seconds(2)-T.Seconds(1))/T.Seconds(2);

    fprintf('Speed advantage over second place: %.1f%%\n',speedGain);

    if speedGain < 5
        fprintf(['The timing difference is <5%%. Treat the windows as computationally\n' ...
                 'equivalent and choose W from the science/SNR/time-resolution tradeoff.\n']);
    end
end

fprintf('\nRemember:\n');
fprintf(['  The final average 2-D difference is essentially independent of W\n' ...
         '  (apart from an incomplete tail / masking details).\n']);
fprintf(['  W mainly changes how many pair-differences are averaged into each\n' ...
         '  DeltaI(q) time point:\n\n']);

for W = WINDOWS_TO_TEST
    fprintf('    W=%-3d -> %d DeltaI(q) curves over the full %d-file run\n', ...
        W,floor(nFiles/(2*W)),nFiles);
end

resultFile = fullfile(benchRoot,'benchmark_radial_windows_result.mat');

save(resultFile,'T','WINDOWS_TO_TEST','NQ','CPU_WORKERS');

fprintf('\nResults saved to:\n%s\n',resultFile);
fprintf('============================================================\n');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function R = run_radial_only_pipeline( ...
    paths,window,ny,nx,nq,nWorkers, ...
    detectorValid,normIdx,qMask,qBins,qFactor)

    nF = numel(paths);
    nOut = floor(nF/(2*window));

    deltaIq = nan(nq,nOut,'single');
    normCell = cell(nOut,1);

    % One reduction only: sum of all normalized pair differences.
    sumDiff2D = zeros(ny,nx,'double');

    if nWorkers > 1

        parfor g = 1:nOut

            [iq,normVals,localDiff] = process_one_window( ...
                paths,g,window,ny,nx,nq, ...
                detectorValid,normIdx,qMask,qBins,qFactor);

            deltaIq(:,g) = iq;
            normCell{g} = normVals;

            sumDiff2D = sumDiff2D + localDiff;
        end

    else

        for g = 1:nOut

            [iq,normVals,localDiff] = process_one_window( ...
                paths,g,window,ny,nx,nq, ...
                detectorValid,normIdx,qMask,qBins,qFactor);

            deltaIq(:,g) = iq;
            normCell{g} = normVals;

            sumDiff2D = sumDiff2D + localDiff;
        end
    end

    normTrace = nan(nF,1);

    for g = 1:nOut

        first = (g-1)*2*window + 1;

        normTrace(first:first+2*window-1) = normCell{g};
    end

    R.avgDiff2D = sumDiff2D/(nOut*window);
    R.deltaIq = deltaIq;
    R.normalization = normTrace;
end


function [iq,normVals,sumDiff2D] = process_one_window( ...
    paths,g,window,ny,nx,nq, ...
    detectorValid,normIdx,qMask,qBins,qFactor)

    first = (g-1)*2*window + 1;

    normVals = nan(2*window,1);

    sumDiff2D = zeros(ny,nx,'double');
    dWindow = zeros(ny,nx,'single');

    for j = 1:window

        i1 = first + 2*j - 2;
        i2 = i1 + 1;

        raw1 = read_tiff_fast(paths(i1));
        raw2 = read_tiff_fast(paths(i2));

        valid1 = isfinite(raw1) & raw1 >= 0;
        valid2 = isfinite(raw2) & raw2 >= 0;

        idx1 = normIdx(valid1(normIdx));
        idx2 = normIdx(valid2(normIdx));

        n1 = sum(double(raw1(idx1)));
        n2 = sum(double(raw2(idx2)));

        if ~(isfinite(n1) && n1 > 0 && isfinite(n2) && n2 > 0)
            error('Invalid normalization scalar in output group %d.',g);
        end

        normVals(2*j-1) = n1;
        normVals(2*j)   = n2;

        good1 = detectorValid & valid1;
        good2 = detectorValid & valid2;

        % Normalize before difference.
        f1 = single(raw1) / single(n1);
        f2 = single(raw2) / single(n2);

        f1(~good1) = NaN;
        f2(~good2) = NaN;

        d = f2 - f1;

        % For the running 2-D average, invalid pixels contribute zero.
        % Static bad pixels remain NaN only if you choose to re-mask the final
        % displayed image afterward.
        dz = d;
        dz(~isfinite(dz)) = 0;

        sumDiff2D = sumDiff2D + double(dz);
        dWindow   = dWindow + dz;
    end

    dWindow = dWindow / single(window);

    iq = binned_mean_cpu(dWindow,qMask,qBins,qFactor,nq);
end


function A = read_tiff_fast(filename)

    t = Tiff(char(filename),'r');
    A = t.read();
    t.close();
end


function out = binned_mean_cpu(I,mask,bins,factor,nBins)

    v = I(mask);

    good = isfinite(v);

    if ~any(good)
        out = nan(nBins,1,'single');
        return
    end

    b = bins(good);

    v = single(v(good)) .* factor(good);

    s = accumarray(b,double(v),[nBins 1],@sum,0);
    c = accumarray(b,1,[nBins 1],@sum,0);

    z = s ./ max(c,1);
    z(c==0) = NaN;

    out = single(z);
end
