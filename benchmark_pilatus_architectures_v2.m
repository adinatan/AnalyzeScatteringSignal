clearvars
close all
clc

%% PILATUS ARCHITECTURE BENCHMARK v2
%
% Stand-alone benchmark for the current Ruett Sept 2026 run.
%
% TEST 1:
%   CPU, 1000 q x 360 phi, time-resolved DeltaI(q,phi,t)
%
% TEST 2:
%   CPU,  500 q x  90 phi, time-resolved DeltaI(q,phi,t)
%
% TEST 3:
%   GPU,  500 q x  90 phi, NO time-resolved DeltaI(q,phi,t)
%   Keeps:
%       DeltaI(q,t)
%       normalization trace
%       sum of all normalized detector images
%       sum of all raw detector images
%       average 2-D difference
%       I(q), I(q,phi) of the all-image sums
%       mean DeltaI(q,phi)
%
% Important efficiency changes in v2:
%   * normalization ROI indexing bug fixed
%   * every TIFF is read exactly once inside a timed pipeline
%   * q and q/phi pixel-bin membership is precomputed once
%   * parfor uses reductions instead of returning detector-sized arrays
%     from every iteration
%   * summed I(q,phi) and mean DeltaI(q,phi) are calculated ONCE at the end
%     from the corresponding accumulated 2-D images (linearity)
%   * geometry/polarization correction is applied only to pixels entering
%     q / q-phi integration
%
% MATLAB R2026a
% Parallel Computing Toolbox required.
%
% RUN THIS FILE DIRECTLY.

%% ========================================================================
% USER SETTINGS
% ========================================================================

inDir = "X:\2026-3\RuettSept26\data\fresh_ag_1\run3_UV\";

benchRoot = ...
    "C:\Users\11idbuser\Downloads\Natan\pilatus_fast_diff_matlab\analysis\benchmark_fast";

window = 3;

% Use the winner from your earlier 1/2/4-worker benchmark.
CPU_WORKERS = 2;

% 180 TIFFs = 30 output windows for window=3.
% Increase to 360 if the timings are too short/noisy.
BENCH_FILES_REQUESTED = 180;

% Untimed warm-up for each architecture.
WARMUP_FILES_REQUESTED = 12;

% geometry/polarization convention
% "divide":   Icorr = I ./ (geometry .* polarization)
% "multiply": Icorr = I .* (geometry .* polarization)
CorrectionMode = "divide";

DELETE_TEMP_OUTPUTS = true;

%% ========================================================================
% INPUT FILES
% ========================================================================

D = dir(fullfile(inDir,"*.tif"));

if isempty(D)
    error('No TIFF files found in %s',inDir);
end

names     = string({D.name})';
paths     = fullfile(string({D.folder})',names);
fileBytes = double([D.bytes])';

% Numeric sort using the LAST integer in each filename.
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
names     = names(ord);
fileBytes = fileBytes(ord);

nFiles = numel(paths);

fprintf('\n============================================================\n');
fprintf('PILATUS ARCHITECTURE BENCHMARK v2\n');
fprintf('============================================================\n');
fprintf('Files found : %d\n',nFiles);
fprintf('Input size  : %.3f TB\n',sum(fileBytes)/1e12);
fprintf('Frame IDs   : %d ... %d\n',frameIDs(1),frameIDs(end));
fprintf('Window      : %d (%d TIFFs/output)\n',window,2*window);
fprintf('CPU workers : %d\n',CPU_WORKERS);

%% ========================================================================
% CALIBRATION / MASK / NORMALIZATION ROI
% ========================================================================

load("qgp.mat");
load("lab6avg_tz400.mat");

qMap         = QGP.qMap;
phiMap       = QGP.phiMap;
geometry     = QGP.geometryCorr;
polarization = QGP.polarizationCorr;

badMask = lab6_avg <= 0;

% Detector pixels that may contribute to detector-space sums/differences.
detectorValid = ~badMask;

% Hard-coded normalization ROI.
normMask = false(size(qMap));
NORM_ROWS = 1:408;
NORM_COLS = 990:1475;
normMask(NORM_ROWS,NORM_COLS) = true;
normMask(badMask) = false;
normIdx = find(normMask);

% Fixed correction factor.
corrMap = geometry .* polarization;

if CorrectionMode == "divide"
    correctionValid = isfinite(corrMap) & corrMap ~= 0;
elseif CorrectionMode == "multiply"
    correctionValid = isfinite(corrMap);
else
    error('CorrectionMode must be "divide" or "multiply".');
end

% Check detector dimensions.
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
% BENCHMARK BLOCKS
% ========================================================================

groupSize = 2*window;

benchN = floor(BENCH_FILES_REQUESTED/groupSize)*groupSize;
warmN  = floor(WARMUP_FILES_REQUESTED/groupSize)*groupSize;

benchN = max(benchN,groupSize);
warmN  = max(warmN,groupSize);

blockN = benchN + warmN;

if nFiles < 3*blockN
    availablePerTest = floor(nFiles/3/groupSize)*groupSize;
    benchN = floor((availablePerTest-warmN)/groupSize)*groupSize;

    if benchN < groupSize
        error('Not enough files for three independent benchmark blocks.');
    end

    blockN = benchN + warmN;
end

starts = round(linspace(1,nFiles-blockN+1,3));
starts = 1 + floor((starts-1)/groupSize)*groupSize;

for k = 2:3
    if starts(k) < starts(k-1)+blockN
        starts(k) = starts(k-1)+blockN;
        starts(k) = 1 + floor((starts(k)-1)/groupSize)*groupSize;
    end
end

if starts(3)+blockN-1 > nFiles
    starts(3) = nFiles-blockN+1;
    starts(3) = 1 + floor((starts(3)-1)/groupSize)*groupSize;
end

fprintf('Timed/test   : %d TIFFs (~%.2f GB)\n', ...
    benchN,median(fileBytes)*benchN/1e9);
fprintf('Warm-up/test : %d TIFFs\n',warmN);

if isfolder(benchRoot)
    try
        rmdir(benchRoot,'s');
    catch
    end
end
mkdir(benchRoot);

%% ========================================================================
% TEST DEFINITIONS
% ========================================================================

testName = [
    "CPU 1000q x 360phi + qphi(t)"
    "CPU  500q x  90phi + qphi(t)"
    "GPU  500q x  90phi, radial(t) only"
];

NQ             = [1000 500 500];
NPHI           = [ 360  90  90];
USE_GPU        = [false false true];
SAVE_QPHI_TIME = [true  true  false];

elapsed      = nan(3,1);
tiffPerSec   = nan(3,1);
inputMBps    = nan(3,1);
projectedMin = nan(3,1);
outputMB     = nan(3,1);

gpuAvailable = false;

try
    gd = gpuDevice;
    gpuAvailable = true;
    fprintf('GPU          : %s, %.1f GB\n',gd.Name,gd.TotalMemory/1e9);
catch
    fprintf('GPU          : unavailable, TEST 3 will be skipped.\n');
end

%% ========================================================================
% RUN BENCHMARK
% ========================================================================

for t = 1:3

    fprintf('\n------------------------------------------------------------\n');
    fprintf('TEST %d/3: %s\n',t,testName(t));
    fprintf('------------------------------------------------------------\n');

    if USE_GPU(t) && ~gpuAvailable
        fprintf('SKIPPED: no usable MATLAB GPU device.\n');
        continue
    end

    nq   = NQ(t);
    nphi = NPHI(t);

    % ------------------------------------------------------------
    % PRECOMPUTE q / qphi bin membership ONCE for this resolution.
    % ------------------------------------------------------------

    qStaticValid = detectorValid & correctionValid & isfinite(qMap);

    qpStaticValid = qStaticValid & isfinite(phiMap);

    qGood = qMap(qStaticValid);
    pGood = phiMap(qpStaticValid);

    qEdges   = linspace(min(qGood),max(qGood),nq+1);
    phiEdges = linspace(min(pGood),max(pGood),nphi+1);

    qBinMap   = discretize(qMap,qEdges);
    phiBinMap = discretize(phiMap,phiEdges);

    qMask = qStaticValid & ~isnan(qBinMap);

    qpMask = qpStaticValid & ...
             ~isnan(qBinMap) & ...
             ~isnan(phiBinMap);

    qBins = double(qBinMap(qMask));

    qpBins = double(qBinMap(qpMask) + ...
                    (phiBinMap(qpMask)-1)*nq);

    if CorrectionMode == "divide"
        qFactor  = 1 ./ corrMap(qMask);
        qpFactor = 1 ./ corrMap(qpMask);
    else
        qFactor  = corrMap(qMask);
        qpFactor = corrMap(qpMask);
    end

    qFactor  = single(qFactor);
    qpFactor = single(qpFactor);

    qCenter   = (qEdges(1:end-1)+qEdges(2:end))/2;
    phiCenter = (phiEdges(1:end-1)+phiEdges(2:end))/2;

    % ------------------------------------------------------------
    % Independent acquisition block for this test.
    % ------------------------------------------------------------

    first = starts(t);

    warmIdx  = first:first+warmN-1;
    timedIdx = first+warmN:first+warmN+benchN-1;

    warmPaths  = paths(warmIdx);
    timedPaths = paths(timedIdx);

    % ------------------------------------------------------------
    % Initialize architecture OUTSIDE timing.
    % ------------------------------------------------------------

    if ~USE_GPU(t)

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

    else

        gd = gpuDevice;
        wait(gd);

    end

    % ------------------------------------------------------------
    % UNTIMED WARM-UP
    % ------------------------------------------------------------

    fprintf('Warm-up...\n');

    if USE_GPU(t)

        run_gpu_pipeline(warmPaths,window,ny,nx,nq,nphi, ...
            detectorValid,normIdx, ...
            qMask,qBins,qFactor, ...
            qpMask,qpBins,qpFactor, ...
            SAVE_QPHI_TIME(t));

        wait(gpuDevice);

    else

        run_cpu_pipeline(warmPaths,window,ny,nx,nq,nphi,CPU_WORKERS, ...
            detectorValid,normIdx, ...
            qMask,qBins,qFactor, ...
            qpMask,qpBins,qpFactor, ...
            SAVE_QPHI_TIME(t));

    end

    % ------------------------------------------------------------
    % TIMED REAL PIPELINE
    % ------------------------------------------------------------

    fprintf('Timing %d TIFFs...\n',numel(timedPaths));
    drawnow

    t0 = tic;

    if USE_GPU(t)

        R = run_gpu_pipeline(timedPaths,window,ny,nx,nq,nphi, ...
            detectorValid,normIdx, ...
            qMask,qBins,qFactor, ...
            qpMask,qpBins,qpFactor, ...
            SAVE_QPHI_TIME(t));

        wait(gpuDevice);

    else

        R = run_cpu_pipeline(timedPaths,window,ny,nx,nq,nphi,CPU_WORKERS, ...
            detectorValid,normIdx, ...
            qMask,qBins,qFactor, ...
            qpMask,qpBins,qpFactor, ...
            SAVE_QPHI_TIME(t));

    end

    % Include a representative uncompressed disk write.
    R.q   = qCenter;
    R.phi = phiCenter;

    testOut = fullfile(benchRoot,sprintf('test_%d.mat',t));

    save(testOut,'R','-v7','-nocompression');

    elapsed(t) = toc(t0);

    info = dir(testOut);
    outputMB(t) = info.bytes/1e6;

    bytesIn = sum(fileBytes(timedIdx));

    tiffPerSec(t)   = numel(timedPaths)/elapsed(t);
    inputMBps(t)    = bytesIn/1e6/elapsed(t);
    projectedMin(t) = nFiles/tiffPerSec(t)/60;

    fprintf('Elapsed           : %.2f s\n',elapsed(t));
    fprintf('Throughput        : %.2f TIFF/s\n',tiffPerSec(t));
    fprintf('Input throughput  : %.1f MB/s\n',inputMBps(t));
    fprintf('Benchmark output  : %.1f MB\n',outputMB(t));
    fprintf('Projected full run: %.1f min\n',projectedMin(t));

    clear R qBinMap phiBinMap qMask qpMask qBins qpBins qFactor qpFactor

    if DELETE_TEMP_OUTPUTS && isfile(testOut)
        delete(testOut);
    end
end

%% ========================================================================
% RESULTS
% ========================================================================

T = table((1:3)',testName,NQ',NPHI',USE_GPU',SAVE_QPHI_TIME', ...
    elapsed,tiffPerSec,inputMBps,outputMB,projectedMin, ...
    'VariableNames', ...
    {'Test','Configuration','Nq','Nphi','GPU','SaveQPhiTime', ...
     'Seconds','TIFF_per_s','Input_MB_per_s','Output_MB', ...
     'Projected_full_min'});

Tvalid = T(isfinite(T.Seconds),:);
Tvalid = sortrows(Tvalid,'Seconds','ascend');

fprintf('\n\n============================================================\n');
fprintf('FINAL RESULT\n');
fprintf('============================================================\n');
disp(Tvalid);

if ~isempty(Tvalid)

    winner = Tvalid(1,:);

    fprintf('\n*** FASTEST: TEST %d - %s ***\n', ...
        winner.Test,winner.Configuration);

    fprintf('%.2f TIFF/s, projected %.1f min for all %d TIFFs.\n', ...
        winner.TIFF_per_s,winner.Projected_full_min,nFiles);

    if height(Tvalid) >= 2

        gain = 100*(Tvalid.Seconds(2)-Tvalid.Seconds(1)) / ...
                   Tvalid.Seconds(2);

        fprintf('Fastest is %.1f%% quicker than second place.\n',gain);

        if gain < 5
            fprintf(['Top two are within 5%%, so treat them as effectively tied\n' ...
                     'unless a longer benchmark separates them.\n']);
        end
    end
end

fprintf('\nHow to interpret the tests:\n');
fprintf('  TEST 1 vs TEST 2: effect of reducing q/phi resolution.\n');
fprintf(['  TEST 2 vs TEST 3: combined gain from GPU processing and not storing\n' ...
         '                     the full time-resolved q/phi cube.\n']);
fprintf(['  If all Input_MB_per_s values are nearly the same, X: I/O is probably\n' ...
         '  the dominant bottleneck and further compute optimization will help less.\n']);

resultFile = fullfile(benchRoot,'benchmark_result.mat');

save(resultFile,'T','Tvalid','window','CPU_WORKERS');

fprintf('\nResults saved to:\n%s\n',resultFile);
fprintf('============================================================\n');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function R = run_cpu_pipeline(paths,window,ny,nx,nq,nphi,nWorkers, ...
    detectorValid,normIdx, ...
    qMask,qBins,qFactor, ...
    qpMask,qpBins,qpFactor, ...
    saveQPhiTime)

    nF = numel(paths);
    nOut = floor(nF/(2*window));

    deltaIq = nan(nq,nOut,'single');

    if saveQPhiTime
        deltaQPhi = nan(nq,nphi,nOut,'single');
    else
        deltaQPhi = [];
    end

    normCell = cell(nOut,1);

    % Reduction variables: workers combine partial sums locally, then only
    % return the reduced arrays to the client.
    sumAll2D    = zeros(ny,nx,'double');
    sumAllRaw2D = zeros(ny,nx,'double');
    sumDiff2D   = zeros(ny,nx,'double');

    if nWorkers > 1

        parfor g = 1:nOut

            [iq,qpt,normVals,localAll,localRaw,localDiff] = ...
                process_one_group_cpu(paths,g,window,ny,nx,nq,nphi, ...
                    detectorValid,normIdx, ...
                    qMask,qBins,qFactor, ...
                    qpMask,qpBins,qpFactor, ...
                    saveQPhiTime);

            deltaIq(:,g) = iq;

            if saveQPhiTime
                deltaQPhi(:,:,g) = qpt;
            end

            normCell{g} = normVals;

            sumAll2D    = sumAll2D    + localAll;
            sumAllRaw2D = sumAllRaw2D + localRaw;
            sumDiff2D   = sumDiff2D   + localDiff;
        end

    else

        for g = 1:nOut

            [iq,qpt,normVals,localAll,localRaw,localDiff] = ...
                process_one_group_cpu(paths,g,window,ny,nx,nq,nphi, ...
                    detectorValid,normIdx, ...
                    qMask,qBins,qFactor, ...
                    qpMask,qpBins,qpFactor, ...
                    saveQPhiTime);

            deltaIq(:,g) = iq;

            if saveQPhiTime
                deltaQPhi(:,:,g) = qpt;
            end

            normCell{g} = normVals;

            sumAll2D    = sumAll2D    + localAll;
            sumAllRaw2D = sumAllRaw2D + localRaw;
            sumDiff2D   = sumDiff2D   + localDiff;
        end
    end

    normTrace = nan(nF,1);

    for g = 1:nOut
        idx0 = (g-1)*2*window + 1;
        normTrace(idx0:idx0+2*window-1) = normCell{g};
    end

    avgDiff2D = sumDiff2D/(nOut*window);

    % These diagnostics are linear, so do them only ONCE after all frames.
    sumAllIq = binned_mean_cpu(sumAll2D,qMask,qBins,qFactor,nq);

    sumAllQPhi = reshape( ...
        binned_mean_cpu(sumAll2D,qpMask,qpBins,qpFactor,nq*nphi), ...
        nq,nphi);

    avgDiffQPhi = reshape( ...
        binned_mean_cpu(avgDiff2D,qpMask,qpBins,qpFactor,nq*nphi), ...
        nq,nphi);

    sumAllRawIq = binned_mean_cpu(sumAllRaw2D,qMask,qBins,qFactor,nq);

    sumAllRawQPhi = reshape( ...
        binned_mean_cpu(sumAllRaw2D,qpMask,qpBins,qpFactor,nq*nphi), ...
        nq,nphi);

    R.deltaIq       = deltaIq;
    R.deltaQPhi     = deltaQPhi;
    R.normalization = normTrace;

    R.sumAll2D      = sumAll2D;
    R.sumAllRaw2D   = sumAllRaw2D;
    R.avgDiff2D     = avgDiff2D;

    R.sumAllIq      = sumAllIq;
    R.sumAllQPhi    = sumAllQPhi;
    R.avgDiffQPhi   = avgDiffQPhi;

    R.sumAllRawIq   = sumAllRawIq;
    R.sumAllRawQPhi = sumAllRawQPhi;
end


function [iq,qpt,normVals,sumAll2D,sumAllRaw2D,sumDiff2D] = ...
    process_one_group_cpu(paths,g,window,ny,nx,nq,nphi, ...
        detectorValid,normIdx, ...
        qMask,qBins,qFactor, ...
        qpMask,qpBins,qpFactor, ...
        saveQPhiTime)

    first = (g-1)*2*window + 1;

    normVals = nan(2*window,1);

    sumAll2D    = zeros(ny,nx,'double');
    sumAllRaw2D = zeros(ny,nx,'double');
    sumDiff2D   = zeros(ny,nx,'double');

    dWindow = zeros(ny,nx,'single');

    for j = 1:window

        i1 = first + 2*j - 2;
        i2 = i1 + 1;

        raw1 = read_tiff_fast(paths(i1));
        raw2 = read_tiff_fast(paths(i2));

        valid1 = isfinite(raw1) & raw1 >= 0;
        valid2 = isfinite(raw2) & raw2 >= 0;

        % CORRECT ROI INDEXING:
        idx1 = normIdx(valid1(normIdx));
        idx2 = normIdx(valid2(normIdx));

        n1 = sum(double(raw1(idx1)));
        n2 = sum(double(raw2(idx2)));

        if ~(isfinite(n1) && n1 > 0 && isfinite(n2) && n2 > 0)
            error('Invalid normalization scalar in group %d.',g);
        end

        normVals(2*j-1) = n1;
        normVals(2*j)   = n2;

        good1 = detectorValid & valid1;
        good2 = detectorValid & valid2;

        % Raw detector sum for later mask construction.
        r1 = double(raw1);
        r2 = double(raw2);

        r1(~good1) = 0;
        r2(~good2) = 0;

        sumAllRaw2D = sumAllRaw2D + r1 + r2;

        % Normalize BEFORE taking the difference.
        f1 = single(raw1) / single(n1);
        f2 = single(raw2) / single(n2);

        f1(~good1) = NaN;
        f2(~good2) = NaN;

        f1z = f1;
        f2z = f2;

        f1z(~isfinite(f1z)) = 0;
        f2z(~isfinite(f2z)) = 0;

        sumAll2D = sumAll2D + double(f1z) + double(f2z);

        d = f2 - f1;

        dz = d;
        dz(~isfinite(dz)) = 0;

        sumDiff2D = sumDiff2D + double(dz);
        dWindow   = dWindow + dz;
    end

    dWindow = dWindow / single(window);

    iq = binned_mean_cpu(dWindow,qMask,qBins,qFactor,nq);

    if saveQPhiTime
        qpt = reshape( ...
            binned_mean_cpu(dWindow,qpMask,qpBins,qpFactor,nq*nphi), ...
            nq,nphi);
    else
        qpt = [];
    end
end


function R = run_gpu_pipeline(paths,window,ny,nx,nq,nphi, ...
    detectorValid,normIdx, ...
    qMask,qBins,qFactor, ...
    qpMask,qpBins,qpFactor, ...
    saveQPhiTime)

    nF = numel(paths);
    nOut = floor(nF/(2*window));

    % Static reciprocal-space lookup information stays on the GPU.
    gQMask   = gpuArray(qMask);
    gQBins   = gpuArray(qBins);
    gQFactor = gpuArray(qFactor);

    gQPMask   = gpuArray(qpMask);
    gQPBins   = gpuArray(qpBins);
    gQPFactor = gpuArray(qpFactor);

    % Single precision is intentional on RTX 4000 for the image pipeline.
    gSumAll2D  = gpuArray.zeros(ny,nx,'single');
    gSumDiff2D = gpuArray.zeros(ny,nx,'single');

    sumAllRaw2D = zeros(ny,nx,'double');

    deltaIq = nan(nq,nOut,'single');

    if saveQPhiTime
        deltaQPhi = nan(nq,nphi,nOut,'single');
    else
        deltaQPhi = [];
    end

    normTrace = nan(nF,1);

    for g = 1:nOut

        first = (g-1)*2*window + 1;

        gDWindow = gpuArray.zeros(ny,nx,'single');

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
                error('Invalid normalization scalar in GPU group %d.',g);
            end

            normTrace(i1) = n1;
            normTrace(i2) = n2;

            good1 = detectorValid & valid1;
            good2 = detectorValid & valid2;

            % Keep raw summed detector image on CPU.
            r1 = double(raw1);
            r2 = double(raw2);

            r1(~good1) = 0;
            r2(~good2) = 0;

            sumAllRaw2D = sumAllRaw2D + r1 + r2;

            % Normalize/mask before GPU transfer.
            f1 = single(raw1) / single(n1);
            f2 = single(raw2) / single(n2);

            f1(~good1) = NaN;
            f2(~good2) = NaN;

            gF1 = gpuArray(f1);
            gF2 = gpuArray(f2);

            gF1z = gF1;
            gF2z = gF2;

            gF1z(~isfinite(gF1z)) = 0;
            gF2z(~isfinite(gF2z)) = 0;

            gSumAll2D = gSumAll2D + gF1z + gF2z;

            gD = gF2 - gF1;
            gD(~isfinite(gD)) = 0;

            gSumDiff2D = gSumDiff2D + gD;
            gDWindow   = gDWindow + gD;
        end

        gDWindow = gDWindow / single(window);

        gIq = binned_mean_gpu(gDWindow,gQMask,gQBins,gQFactor,nq);
        deltaIq(:,g) = gather(gIq);

        if saveQPhiTime
            gQP = binned_mean_gpu( ...
                gDWindow,gQPMask,gQPBins,gQPFactor,nq*nphi);

            deltaQPhi(:,:,g) = reshape(gather(gQP),nq,nphi);
        end
    end

    gAvgDiff2D = gSumDiff2D / single(nOut*window);

    % Calculate diagnostic reciprocal-space maps only ONCE at the end.
    gSumAllIq = binned_mean_gpu( ...
        gSumAll2D,gQMask,gQBins,gQFactor,nq);

    gSumAllQPhi = binned_mean_gpu( ...
        gSumAll2D,gQPMask,gQPBins,gQPFactor,nq*nphi);

    gAvgDiffQPhi = binned_mean_gpu( ...
        gAvgDiff2D,gQPMask,gQPBins,gQPFactor,nq*nphi);

    % Raw sum is accumulated on CPU. Transfer ONCE for its diagnostics.
    gRaw = gpuArray(single(sumAllRaw2D));

    gSumAllRawIq = binned_mean_gpu( ...
        gRaw,gQMask,gQBins,gQFactor,nq);

    gSumAllRawQPhi = binned_mean_gpu( ...
        gRaw,gQPMask,gQPBins,gQPFactor,nq*nphi);

    R.deltaIq       = deltaIq;
    R.deltaQPhi     = deltaQPhi;
    R.normalization = normTrace;

    R.sumAll2D      = double(gather(gSumAll2D));
    R.sumAllRaw2D   = sumAllRaw2D;
    R.avgDiff2D     = double(gather(gAvgDiff2D));

    R.sumAllIq      = gather(gSumAllIq);
    R.sumAllQPhi    = reshape(gather(gSumAllQPhi),nq,nphi);
    R.avgDiffQPhi   = reshape(gather(gAvgDiffQPhi),nq,nphi);

    R.sumAllRawIq   = gather(gSumAllRawIq);
    R.sumAllRawQPhi = reshape(gather(gSumAllRawQPhi),nq,nphi);
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


function out = binned_mean_gpu(I,gMask,gBins,gFactor,nBins)

    v = I(gMask);

    good = isfinite(v);

    b = gBins(good);
    v = single(v(good)) .* gFactor(good);

    % accumarray supports gpuArray input for @sum.
    s = accumarray(b,v,[nBins 1],@sum,single(0));
    c = accumarray(b,ones(size(v),'single'),[nBins 1],@sum,single(0));

    out = s ./ max(c,single(1));
    out(c==0) = NaN;
end
