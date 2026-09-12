clear all
close all
clc

%% PILATUS PERFORMANCE BENCHMARK - CPU resolution vs GPU architecture
%
% This is a STAND-ALONE benchmark for your current experiment.
%
% It tests three architectures on different short blocks of the acquisition:
%
%   TEST 1  CPU, 1000 q x 360 phi, saves DeltaI(q,phi,t)
%           -> close to your current full-resolution workflow
%
%   TEST 2  CPU,  500 q x  90 phi, saves DeltaI(q,phi,t)
%           -> tests whether coarser q/phi binning materially helps
%
%   TEST 3  GPU,  500 q x  90 phi, DOES NOT save DeltaI(q,phi,t)
%           -> recommended architecture:
%              keep DeltaI(q,t), one summed I(q,phi), one mean DeltaI(q,phi),
%              normalization trace, summed 2D image, mean 2D difference
%
% Important:
% - q and q/phi bin memberships are precomputed ONCE.
% - TIFFs are read only once per benchmark block.
% - normalization is measured before subtraction.
% - bad pixels are masked before normalization/difference.
% - geometry and polarization corrections are applied after the normalized
%   windowed difference and before reciprocal-space integration.
% - each test uses a DIFFERENT acquisition block to reduce cache bias.
% - process-pool/GPU initialization and warm-up are NOT timed.
%
% MATLAB: R2026a
% Parallel Computing Toolbox required for TEST 1/2 if CPU_WORKERS > 1.
% TEST 3 requires a supported NVIDIA GPU.

%% ========================================================================
% USER SETTINGS
% ========================================================================

inDir = "X:\2026-3\RuettSept26\data\fresh_ag_1\run3_UV\";

% Only small temporary benchmark files are written here.
benchRoot = "C:\Users\11idbuser\Downloads\Natan\pilatus_fast_diff_matlab\analysis\benchmark_fast";

window = 3;

% CPU worker count to use in TEST 1 and TEST 2.
% If your earlier worker benchmark found another winner, change this ONE line.
CPU_WORKERS = 2;

% 180 TIFFs = 30 output windows when window=3.
% Increase to 360 if TEST 1/2 are too short/noisy.
BENCH_FILES_REQUESTED = 180;

% Small untimed warm-up before each timed test.
WARMUP_FILES_REQUESTED = 12;

% Geometry/polarization convention:
% "divide"   -> corrected = image ./ (geometry .* polarization)
% "multiply" -> corrected = image .* (geometry .* polarization)
CorrectionMode = "divide";

% Remove temporary benchmark outputs after measuring them.
DELETE_TEMP_OUTPUTS = true;

%% ========================================================================
% LOAD FILE LIST - numeric sort from filename
% ========================================================================

D = dir(fullfile(inDir,"*.tif"));
if isempty(D)
    error('No TIFF files found in %s',inDir);
end

names = string({D.name})';
paths = fullfile(string({D.folder})',names);
fileBytes = double([D.bytes])';

% Numeric key = last integer appearing in each filename.
ids = nan(numel(names),1);
for k = 1:numel(names)
    tok = regexp(names(k),'\d+','match');
    if isempty(tok)
        error('Could not find a numeric frame ID in filename: %s',names(k));
    end
    ids(k) = str2double(tok{end});
end

[ids,ord] = sort(ids);
paths = paths(ord);
names = names(ord);
fileBytes = fileBytes(ord);

nFiles = numel(paths);

fprintf('\n============================================================\n');
fprintf('PILATUS ARCHITECTURE BENCHMARK\n');
fprintf('============================================================\n');
fprintf('Files found : %d\n',nFiles);
fprintf('Input size  : %.3f TB\n',sum(fileBytes)/1e12);
fprintf('Frame IDs   : %d ... %d\n',ids(1),ids(end));
fprintf('Window      : %d (%d TIFFs/output)\n',window,2*window);
fprintf('CPU workers : %d\n',CPU_WORKERS);

%% ========================================================================
% CALIBRATION / MASK / NORMALIZATION
% ========================================================================

load("qgp.mat");
load("lab6avg_tz400.mat");

qMap         = QGP.qMap;
phiMap       = QGP.phiMap;
geometry     = QGP.geometryCorr;
polarization = QGP.polarizationCorr;

badMask = lab6_avg <= 0;

normMask = false(size(qMap));
NORM_ROWS = 1:408;
NORM_COLS = 990:1475;
normMask(NORM_ROWS,NORM_COLS) = true;
normMask(badMask) = false;

% Static reciprocal-space-valid pixels.
staticValid = ~badMask & isfinite(qMap) & isfinite(phiMap) & ...
              isfinite(geometry) & isfinite(polarization);

if CorrectionMode == "divide"
    staticValid = staticValid & geometry ~= 0 & polarization ~= 0;
    corrMap = geometry .* polarization;
elseif CorrectionMode == "multiply"
    corrMap = geometry .* polarization;
else
    error('CorrectionMode must be "divide" or "multiply".');
end

% Check detector dimensions against first TIFF.
A0 = read_tiff_fast(paths(1));
if ~isequal(size(A0),size(qMap))
    error('TIFF is %dx%d but qMap is %dx%d.', ...
        size(A0,1),size(A0,2),size(qMap,1),size(qMap,2));
end
clear A0

normIdx = find(normMask);

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
    benchN = availablePerTest - warmN;
    benchN = floor(benchN/groupSize)*groupSize;
    if benchN < groupSize
        error('Not enough files for three independent benchmark blocks.');
    end
    blockN = benchN + warmN;
end

% Spread three blocks through the run.
starts = round(linspace(1,nFiles-blockN+1,3));
starts = 1 + floor((starts-1)/groupSize)*groupSize;

% Make sure blocks do not overlap.
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

fprintf('Timed TIFFs/test: %d (~%.2f GB)\n', ...
    benchN,median(fileBytes)*benchN/1e9);
fprintf('Warm-up/test    : %d TIFFs\n',warmN);

if isfolder(benchRoot)
    try
        rmdir(benchRoot,'s');
    catch
    end
end
mkdir(benchRoot);

%% ========================================================================
% DEFINE THREE TESTS
% ========================================================================

testName = [
    "CPU 1000q x 360phi + qphi(t)"
    "CPU  500q x  90phi + qphi(t)"
    "GPU  500q x  90phi, no qphi(t)"
];

NQ   = [1000 500 500];
NPHI = [ 360  90  90];
USE_GPU = [false false true];
SAVE_QPHI_T = [true true false];

elapsed = nan(3,1);
tiffPerSec = nan(3,1);
inputMBps = nan(3,1);
projectedMin = nan(3,1);
outputMB = nan(3,1);
gpuAvailable = false;

try
    gd = gpuDevice;
    gpuAvailable = true;
    fprintf('GPU          : %s, %.1f GB\n',gd.Name,gd.TotalMemory/1e9);
catch
    fprintf('GPU          : unavailable - TEST 3 will be skipped.\n');
end

%% ========================================================================
% RUN TESTS
% ========================================================================

for t = 1:3

    fprintf('\n------------------------------------------------------------\n');
    fprintf('TEST %d/3: %s\n',t,testName(t));
    fprintf('------------------------------------------------------------\n');

    if USE_GPU(t) && ~gpuAvailable
        fprintf('SKIPPED: no usable MATLAB GPU device.\n');
        continue
    end

    nq = NQ(t);
    nphi = NPHI(t);

    qGood = qMap(staticValid);
    phiGood = phiMap(staticValid);

    qEdges = linspace(min(qGood),max(qGood),nq+1);
    phiEdges = linspace(min(phiGood),max(phiGood),nphi+1);

    % PRECOMPUTE bin memberships ONCE.
    qBin = discretize(qMap,qEdges);
    phiBin = discretize(phiMap,phiEdges);

    qphiBin = nan(size(qBin));
    okBin = staticValid & ~isnan(qBin) & ~isnan(phiBin);
    qphiBin(okBin) = qBin(okBin) + (phiBin(okBin)-1)*nq;

    qIdxStatic = find(staticValid & ~isnan(qBin));
    qBinStatic = uint32(qBin(qIdxStatic));

    qpIdxStatic = find(okBin);
    qpBinStatic = uint32(qphiBin(qpIdxStatic));

    qCenter = (qEdges(1:end-1)+qEdges(2:end))/2;
    phiCenter = (phiEdges(1:end-1)+phiEdges(2:end))/2;

    first = starts(t);
    warmIdx  = first:first+warmN-1;
    timedIdx = first+warmN:first+warmN+benchN-1;

    warmPaths = paths(warmIdx);
    timedPaths = paths(timedIdx);

    % Pool setup outside timing.
    if ~USE_GPU(t)
        p = gcp('nocreate');
        if CPU_WORKERS <= 1
            if ~isempty(p), delete(p); end
        else
            if isempty(p) || p.NumWorkers ~= CPU_WORKERS
                if ~isempty(p), delete(p); end
                fprintf('Starting %d-process pool (not timed)...\n',CPU_WORKERS);
                parpool("Processes",CPU_WORKERS);
            end
        end
    else
        % GPU warm initialization outside timing.
        gd = gpuDevice;
        wait(gd);
    end

    % ---------------- untimed warm-up ----------------
    fprintf('Warm-up...\n');
    if USE_GPU(t)
        run_gpu_pipeline(warmPaths,window,nq,nphi, ...
            normIdx,staticValid,corrMap,CorrectionMode, ...
            qIdxStatic,qBinStatic,qpIdxStatic,qpBinStatic,false);
        wait(gpuDevice);
    else
        run_cpu_pipeline(warmPaths,window,nq,nphi,CPU_WORKERS, ...
            normIdx,staticValid,corrMap,CorrectionMode, ...
            qIdxStatic,qBinStatic,qpIdxStatic,qpBinStatic,false);
    end

    % ---------------- timed pipeline ----------------
    fprintf('Timing %d TIFFs...\n',numel(timedPaths));
    drawnow

    t0 = tic;

    if USE_GPU(t)
        R = run_gpu_pipeline(timedPaths,window,nq,nphi, ...
            normIdx,staticValid,corrMap,CorrectionMode, ...
            qIdxStatic,qBinStatic,qpIdxStatic,qpBinStatic,SAVE_QPHI_T(t));
        wait(gpuDevice);
    else
        R = run_cpu_pipeline(timedPaths,window,nq,nphi,CPU_WORKERS, ...
            normIdx,staticValid,corrMap,CorrectionMode, ...
            qIdxStatic,qBinStatic,qpIdxStatic,qpBinStatic,SAVE_QPHI_T(t));
    end

    % Include realistic output write in timing.
    testOut = fullfile(benchRoot,sprintf('test_%d.mat',t));

    R.q = qCenter;
    R.phi = phiCenter;
    save(testOut,'R','-v7','-nocompression');

    elapsed(t) = toc(t0);

    info = dir(testOut);
    outputMB(t) = info.bytes/1e6;

    bytesIn = sum(fileBytes(timedIdx));
    tiffPerSec(t) = numel(timedPaths)/elapsed(t);
    inputMBps(t) = bytesIn/1e6/elapsed(t);
    projectedMin(t) = nFiles/tiffPerSec(t)/60;

    fprintf('Elapsed           : %.2f s\n',elapsed(t));
    fprintf('Throughput        : %.2f TIFF/s\n',tiffPerSec(t));
    fprintf('Input throughput  : %.1f MB/s\n',inputMBps(t));
    fprintf('Benchmark output  : %.1f MB\n',outputMB(t));
    fprintf('Projected full run: %.1f min\n',projectedMin(t));

    clear R qBin phiBin qphiBin qIdxStatic qBinStatic qpIdxStatic qpBinStatic

    if DELETE_TEMP_OUTPUTS && isfile(testOut)
        delete(testOut);
    end
end

%% ========================================================================
% RESULT
% ========================================================================

T = table((1:3)',testName,NQ',NPHI',USE_GPU',SAVE_QPHI_T', ...
    elapsed,tiffPerSec,inputMBps,outputMB,projectedMin, ...
    'VariableNames',{'Test','Configuration','Nq','Nphi','GPU','SaveQPhiTime', ...
    'Seconds','TIFF_per_s','Input_MB_per_s','Output_MB','Projected_full_min'});

Tvalid = T(isfinite(T.Seconds),:);
Tvalid = sortrows(Tvalid,'Seconds','ascend');

fprintf('\n\n============================================================\n');
fprintf('FINAL RESULT\n');
fprintf('============================================================\n');
disp(Tvalid);

if ~isempty(Tvalid)
    winner = Tvalid(1,:);
    fprintf('\n*** FASTEST: TEST %d - %s ***\n',winner.Test,winner.Configuration);
    fprintf('%.2f TIFF/s, projected %.1f min for all %d TIFFs.\n', ...
        winner.TIFF_per_s,winner.Projected_full_min,nFiles);

    if height(Tvalid) >= 2
        gain = 100*(Tvalid.Seconds(2)-Tvalid.Seconds(1))/Tvalid.Seconds(2);
        fprintf('Fastest is %.1f%% quicker than second place in this test.\n',gain);
    end
end

fprintf('\nInterpretation:\n');
fprintf(' TEST 1 vs TEST 2 = benefit from reducing q/phi resolution.\n');
fprintf(' TEST 2 vs TEST 3 = benefit from GPU + not storing qphi(t).\n');
fprintf([' If TEST 3 wins, use GPU for the reciprocal-space accumulation and keep only\n' ...
         ' DeltaI(q,t), summed I(q,phi), and mean DeltaI(q,phi) during production.\n']);

resultFile = fullfile(benchRoot,'benchmark_result.mat');
save(resultFile,'T','Tvalid','window','CPU_WORKERS');
fprintf('\nResults saved to:\n%s\n',resultFile);
fprintf('============================================================\n');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function R = run_cpu_pipeline(paths,window,nq,nphi,nWorkers, ...
    normIdx,staticValid,corrMap,CorrectionMode, ...
    qIdxStatic,qBinStatic,qpIdxStatic,qpBinStatic,saveQPhiT)

    nF = numel(paths);
    groupSize = 2*window;
    nOut = floor(nF/groupSize);

    A = read_tiff_fast(paths(1));
    [ny,nx] = size(A);
    clear A

    deltaIq = nan(nq,nOut,'single');
    if saveQPhiT
        deltaQPhi = nan(nq,nphi,nOut,'single');
    else
        deltaQPhi = [];
    end

    normTrace = nan(nF,1);

    % Per-window outputs for reduction after parfor. This avoids transmitting
    % detector-sized differences back to the client.
    sumAllLocal = cell(nOut,1);
    sumDiffLocal = cell(nOut,1);
    sumQPhiAllLocal = cell(nOut,1);
    sumQPhiDiffLocal = cell(nOut,1);
    normLocal = cell(nOut,1);

    if nWorkers <= 1
        loopType = 0;
    else
        loopType = 1;
    end

    if loopType == 1
        parfor g = 1:nOut
            [deltaIq(:,g),qpt,sumAllLocal{g},sumDiffLocal{g}, ...
                sumQPhiAllLocal{g},sumQPhiDiffLocal{g},normLocal{g}] = ...
                process_one_group_cpu(paths,g,window,nq,nphi,normIdx, ...
                staticValid,corrMap,CorrectionMode,qIdxStatic,qBinStatic, ...
                qpIdxStatic,qpBinStatic,saveQPhiT);
            if saveQPhiT
                deltaQPhi(:,:,g) = qpt;
            end
        end
    else
        for g = 1:nOut
            [deltaIq(:,g),qpt,sumAllLocal{g},sumDiffLocal{g}, ...
                sumQPhiAllLocal{g},sumQPhiDiffLocal{g},normLocal{g}] = ...
                process_one_group_cpu(paths,g,window,nq,nphi,normIdx, ...
                staticValid,corrMap,CorrectionMode,qIdxStatic,qBinStatic, ...
                qpIdxStatic,qpBinStatic,saveQPhiT);
            if saveQPhiT
                deltaQPhi(:,:,g) = qpt;
            end
        end
    end

    sumAll2D = zeros(ny,nx,'double');
    sumDiff2D = zeros(ny,nx,'double');
    sumAllQPhi = zeros(nq,nphi,'double');
    sumDiffQPhi = zeros(nq,nphi,'double');

    for g = 1:nOut
        sumAll2D = sumAll2D + sumAllLocal{g};
        sumDiff2D = sumDiff2D + sumDiffLocal{g};
        sumAllQPhi = sumAllQPhi + sumQPhiAllLocal{g};
        sumDiffQPhi = sumDiffQPhi + sumQPhiDiffLocal{g};

        idx0 = (g-1)*2*window + 1;
        normTrace(idx0:idx0+2*window-1) = normLocal{g};
    end

    R.deltaIq = deltaIq;
    R.deltaQPhi = deltaQPhi;
    R.normalization = normTrace;
    R.sumAll2D = sumAll2D;
    R.avgDiff2D = sumDiff2D/(nOut*window);
    R.sumAllQPhi = single(sumAllQPhi);
    R.avgDiffQPhi = single(sumDiffQPhi/(nOut*window));
end

function [iq,qpt,sumAll2D,sumDiff2D,sumAllQPhi,sumDiffQPhi,normVals] = ...
    process_one_group_cpu(paths,g,window,nq,nphi,normIdx,staticValid,corrMap, ...
    CorrectionMode,qIdxStatic,qBinStatic,qpIdxStatic,qpBinStatic,saveQPhiT)

    first = (g-1)*2*window + 1;
    A0 = read_tiff_fast(paths(first));
    [ny,nx] = size(A0);
    clear A0

    sumAll2D = zeros(ny,nx,'double');
    sumDiff2D = zeros(ny,nx,'double');
    sumAllQPhi = zeros(nq,nphi,'double');
    sumDiffQPhi = zeros(nq,nphi,'double');
    normVals = nan(2*window,1);

    dWindow = zeros(ny,nx,'single');

    for j = 1:window
        i1 = first + 2*j - 2;
        i2 = i1 + 1;

        raw1 = read_tiff_fast(paths(i1));
        raw2 = read_tiff_fast(paths(i2));

        valid1 = ~isnan(raw1) & raw1 >= 0;
        valid2 = ~isnan(raw2) & raw2 >= 0;

        n1 = sum(double(raw1(normIdx & valid1(:))));
        n2 = sum(double(raw2(normIdx & valid2(:))));

        if ~(isfinite(n1) && n1 > 0 && isfinite(n2) && n2 > 0)
            error('Invalid normalization scalar.');
        end

        normVals(2*j-1) = n1;
        normVals(2*j) = n2;

        f1 = single(raw1) / single(n1);
        f2 = single(raw2) / single(n2);

        f1(~staticValid | ~valid1) = NaN;
        f2(~staticValid | ~valid2) = NaN;

        sumAll2D = sumAll2D + double(fillmissing_zero(f1)) + double(fillmissing_zero(f2));

        d = f2 - f1;
        dWindow = dWindow + fillmissing_zero(d);

        % Accumulate all-image qphi diagnostic after correction.
        c1 = apply_corr_cpu(f1,corrMap,CorrectionMode);
        c2 = apply_corr_cpu(f2,corrMap,CorrectionMode);
        sumAllQPhi = sumAllQPhi + ...
            qphi_sum_cpu(c1,nq,nphi,qpIdxStatic,qpBinStatic) + ...
            qphi_sum_cpu(c2,nq,nphi,qpIdxStatic,qpBinStatic);

        raw1 = []; raw2 = []; f1 = []; f2 = [];
    end

    dWindow = dWindow / single(window);
    sumDiff2D = double(dWindow) * window;

    dCorr = apply_corr_cpu(dWindow,corrMap,CorrectionMode);

    iq = q_mean_cpu(dCorr,nq,qIdxStatic,qBinStatic);
    qpt = [];
    if saveQPhiT
        qpt = qphi_mean_cpu(dCorr,nq,nphi,qpIdxStatic,qpBinStatic);
    end

    % mean-difference diagnostic accumulated as SUM of pair differences
    sumDiffQPhi = qphi_sum_cpu(dCorr,nq,nphi,qpIdxStatic,qpBinStatic) * window;
end

function R = run_gpu_pipeline(paths,window,nq,nphi, ...
    normIdx,staticValid,corrMap,CorrectionMode, ...
    qIdxStatic,qBinStatic,qpIdxStatic,qpBinStatic,saveQPhiT)

    nF = numel(paths);
    groupSize = 2*window;
    nOut = floor(nF/groupSize);

    A = read_tiff_fast(paths(1));
    [ny,nx] = size(A);
    clear A

    % Static arrays stay on GPU.
    gStaticValid = gpuArray(staticValid);
    gCorr = gpuArray(single(corrMap));

    gQIdx = gpuArray(uint32(qIdxStatic));
    gQBin = gpuArray(uint32(qBinStatic));

    gQPIdx = gpuArray(uint32(qpIdxStatic));
    gQPBin = gpuArray(uint32(qpBinStatic));

    gSumAll2D = gpuArray.zeros(ny,nx,'double');
    gSumDiff2D = gpuArray.zeros(ny,nx,'double');
    gSumAllQPhi = gpuArray.zeros(nq*nphi,1,'double');
    gSumDiffQPhi = gpuArray.zeros(nq*nphi,1,'double');

    deltaIq = nan(nq,nOut,'single');

    if saveQPhiT
        deltaQPhi = nan(nq,nphi,nOut,'single');
    else
        deltaQPhi = [];
    end

    normTrace = nan(nF,1);

    for g = 1:nOut
        first = (g-1)*groupSize + 1;
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
                error('Invalid normalization scalar.');
            end

            normTrace(i1) = n1;
            normTrace(i2) = n2;

            gF1 = gpuArray(single(raw1)) / single(n1);
            gF2 = gpuArray(single(raw2)) / single(n2);

            gV1 = gpuArray(valid1) & gStaticValid;
            gV2 = gpuArray(valid2) & gStaticValid;

            gF1(~gV1) = NaN;
            gF2(~gV2) = NaN;

            gF1z = gF1; gF1z(~isfinite(gF1z)) = 0;
            gF2z = gF2; gF2z(~isfinite(gF2z)) = 0;

            gSumAll2D = gSumAll2D + double(gF1z) + double(gF2z);

            gD = gF2 - gF1;
            gD(~isfinite(gD)) = 0;
            gDWindow = gDWindow + gD;

            gC1 = apply_corr_gpu(gF1,gCorr,CorrectionMode);
            gC2 = apply_corr_gpu(gF2,gCorr,CorrectionMode);

            gSumAllQPhi = gSumAllQPhi + ...
                bin_sum_gpu(gC1,gQPIdx,gQPBin,nq*nphi) + ...
                bin_sum_gpu(gC2,gQPIdx,gQPBin,nq*nphi);
        end

        gDWindow = gDWindow / single(window);
        gSumDiff2D = gSumDiff2D + double(gDWindow) * window;

        gDCorr = apply_corr_gpu(gDWindow,gCorr,CorrectionMode);

        giq = bin_mean_gpu(gDCorr,gQIdx,gQBin,nq);
        deltaIq(:,g) = gather(single(giq));

        if saveQPhiT
            gqp = bin_mean_gpu(gDCorr,gQPIdx,gQPBin,nq*nphi);
            deltaQPhi(:,:,g) = reshape(gather(single(gqp)),nq,nphi);
        end

        gSumDiffQPhi = gSumDiffQPhi + ...
            bin_sum_gpu(gDCorr,gQPIdx,gQPBin,nq*nphi) * window;
    end

    R.deltaIq = deltaIq;
    R.deltaQPhi = deltaQPhi;
    R.normalization = normTrace;
    R.sumAll2D = gather(gSumAll2D);
    R.avgDiff2D = gather(gSumDiff2D)/(nOut*window);
    R.sumAllQPhi = reshape(single(gather(gSumAllQPhi)),nq,nphi);
    R.avgDiffQPhi = reshape(single(gather(gSumDiffQPhi))/(nOut*window),nq,nphi);
end

function A = read_tiff_fast(filename)
    % Tiff class avoids some general imread dispatch overhead for simple TIFFs.
    t = Tiff(char(filename),'r');
    A = t.read();
    t.close();
end

function Y = fillmissing_zero(X)
    Y = X;
    Y(~isfinite(Y)) = 0;
end

function C = apply_corr_cpu(I,corrMap,mode)
    if mode == "divide"
        C = I ./ single(corrMap);
    else
        C = I .* single(corrMap);
    end
end

function C = apply_corr_gpu(I,gCorr,mode)
    if mode == "divide"
        C = I ./ gCorr;
    else
        C = I .* gCorr;
    end
end

function iq = q_mean_cpu(I,nq,pixIdx,binIdx)
    v = I(pixIdx);
    good = isfinite(v);
    b = double(binIdx(good));
    v = double(v(good));
    s = accumarray(b,v,[nq 1],@sum,0);
    c = accumarray(b,1,[nq 1],@sum,0);
    iq = single(s ./ max(c,1));
    iq(c==0) = NaN;
end

function out = qphi_mean_cpu(I,nq,nphi,pixIdx,binIdx)
    v = I(pixIdx);
    good = isfinite(v);
    b = double(binIdx(good));
    v = double(v(good));
    nb = nq*nphi;
    s = accumarray(b,v,[nb 1],@sum,0);
    c = accumarray(b,1,[nb 1],@sum,0);
    z = s ./ max(c,1);
    z(c==0) = NaN;
    out = reshape(single(z),nq,nphi);
end

function out = qphi_sum_cpu(I,nq,nphi,pixIdx,binIdx)
    v = I(pixIdx);
    good = isfinite(v);
    b = double(binIdx(good));
    v = double(v(good));
    z = accumarray(b,v,[nq*nphi 1],@sum,0);
    out = reshape(z,nq,nphi);
end

function s = bin_sum_gpu(I,gPixIdx,gBinIdx,nBins)
    v = I(gPixIdx);
    good = isfinite(v);
    b = double(gBinIdx(good));
    v = double(v(good));
    s = accumarray(b,v,[nBins 1],@sum,0);
end

function m = bin_mean_gpu(I,gPixIdx,gBinIdx,nBins)
    v = I(gPixIdx);
    good = isfinite(v);
    b = double(gBinIdx(good));
    v = single(v(good));

    s = accumarray(b,v,[nBins 1],@sum,single(0));
    c = accumarray(b,single(1),[nBins 1],@sum,single(0));

    m = s ./ max(c,single(1));
    m(c==0) = NaN;
end
