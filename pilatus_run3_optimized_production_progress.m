clearvars
close all
clc

%% PILATUS PRODUCTION ANALYSIS - OPTIMIZED FROM FULL BENCHMARK
%
% Benchmarked winner on this workstation / storage:
%   CPU
%   Tiff reader
%   4 process workers
%   sparse radial projection
%   window = 3
%   Nq = 1000
%   process directly from X:
%
% Production outputs:
%   R.avgDiff2D                  one final detector-space average difference
%   R.deltaIq                    [Nq x Nwindows] radial difference data
%   R.q                          q-bin centers
%   R.normalization.value        raw normalization ROI sum for every TIFF
%   R.normalization.scale        normalization multiplier for every TIFF
%   R.normalization.frameID      source frame IDs
%   R.windowFrameIDStart/End     source-frame range for each DeltaI(q) column
%
% Processing:
%
%   raw TIFF
%      -> permanent detector mask + invalid-pixel mask
%      -> normalization scalar from hard-coded ROI
%      -> normalize each frame
%      -> I_even - I_odd
%      -> average WINDOW pair differences
%      -> geometry/polarization correction
%      -> sparse radial projection
%
% Dynamic invalid pixels are handled with a per-pixel validity count rather
% than silently treated as zeros.
%
% MATLAB R2026a + Parallel Computing Toolbox.

%% ========================================================================
% USER SETTINGS
% ========================================================================

inDir = "X:\2026-3\RuettSept26\data\fresh_ag_1\run3_UV\";

outDir = ...
    "C:\Users\11idbuser\Downloads\Natan\pilatus_fast_diff_matlab\analysis\run003_optimized";

% Science parameter.
% Benchmark winner was WINDOW = 3.
WINDOW = 3;

% Benchmark winner was 1000.
NQ = 1000;

% Benchmark winner was 4 process workers.
N_WORKERS = 4;

% Normalization:
% normalized frame = raw frame * NORM_REFERENCE / normalization_ROI_sum
%
% Keep 1 for the same convention used in the benchmark.
NORM_REFERENCE = 1;

% Geometry / polarization convention used in the benchmark.
% "divide"   -> Icorr = I ./ (geometry .* polarization)
% "multiply" -> Icorr = I .* (geometry .* polarization)
CORRECTION_MODE = "divide";

% Output filename.
outputFile = fullfile(outDir, ...
    sprintf("run003_W%d_Nq%d_optimized.mat",WINDOW,NQ));

%% ========================================================================
% INPUT FILE LIST - NUMERIC FRAME ORDER
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
        error('Could not find a numeric frame ID in filename: %s',names(k));
    end

    frameIDs(k) = str2double(tok{end});
end

[frameIDs,ord] = sort(frameIDs);

paths     = paths(ord);
names     = names(ord);
fileBytes = fileBytes(ord);

nFiles = numel(paths);

% Pairing depends on exact acquisition order.  Fail rather than silently
% shifting the even/odd pairing if a TIFF is missing or duplicated.
if any(diff(frameIDs) ~= 1)
    bad = find(diff(frameIDs) ~= 1,1,'first');

    error(['Frame IDs are not consecutive near %d -> %d. ' ...
           'Fix missing/duplicate TIFFs before processing.'], ...
           frameIDs(bad),frameIDs(bad+1));
end

fprintf('\n============================================================\n');
fprintf('PILATUS OPTIMIZED PRODUCTION ANALYSIS\n');
fprintf('============================================================\n');
fprintf('Input files     : %d\n',nFiles);
fprintf('Input size      : %.3f TB\n',sum(fileBytes)/1e12);
fprintf('Frame IDs       : %d ... %d\n',frameIDs(1),frameIDs(end));
fprintf('Window          : %d (%d TIFFs/output)\n',WINDOW,2*WINDOW);
fprintf('q bins          : %d\n',NQ);
fprintf('Process workers : %d\n',N_WORKERS);

%% ========================================================================
% LOAD CALIBRATION
% ========================================================================

load("qgp.mat");
load("lab6avg_tz400.mat");

qMap         = QGP.qMap;
geometry     = QGP.geometryCorr;
polarization = QGP.polarizationCorr;

%% ========================================================================
% STATIC DETECTOR MASK
% ========================================================================

% TRUE means permanently invalid.
badMask = lab6_avg <= 0;

detectorValid = ~badMask;

%% ========================================================================
% HARD-CODED NORMALIZATION ROI
% ========================================================================

normMask = false(size(qMap));

NORM_ROWS = 1:408;
NORM_COLS = 990:1475;

normMask(NORM_ROWS,NORM_COLS) = true;

% Never use permanently bad detector pixels for normalization.
normMask(badMask) = false;

normIdx = find(normMask);

%% ========================================================================
% CHECK DETECTOR / CALIBRATION DIMENSIONS
% ========================================================================

A0 = read_tiff_fast(paths(1));

[ny,nx] = size(A0);

if ~isequal(size(A0),size(qMap),size(geometry),size(polarization),size(badMask))
    error(['TIFF, qMap, geometry, polarization and badMask must all have ' ...
           'exactly the same detector dimensions.']);
end

clear A0

fprintf('Detector         : %d x %d\n',ny,nx);
fprintf('Normalization ROI: %d good pixels\n',numel(normIdx));

%% ========================================================================
% PRECOMPUTE OPTIMIZED SPARSE RADIAL OPERATOR
% ========================================================================

corrMap = geometry .* polarization;

if CORRECTION_MODE == "divide"

    qStaticValid = detectorValid & ...
                   isfinite(qMap) & ...
                   isfinite(corrMap) & ...
                   corrMap ~= 0;

elseif CORRECTION_MODE == "multiply"

    qStaticValid = detectorValid & ...
                   isfinite(qMap) & ...
                   isfinite(corrMap);

else

    error('CORRECTION_MODE must be "divide" or "multiply".');
end

qGood = qMap(qStaticValid);

qEdges  = linspace(min(qGood),max(qGood),NQ+1);
qCenter = (qEdges(1:end-1) + qEdges(2:end))/2;

qBinMap = discretize(qMap,qEdges);

qMask = qStaticValid & ~isnan(qBinMap);

qBins = double(qBinMap(qMask));

if CORRECTION_MODE == "divide"
    qFactor = single(1 ./ corrMap(qMask));
else
    qFactor = single(corrMap(qMask));
end

nRadialPixels = numel(qBins);

% Sparse projection matrix: one non-zero per valid detector pixel.
% Rows = q bins, columns = valid detector pixels.
Srad = sparse( ...
    qBins, ...
    (1:nRadialPixels)', ...
    ones(nRadialPixels,1), ...
    NQ, ...
    nRadialPixels);

fprintf('Radial pixels    : %d\n',nRadialPixels);

% Operator packaged once and broadcast once to each process worker.
radialOp = struct;
radialOp.mask    = qMask;
radialOp.factor  = qFactor;
radialOp.S       = Srad;
radialOp.Nq      = NQ;
radialOp.qCenter = qCenter;

clear qBinMap qBins qGood Srad

%% ========================================================================
% COMPLETE WINDOWS
% ========================================================================

groupSize = 2*WINDOW;

nOut = floor(nFiles/groupSize);

nUsed = nOut*groupSize;
nTail = nFiles-nUsed;

if nOut < 1
    error('Not enough TIFFs for WINDOW=%d.',WINDOW);
end

if nTail > 0
    warning('%d TIFF(s) at the end do not make a complete window and will be ignored.',nTail);
end

usedPaths    = paths(1:nUsed);
usedFrameIDs = frameIDs(1:nUsed);

fprintf('Output curves    : %d\n',nOut);
fprintf('TIFFs used       : %d\n',nUsed);

%% ========================================================================
% START BENCHMARKED 4-PROCESS POOL
% ========================================================================

p = gcp('nocreate');

if isempty(p) || p.NumWorkers ~= N_WORKERS

    if ~isempty(p)
        delete(p);
    end

    fprintf('Starting %d-process pool...\n',N_WORKERS);

    parpool("Processes",N_WORKERS);
end

%% ========================================================================
% PREALLOCATE OUTPUTS
% ========================================================================

deltaIq = nan(NQ,nOut,'single');

% One small normalization vector per output group.  Using cells avoids
% workers writing overlapping sections of one shared vector.
normCell = cell(nOut,1);

% parfor reduction variables.
%
% diffSumAll:
%   sum of every valid normalized pair difference
%
% diffCountAll:
%   number of valid pair differences contributing to each detector pixel
%
% Final:
%   avgDiff2D = diffSumAll ./ diffCountAll
diffSumAll   = zeros(ny,nx,'double');
diffCountAll = zeros(ny,nx,'double');

%% ========================================================================
% ONE-PASS PARALLEL ANALYSIS
% ========================================================================

fprintf('\nProcessing...\n');

% Text progress / ETA for parfor
progressQ = parallel.pool.DataQueue;

ticRun = tic;

progress_text("reset",nOut,nUsed,groupSize,ticRun);

afterEach(progressQ, ...
    @(msg) progress_text(msg,nOut,nUsed,groupSize,ticRun));

parfor g = 1:nOut

    first = (g-1)*groupSize + 1;

    normVals = nan(groupSize,1);

    % Detector-space running sums for this output window.
    windowSum   = zeros(ny,nx,'single');
    windowCount = zeros(ny,nx,'single');

    % Reduction contributions for the final all-run average difference.
    localDiffSum   = zeros(ny,nx,'double');
    localDiffCount = zeros(ny,nx,'double');

    for j = 1:WINDOW

        i1 = first + 2*j - 2;
        i2 = i1 + 1;

        % -----------------------------------------------------------------
        % READ PAIR
        % -----------------------------------------------------------------

        raw1 = read_tiff_fast(usedPaths(i1));
        raw2 = read_tiff_fast(usedPaths(i2));

        % Pilatus invalid pixels may be negative.
        valid1 = isfinite(raw1) & raw1 >= 0;
        valid2 = isfinite(raw2) & raw2 >= 0;

        % -----------------------------------------------------------------
        % NORMALIZATION SCALAR - BEFORE DIFFERENCE
        % -----------------------------------------------------------------

        idx1 = normIdx(valid1(normIdx));
        idx2 = normIdx(valid2(normIdx));

        n1 = sum(double(raw1(idx1)));
        n2 = sum(double(raw2(idx2)));

        if ~(isfinite(n1) && n1 > 0 && isfinite(n2) && n2 > 0)
            error('Invalid normalization scalar in output group %d.',g);
        end

        normVals(2*j-1) = n1;
        normVals(2*j)   = n2;

        % Pixel contributes only when valid in BOTH members of this pair.
        pairValid = detectorValid & valid1 & valid2;

        % -----------------------------------------------------------------
        % NORMALIZED PAIR DIFFERENCE
        %
        % Avoid creating two full normalized detector arrays.  The single
        % expression reduces memory traffic and temporary allocations.
        % -----------------------------------------------------------------

        d = ...
            single(raw2) * single(NORM_REFERENCE/n2) - ...
            single(raw1) * single(NORM_REFERENCE/n1);

        d(~pairValid) = 0;

        % Final detector-space average difference.
        localDiffSum   = localDiffSum   + double(d);
        localDiffCount = localDiffCount + double(pairValid);

        % Windowed detector difference for radial integration.
        windowSum   = windowSum   + d;
        windowCount = windowCount + single(pairValid);
    end

    % ---------------------------------------------------------------------
    % WINDOW AVERAGE IN DETECTOR SPACE
    %
    % Per-pixel count is used instead of blindly dividing by WINDOW so a
    % transient invalid pixel does not bias the average toward zero.
    % ---------------------------------------------------------------------

    dWindow = windowSum ./ max(windowCount,single(1));

    dWindow(windowCount == 0) = NaN;

    % ---------------------------------------------------------------------
    % OPTIMIZED SPARSE RADIAL INTEGRATION
    % ---------------------------------------------------------------------

    deltaIq(:,g) = radial_sparse_mean(dWindow,radialOp);

    normCell{g} = normVals;

    % parfor reductions
    diffSumAll   = diffSumAll   + localDiffSum;
    diffCountAll = diffCountAll + localDiffCount;

    send(progressQ,1);
end

elapsedSeconds = toc(ticRun);

fprintf('\n');

%% ========================================================================
% ASSEMBLE NORMALIZATION TRACE
% ========================================================================

normalizationValue = nan(nUsed,1);

for g = 1:nOut

    first = (g-1)*groupSize + 1;

    normalizationValue(first:first+groupSize-1) = normCell{g};
end

normalizationScale = NORM_REFERENCE ./ normalizationValue;

%% ========================================================================
% FINAL 2-D AVERAGE DIFFERENCE
% ========================================================================

avgDiff2D = diffSumAll ./ max(diffCountAll,1);

avgDiff2D(diffCountAll == 0) = NaN;

% Permanently masked pixels should remain visibly masked.
avgDiff2D(badMask) = NaN;

%% ========================================================================
% SOURCE FRAME RANGE FOR EACH DeltaI(q) COLUMN
% ========================================================================

windowFrameIDStart = nan(nOut,1);
windowFrameIDEnd   = nan(nOut,1);

for g = 1:nOut

    first = (g-1)*groupSize + 1;
    last  = first + groupSize - 1;

    windowFrameIDStart(g) = usedFrameIDs(first);
    windowFrameIDEnd(g)   = usedFrameIDs(last);
end

%% ========================================================================
% RESULT STRUCTURE
% ========================================================================

R = struct;

R.q = qCenter(:);

R.deltaIq = deltaIq;

R.avgDiff2D = avgDiff2D;

R.avgDiffValidCount = diffCountAll;

R.normalization = struct;
R.normalization.frameID = usedFrameIDs;
R.normalization.value   = normalizationValue;
R.normalization.scale   = normalizationScale;
R.normalization.roiRows = NORM_ROWS;
R.normalization.roiCols = NORM_COLS;

R.window = WINDOW;
R.Nq = NQ;

R.windowFrameIDStart = windowFrameIDStart;
R.windowFrameIDEnd   = windowFrameIDEnd;

R.inputDirectory = inDir;
R.inputFiles = names(1:nUsed);
R.inputFrameIDs = usedFrameIDs;

R.nInputFiles = nFiles;
R.nUsedFiles = nUsed;
R.nDiscardedTailFiles = nTail;

R.reader = "Tiff";
R.workers = N_WORKERS;
R.radialMethod = "sparse";
R.correctionMode = CORRECTION_MODE;
R.normReference = NORM_REFERENCE;

R.elapsedSeconds = elapsedSeconds;
R.TIFF_per_second = nUsed/elapsedSeconds;
R.input_MB_per_second = sum(fileBytes(1:nUsed))/1e6/elapsedSeconds;

%% ========================================================================
% SAVE
% ========================================================================

if ~isfolder(outDir)
    mkdir(outDir);
end

fprintf('\nSaving result...\n');

ticSave = tic;

% Result is comfortably below the v7 per-variable size limit for this run.
% Uncompressed output minimizes CPU overhead.
save(outputFile,'R','-v7','-nocompression');

saveSeconds = toc(ticSave);

fprintf('\n============================================================\n');
fprintf('DONE\n');
fprintf('============================================================\n');
fprintf('Processing time : %.1f s = %.2f min\n', ...
    elapsedSeconds,elapsedSeconds/60);
fprintf('Throughput      : %.2f TIFF/s\n',R.TIFF_per_second);
fprintf('Input rate      : %.1f MB/s\n',R.input_MB_per_second);
fprintf('Save time       : %.2f s\n',saveSeconds);
fprintf('Output          : %s\n',outputFile);

%% ========================================================================
% QUICK DIAGNOSTIC PLOTS
% ========================================================================

% 1) Normalization monitor
figure('Name','Normalization');

plot( ...
    R.normalization.frameID, ...
    R.normalization.value, ...
    'LineWidth',1);

grid on

xlabel('Frame ID');
ylabel('Normalization ROI sum');
title('Per-frame normalization');

% 2) Final 2-D average normalized difference
figure('Name','Average 2-D difference');

imagesc(R.avgDiff2D);

axis image
colorbar

title(sprintf( ...
    'Average normalized I_{even}-I_{odd}, W=%d',WINDOW));

xlabel('Detector column');
ylabel('Detector row');

% 3) Radial DeltaI(q,t/window index)
figure('Name','Delta I(q)');

imagesc( ...
    1:size(R.deltaIq,2), ...
    R.q, ...
    R.deltaIq);

axis xy
colorbar

xlabel('Output window index');
ylabel('q');

title(sprintf('\\DeltaI(q), W=%d',WINDOW));

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function A = read_tiff_fast(filename)
% Fast reader selected by the benchmark.

    t = Tiff(char(filename),'r');

    A = t.read();

    t.close();
end


function iq = radial_sparse_mean(I,op)
% Corrected radial mean using the benchmark-winning sparse projection.
%
% op.S:
%   [Nq x Nvalidpixels], one nonzero per valid detector pixel.
%
% Geometry / polarization correction is already stored in op.factor.

    v = I(op.mask);

    good = isfinite(v);

    % Corrected detector values.
    values = double(v) .* double(op.factor);

    % Invalid dynamic pixels must contribute neither signal nor denominator.
    values(~good) = 0;

    sums = op.S * values;

    counts = op.S * double(good);

    z = sums ./ max(counts,1);

    z(counts == 0) = NaN;

    iq = single(z);
end


function progress_text(msg,nOut,nUsed,groupSize,startToken)
% Single-line command-window progress display for parfor.

    persistent nDone lastPrintedPercent lastPrintedTime

    if isstring(msg) || ischar(msg)
        if string(msg) == "reset"
            nDone = 0;
            lastPrintedPercent = -Inf;
            lastPrintedTime = -Inf;

            fprintf(['  0.00%% | 0/%d windows | 0/%d TIFFs | ' ...
                     'elapsed 00:00 | ETA estimating...'], ...
                     nOut,nUsed);
            return
        end
    end

    if isempty(nDone)
        nDone = 0;
        lastPrintedPercent = -Inf;
        lastPrintedTime = -Inf;
    end

    nDone = min(nDone + double(msg), nOut);

    elapsed = toc(startToken);
    pct = 100*nDone/nOut;

    % Avoid excessive command-window printing.
    shouldPrint = ...
        (pct-lastPrintedPercent >= 1) || ...
        (elapsed-lastPrintedTime >= 5) || ...
        (nDone == nOut);

    if ~shouldPrint
        return
    end

    doneTIFF = min(nDone*groupSize,nUsed);

    if nDone > 0
        etaSeconds = elapsed/nDone * (nOut-nDone);
        etaText = duration_text(etaSeconds);
    else
        etaText = "estimating...";
    end

    elapsedText = duration_text(elapsed);

    fprintf(['\r%6.2f%% | %d/%d windows | %d/%d TIFFs | ' ...
             'elapsed %s | ETA %s          '], ...
             pct,nDone,nOut,doneTIFF,nUsed, ...
             char(elapsedText),char(etaText));

    lastPrintedPercent = pct;
    lastPrintedTime = elapsed;

    if nDone == nOut
        fprintf('\n');
    end
end


function txt = duration_text(secondsValue)
% Compact hh:mm:ss / mm:ss formatter.

    secondsValue = max(0,round(double(secondsValue)));

    h = floor(secondsValue/3600);
    m = floor(mod(secondsValue,3600)/60);
    s = mod(secondsValue,60);

    if h > 0
        txt = string(sprintf('%02d:%02d:%02d',h,m,s));
    else
        txt = string(sprintf('%02d:%02d',m,s));
    end
end
