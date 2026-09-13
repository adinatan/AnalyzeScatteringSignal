clearvars
close all
clc

%% PILATUS OPTIMIZED PRODUCTION ANALYSIS
% Configurable file pairing + subset selection + progress/ETA
%
% Benchmarked architecture retained:
%   CPU
%   Tiff reader
%   4 process workers
%   sparse radial projection
%   Nq = 1000
%   direct processing from X:
%
% NEW: arbitrary regular subtraction pattern
%
%   pair k:
%       I(1 + PAIR_OFFSET + (k-1)*PAIR_STRIDE)
%       -
%       I(1               + (k-1)*PAIR_STRIDE)
%
% Examples, relative to the SELECTED file subset:
%
%   PAIR_OFFSET = 1; PAIR_STRIDE = 2;
%       -> 2-1, 4-3, 6-5, ...
%
%   PAIR_OFFSET = 2; PAIR_STRIDE = 2;
%       -> 3-1, 5-3, 7-5, ...
%
%   PAIR_OFFSET = 2; PAIR_STRIDE = 1;
%       -> 3-1, 4-2, 5-3, ...
%
% WINDOW is the number of consecutive pair differences averaged before
% producing one DeltaI(q) column.
%
% Overlapping source frames within one WINDOW are read and normalized only
% once, then reused for all pair differences in that window.

%% ========================================================================
% USER SETTINGS
% ========================================================================

inDir = "X:\2026-3\RuettSept26\data\fresh_ag_1\run3_UV\";

outDir = ...
    "C:\Users\11idbuser\Downloads\Natan\pilatus_fast_diff_matlab\analysis\run003_optimized";

% ------------------------------------------------------------------------
% FILE SUBSET
%
% Positions in the numerically sorted TIFF list.
%
% Examples:
%   FILE_START = 1;     FILE_END = 1000;
%   FILE_START = 1001;  FILE_END = 5000;
%   FILE_START = 1;     FILE_END = inf;
%
% Pair numbering below is RELATIVE TO THIS SELECTED SUBSET.
% ------------------------------------------------------------------------
FILE_START = 1;
FILE_END   = inf;

% ------------------------------------------------------------------------
% SUBTRACTION PATTERN
%
% result = later file - earlier file
%
% Standard:
%   2-1, 4-3, 6-5, ...
%       PAIR_OFFSET = 1;
%       PAIR_STRIDE = 2;
%
% Requested alternative:
%   3-1, 5-3, 7-5, ...
%       PAIR_OFFSET = 2;
%       PAIR_STRIDE = 2;
% ------------------------------------------------------------------------
PAIR_OFFSET = 2;
PAIR_STRIDE = 2;

% Number of pair differences averaged into each DeltaI(q) output.
WINDOW = 3;

% Benchmark winner.
NQ = 1000;

% Benchmark winner.
N_WORKERS = 4;

% normalized frame = raw frame * NORM_REFERENCE / normalization_ROI_sum
NORM_REFERENCE = 1;

% "divide"   -> Icorr = I ./ (geometry .* polarization)
% "multiply" -> Icorr = I .* (geometry .* polarization)
CORRECTION_MODE = "divide";

% Output filename encodes the pairing pattern.
outputFile = fullfile(outDir, ...
    sprintf("run003_offset%d_stride%d_W%d_Nq%d_optimized.mat", ...
    PAIR_OFFSET,PAIR_STRIDE,WINDOW,NQ));

%% ========================================================================
% VALIDATE USER SETTINGS
% ========================================================================

mustBePositiveInteger(FILE_START,'FILE_START');

if ~(isscalar(FILE_END) && (isinf(FILE_END) || ...
        (FILE_END >= FILE_START && FILE_END == floor(FILE_END))))
    error('FILE_END must be an integer >= FILE_START, or inf.');
end

mustBePositiveInteger(PAIR_OFFSET,'PAIR_OFFSET');
mustBePositiveInteger(PAIR_STRIDE,'PAIR_STRIDE');
mustBePositiveInteger(WINDOW,'WINDOW');
mustBePositiveInteger(NQ,'NQ');
mustBePositiveInteger(N_WORKERS,'N_WORKERS');

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

nFilesTotal = numel(paths);

% Acquisition order must be intact.
if any(diff(frameIDs) ~= 1)

    bad = find(diff(frameIDs) ~= 1,1,'first');

    error(['Frame IDs are not consecutive near %d -> %d. ' ...
           'Fix missing/duplicate TIFFs before processing.'], ...
           frameIDs(bad),frameIDs(bad+1));
end

%% ========================================================================
% APPLY FILE SUBSET
% ========================================================================

if isinf(FILE_END)
    FILE_END = nFilesTotal;
end

if FILE_END > nFilesTotal
    error('FILE_END=%d exceeds the %d available TIFFs.', ...
        FILE_END,nFilesTotal);
end

selectionGlobalPosition = (FILE_START:FILE_END)';

paths     = paths(selectionGlobalPosition);
names     = names(selectionGlobalPosition);
fileBytes = fileBytes(selectionGlobalPosition);
frameIDs  = frameIDs(selectionGlobalPosition);

nFiles = numel(paths);

if 1 + PAIR_OFFSET > nFiles
    error(['The selected subset has only %d TIFFs, but PAIR_OFFSET=%d ' ...
           'requires at least %d.'], ...
           nFiles,PAIR_OFFSET,1+PAIR_OFFSET);
end

%% ========================================================================
% BUILD PAIR LIST
%
% Pair positions are relative to the selected subset:
%
%   A_k = 1 + (k-1)*PAIR_STRIDE
%   B_k = A_k + PAIR_OFFSET
%   D_k = I(B_k) - I(A_k)
% ========================================================================

nPairsAvailable = ...
    floor((nFiles - 1 - PAIR_OFFSET)/PAIR_STRIDE) + 1;

pairA = 1 + (0:nPairsAvailable-1)'*PAIR_STRIDE;
pairB = pairA + PAIR_OFFSET;

% Only complete WINDOW groups are processed.
nOut = floor(nPairsAvailable/WINDOW);

if nOut < 1
    error(['Only %d subtraction pairs are available, not enough for ' ...
           'WINDOW=%d.'],nPairsAvailable,WINDOW);
end

nPairsUsed = nOut*WINDOW;
nTailPairs = nPairsAvailable-nPairsUsed;

pairAUsed = pairA(1:nPairsUsed);
pairBUsed = pairB(1:nPairsUsed);

uniqueFilesUsed = unique([pairAUsed; pairBUsed]);

fprintf('\n============================================================\n');
fprintf('PILATUS OPTIMIZED PRODUCTION ANALYSIS\n');
fprintf('============================================================\n');
fprintf('Available TIFFs   : %d\n',nFilesTotal);
fprintf('Selected files    : %d ... %d (%d TIFFs)\n', ...
    FILE_START,FILE_END,nFiles);
fprintf('Selected size     : %.3f GB\n',sum(fileBytes)/1e9);
fprintf('Selected frame IDs: %d ... %d\n',frameIDs(1),frameIDs(end));
fprintf('Pair pattern      : I(n+%d) - I(n), step n by %d\n', ...
    PAIR_OFFSET,PAIR_STRIDE);
fprintf('First pairs       : ');

nShow = min(4,nPairsAvailable);

for k = 1:nShow
    if k > 1
        fprintf(', ');
    end
    fprintf('%d-%d',pairB(k),pairA(k));
end

if nPairsAvailable > nShow
    fprintf(', ...');
end

fprintf('\n');
fprintf('Pair differences  : %d available, %d used\n', ...
    nPairsAvailable,nPairsUsed);
fprintf('Window            : %d pairs/output\n',WINDOW);
fprintf('Output curves     : %d\n',nOut);
fprintf('Unique TIFFs used : %d\n',numel(uniqueFilesUsed));
fprintf('q bins            : %d\n',NQ);
fprintf('Process workers   : %d\n',N_WORKERS);

if nTailPairs > 0
    warning('%d final pair(s) do not complete a WINDOW and will be ignored.', ...
        nTailPairs);
end

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

badMask = lab6_avg <= 0;

detectorValid = ~badMask;

%% ========================================================================
% HARD-CODED NORMALIZATION ROI
% ========================================================================

normMask = false(size(qMap));

NORM_ROWS = 1:408;
NORM_COLS = 990:1475;

normMask(NORM_ROWS,NORM_COLS) = true;
normMask(badMask) = false;

normIdx = find(normMask);

%% ========================================================================
% CHECK DETECTOR / CALIBRATION DIMENSIONS
% ========================================================================

A0 = read_tiff_fast(paths(1));

[ny,nx] = size(A0);

if ~isequal(size(A0),size(qMap),size(geometry), ...
        size(polarization),size(badMask))

    error(['TIFF, qMap, geometry, polarization and badMask must all ' ...
           'have exactly the same detector dimensions.']);
end

clear A0

fprintf('Detector          : %d x %d\n',ny,nx);
fprintf('Normalization ROI : %d good pixels\n',numel(normIdx));

%% ========================================================================
% PRECOMPUTE BENCHMARK-WINNING SPARSE RADIAL OPERATOR
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
qCenter = (qEdges(1:end-1)+qEdges(2:end))/2;

qBinMap = discretize(qMap,qEdges);

qMask = qStaticValid & ~isnan(qBinMap);

qBins = double(qBinMap(qMask));

if CORRECTION_MODE == "divide"
    qFactor = single(1 ./ corrMap(qMask));
else
    qFactor = single(corrMap(qMask));
end

nRadialPixels = numel(qBins);

Srad = sparse( ...
    qBins, ...
    (1:nRadialPixels)', ...
    ones(nRadialPixels,1), ...
    NQ, ...
    nRadialPixels);

radialOp = struct;
radialOp.mask    = qMask;
radialOp.factor  = qFactor;
radialOp.S       = Srad;
radialOp.Nq      = NQ;
radialOp.qCenter = qCenter;

clear qBinMap qBins qGood Srad

fprintf('Radial pixels     : %d\n',nRadialPixels);

%% ========================================================================
% PRECOMPUTE EXACT FILE READ COUNT / EXPECTED INPUT BYTES
%
% Overlapping frames INSIDE one WINDOW are read only once.
% A frame shared by two different parallel windows may be read once by each
% window, because the windows are independent parfor jobs.
% ========================================================================

groupReadCount = zeros(nOut,1);
groupReadBytes = zeros(nOut,1);

for g = 1:nOut

    p = (g-1)*WINDOW + (1:WINDOW);

    groupFiles = unique([pairAUsed(p); pairBUsed(p)]);

    groupReadCount(g) = numel(groupFiles);
    groupReadBytes(g) = sum(fileBytes(groupFiles));
end

nTIFFReadsExpected = sum(groupReadCount);
inputBytesExpected = sum(groupReadBytes);

fprintf('Expected TIFF reads: %d (%.3f GB)\n', ...
    nTIFFReadsExpected,inputBytesExpected/1e9);

%% ========================================================================
% START BENCHMARKED PROCESS POOL
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

% Each parallel window may normalize a small set of source TIFFs.
% Store the relative subset positions and normalization values for later
% assembly on the MATLAB client.
normFileCell = cell(nOut,1);
normValueCell = cell(nOut,1);

% Final average over EVERY pair difference.
diffSumAll   = zeros(ny,nx,'double');
diffCountAll = zeros(ny,nx,'double');

%% ========================================================================
% ONE-PASS PARALLEL ANALYSIS
% ========================================================================

fprintf('\nProcessing...\n');

progressQ = parallel.pool.DataQueue;

ticRun = tic;

progress_text("reset",nOut,nPairsUsed,WINDOW,ticRun);

afterEach(progressQ, ...
    @(msg) progress_text(msg,nOut,nPairsUsed,WINDOW,ticRun));

parfor g = 1:nOut

    % Pair numbers belonging to this output window.
    p = (g-1)*WINDOW + (1:WINDOW);

    aIdx = pairAUsed(p);
    bIdx = pairBUsed(p);

    % ---------------------------------------------------------------
    % READ EVERY UNIQUE TIFF IN THIS WINDOW ONLY ONCE.
    %
    % Example for 3-1,5-3,7-5:
    %   source positions = [1 3 5 7]
    % rather than reading 6 source TIFFs.
    % ---------------------------------------------------------------

    groupFiles = unique([aIdx; bIdx]);

    nGroupFiles = numel(groupFiles);

    normalizedFrames = cell(nGroupFiles,1);
    groupNorm = nan(nGroupFiles,1);

    for r = 1:nGroupFiles

        fileIdx = groupFiles(r);

        raw = read_tiff_fast(paths(fileIdx));

        valid = isfinite(raw) & raw >= 0;

        idxNorm = normIdx(valid(normIdx));

        nrm = sum(double(raw(idxNorm)));

        if ~(isfinite(nrm) && nrm > 0)
            error(['Invalid normalization scalar in output window %d, ' ...
                   'selected file position %d.'],g,fileIdx);
        end

        f = single(raw) * single(NORM_REFERENCE/nrm);

        % Apply permanent + dynamic mask before any subtraction.
        f(~(detectorValid & valid)) = NaN;

        normalizedFrames{r} = f;
        groupNorm(r) = nrm;
    end

    % Map each pair endpoint to its cached frame.
    [~,aLoc] = ismember(aIdx,groupFiles);
    [~,bLoc] = ismember(bIdx,groupFiles);

    % Detector-space accumulation for this output window.
    windowSum   = zeros(ny,nx,'single');
    windowCount = zeros(ny,nx,'single');

    % Reduction contribution to all-run mean difference.
    localDiffSum   = zeros(ny,nx,'double');
    localDiffCount = zeros(ny,nx,'double');

    for j = 1:WINDOW

        f1 = normalizedFrames{aLoc(j)};
        f2 = normalizedFrames{bLoc(j)};

        pairValid = isfinite(f1) & isfinite(f2);

        d = f2-f1;

        d(~pairValid) = 0;

        % Average over every pair for final 2-D result.
        localDiffSum   = localDiffSum   + double(d);
        localDiffCount = localDiffCount + double(pairValid);

        % Average WINDOW pairs before radial integration.
        windowSum   = windowSum   + d;
        windowCount = windowCount + single(pairValid);
    end

    dWindow = windowSum ./ max(windowCount,single(1));
    dWindow(windowCount == 0) = NaN;

    deltaIq(:,g) = radial_sparse_mean(dWindow,radialOp);

    normFileCell{g}  = groupFiles;
    normValueCell{g} = groupNorm;

    diffSumAll   = diffSumAll   + localDiffSum;
    diffCountAll = diffCountAll + localDiffCount;

    send(progressQ,1);
end

elapsedSeconds = toc(ticRun);

fprintf('\n');

%% ========================================================================
% ASSEMBLE NORMALIZATION TRACE
%
% Files that are inside the selected subset but never used by the requested
% pairing pattern remain NaN.
% ========================================================================

normalizationValue = nan(nFiles,1);
normalizationUsed  = false(nFiles,1);

for g = 1:nOut

    idx = normFileCell{g};
    val = normValueCell{g};

    normalizationValue(idx) = val;
    normalizationUsed(idx) = true;
end

normalizationScale = nan(nFiles,1);

normalizationScale(normalizationUsed) = ...
    NORM_REFERENCE ./ normalizationValue(normalizationUsed);

%% ========================================================================
% FINAL 2-D AVERAGE DIFFERENCE
% ========================================================================

avgDiff2D = diffSumAll ./ max(diffCountAll,1);

avgDiff2D(diffCountAll == 0) = NaN;
avgDiff2D(badMask) = NaN;

%% ========================================================================
% PAIR / WINDOW METADATA
% ========================================================================

pairAFrameID = frameIDs(pairAUsed);
pairBFrameID = frameIDs(pairBUsed);

pairAGlobalFilePosition = selectionGlobalPosition(pairAUsed);
pairBGlobalFilePosition = selectionGlobalPosition(pairBUsed);

windowFrameIDStart = nan(nOut,1);
windowFrameIDEnd   = nan(nOut,1);

windowPairStart = nan(nOut,1);
windowPairEnd   = nan(nOut,1);

for g = 1:nOut

    p = (g-1)*WINDOW + (1:WINDOW);

    allFrames = [pairAFrameID(p); pairBFrameID(p)];

    windowFrameIDStart(g) = min(allFrames);
    windowFrameIDEnd(g)   = max(allFrames);

    windowPairStart(g) = p(1);
    windowPairEnd(g)   = p(end);
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
R.normalization.filePositionInSubset = (1:nFiles)';
R.normalization.globalFilePosition = selectionGlobalPosition;
R.normalization.frameID = frameIDs;
R.normalization.value = normalizationValue;
R.normalization.scale = normalizationScale;
R.normalization.used = normalizationUsed;
R.normalization.roiRows = NORM_ROWS;
R.normalization.roiCols = NORM_COLS;

R.pairing = struct;
R.pairing.offset = PAIR_OFFSET;
R.pairing.stride = PAIR_STRIDE;
R.pairing.expression = sprintf('I(n+%d)-I(n), n advances by %d', ...
    PAIR_OFFSET,PAIR_STRIDE);

R.pairing.A_positionInSubset = pairAUsed;
R.pairing.B_positionInSubset = pairBUsed;

R.pairing.A_globalFilePosition = pairAGlobalFilePosition;
R.pairing.B_globalFilePosition = pairBGlobalFilePosition;

R.pairing.A_frameID = pairAFrameID;
R.pairing.B_frameID = pairBFrameID;

R.window = WINDOW;
R.Nq = NQ;

R.windowPairStart = windowPairStart;
R.windowPairEnd = windowPairEnd;
R.windowFrameIDStart = windowFrameIDStart;
R.windowFrameIDEnd = windowFrameIDEnd;

R.inputDirectory = inDir;
R.inputFiles = names;
R.inputFrameIDs = frameIDs;

R.nAvailableFiles = nFilesTotal;
R.selectedFileStart = FILE_START;
R.selectedFileEnd = FILE_END;
R.nSelectedFiles = nFiles;

R.nPairsAvailable = nPairsAvailable;
R.nPairsUsed = nPairsUsed;
R.nDiscardedTailPairs = nTailPairs;
R.nUniqueFilesUsed = numel(uniqueFilesUsed);

R.expectedTIFFReads = nTIFFReadsExpected;
R.expectedInputBytesRead = inputBytesExpected;

R.reader = "Tiff";
R.workers = N_WORKERS;
R.radialMethod = "sparse";
R.correctionMode = CORRECTION_MODE;
R.normReference = NORM_REFERENCE;

R.elapsedSeconds = elapsedSeconds;

R.pairs_per_second = nPairsUsed/elapsedSeconds;
R.TIFF_reads_per_second = nTIFFReadsExpected/elapsedSeconds;
R.input_MB_per_second = inputBytesExpected/1e6/elapsedSeconds;

%% ========================================================================
% SAVE
% ========================================================================

if ~isfolder(outDir)
    mkdir(outDir);
end

fprintf('\nSaving result...\n');

ticSave = tic;

save(outputFile,'R','-v7','-nocompression');

saveSeconds = toc(ticSave);

fprintf('\n============================================================\n');
fprintf('DONE\n');
fprintf('============================================================\n');
fprintf('Processing time : %.1f s = %.2f min\n', ...
    elapsedSeconds,elapsedSeconds/60);
fprintf('Pair rate       : %.2f pair differences/s\n',R.pairs_per_second);
fprintf('TIFF read rate  : %.2f reads/s\n',R.TIFF_reads_per_second);
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
    '.-', ...
    'LineWidth',1);

grid on

xlabel('Frame ID');
ylabel('Normalization ROI sum');

title(sprintf( ...
    'Normalization, offset=%d, stride=%d', ...
    PAIR_OFFSET,PAIR_STRIDE));

% 2) Final 2-D average normalized difference
figure('Name','Average 2-D difference');

imagesc(R.avgDiff2D);

axis image
colorbar

title(sprintf( ...
    'Average normalized difference, offset=%d, stride=%d, W=%d', ...
    PAIR_OFFSET,PAIR_STRIDE,WINDOW));

xlabel('Detector column');
ylabel('Detector row');

% 3) Radial DeltaI(q,window)
figure('Name','Delta I(q)');

imagesc( ...
    1:size(R.deltaIq,2), ...
    R.q, ...
    R.deltaIq);

axis xy
colorbar

xlabel('Output window index');
ylabel('q');

title(sprintf( ...
    '\\DeltaI(q), offset=%d, stride=%d, W=%d', ...
    PAIR_OFFSET,PAIR_STRIDE,WINDOW));

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function A = read_tiff_fast(filename)
% Fast reader selected by benchmark.

    t = Tiff(char(filename),'r');
    A = t.read();
    t.close();
end


function iq = radial_sparse_mean(I,op)
% Geometry/polarization-corrected radial mean using sparse projection.

    v = I(op.mask);

    good = isfinite(v);

    values = double(v) .* double(op.factor);

    values(~good) = 0;

    sums = op.S * values;
    counts = op.S * double(good);

    z = sums ./ max(counts,1);

    z(counts == 0) = NaN;

    iq = single(z);
end


function progress_text(msg,nOut,nPairsUsed,window,startToken)
% Single-line command-window progress / ETA for parfor.

    persistent nDone lastPrintedPercent lastPrintedTime

    if isstring(msg) || ischar(msg)

        if string(msg) == "reset"

            nDone = 0;
            lastPrintedPercent = -Inf;
            lastPrintedTime = -Inf;

            fprintf(['  0.00%% | 0/%d windows | 0/%d pairs | ' ...
                     'elapsed 00:00 | ETA estimating...'], ...
                     nOut,nPairsUsed);

            return
        end
    end

    if isempty(nDone)

        nDone = 0;
        lastPrintedPercent = -Inf;
        lastPrintedTime = -Inf;
    end

    nDone = min(nDone+double(msg),nOut);

    elapsed = toc(startToken);

    pct = 100*nDone/nOut;

    shouldPrint = ...
        (pct-lastPrintedPercent >= 1) || ...
        (elapsed-lastPrintedTime >= 5) || ...
        (nDone == nOut);

    if ~shouldPrint
        return
    end

    donePairs = min(nDone*window,nPairsUsed);

    etaSeconds = elapsed/max(nDone,1)*(nOut-nDone);

    elapsedText = duration_text(elapsed);
    etaText = duration_text(etaSeconds);

    fprintf(['\r%6.2f%% | %d/%d windows | %d/%d pairs | ' ...
             'elapsed %s | ETA %s          '], ...
             pct,nDone,nOut,donePairs,nPairsUsed, ...
             char(elapsedText),char(etaText));

    lastPrintedPercent = pct;
    lastPrintedTime = elapsed;

    if nDone == nOut
        fprintf('\n');
    end
end


function txt = duration_text(secondsValue)

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


function mustBePositiveInteger(x,name)

    if ~(isscalar(x) && isfinite(x) && x >= 1 && x == floor(x))
        error('%s must be a positive integer.',name);
    end
end
