clearvars
close all
clc

%% PILATUS OPTIMIZED SLIDING-TIME ANALYSIS
%
% Process TIFFs ONCE, then construct overlapping time windows such as:
%   0.00-1.00 h, 0.25-1.25 h, 0.50-1.50 h, ...
%
% The expensive TIFF read / normalization / subtraction is not repeated.
% Pair-resolved radial curves are saved once. Detector-space differences
% are accumulated into atomic chunks defined by the sliding-window
% boundaries, then reused to build every overlapping hourly average.
%
% Benchmarked architecture retained:
%   Tiff reader, 4 process workers, sparse radial projection, Nq=1000.

%% ========================================================================
% USER SETTINGS
% ========================================================================

inDir = "X:\2026-3\RuettSept26\data\GC\run12_UV\";

outDir = ...
    "C:\Users\11idbuser\Downloads\Natan\pilatus_fast_diff_matlab\analysis\run012_sliding";

% Optional overall subset. Time zero is the first selected TIFF.
FILE_START = 1;
FILE_END   = inf;

% Acquisition / sliding-window settings.
REP_RATE_HZ        = 1;
SLIDING_WINDOW_MIN = 60;
SLIDING_STEP_MIN   = 15;

% Pair subtraction pattern.
% 1,2 -> 2-1,4-3,6-5,...
% 2,2 -> 3-1,5-3,7-5,...
PAIR_OFFSET = 1;
PAIR_STRIDE = 2;

% Analysis settings.
NQ = 1000;
N_WORKERS = 4;
NORM_REFERENCE = 1;
CORRECTION_MODE = "divide";  % "divide" or "multiply"

% Preserve your previous radial summary statistic.
RADIAL_TRIM_PERCENT = 37;
RADIAL_SMOOTH_SPAN  = 5;

% Plot settings.
DIFF_CLIM = [-1e-6 1e-6];
MAKE_PLOTS        = true;
SAVE_PLOTS        = true;
KEEP_FIGURES_OPEN = false;
PLOT_RESOLUTION   = 150;

outputFile = fullfile(outDir, ...
    sprintf("run012_sliding_%dmin_step%dmin_offset%d_stride%d_Nq%d.mat", ...
    SLIDING_WINDOW_MIN,SLIDING_STEP_MIN,PAIR_OFFSET,PAIR_STRIDE,NQ));

plotDir = fullfile(outDir,"sliding_plots");

%% ========================================================================
% VALIDATE SETTINGS
% ========================================================================

mustBePositiveInteger(FILE_START,'FILE_START');
mustBePositiveInteger(PAIR_OFFSET,'PAIR_OFFSET');
mustBePositiveInteger(PAIR_STRIDE,'PAIR_STRIDE');
mustBePositiveInteger(NQ,'NQ');
mustBePositiveInteger(N_WORKERS,'N_WORKERS');

if ~(isscalar(REP_RATE_HZ) && isfinite(REP_RATE_HZ) && REP_RATE_HZ > 0)
    error('REP_RATE_HZ must be positive.');
end

if ~(isscalar(SLIDING_WINDOW_MIN) && SLIDING_WINDOW_MIN > 0)
    error('SLIDING_WINDOW_MIN must be positive.');
end

if ~(isscalar(SLIDING_STEP_MIN) && SLIDING_STEP_MIN > 0)
    error('SLIDING_STEP_MIN must be positive.');
end

if ~(isscalar(RADIAL_TRIM_PERCENT) && ...
        RADIAL_TRIM_PERCENT >= 0 && RADIAL_TRIM_PERCENT < 100)
    error('RADIAL_TRIM_PERCENT must be in [0,100).');
end

filesPerWindowExact = REP_RATE_HZ * SLIDING_WINDOW_MIN * 60;
filesPerStepExact   = REP_RATE_HZ * SLIDING_STEP_MIN   * 60;

filesPerWindow = round(filesPerWindowExact);
filesPerStep   = round(filesPerStepExact);

if abs(filesPerWindow-filesPerWindowExact) > 1e-9
    error('REP_RATE_HZ * SLIDING_WINDOW_MIN * 60 must be an integer.');
end

if abs(filesPerStep-filesPerStepExact) > 1e-9
    error('REP_RATE_HZ * SLIDING_STEP_MIN * 60 must be an integer.');
end

%% ========================================================================
% INPUT FILE LIST - NUMERIC ORDER
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

if any(diff(frameIDs) ~= 1)
    bad = find(diff(frameIDs) ~= 1,1,'first');
    error(['Frame IDs are not consecutive near %d -> %d. ' ...
           'Fix missing/duplicate TIFFs before processing.'], ...
           frameIDs(bad),frameIDs(bad+1));
end

%% ========================================================================
% APPLY OPTIONAL FILE SUBSET
% ========================================================================

if isinf(FILE_END)
    FILE_END = nFilesTotal;
end

if ~(isscalar(FILE_END) && FILE_END >= FILE_START && FILE_END == floor(FILE_END))
    error('FILE_END must be an integer >= FILE_START, or inf.');
end

if FILE_END > nFilesTotal
    error('FILE_END=%d exceeds the %d available TIFFs.',FILE_END,nFilesTotal);
end

selectionGlobalPosition = (FILE_START:FILE_END)';

paths     = paths(selectionGlobalPosition);
names     = names(selectionGlobalPosition);
fileBytes = fileBytes(selectionGlobalPosition);
frameIDs  = frameIDs(selectionGlobalPosition);

nFiles = numel(paths);

if nFiles < filesPerWindow
    error(['Selected subset contains only %d TIFFs (%.3f h), but one ' ...
           '%.2f min sliding window requires %d TIFFs.'], ...
           nFiles,nFiles/REP_RATE_HZ/3600,SLIDING_WINDOW_MIN,filesPerWindow);
end

%% ========================================================================
% DEFINE SLIDING SOURCE-FILE WINDOWS
% ========================================================================

slideStartFile = (1:filesPerStep:(nFiles-filesPerWindow+1))';
slideEndFile   = slideStartFile + filesPerWindow - 1;

nSlides = numel(slideStartFile);

slideStartHourRelative = (slideStartFile-1)/REP_RATE_HZ/3600;
slideEndHourRelative   = slideStartHourRelative + SLIDING_WINDOW_MIN/60;

slideStartHourAcquisition = ...
    (selectionGlobalPosition(slideStartFile)-1)/REP_RATE_HZ/3600;
slideEndHourAcquisition = ...
    slideStartHourAcquisition + SLIDING_WINDOW_MIN/60;

%% ========================================================================
% GLOBAL PAIR LIST
% ========================================================================

nPairsAvailable = floor((nFiles-1-PAIR_OFFSET)/PAIR_STRIDE) + 1;

if nPairsAvailable < 1
    error('No requested subtraction pairs fit inside the selected file range.');
end

pairA = 1 + (0:nPairsAvailable-1)'*PAIR_STRIDE;
pairB = pairA + PAIR_OFFSET;

slideFirstPair = nan(nSlides,1);
slideLastPair  = nan(nSlides,1);

for s = 1:nSlides
    k1 = find(pairA >= slideStartFile(s),1,'first');
    k2 = find(pairB <= slideEndFile(s),1,'last');

    if isempty(k1) || isempty(k2) || k2 < k1
        error('Sliding window %d contains no complete subtraction pairs.',s);
    end

    slideFirstPair(s) = k1;
    slideLastPair(s)  = k2;
end

pairProcessFirst = min(slideFirstPair);
pairProcessLast  = max(slideLastPair);
nPairsProcessed  = pairProcessLast-pairProcessFirst+1;

%% ========================================================================
% BUILD ATOMIC PAIR CHUNKS
%
% Every sliding-window start/end becomes a chunk boundary. Therefore each
% 1-hour result can later be assembled exactly by adding whole chunks.
% ========================================================================

cuts = unique([ ...
    pairProcessFirst; ...
    pairProcessLast+1; ...
    slideFirstPair; ...
    slideLastPair+1]);

cuts = sort(cuts);

chunkFirstPair = cuts(1:end-1);
chunkLastPair  = cuts(2:end)-1;

validChunk = chunkLastPair >= chunkFirstPair;
chunkFirstPair = chunkFirstPair(validChunk);
chunkLastPair  = chunkLastPair(validChunk);

nChunks = numel(chunkFirstPair);
chunkPairCount = chunkLastPair-chunkFirstPair+1;

if max(chunkPairCount) <= intmax('uint16')
    COUNT_CLASS = "uint16";
else
    COUNT_CLASS = "uint32";
end

slideChunkIndices = cell(nSlides,1);

for s = 1:nSlides
    use = find( ...
        chunkFirstPair >= slideFirstPair(s) & ...
        chunkLastPair  <= slideLastPair(s));

    if isempty(use) || ...
            chunkFirstPair(use(1)) ~= slideFirstPair(s) || ...
            chunkLastPair(use(end)) ~= slideLastPair(s)
        error('Internal chunk construction error for sliding window %d.',s);
    end

    slideChunkIndices{s} = use;
end

%% ========================================================================
% REPORT SCHEDULE
% ========================================================================

fprintf('\n============================================================\n');
fprintf('PILATUS OPTIMIZED SLIDING-TIME ANALYSIS\n');
fprintf('============================================================\n');
fprintf('Available TIFFs      : %d\n',nFilesTotal);
fprintf('Selected TIFFs       : %d ... %d (%d files)\n', ...
    FILE_START,FILE_END,nFiles);
fprintf('Selected duration    : %.3f h\n',nFiles/REP_RATE_HZ/3600);
fprintf('Repetition rate      : %.6g Hz\n',REP_RATE_HZ);
fprintf('Sliding window       : %.2f min = %d TIFFs\n', ...
    SLIDING_WINDOW_MIN,filesPerWindow);
fprintf('Sliding step         : %.2f min = %d TIFFs\n', ...
    SLIDING_STEP_MIN,filesPerStep);
fprintf('Sliding outputs      : %d\n',nSlides);
fprintf('First interval       : %.2f - %.2f h\n', ...
    slideStartHourRelative(1),slideEndHourRelative(1));
fprintf('Last interval        : %.2f - %.2f h\n', ...
    slideStartHourRelative(end),slideEndHourRelative(end));
fprintf('Pair pattern         : I(n+%d)-I(n), n step %d\n', ...
    PAIR_OFFSET,PAIR_STRIDE);
fprintf('Pairs processed once : %d\n',nPairsProcessed);
fprintf('Atomic 2-D chunks    : %d\n',nChunks);
fprintf('q bins               : %d\n',NQ);
fprintf('Process workers      : %d\n',N_WORKERS);

fprintf('\nSliding intervals:\n');
for s = 1:nSlides
    fprintf('  %2d: %5.2f - %5.2f h | files %d:%d | pairs %d:%d\n', ...
        s,slideStartHourRelative(s),slideEndHourRelative(s), ...
        slideStartFile(s),slideEndFile(s), ...
        slideFirstPair(s),slideLastPair(s));
end

%% ========================================================================
% LOAD CALIBRATION / MASK
% ========================================================================

load("qgp.mat");
load("lab6avg_tz400.mat");

qMap         = QGP.qMap;
geometry     = QGP.geometryCorr;
polarization = QGP.polarizationCorr;

badMask = lab6_avg <= 0;
detectorValid = ~badMask;

%% Normalization ROI
normMask = false(size(qMap));
NORM_ROWS = 1:408;
NORM_COLS = 990:1475;
normMask(NORM_ROWS,NORM_COLS) = true;
normMask(badMask) = false;
normIdx = find(normMask);

%% Check dimensions
A0 = read_tiff_fast(paths(1));
[ny,nx] = size(A0);

if ~isequal(size(A0),size(qMap),size(geometry),size(polarization),size(badMask))
    error(['TIFF, qMap, geometry, polarization and badMask must all ' ...
           'have exactly the same detector dimensions.']);
end
clear A0

fprintf('\nDetector             : %d x %d\n',ny,nx);
fprintf('Normalization ROI    : %d good pixels\n',numel(normIdx));

%% ========================================================================
% PRECOMPUTE SPARSE RADIAL OPERATOR
% ========================================================================

corrMap = geometry .* polarization;

if CORRECTION_MODE == "divide"
    qStaticValid = detectorValid & isfinite(qMap) & ...
                   isfinite(corrMap) & corrMap ~= 0;
elseif CORRECTION_MODE == "multiply"
    qStaticValid = detectorValid & isfinite(qMap) & isfinite(corrMap);
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
    NQ,nRadialPixels);

radialOp = struct;
radialOp.mask    = qMask;
radialOp.factor  = qFactor;
radialOp.S       = Srad;
radialOp.Nq      = NQ;
radialOp.qCenter = qCenter;

clear qBinMap qBins qGood Srad

fprintf('Radial pixels        : %d\n',nRadialPixels);

%% ========================================================================
% START PROCESS POOL
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
% PREALLOCATE CHUNK OUTPUTS
% ========================================================================

chunkDiffSum   = cell(nChunks,1);
chunkDiffCount = cell(nChunks,1);
chunkIq        = cell(nChunks,1);
chunkNormFile  = cell(nChunks,1);
chunkNormValue = cell(nChunks,1);

%% ========================================================================
% ONE PASS THROUGH TIFF DATA
% ========================================================================

fprintf('\nProcessing TIFFs once...\n');
progressQ = parallel.pool.DataQueue;
ticRun = tic;
progress_pairs("reset",nPairsProcessed,ticRun);
afterEach(progressQ,@(nDone) progress_pairs(nDone,nPairsProcessed,ticRun));

parfor c = 1:nChunks

    kFirst = chunkFirstPair(c);
    kLast  = chunkLastPair(c);

    [localSum,localCount,localIq,normFiles,normValues] = ...
        process_pair_chunk( ...
            paths,pairA,pairB,kFirst,kLast, ...
            detectorValid,normIdx,NORM_REFERENCE,radialOp, ...
            ny,nx,COUNT_CLASS);

    chunkDiffSum{c}   = localSum;
    chunkDiffCount{c} = localCount;
    chunkIq{c}        = localIq;
    chunkNormFile{c}  = normFiles;
    chunkNormValue{c} = normValues;

    send(progressQ,kLast-kFirst+1);
end

elapsedSeconds = toc(ticRun);
fprintf('\n');

%% ========================================================================
% ASSEMBLE PAIR-RESOLVED DeltaI(q) ONCE
% ========================================================================

deltaIq = nan(NQ,nPairsAvailable,'single');

for c = 1:nChunks
    kk = chunkFirstPair(c):chunkLastPair(c);
    deltaIq(:,kk) = chunkIq{c};
end

%% ========================================================================
% ASSEMBLE NORMALIZATION TRACE
% ========================================================================

normalizationValue = nan(nFiles,1);
normalizationUsed  = false(nFiles,1);

for c = 1:nChunks
    idx = chunkNormFile{c};
    val = chunkNormValue{c};
    normalizationValue(idx) = val;
    normalizationUsed(idx) = true;
end

normalizationScale = nan(nFiles,1);
normalizationScale(normalizationUsed) = ...
    NORM_REFERENCE ./ normalizationValue(normalizationUsed);

%% ========================================================================
% BUILD ALL SLIDING RESULTS FROM CACHED DATA ONLY
% NO TIFF ACCESS BELOW THIS LINE.
% ========================================================================

fprintf('Building %d overlapping time windows from cached chunks...\n',nSlides);

slidingAvgDiff2D       = nan(ny,nx,nSlides,'single');
slidingRadialTrimMean  = nan(NQ,nSlides,'single');
slidingRadialSmooth    = nan(NQ,nSlides,'single');
slidingPairCount       = nan(nSlides,1);

for s = 1:nSlides

    useChunks = slideChunkIndices{s};

    sum2D   = zeros(ny,nx,'single');
    count2D = zeros(ny,nx,'uint32');

    for cc = useChunks(:)'
        sum2D = sum2D + chunkDiffSum{cc};
        count2D = count2D + uint32(chunkDiffCount{cc});
    end

    avg2D = sum2D ./ max(single(count2D),single(1));
    avg2D(count2D == 0) = NaN;
    avg2D(badMask) = NaN;

    slidingAvgDiff2D(:,:,s) = avg2D;

    k1 = slideFirstPair(s);
    k2 = slideLastPair(s);
    slidingPairCount(s) = k2-k1+1;

    % Same robust radial statistic used in your previous final plot.
    radialTrim = trimmean( ...
        double(deltaIq(:,k1:k2)), ...
        RADIAL_TRIM_PERCENT,2);

    slidingRadialTrimMean(:,s) = single(radialTrim);
    slidingRadialSmooth(:,s) = single( ...
        smoothdata(radialTrim,'movmean',RADIAL_SMOOTH_SPAN));
end

%% ========================================================================
% PAIR METADATA
% ========================================================================

pairAFrameID = frameIDs(pairA);
pairBFrameID = frameIDs(pairB);
pairAGlobalFilePosition = selectionGlobalPosition(pairA);
pairBGlobalFilePosition = selectionGlobalPosition(pairB);

pairMidTimeHourRelative = ...
    (((pairA-1)+(pairB-1))/2)/REP_RATE_HZ/3600;

%% ========================================================================
% RESULT STRUCTURE
% ========================================================================

R = struct;
R.q = qCenter(:);
R.deltaIq = deltaIq;

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
R.pairing.A_positionInSubset = pairA;
R.pairing.B_positionInSubset = pairB;
R.pairing.A_globalFilePosition = pairAGlobalFilePosition;
R.pairing.B_globalFilePosition = pairBGlobalFilePosition;
R.pairing.A_frameID = pairAFrameID;
R.pairing.B_frameID = pairBFrameID;
R.pairing.midTimeHourRelative = pairMidTimeHourRelative;

R.sliding = struct;
R.sliding.windowMinutes = SLIDING_WINDOW_MIN;
R.sliding.stepMinutes = SLIDING_STEP_MIN;
R.sliding.repRateHz = REP_RATE_HZ;
R.sliding.startFileInSubset = slideStartFile;
R.sliding.endFileInSubset = slideEndFile;
R.sliding.startGlobalFilePosition = selectionGlobalPosition(slideStartFile);
R.sliding.endGlobalFilePosition = selectionGlobalPosition(slideEndFile);
R.sliding.startHourRelative = slideStartHourRelative;
R.sliding.endHourRelative = slideEndHourRelative;
R.sliding.startHourAcquisition = slideStartHourAcquisition;
R.sliding.endHourAcquisition = slideEndHourAcquisition;
R.sliding.firstPair = slideFirstPair;
R.sliding.lastPair = slideLastPair;
R.sliding.pairCount = slidingPairCount;
R.sliding.avgDiff2D = slidingAvgDiff2D;
R.sliding.radialTrimMean = slidingRadialTrimMean;
R.sliding.radialSmooth = slidingRadialSmooth;
R.sliding.chunkFirstPair = chunkFirstPair;
R.sliding.chunkLastPair = chunkLastPair;

R.inputDirectory = inDir;
R.inputFiles = names;
R.inputFrameIDs = frameIDs;
R.nAvailableFiles = nFilesTotal;
R.selectedFileStart = FILE_START;
R.selectedFileEnd = FILE_END;
R.nSelectedFiles = nFiles;
R.Nq = NQ;
R.reader = "Tiff";
R.workers = N_WORKERS;
R.radialMethod = "sparse";
R.correctionMode = CORRECTION_MODE;
R.normReference = NORM_REFERENCE;
R.radialTrimPercent = RADIAL_TRIM_PERCENT;
R.radialSmoothSpan = RADIAL_SMOOTH_SPAN;
R.processingSeconds = elapsedSeconds;

%% ========================================================================
% SAVE ALL DATA
% ========================================================================

if ~isfolder(outDir)
    mkdir(outDir);
end
if SAVE_PLOTS && ~isfolder(plotDir)
    mkdir(plotDir);
end

fprintf('\nSaving sliding-window data...\n');
ticSave = tic;
save(outputFile,'R','-v7.3');
saveSeconds = toc(ticSave);

fprintf('\n============================================================\n');
fprintf('DATA COMPLETE\n');
fprintf('============================================================\n');
fprintf('TIFF processing : %.1f s = %.2f min\n',elapsedSeconds,elapsedSeconds/60);
fprintf('Save time       : %.2f s\n',saveSeconds);
fprintf('Sliding windows : %d\n',nSlides);
fprintf('Output          : %s\n',outputFile);

%% ========================================================================
% GENERATE ALL SLIDING-WINDOW PLOTS AT THE VERY END
% ========================================================================

if MAKE_PLOTS

    fprintf('\nGenerating %d final sliding-window plots...\n',nSlides);

    for s = 1:nSlides

        t0 = slideStartHourRelative(s);
        t1 = slideEndHourRelative(s);

        fig = figure( ...
            'Name',sprintf('%.2f - %.2f h',t0,t1), ...
            'Color','w');

        tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

        ax1 = nexttile(tl,1);
        imagesc(ax1,R.sliding.avgDiff2D(:,:,s));
        axis(ax1,'image');
        clim(ax1,DIFF_CLIM);
        colorbar(ax1);
        colormap(ax1,jet(256));
        xlabel(ax1,'Detector column');
        ylabel(ax1,'Detector row');
        title(ax1,sprintf('<\Delta I(x,y)>, %.2f-%.2f h',t0,t1));

        ax2 = nexttile(tl,2);
        plot(ax2,R.q,R.sliding.radialSmooth(:,s),'LineWidth',1.2);
        grid(ax2,'on');
        xlabel(ax2,'Q');
        ylabel(ax2,'Trimmed mean \DeltaI(q)');
        title(ax2,sprintf('%.0f%% trimmed mean, %.2f-%.2f h', ...
            RADIAL_TRIM_PERCENT,t0,t1));

        title(tl,sprintf( ...
            'Sliding window %d/%d | %.2f - %.2f h | %d pairs', ...
            s,nSlides,t0,t1,R.sliding.pairCount(s)));

        if SAVE_PLOTS
            plotName = sprintf('window_%02d_%05.2fh_to_%05.2fh.png',s,t0,t1);
            exportgraphics(fig,fullfile(plotDir,plotName), ...
                'Resolution',PLOT_RESOLUTION);
        end

        if ~KEEP_FIGURES_OPEN
            close(fig);
        end
    end

    fprintf('Plots saved to: %s\n',plotDir);
end

fprintf('============================================================\n');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function [diffSum,diffCount,iq,normFiles,normValues] = ...
    process_pair_chunk( ...
        paths,pairA,pairB,kFirst,kLast, ...
        detectorValid,normIdx,normReference,radialOp, ...
        ny,nx,countClass)

    kk = kFirst:kLast;
    nPair = numel(kk);

    diffSum = zeros(ny,nx,'single');
    diffCount = zeros(ny,nx,char(countClass));
    iq = nan(radialOp.Nq,nPair,'single');

    % Small rolling cache for overlapping patterns such as 3-1,5-3,7-5.
    cacheIdx = zeros(0,1);
    cacheFrames = cell(0,1);

    normFiles = zeros(0,1);
    normValues = zeros(0,1);

    for ii = 1:nPair

        k = kk(ii);
        aIdx = pairA(k);
        bIdx = pairB(k);

        locA = find(cacheIdx == aIdx,1);
        if isempty(locA)
            [fA,nA] = load_normalized_frame( ...
                paths(aIdx),detectorValid,normIdx,normReference);
            cacheIdx(end+1,1) = aIdx;
            cacheFrames{end+1,1} = fA;
            normFiles(end+1,1) = aIdx;
            normValues(end+1,1) = nA;
        else
            fA = cacheFrames{locA};
        end

        locB = find(cacheIdx == bIdx,1);
        if isempty(locB)
            [fB,nB] = load_normalized_frame( ...
                paths(bIdx),detectorValid,normIdx,normReference);
            cacheIdx(end+1,1) = bIdx;
            cacheFrames{end+1,1} = fB;
            normFiles(end+1,1) = bIdx;
            normValues(end+1,1) = nB;
        else
            fB = cacheFrames{locB};
        end

        pairValid = isfinite(fA) & isfinite(fB);
        d = fB-fA;
        d(~pairValid) = 0;

        diffSum = diffSum + d;
        diffCount = diffCount + cast(pairValid,char(countClass));

        dForRadial = d;
        dForRadial(~pairValid) = NaN;
        iq(:,ii) = radial_sparse_mean(dForRadial,radialOp);

        % Anything below the next A index can never be needed again.
        if ii < nPair
            nextA = pairA(kk(ii+1));
            keep = cacheIdx >= nextA;
            cacheIdx = cacheIdx(keep);
            cacheFrames = cacheFrames(keep);
        end
    end

    [normFiles,ia] = unique(normFiles,'stable');
    normValues = normValues(ia);
end


function [f,nrm] = load_normalized_frame( ...
    filename,detectorValid,normIdx,normReference)

    raw = read_tiff_fast(filename);
    valid = isfinite(raw) & raw >= 0;
    idxNorm = normIdx(valid(normIdx));
    nrm = sum(double(raw(idxNorm)));

    if ~(isfinite(nrm) && nrm > 0)
        error('Invalid normalization scalar for %s',filename);
    end

    f = single(raw) * single(normReference/nrm);
    f(~(detectorValid & valid)) = NaN;
end


function A = read_tiff_fast(filename)
    t = Tiff(char(filename),'r');
    A = t.read();
    t.close();
end


function iq = radial_sparse_mean(I,op)
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


function progress_pairs(msg,nPairs,startToken)

    persistent nDone lastPrintedPercent lastPrintedTime

    if isstring(msg) || ischar(msg)
        if string(msg) == "reset"
            nDone = 0;
            lastPrintedPercent = -Inf;
            lastPrintedTime = -Inf;
            fprintf(['  0.00%% | 0/%d pairs | elapsed 00:00 | ' ...
                     'ETA estimating...'],nPairs);
            return
        end
    end

    if isempty(nDone)
        nDone = 0;
        lastPrintedPercent = -Inf;
        lastPrintedTime = -Inf;
    end

    nDone = min(nDone+double(msg),nPairs);
    elapsed = toc(startToken);
    pct = 100*nDone/nPairs;

    shouldPrint = ...
        (pct-lastPrintedPercent >= 1) || ...
        (elapsed-lastPrintedTime >= 5) || ...
        (nDone == nPairs);

    if ~shouldPrint
        return
    end

    etaSeconds = elapsed/max(nDone,1)*(nPairs-nDone);

    fprintf(['\r%6.2f%% | %d/%d pairs | elapsed %s | ETA %s          '], ...
        pct,nDone,nPairs, ...
        char(duration_text(elapsed)),char(duration_text(etaSeconds)));

    lastPrintedPercent = pct;
    lastPrintedTime = elapsed;

    if nDone == nPairs
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
