clear all
close all
clc

%% QUICK BENCHMARK: choose the fastest configuration for the real PILATUS pipeline
%
% Tests the SAME operations used in the full analysis:
%   raw TIFF read
%   static bad-pixel mask
%   per-frame normalization ROI
%   normalized even-odd pair difference
%   window averaging
%   q integration
%   q/phi integration
%   all-image accumulation
%   q/phi output writing
%
% CONFIGURATIONS TESTED:
%   1) 1 process / Tiff reader
%   2) 2 processes / Tiff reader
%   3) 4 processes / Tiff reader
%
% Pool startup and benchmark-file deletion are NOT included in the timing.
%
% REQUIREMENTS:
%   pilatus_filelist.m
%   pilatus_stream_analysis.m
%
% Put this script in the same folder as those functions, then RUN it.

%% USER SETTINGS copied from your current analysis
inDir  = "X:\2026-3\RuettSept26\data\fresh_ag_1\run3_UV\";
outDir = "C:\Users\11idbuser\Downloads\Natan\pilatus_fast_diff_matlab\analysis\run003_w3";

window = 3;

% Quick benchmark size. Must be a multiple of 2*window after adjustment.
% 180 TIFFs is about 1.8 GB if each TIFF is about 10 MB.
BENCH_FILES_REQUESTED = 180;

% Small untimed warm-up for MATLAB/TIFF/network/filesystem initialization.
WARMUP_FILES_REQUESTED = 12;

% Delete large temporary q/phi benchmark outputs after each test.
DELETE_TEMP_OUTPUTS = true;

%% Find acquisition
[files,frameIDs,fileBytes] = pilatus_filelist(inDir,"*.tif");

nFiles = numel(files);

fprintf('\n============================================================\n');
fprintf('PILATUS QUICK BENCHMARK\n');
fprintf('============================================================\n');
fprintf('Input files : %d\n',nFiles);
fprintf('Input size  : %.3f TB\n',sum(fileBytes)/1e12);
fprintf('Window      : %d  (%d TIFFs per output window)\n',window,2*window);
fprintf('First IDs   : %d, %d\n',frameIDs(1),frameIDs(2));

%% Load detector calibration
load("qgp.mat");
load("lab6avg_tz400.mat");

qMap         = QGP.qMap;
phiMap       = QGP.phiMap;
geometry     = QGP.geometryCorr;
polarization = QGP.polarizationCorr;

%% Static detector mask
badMask = lab6_avg <= 0;

%% Hard-coded normalization ROI
normMask = false(size(qMap));

NORM_ROWS = 1:408;
NORM_COLS = 990:1475;

normMask(NORM_ROWS,NORM_COLS) = true;
normMask(badMask) = false;

%% q / phi bins
qGood = qMap(isfinite(qMap) & ~badMask);
qEdges = linspace(min(qGood),max(qGood),1001);

phiGood = phiMap(isfinite(phiMap) & ~badMask);
phiEdges = linspace(min(phiGood),max(phiGood),361);

%% Prepare fair benchmark blocks
groupSize = 2*window;

benchFiles = floor(BENCH_FILES_REQUESTED/groupSize)*groupSize;
warmFiles  = floor(WARMUP_FILES_REQUESTED/groupSize)*groupSize;

benchFiles = max(benchFiles,groupSize);
warmFiles  = max(warmFiles,groupSize);

blockLength = warmFiles + benchFiles;

if nFiles < 3*blockLength
    availablePerTest = floor(nFiles/3/groupSize)*groupSize;
    benchFiles = availablePerTest - warmFiles;
    benchFiles = floor(benchFiles/groupSize)*groupSize;

    if benchFiles < groupSize
        error('Not enough TIFFs to benchmark three independent configurations.');
    end

    blockLength = warmFiles + benchFiles;
end

% Three non-overlapping blocks spread through the acquisition.
maxStart = nFiles - blockLength + 1;

starts = round(linspace(1,maxStart,3));
starts = 1 + floor((starts-1)/groupSize)*groupSize;

for k = 2:3
    if starts(k) < starts(k-1) + blockLength
        starts(k) = starts(k-1) + blockLength;
        starts(k) = 1 + floor((starts(k)-1)/groupSize)*groupSize;
    end
end

if starts(3) + blockLength - 1 > nFiles
    starts(3) = nFiles - blockLength + 1;
    starts(3) = 1 + floor((starts(3)-1)/groupSize)*groupSize;
end

fprintf('\nQuick timing uses %d TIFFs/configuration + %d warm-up TIFFs.\n', ...
    benchFiles,warmFiles);
fprintf('Approx timed input/configuration: %.2f GB\n', ...
    median(fileBytes)*benchFiles/1e9);

%% Benchmark configurations
workersToTest = [1 2 4];
readerToTest  = ["tiff" "tiff" "tiff"];

nCfg = numel(workersToTest);

elapsed_s     = nan(nCfg,1);
filesPerSec   = nan(nCfg,1);
inputMBperSec = nan(nCfg,1);
estFull_min   = nan(nCfg,1);

benchRoot = fullfile(fileparts(outDir), ...
    "__pilatus_benchmark_" + string(datetime('now','Format','yyyyMMdd_HHmmss')));
mkdir(benchRoot);

fprintf('\nBenchmark temporary folder:\n%s\n',benchRoot);

for c = 1:nCfg

    workers = workersToTest(c);
    reader  = readerToTest(c);

    if workers == 1
        workerLabel = "worker";
    else
        workerLabel = "workers";
    end

    fprintf('\n------------------------------------------------------------\n');
    fprintf('TEST %d/%d : %d %s, %s reader\n', ...
        c,nCfg,workers,workerLabel,upper(reader));
    fprintf('------------------------------------------------------------\n');

    % Start the requested process pool BEFORE timing.
    p = gcp('nocreate');

    if workers == 1
        if ~isempty(p)
            delete(p);
        end
    else
        if isempty(p) || p.NumWorkers ~= workers
            if ~isempty(p)
                delete(p);
            end
            fprintf('Starting %d-process pool (not timed)...\n',workers);
            parpool("Processes",workers);
        end
    end

    first = starts(c);

    warmIdx  = first : first + warmFiles - 1;
    timedIdx = first + warmFiles : first + warmFiles + benchFiles - 1;

    warmFilesList  = files(warmIdx);
    timedFilesList = files(timedIdx);

    warmOut = fullfile(benchRoot,sprintf('cfg%d_warm',c));
    testOut = fullfile(benchRoot,sprintf('cfg%d_test',c));

    % -------------------- untimed warm-up --------------------
    fprintf('Warm-up: %d TIFFs...\n',numel(warmFilesList));

    pilatus_stream_analysis(warmFilesList,warmOut,window, ...
        Workers=workers, ...
        Reader=reader, ...
        ChunkFrames=64, ...
        Accumulator="double", ...
        OutputClass="single", ...
        BadMask=badMask, ...
        MaskNegative=true, ...
        NormMask=normMask, ...
        NormReference=1, ...
        QMap=qMap, ...
        PhiMap=phiMap, ...
        Geometry=geometry, ...
        Polarization=polarization, ...
        CorrectionMode="divide", ...
        QEdges=qEdges, ...
        PhiEdges=phiEdges, ...
        SaveQPhi=true, ...
        SaveDiff2D=false, ...
        SaveRawSum=true);

    if DELETE_TEMP_OUTPUTS && isfolder(warmOut)
        rmdir(warmOut,'s');
    end

    % -------------------- timed real pipeline --------------------
    timedBytes = sum(fileBytes(timedIdx));

    fprintf('TIMING: %d TIFFs, %.2f GB...\n', ...
        numel(timedFilesList),timedBytes/1e9);

    drawnow;
    t0 = tic;

    pilatus_stream_analysis(timedFilesList,testOut,window, ...
        Workers=workers, ...
        Reader=reader, ...
        ChunkFrames=64, ...
        Accumulator="double", ...
        OutputClass="single", ...
        BadMask=badMask, ...
        MaskNegative=true, ...
        NormMask=normMask, ...
        NormReference=1, ...
        QMap=qMap, ...
        PhiMap=phiMap, ...
        Geometry=geometry, ...
        Polarization=polarization, ...
        CorrectionMode="divide", ...
        QEdges=qEdges, ...
        PhiEdges=phiEdges, ...
        SaveQPhi=true, ...
        SaveDiff2D=false, ...
        SaveRawSum=true);

    elapsed_s(c) = toc(t0);

    filesPerSec(c)   = numel(timedFilesList)/elapsed_s(c);
    inputMBperSec(c) = (timedBytes/1e6)/elapsed_s(c);
    estFull_min(c)   = (nFiles/filesPerSec(c))/60;

    fprintf('\nRESULT: %.2f s\n',elapsed_s(c));
    fprintf('        %.2f TIFF/s\n',filesPerSec(c));
    fprintf('        %.1f input MB/s\n',inputMBperSec(c));
    fprintf('        estimated %.1f min for all %d TIFFs\n', ...
        estFull_min(c),nFiles);

    % Deletion is deliberately OUTSIDE the timed section.
    if DELETE_TEMP_OUTPUTS && isfolder(testOut)
        fprintf('Deleting temporary benchmark output (not timed)...\n');
        rmdir(testOut,'s');
    end
end

%% Results / recommendation
T = table((1:nCfg)',workersToTest(:),readerToTest(:), ...
    elapsed_s,filesPerSec,inputMBperSec,estFull_min, ...
    'VariableNames', ...
    {'Test','Workers','Reader','Seconds','TIFF_per_s','Input_MB_per_s','Estimated_full_min'});

T = sortrows(T,'Seconds','ascend');

fprintf('\n\n============================================================\n');
fprintf('FINAL BENCHMARK RESULT\n');
fprintf('============================================================\n');
disp(T);

bestWorkers = T.Workers(1);
bestReader  = T.Reader(1);

if height(T) >= 2
    marginPct = 100*(T.Seconds(2)-T.Seconds(1))/T.Seconds(1);
else
    marginPct = NaN;
end

fprintf('\n*** WINNER: Workers=%d, Reader="%s" ***\n',bestWorkers,bestReader);
fprintf('Measured throughput: %.2f TIFF/s (%.1f input MB/s)\n', ...
    T.TIFF_per_s(1),T.Input_MB_per_s(1));
fprintf('Projected full-run time: %.1f minutes for %d TIFFs\n', ...
    T.Estimated_full_min(1),nFiles);

if isfinite(marginPct)
    fprintf('Winner is %.1f%% faster than second place in this quick test.\n',marginPct);

    if marginPct < 5
        fprintf(['NOTE: The top two are within 5%%. Treat that as a practical tie.\n' ...
                 'Use the LOWER worker count of those two for the long run unless\n' ...
                 'a longer test separates them clearly.\n']);
    elseif marginPct < 10
        fprintf(['NOTE: The top two are fairly close (<10%%). Storage/network\n' ...
                 'variability may swap their order on a repeated short test.\n']);
    end
end

fprintf('\nUse this in the real call:\n');
fprintf('    Workers=%d, ...\n',bestWorkers);
fprintf('    Reader="%s", ...\n',bestReader);

resultsFile = fullfile(fileparts(outDir), ...
    "pilatus_benchmark_result_" + string(datetime('now','Format','yyyyMMdd_HHmmss')) + ".mat");

save(resultsFile,'T','bestWorkers','bestReader','window','benchFiles','warmFiles');
fprintf('\nTiming table saved to:\n%s\n',resultsFile);

if DELETE_TEMP_OUTPUTS && isfolder(benchRoot)
    try
        rmdir(benchRoot,'s');
    catch
    end
end

fprintf('============================================================\n');
