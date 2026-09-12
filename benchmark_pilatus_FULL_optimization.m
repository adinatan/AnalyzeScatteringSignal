clearvars
close all
clc

%% PILATUS FULL OPTIMIZATION BENCHMARK - RADIAL-ONLY PRODUCTION GOAL
%
% This benchmark searches the meaningful performance space for the workflow
% you said you actually want:
%
%   TIFF -> mask -> per-frame normalization -> even/odd pair difference
%        -> window average -> DeltaI(q,t)
%
% and one final:
%   avgDiff2D = average normalized detector-space pair difference
%
% It does NOT calculate q/phi maps because you explicitly asked to benchmark
% the case where you only want the full 2-D average difference and radial data.
%
% WHAT IT TESTS
% -------------------------------------------------------------------------
% 0) Radial reduction kernel:
%       accumarray vs precomputed sparse-matrix binning
%
% 1) FULL I/O / CPU parallel search:
%       readers = Tiff class, imread
%       workers = 1, 2, 3, 4, 6, 8 process workers
%
% 2) Science/output search on the best THREE I/O configurations:
%       window = 1, 3, 10
%       Nq     = 500, 1000, 1500
%
% 3) GPU search:
%       readers       = Tiff class, imread
%       GPU transfer  = frame-by-frame vs whole-window batch transfer
%       window        = 1, 3, 10
%       Nq            = 500, 1000, 1500
%
% 4) Final validation:
%       the best 2 CPU + best 2 GPU configurations are re-run on a much
%       larger, previously unused acquisition block.
%
% 5) Storage test:
%       the final winner is run directly from X:, then the same size block
%       is copied to C: and run locally.  The script reports direct-X time
%       versus estimated copy-to-local + local-processing time.
%
% IMPORTANT
% -------------------------------------------------------------------------
% - Pool startup, q-bin construction, sparse-operator construction and GPU
%   initialization are excluded from timed processing.
% - Timed quick candidates use different acquisition blocks, so we do not
%   repeatedly benchmark the same cached TIFFs.
% - The final validation uses larger fresh blocks.
% - Normalization sums are calculated in double.
% - Image arithmetic is single precision for speed.
% - The final avgDiff2D accumulators use double precision.
% - Dynamic invalid pixels are NOT silently counted as zero in the radial
%   denominator.  Pair-validity is tracked explicitly.
%
% MATLAB R2026a + Parallel Computing Toolbox.
%
% RUN THIS FILE DIRECTLY.

%% ========================================================================
% USER SETTINGS
% ========================================================================

inDir = "X:\2026-3\RuettSept26\data\fresh_ag_1\run3_UV\";

benchRoot = ...
    "C:\Users\11idbuser\Downloads\Natan\pilatus_fast_diff_matlab\analysis\full_optimization_benchmark";

% Candidate space
READERS         = ["tiff","imread"];
WORKERS_TO_TEST = [1 2 3 4 6 8];
WINDOWS_TO_TEST = [1 3 10];
NQ_TO_TEST      = [500 1000 1500];
GPU_TRANSFER    = ["frame","batch"];

% Broad search: must be divisible by 2*LCM(WINDOWS_TO_TEST)=20.
QUICK_FILES = 120;

% Finalists: also divisible by 20.
FINAL_FILES = 600;

% Local staging test block.
LOCAL_FILES = 300;

% Number of best CPU reader/worker combinations carried into the full
% window x Nq search.
TOP_IO_CONFIGS = 3;

% Number of CPU/GPU finalists for the larger validation.
TOP_CPU_FINALISTS = 2;
TOP_GPU_FINALISTS = 2;

% Radial correction convention.
% "divide"   : Icorr = I ./ (geometry .* polarization)
% "multiply" : Icorr = I .* (geometry .* polarization)
CorrectionMode = "divide";

% The benchmark writes a realistic small result file.  Delete candidate
% outputs immediately after measuring them.
DELETE_CANDIDATE_OUTPUTS = true;

% Run the local C: staging comparison.
RUN_LOCAL_STAGING_TEST = true;

% Location for the temporary local TIFF copy.
localStageDir = fullfile(benchRoot,"local_tiff_stage");

rng(1);  % reproducible ordering when needed

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
names     = names(ord);
fileBytes = fileBytes(ord);

nFiles = numel(paths);

fprintf('\n============================================================\n');
fprintf('PILATUS FULL OPTIMIZATION BENCHMARK\n');
fprintf('============================================================\n');
fprintf('Files found : %d\n',nFiles);
fprintf('Input size  : %.3f TB\n',sum(fileBytes)/1e12);
fprintf('Frame IDs   : %d ... %d\n',frameIDs(1),frameIDs(end));
fprintf('Readers     : %s\n',strjoin(READERS,", "));
fprintf('Workers     : %s\n',mat2str(WORKERS_TO_TEST));
fprintf('Windows     : %s\n',mat2str(WINDOWS_TO_TEST));
fprintf('q bins      : %s\n',mat2str(NQ_TO_TEST));
fprintf('GPU modes   : %s\n',strjoin(GPU_TRANSFER,", "));

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

corrMap = geometry .* polarization;

if CorrectionMode == "divide"

    qStaticValid = detectorValid & ...
                   isfinite(qMap) & ...
                   isfinite(corrMap) & ...
                   corrMap ~= 0;

elseif CorrectionMode == "multiply"

    qStaticValid = detectorValid & ...
                   isfinite(qMap) & ...
                   isfinite(corrMap);

else

    error('CorrectionMode must be "divide" or "multiply".');
end

A0 = read_tiff_any(paths(1),"tiff");

[ny,nx] = size(A0);

if ~isequal(size(A0),size(qMap))
    error('TIFF is %dx%d but qMap is %dx%d.', ...
        ny,nx,size(qMap,1),size(qMap,2));
end

clear A0

fprintf('Detector     : %d x %d\n',ny,nx);
fprintf('Norm ROI     : %d good pixels\n',numel(normIdx));

%% ========================================================================
% OUTPUT DIRECTORY
% ========================================================================

if isfolder(benchRoot)
    try
        rmdir(benchRoot,'s');
    catch
    end
end

mkdir(benchRoot);

%% ========================================================================
% GPU AVAILABILITY
% ========================================================================

gpuAvailable = false;

try
    gd = gpuDevice;
    gpuAvailable = true;

    fprintf('GPU          : %s, %.1f GB\n', ...
        gd.Name,gd.TotalMemory/1e9);

catch ME
    fprintf('GPU          : unavailable (%s)\n',ME.message);
end

%% ========================================================================
% CHECK BENCHMARK BLOCK SIZES
% ========================================================================

commonGroup = 2;

for W = WINDOWS_TO_TEST
    commonGroup = lcm(commonGroup,2*W);
end

QUICK_FILES = floor(QUICK_FILES/commonGroup)*commonGroup;
FINAL_FILES = floor(FINAL_FILES/commonGroup)*commonGroup;
LOCAL_FILES = floor(LOCAL_FILES/commonGroup)*commonGroup;

if QUICK_FILES < commonGroup
    QUICK_FILES = commonGroup;
end

if FINAL_FILES < commonGroup
    FINAL_FILES = commonGroup;
end

if LOCAL_FILES < commonGroup
    LOCAL_FILES = commonGroup;
end

fprintf('Quick block  : %d TIFFs (~%.2f GB)\n', ...
    QUICK_FILES,median(fileBytes)*QUICK_FILES/1e9);

fprintf('Final block  : %d TIFFs (~%.2f GB)\n', ...
    FINAL_FILES,median(fileBytes)*FINAL_FILES/1e9);

%% ========================================================================
% STAGE 0 - RADIAL REDUCTION KERNEL: accumarray vs sparse
% ========================================================================

fprintf('\n============================================================\n');
fprintf('STAGE 0: RADIAL REDUCTION KERNEL\n');
fprintf('============================================================\n');

kernelWindow = 3;
kernelNq = 1000;

% Make one real normalized window in memory.
kernelPaths = paths(1:2*kernelWindow);

dExample = make_one_window_cpu( ...
    kernelPaths,kernelWindow,"tiff", ...
    detectorValid,normIdx,ny,nx);

opAccum = prepare_q_operator( ...
    qMap,corrMap,qStaticValid,CorrectionMode,kernelNq,"accumarray");

opSparse = prepare_q_operator( ...
    qMap,corrMap,qStaticValid,CorrectionMode,kernelNq,"sparse");

% Warm each method.
tmp = radial_cpu(dExample,opAccum);
tmp = radial_cpu(dExample,opSparse);
clear tmp

N_KERNEL_REPEAT = 20;

tic
for k = 1:N_KERNEL_REPEAT
    tmp = radial_cpu(dExample,opAccum);
end
tAccum = toc;
clear tmp

tic
for k = 1:N_KERNEL_REPEAT
    tmp = radial_cpu(dExample,opSparse);
end
tSparse = toc;
clear tmp dExample

fprintf('accumarray: %.4f s for %d reductions\n',tAccum,N_KERNEL_REPEAT);
fprintf('sparse    : %.4f s for %d reductions\n',tSparse,N_KERNEL_REPEAT);

if tSparse < tAccum
    CPU_RADIAL_METHOD = "sparse";
else
    CPU_RADIAL_METHOD = "accumarray";
end

fprintf('Selected CPU radial kernel: %s\n',upper(CPU_RADIAL_METHOD));

%% ========================================================================
% BENCHMARK RESULT STORAGE
% ========================================================================

emptyResult = struct( ...
    'Stage',"", ...
    'Platform',"", ...
    'Reader',"", ...
    'Workers',NaN, ...
    'Window',NaN, ...
    'Nq',NaN, ...
    'RadialMethod',"", ...
    'GPUTransfer',"", ...
    'Seconds',NaN, ...
    'TIFF_per_s',NaN, ...
    'Input_MB_per_s',NaN, ...
    'Output_MB',NaN, ...
    'Projected_full_min',NaN, ...
    'Success',false, ...
    'Message',"");

results = repmat(emptyResult,0,1);

% Each timed quick candidate gets a fresh contiguous source block.
nextQuickStart = 1;

% Warm-up block is deliberately not part of timing.
warmPaths = paths(1:min(60,nFiles));

%% ========================================================================
% STAGE 1 - FULL READER x WORKER SEARCH
% ========================================================================

fprintf('\n============================================================\n');
fprintf('STAGE 1: READER x PROCESS-WORKER SEARCH\n');
fprintf('Representative science setting: W=3, Nq=1000\n');
fprintf('============================================================\n');

stage1Start = numel(results)+1;

for workers = WORKERS_TO_TEST

    for reader = READERS

        cfg = make_cfg( ...
            "CPU",reader,workers,3,1000, ...
            CPU_RADIAL_METHOD,"");

        fprintf('\nCPU reader=%s, workers=%d\n',upper(reader),workers);

        [poolOK,poolMsg] = ensure_process_pool(workers);

        if ~poolOK

            r = emptyResult;
            r.Stage = "IO";
            r.Platform = "CPU";
            r.Reader = reader;
            r.Workers = workers;
            r.Window = 3;
            r.Nq = 1000;
            r.RadialMethod = CPU_RADIAL_METHOD;
            r.Success = false;
            r.Message = poolMsg;

            results(end+1) = r; %#ok<SAGROW>

            fprintf('SKIPPED: %s\n',poolMsg);
            continue
        end

        op = prepare_q_operator( ...
            qMap,corrMap,qStaticValid,CorrectionMode, ...
            cfg.Nq,cfg.RadialMethod);

        % Warm this reader / worker architecture outside timing.
        try
            run_cpu_radial( ...
                warmPaths(1:floor(numel(warmPaths)/(2*cfg.Window))*2*cfg.Window), ...
                cfg,detectorValid,normIdx,ny,nx,op);
        catch ME
            fprintf('Warm-up warning: %s\n',ME.message);
        end

        [blockPaths,blockBytes,nextQuickStart] = ...
            take_fresh_block(paths,fileBytes,nextQuickStart,QUICK_FILES);

        outFile = fullfile(benchRoot, ...
            sprintf('stage1_%s_w%d.mat',reader,workers));

        r = benchmark_one( ...
            blockPaths,blockBytes,cfg, ...
            detectorValid,normIdx,ny,nx, ...
            qMap,corrMap,qStaticValid,CorrectionMode, ...
            outFile);

        r.Stage = "IO";
        results(end+1) = r; %#ok<SAGROW>

        print_result(r);

        if DELETE_CANDIDATE_OUTPUTS && isfile(outFile)
            delete(outFile);
        end
    end
end

stage1End = numel(results);

T1 = struct2table(results(stage1Start:stage1End));
T1good = T1(T1.Success,:);

T1good = sortrows(T1good,'Seconds','ascend');

fprintf('\nSTAGE 1 RANKING\n');
disp(T1good(:, ...
    {'Reader','Workers','Seconds','TIFF_per_s','Input_MB_per_s','Projected_full_min'}));

if isempty(T1good)
    error('No CPU reader/worker configuration completed successfully.');
end

nCarry = min(TOP_IO_CONFIGS,height(T1good));

topIO = T1good(1:nCarry,:);

fprintf('Carrying the top %d I/O configurations into Stage 2.\n',nCarry);

%% ========================================================================
% STAGE 2 - WINDOW x Nq ON TOP I/O CONFIGURATIONS
% ========================================================================

fprintf('\n============================================================\n');
fprintf('STAGE 2: WINDOW x q-RESOLUTION SEARCH\n');
fprintf('============================================================\n');

stage2Start = numel(results)+1;

for ii = 1:height(topIO)

    reader  = topIO.Reader(ii);
    workers = topIO.Workers(ii);

    [poolOK,poolMsg] = ensure_process_pool(workers);

    if ~poolOK
        fprintf('Skipping top-I/O config because pool failed: %s\n',poolMsg);
        continue
    end

    for W = WINDOWS_TO_TEST

        for nq = NQ_TO_TEST

            cfg = make_cfg( ...
                "CPU",reader,workers,W,nq, ...
                CPU_RADIAL_METHOD,"");

            fprintf('\nCPU %s, workers=%d, W=%d, Nq=%d\n', ...
                upper(reader),workers,W,nq);

            [blockPaths,blockBytes,nextQuickStart] = ...
                take_fresh_block(paths,fileBytes,nextQuickStart,QUICK_FILES);

            outFile = fullfile(benchRoot, ...
                sprintf('stage2_%s_wk%d_W%d_q%d.mat', ...
                reader,workers,W,nq));

            r = benchmark_one( ...
                blockPaths,blockBytes,cfg, ...
                detectorValid,normIdx,ny,nx, ...
                qMap,corrMap,qStaticValid,CorrectionMode, ...
                outFile);

            r.Stage = "CPU_GRID";
            results(end+1) = r; %#ok<SAGROW>

            print_result(r);

            if DELETE_CANDIDATE_OUTPUTS && isfile(outFile)
                delete(outFile);
            end
        end
    end
end

stage2End = numel(results);

T2 = struct2table(results(stage2Start:stage2End));
T2good = T2(T2.Success,:);

T2good = sortrows(T2good,'Seconds','ascend');

fprintf('\nSTAGE 2 CPU RANKING\n');

disp(T2good(:, ...
    {'Reader','Workers','Window','Nq','Seconds','TIFF_per_s', ...
     'Input_MB_per_s','Projected_full_min'}));

%% ========================================================================
% STAGE 3 - GPU SEARCH: reader x transfer-mode x window x Nq
% ========================================================================

fprintf('\n============================================================\n');
fprintf('STAGE 3: GPU SEARCH\n');
fprintf('============================================================\n');

stage3Start = numel(results)+1;

if gpuAvailable

    % GPU uses one MATLAB process.  Close process pool so it does not occupy
    % CPU/RAM while benchmarking the GPU architecture.
    p = gcp('nocreate');

    if ~isempty(p)
        delete(p);
    end

    for reader = READERS

        for gpuMode = GPU_TRANSFER

            % One architecture warm-up outside timing.
            warmCfg = make_cfg( ...
                "GPU",reader,0,3,1000, ...
                "gpu_accumarray",gpuMode);

            warmOp = prepare_q_operator( ...
                qMap,corrMap,qStaticValid,CorrectionMode, ...
                warmCfg.Nq,"accumarray");

            try
                run_gpu_radial( ...
                    warmPaths(1:floor(numel(warmPaths)/(2*warmCfg.Window))*2*warmCfg.Window), ...
                    warmCfg,detectorValid,normIdx,ny,nx,warmOp);

                wait(gpuDevice);

            catch ME
                fprintf('GPU warm-up %s/%s failed: %s\n', ...
                    reader,gpuMode,ME.message);
            end

            for W = WINDOWS_TO_TEST

                for nq = NQ_TO_TEST

                    cfg = make_cfg( ...
                        "GPU",reader,0,W,nq, ...
                        "gpu_accumarray",gpuMode);

                    fprintf('\nGPU reader=%s, transfer=%s, W=%d, Nq=%d\n', ...
                        upper(reader),upper(gpuMode),W,nq);

                    [blockPaths,blockBytes,nextQuickStart] = ...
                        take_fresh_block(paths,fileBytes,nextQuickStart,QUICK_FILES);

                    outFile = fullfile(benchRoot, ...
                        sprintf('stage3_gpu_%s_%s_W%d_q%d.mat', ...
                        reader,gpuMode,W,nq));

                    r = benchmark_one( ...
                        blockPaths,blockBytes,cfg, ...
                        detectorValid,normIdx,ny,nx, ...
                        qMap,corrMap,qStaticValid,CorrectionMode, ...
                        outFile);

                    r.Stage = "GPU_GRID";
                    results(end+1) = r; %#ok<SAGROW>

                    print_result(r);

                    if DELETE_CANDIDATE_OUTPUTS && isfile(outFile)
                        delete(outFile);
                    end
                end
            end
        end
    end

else

    fprintf('GPU stage skipped because no usable GPU was found.\n');

end

stage3End = numel(results);

if stage3End >= stage3Start

    T3 = struct2table(results(stage3Start:stage3End));
    T3good = T3(T3.Success,:);

    if ~isempty(T3good)

        T3good = sortrows(T3good,'Seconds','ascend');

        fprintf('\nSTAGE 3 GPU RANKING\n');

        disp(T3good(:, ...
            {'Reader','GPUTransfer','Window','Nq','Seconds','TIFF_per_s', ...
             'Input_MB_per_s','Projected_full_min'}));

    else
        T3good = table();
    end

else

    T3good = table();

end

%% ========================================================================
% STAGE 4 - LARGE-BLOCK FINAL VALIDATION
% ========================================================================

fprintf('\n============================================================\n');
fprintf('STAGE 4: LARGE-BLOCK FINAL VALIDATION\n');
fprintf('============================================================\n');

% Top CPU candidates.
T2good = sortrows(T2good,'Seconds','ascend');

nCPUFinal = min(TOP_CPU_FINALISTS,height(T2good));
cpuFinal = T2good(1:nCPUFinal,:);

% Top GPU candidates.
if ~isempty(T3good)
    nGPUFinal = min(TOP_GPU_FINALISTS,height(T3good));
    gpuFinal = T3good(1:nGPUFinal,:);
else
    nGPUFinal = 0;
    gpuFinal = table();
end

finalResults = repmat(emptyResult,0,1);

% CPU finalists
for k = 1:nCPUFinal

    cfg = cfg_from_table_row(cpuFinal(k,:));

    fprintf('\nFINAL CPU: %s workers=%d W=%d Nq=%d\n', ...
        upper(cfg.Reader),cfg.Workers,cfg.Window,cfg.Nq);

    [poolOK,poolMsg] = ensure_process_pool(cfg.Workers);

    if ~poolOK
        fprintf('SKIPPED: %s\n',poolMsg);
        continue
    end

    [blockPaths,blockBytes,nextQuickStart] = ...
        take_fresh_block(paths,fileBytes,nextQuickStart,FINAL_FILES);

    outFile = fullfile(benchRoot,sprintf('final_cpu_%d.mat',k));

    r = benchmark_one( ...
        blockPaths,blockBytes,cfg, ...
        detectorValid,normIdx,ny,nx, ...
        qMap,corrMap,qStaticValid,CorrectionMode, ...
        outFile);

    r.Stage = "FINAL";
    finalResults(end+1) = r; %#ok<SAGROW>

    print_result(r);

    if DELETE_CANDIDATE_OUTPUTS && isfile(outFile)
        delete(outFile);
    end
end

% GPU finalists
if nGPUFinal > 0

    p = gcp('nocreate');
    if ~isempty(p)
        delete(p);
    end

    for k = 1:nGPUFinal

        cfg = cfg_from_table_row(gpuFinal(k,:));

        fprintf('\nFINAL GPU: %s %s W=%d Nq=%d\n', ...
            upper(cfg.Reader),upper(cfg.GPUTransfer),cfg.Window,cfg.Nq);

        [blockPaths,blockBytes,nextQuickStart] = ...
            take_fresh_block(paths,fileBytes,nextQuickStart,FINAL_FILES);

        outFile = fullfile(benchRoot,sprintf('final_gpu_%d.mat',k));

        r = benchmark_one( ...
            blockPaths,blockBytes,cfg, ...
            detectorValid,normIdx,ny,nx, ...
            qMap,corrMap,qStaticValid,CorrectionMode, ...
            outFile);

        r.Stage = "FINAL";
        finalResults(end+1) = r; %#ok<SAGROW>

        print_result(r);

        if DELETE_CANDIDATE_OUTPUTS && isfile(outFile)
            delete(outFile);
        end
    end
end

TF = struct2table(finalResults);
TFgood = TF(TF.Success,:);

TFgood = sortrows(TFgood,'Seconds','ascend');

fprintf('\nFINAL VALIDATION RANKING\n');

disp(TFgood(:, ...
    {'Platform','Reader','Workers','Window','Nq','RadialMethod', ...
     'GPUTransfer','Seconds','TIFF_per_s','Input_MB_per_s', ...
     'Projected_full_min'}));

if isempty(TFgood)
    error('No final configuration completed successfully.');
end

winner = TFgood(1,:);

fprintf('\n============================================================\n');
fprintf('WINNING PRODUCTION CONFIGURATION\n');
fprintf('============================================================\n');

disp(winner(:, ...
    {'Platform','Reader','Workers','Window','Nq','RadialMethod', ...
     'GPUTransfer','TIFF_per_s','Input_MB_per_s','Projected_full_min'}));

%% ========================================================================
% STAGE 5 - X: DIRECT VS LOCAL C: STAGING
% ========================================================================

localStageSummary = table();

if RUN_LOCAL_STAGING_TEST

    fprintf('\n============================================================\n');
    fprintf('STAGE 5: X: DIRECT vs LOCAL C: STAGING\n');
    fprintf('============================================================\n');

    cfgWin = cfg_from_table_row(winner);

    if cfgWin.Platform == "CPU"

        [poolOK,poolMsg] = ensure_process_pool(cfgWin.Workers);

        if ~poolOK
            fprintf('Local staging test skipped: %s\n',poolMsg);
        end
    else

        p = gcp('nocreate');
        if ~isempty(p)
            delete(p);
        end
    end

    [stagePaths,stageBytes,nextQuickStart] = ...
        take_fresh_block(paths,fileBytes,nextQuickStart,LOCAL_FILES);

    % Direct X:
    directOut = fullfile(benchRoot,'localtest_directX.mat');

    rDirect = benchmark_one( ...
        stagePaths,stageBytes,cfgWin, ...
        detectorValid,normIdx,ny,nx, ...
        qMap,corrMap,qStaticValid,CorrectionMode, ...
        directOut);

    if DELETE_CANDIDATE_OUTPUTS && isfile(directOut)
        delete(directOut);
    end

    % Copy same block to local C:
    if isfolder(localStageDir)
        rmdir(localStageDir,'s');
    end
    mkdir(localStageDir);

    fprintf('Copying %d TIFFs to local staging folder...\n',LOCAL_FILES);

    copyTic = tic;

    localPaths = strings(LOCAL_FILES,1);

    for k = 1:LOCAL_FILES

        [~,base,ext] = fileparts(stagePaths(k));

        dst = fullfile(localStageDir,base+ext);

        ok = copyfile(char(stagePaths(k)),char(dst));

        if ~ok
            error('Failed to copy %s',stagePaths(k));
        end

        localPaths(k) = string(dst);
    end

    copySeconds = toc(copyTic);

    localOut = fullfile(benchRoot,'localtest_localC.mat');

    rLocal = benchmark_one( ...
        localPaths,stageBytes,cfgWin, ...
        detectorValid,normIdx,ny,nx, ...
        qMap,corrMap,qStaticValid,CorrectionMode, ...
        localOut);

    if DELETE_CANDIDATE_OUTPUTS && isfile(localOut)
        delete(localOut);
    end

    stageGB = sum(stageBytes)/1e9;
    fullGB  = sum(fileBytes)/1e9;

    copyGBps = stageGB/copySeconds;

    projectedCopyMin = (fullGB/copyGBps)/60;

    projectedLocalProcessMin = ...
        (nFiles/rLocal.TIFF_per_s)/60;

    projectedStagePlusProcessMin = ...
        projectedCopyMin + projectedLocalProcessMin;

    projectedDirectMin = ...
        (nFiles/rDirect.TIFF_per_s)/60;

    localStageSummary = table( ...
        rDirect.TIFF_per_s, ...
        rLocal.TIFF_per_s, ...
        copySeconds, ...
        copyGBps, ...
        projectedDirectMin, ...
        projectedCopyMin, ...
        projectedLocalProcessMin, ...
        projectedStagePlusProcessMin, ...
        'VariableNames', ...
        {'DirectX_TIFF_per_s','LocalC_TIFF_per_s','Copy_seconds', ...
         'Copy_GB_per_s','Projected_DirectX_min','Projected_Copy_min', ...
         'Projected_LocalProcess_min','Projected_CopyPlusProcess_min'});

    disp(localStageSummary);

    if projectedStagePlusProcessMin < projectedDirectMin

        fprintf(['LOCAL STAGING WINS including measured copy time.\n' ...
                 'For long runs, copy to C: first, then process locally.\n']);

    else

        fprintf(['DIRECT X: PROCESSING WINS when copy time is included.\n' ...
                 'There is no reason to stage the entire dataset locally for speed.\n']);
    end

    try
        rmdir(localStageDir,'s');
    catch
    end
end

%% ========================================================================
% SAVE ALL BENCHMARK RESULTS
% ========================================================================

Tall = struct2table(results);

resultsFile = fullfile(benchRoot,'FULL_BENCHMARK_RESULTS.mat');

save(resultsFile, ...
    'Tall','TFgood','winner','localStageSummary', ...
    'CPU_RADIAL_METHOD','READERS','WORKERS_TO_TEST', ...
    'WINDOWS_TO_TEST','NQ_TO_TEST','GPU_TRANSFER', ...
    'QUICK_FILES','FINAL_FILES','LOCAL_FILES');

fprintf('\n============================================================\n');
fprintf('BENCHMARK COMPLETE\n');
fprintf('============================================================\n');

fprintf('Full results saved to:\n%s\n',resultsFile);

fprintf('\nFor production use the FINAL validation winner, not merely the\n');
fprintf('fastest 120-file quick candidate.\n');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function cfg = make_cfg(platform,reader,workers,window,nq,radialMethod,gpuTransfer)

    cfg = struct;

    cfg.Platform     = string(platform);
    cfg.Reader       = string(reader);
    cfg.Workers      = workers;
    cfg.Window       = window;
    cfg.Nq           = nq;
    cfg.RadialMethod = string(radialMethod);
    cfg.GPUTransfer  = string(gpuTransfer);
end


function cfg = cfg_from_table_row(Trow)

    cfg = struct;

    cfg.Platform     = string(Trow.Platform);
    cfg.Reader       = string(Trow.Reader);
    cfg.Workers      = Trow.Workers;
    cfg.Window       = Trow.Window;
    cfg.Nq           = Trow.Nq;
    cfg.RadialMethod = string(Trow.RadialMethod);
    cfg.GPUTransfer  = string(Trow.GPUTransfer);
end


function [ok,msg] = ensure_process_pool(workers)

    ok = true;
    msg = "";

    try

        p = gcp('nocreate');

        if workers <= 1

            if ~isempty(p)
                delete(p);
            end

            return
        end

        if isempty(p) || p.NumWorkers ~= workers

            if ~isempty(p)
                delete(p);
            end

            parpool("Processes",workers);
        end

    catch ME

        ok = false;
        msg = string(ME.message);

    end
end


function [blockPaths,blockBytes,nextStart] = ...
    take_fresh_block(paths,fileBytes,startIdx,nBlock)

    nFiles = numel(paths);

    if startIdx+nBlock-1 > nFiles
        error(['Benchmark has exhausted fresh acquisition blocks. ' ...
               'Reduce candidate lists or block sizes.']);
    end

    idx = startIdx:startIdx+nBlock-1;

    blockPaths = paths(idx);
    blockBytes = fileBytes(idx);

    nextStart = startIdx+nBlock;
end


function r = benchmark_one( ...
    blockPaths,blockBytes,cfg, ...
    detectorValid,normIdx,ny,nx, ...
    qMap,corrMap,qStaticValid,CorrectionMode, ...
    outFile)

    r = struct( ...
        'Stage',"", ...
        'Platform',cfg.Platform, ...
        'Reader',cfg.Reader, ...
        'Workers',cfg.Workers, ...
        'Window',cfg.Window, ...
        'Nq',cfg.Nq, ...
        'RadialMethod',cfg.RadialMethod, ...
        'GPUTransfer',cfg.GPUTransfer, ...
        'Seconds',NaN, ...
        'TIFF_per_s',NaN, ...
        'Input_MB_per_s',NaN, ...
        'Output_MB',NaN, ...
        'Projected_full_min',NaN, ...
        'Success',false, ...
        'Message',"");

    try

        % q operator construction is a one-time startup operation and is
        % deliberately excluded from the timed processing.
        if cfg.Platform == "CPU"

            op = prepare_q_operator( ...
                qMap,corrMap,qStaticValid,CorrectionMode, ...
                cfg.Nq,cfg.RadialMethod);

        else

            op = prepare_q_operator( ...
                qMap,corrMap,qStaticValid,CorrectionMode, ...
                cfg.Nq,"accumarray");

            gd = gpuDevice;
            wait(gd);
        end

        drawnow

        t0 = tic;

        if cfg.Platform == "CPU"

            R = run_cpu_radial( ...
                blockPaths,cfg,detectorValid,normIdx,ny,nx,op);

        else

            R = run_gpu_radial( ...
                blockPaths,cfg,detectorValid,normIdx,ny,nx,op);

            wait(gpuDevice);
        end

        R.q = op.qCenter;
        R.window = cfg.Window;

        save(outFile,'R','-v7','-nocompression');

        sec = toc(t0);

        info = dir(outFile);

        bytesIn = sum(blockBytes);
        nF = numel(blockPaths);

        r.Seconds = sec;
        r.TIFF_per_s = nF/sec;
        r.Input_MB_per_s = bytesIn/1e6/sec;
        r.Output_MB = info.bytes/1e6;

        % Filled against the full 14.5k run by caller's rate.
        % nF/sec is the only necessary quantity.
        % The main script uses this field directly for ranking.
        %
        % Since benchmark_one does not know the global nFiles, infer from
        % caller later is inconvenient; use NaN here then main print still
        % has throughput.  The next line uses the known 14520-style scale
        % only if passed through global?  Avoid globals; compute from rate
        % in caller's result postprocessing is not available.
        %
        % Instead use the source directory total via evalin safely.
        fullN = evalin('base','nFiles');

        r.Projected_full_min = fullN/r.TIFF_per_s/60;

        r.Success = true;

    catch ME

        r.Message = string(getReport(ME,'basic','hyperlinks','off'));

    end
end


function print_result(r)

    if r.Success

        fprintf('%.2f s | %.2f TIFF/s | %.1f MB/s | projected %.1f min\n', ...
            r.Seconds,r.TIFF_per_s,r.Input_MB_per_s,r.Projected_full_min);

    else

        fprintf('FAILED: %s\n',r.Message);
    end
end


function op = prepare_q_operator( ...
    qMap,corrMap,qStaticValid,CorrectionMode,nq,method)

    qGood = qMap(qStaticValid);

    qEdges = linspace(min(qGood),max(qGood),nq+1);

    qBinMap = discretize(qMap,qEdges);

    qMask = qStaticValid & ~isnan(qBinMap);

    bins = double(qBinMap(qMask));

    if CorrectionMode == "divide"
        factor = 1 ./ corrMap(qMask);
    else
        factor = corrMap(qMask);
    end

    op = struct;

    op.mask = qMask;
    op.bins = bins;
    op.factor = single(factor);
    op.Nq = nq;
    op.qCenter = (qEdges(1:end-1)+qEdges(2:end))/2;
    op.method = string(method);

    if op.method == "sparse"

        nPix = numel(bins);

        % One nonzero per valid detector pixel.
        op.S = sparse( ...
            bins, ...
            (1:nPix)', ...
            ones(nPix,1), ...
            nq,nPix);

    else

        op.S = [];
    end
end


function R = run_cpu_radial( ...
    paths,cfg,detectorValid,normIdx,ny,nx,op)

    nF = numel(paths);
    W = cfg.Window;

    nOut = floor(nF/(2*W));

    deltaIq = nan(cfg.Nq,nOut,'single');
    normCell = cell(nOut,1);

    diffSumAll   = zeros(ny,nx,'double');
    diffCountAll = zeros(ny,nx,'double');

    if cfg.Workers > 1

        parfor g = 1:nOut

            [iq,normVals,localDiffSum,localDiffCount] = ...
                process_group_cpu( ...
                    paths,g,W,cfg.Reader, ...
                    detectorValid,normIdx,ny,nx,op);

            deltaIq(:,g) = iq;
            normCell{g} = normVals;

            diffSumAll   = diffSumAll   + localDiffSum;
            diffCountAll = diffCountAll + localDiffCount;
        end

    else

        for g = 1:nOut

            [iq,normVals,localDiffSum,localDiffCount] = ...
                process_group_cpu( ...
                    paths,g,W,cfg.Reader, ...
                    detectorValid,normIdx,ny,nx,op);

            deltaIq(:,g) = iq;
            normCell{g} = normVals;

            diffSumAll   = diffSumAll   + localDiffSum;
            diffCountAll = diffCountAll + localDiffCount;
        end
    end

    normTrace = nan(nF,1);

    for g = 1:nOut

        first = (g-1)*2*W+1;

        normTrace(first:first+2*W-1) = normCell{g};
    end

    avgDiff2D = diffSumAll ./ max(diffCountAll,1);
    avgDiff2D(diffCountAll==0) = NaN;

    R = struct;

    R.avgDiff2D = avgDiff2D;
    R.deltaIq = deltaIq;
    R.normalization = normTrace;
end


function [iq,normVals,diffSum,diffCount] = process_group_cpu( ...
    paths,g,W,reader,detectorValid,normIdx,ny,nx,op)

    first = (g-1)*2*W+1;

    normVals = nan(2*W,1);

    diffSum   = zeros(ny,nx,'double');
    diffCount = zeros(ny,nx,'double');

    windowSum   = zeros(ny,nx,'single');
    windowCount = zeros(ny,nx,'single');

    for j = 1:W

        i1 = first+2*j-2;
        i2 = i1+1;

        raw1 = read_tiff_any(paths(i1),reader);
        raw2 = read_tiff_any(paths(i2),reader);

        valid1 = isfinite(raw1) & raw1 >= 0;
        valid2 = isfinite(raw2) & raw2 >= 0;

        idx1 = normIdx(valid1(normIdx));
        idx2 = normIdx(valid2(normIdx));

        n1 = sum(double(raw1(idx1)));
        n2 = sum(double(raw2(idx2)));

        if ~(isfinite(n1) && n1>0 && isfinite(n2) && n2>0)
            error('Invalid normalization scalar in output group %d.',g);
        end

        normVals(2*j-1) = n1;
        normVals(2*j)   = n2;

        pairValid = detectorValid & valid1 & valid2;

        f1 = single(raw1) / single(n1);
        f2 = single(raw2) / single(n2);

        d = f2-f1;

        d(~pairValid) = 0;

        diffSum   = diffSum   + double(d);
        diffCount = diffCount + double(pairValid);

        windowSum   = windowSum   + d;
        windowCount = windowCount + single(pairValid);
    end

    dWindow = windowSum ./ max(windowCount,single(1));
    dWindow(windowCount==0) = NaN;

    iq = radial_cpu(dWindow,op);
end


function R = run_gpu_radial( ...
    paths,cfg,detectorValid,normIdx,ny,nx,op)

    nF = numel(paths);
    W = cfg.Window;

    nOut = floor(nF/(2*W));

    deltaIq = nan(cfg.Nq,nOut,'single');
    normTrace = nan(nF,1);

    gDiffSumAll   = gpuArray.zeros(ny,nx,'single');
    gDiffCountAll = gpuArray.zeros(ny,nx,'single');

    gMask   = gpuArray(op.mask);
    gBins   = gpuArray(op.bins);
    gFactor = gpuArray(op.factor);

    for g = 1:nOut

        first = (g-1)*2*W+1;

        if cfg.GPUTransfer == "batch"

            stack = nan(ny,nx,2*W,'single');

            for k = 1:2*W

                idx = first+k-1;

                raw = read_tiff_any(paths(idx),cfg.Reader);

                valid = isfinite(raw) & raw >= 0;

                normGood = normIdx(valid(normIdx));

                nrm = sum(double(raw(normGood)));

                if ~(isfinite(nrm) && nrm>0)
                    error('Invalid normalization scalar in GPU group %d.',g);
                end

                normTrace(idx) = nrm;

                frame = single(raw)/single(nrm);

                frame(~(detectorValid & valid)) = NaN;

                stack(:,:,k) = frame;
            end

            gStack = gpuArray(stack);

            gOdd  = gStack(:,:,1:2:end);
            gEven = gStack(:,:,2:2:end);

            gPair = gEven-gOdd;

            gValidPair = isfinite(gPair);

            gPair(~gValidPair) = 0;

            gLocalSum   = sum(gPair,3);
            gLocalCount = sum(single(gValidPair),3);

            gDiffSumAll   = gDiffSumAll   + gLocalSum;
            gDiffCountAll = gDiffCountAll + gLocalCount;

            gDWindow = gLocalSum ./ max(gLocalCount,single(1));
            gDWindow(gLocalCount==0) = NaN;

        else

            gLocalSum   = gpuArray.zeros(ny,nx,'single');
            gLocalCount = gpuArray.zeros(ny,nx,'single');

            for j = 1:W

                i1 = first+2*j-2;
                i2 = i1+1;

                raw1 = read_tiff_any(paths(i1),cfg.Reader);
                raw2 = read_tiff_any(paths(i2),cfg.Reader);

                valid1 = isfinite(raw1) & raw1 >= 0;
                valid2 = isfinite(raw2) & raw2 >= 0;

                idx1 = normIdx(valid1(normIdx));
                idx2 = normIdx(valid2(normIdx));

                n1 = sum(double(raw1(idx1)));
                n2 = sum(double(raw2(idx2)));

                if ~(isfinite(n1) && n1>0 && isfinite(n2) && n2>0)
                    error('Invalid normalization scalar in GPU group %d.',g);
                end

                normTrace(i1) = n1;
                normTrace(i2) = n2;

                f1 = single(raw1)/single(n1);
                f2 = single(raw2)/single(n2);

                pairValid = detectorValid & valid1 & valid2;

                gD = gpuArray(f2-f1);
                gV = gpuArray(pairValid);

                gD(~gV) = 0;

                gLocalSum   = gLocalSum   + gD;
                gLocalCount = gLocalCount + single(gV);
            end

            gDiffSumAll   = gDiffSumAll   + gLocalSum;
            gDiffCountAll = gDiffCountAll + gLocalCount;

            gDWindow = gLocalSum ./ max(gLocalCount,single(1));
            gDWindow(gLocalCount==0) = NaN;
        end

        gIq = radial_gpu(gDWindow,gMask,gBins,gFactor,cfg.Nq);

        deltaIq(:,g) = gather(gIq);
    end

    gAvg = gDiffSumAll ./ max(gDiffCountAll,single(1));
    gAvg(gDiffCountAll==0) = NaN;

    R = struct;

    R.avgDiff2D = double(gather(gAvg));
    R.deltaIq = deltaIq;
    R.normalization = normTrace;
end


function out = radial_cpu(I,op)

    v = I(op.mask);

    good = isfinite(v);

    if op.method == "sparse"

        values = double(v) .* double(op.factor);

        values(~good) = 0;

        sums = op.S * values;
        counts = op.S * double(good);

    else

        b = op.bins(good);

        values = double(v(good)) .* double(op.factor(good));

        sums = accumarray( ...
            b,values,[op.Nq 1],@sum,0);

        counts = accumarray( ...
            b,1,[op.Nq 1],@sum,0);
    end

    z = sums ./ max(counts,1);
    z(counts==0) = NaN;

    out = single(z);
end


function out = radial_gpu(I,gMask,gBins,gFactor,nq)

    v = I(gMask);

    good = isfinite(v);

    b = gBins(good);

    values = single(v(good)) .* gFactor(good);

    sums = accumarray( ...
        b,values,[nq 1],@sum,single(0));

    counts = accumarray( ...
        b,ones(size(values),'single'), ...
        [nq 1],@sum,single(0));

    out = sums ./ max(counts,single(1));
    out(counts==0) = NaN;
end


function dWindow = make_one_window_cpu( ...
    paths,W,reader,detectorValid,normIdx,ny,nx)

    windowSum   = zeros(ny,nx,'single');
    windowCount = zeros(ny,nx,'single');

    for j = 1:W

        raw1 = read_tiff_any(paths(2*j-1),reader);
        raw2 = read_tiff_any(paths(2*j),reader);

        valid1 = isfinite(raw1) & raw1>=0;
        valid2 = isfinite(raw2) & raw2>=0;

        idx1 = normIdx(valid1(normIdx));
        idx2 = normIdx(valid2(normIdx));

        n1 = sum(double(raw1(idx1)));
        n2 = sum(double(raw2(idx2)));

        pairValid = detectorValid & valid1 & valid2;

        d = single(raw2)/single(n2) - ...
            single(raw1)/single(n1);

        d(~pairValid) = 0;

        windowSum   = windowSum+d;
        windowCount = windowCount+single(pairValid);
    end

    dWindow = windowSum ./ max(windowCount,single(1));
    dWindow(windowCount==0) = NaN;
end


function A = read_tiff_any(filename,reader)

    if reader == "tiff"

        t = Tiff(char(filename),'r');
        A = t.read();
        t.close();

    elseif reader == "imread"

        A = imread(char(filename));

    else

        error('Unknown TIFF reader: %s',reader);
    end
end
