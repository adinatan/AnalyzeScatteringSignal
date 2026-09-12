clear all
close all
clc

%% PILATUS one-pass normalized difference + q/phi analysis
% EDIT these paths and the hard-coded normalization/mask/calibration sections.
% EDIT THESE PATHS. Input must contain ONE numbered acquisition series.
inDir  = "X:\2026-3\RuettSept26\data\fresh_ag_1\run3_UV\";
%outDir = "C:\Users\11idbuser\Downloads\Natan\pilatus_fast_diff_matlab\analysis\run001_w3";  % New or empty folder, not the input folder.
outDir = "C:\Users\11idbuser\Downloads\Natan\pilatus_fast_diff_matlab\analysis\run003_w3"
% 1 -> normalized I2-I1, I4-I3, ...
% 3 -> mean of three normalized pair differences per output, etc.
window = 3;

[files,frameIDs,fileBytes] = pilatus_filelist(inDir,"*.tif");
fprintf('Input: %d files, %.3f TB\n',numel(files),sum(fileBytes)/1e12);
fprintf('First two file IDs: %d, %d\n',frameIDs(1),frameIDs(2));

%% Load your detector calibration arrays
% Replace this with your actual calibration file/variables.
% All maps must have EXACTLY the same [rows,cols] size as one PILATUS TIFF.
load("qgp.mat");%, ...
load("lab6avg_tz400.mat");

%    "qMap","phiMap","geometry","polarization");
qMap=QGP.qMap;
phiMap=QGP.phiMap;
geometry=QGP.geometryCorr;

polarization=QGP.polarizationCorr;
%% Hard-coded detector mask
% TRUE = invalid pixel BEFORE normalization/difference.
% Use your existing static mask here. If you do not yet have one, start with
% false(size(qMap)) and later replace it with your detector/module/hot-pixel mask.
badMask =  lab6_avg<=0;
 

%% Hard-coded normalization ROI
% TRUE = pixel contributes to the per-frame normalization sum.
% Put your TOP OUTERMOST TILE / chosen monitor-like detector region here.
% Replace these row/column limits with the actual PILATUS pixel range.
normMask = false(size(qMap));
NORM_ROWS = 1:408;       % 
NORM_COLS = 990:1475;       %  
normMask(NORM_ROWS,NORM_COLS) = true;

% Do not normalize from permanently masked pixels.
normMask(badMask) = false;

%% q and phi binning
% Prefer explicit edges from your existing calibration when you have them.
% These examples give 1000 q bins and 360 phi bins over the map ranges.
qGood = qMap(isfinite(qMap) & ~badMask);
qEdges = linspace(min(qGood),max(qGood),1001);
phiGood = phiMap(isfinite(phiMap) & ~badMask);
phiEdges = linspace(min(phiGood),max(phiGood),361);

%% Process once
% Normalization is raw ROI sum -> frame / ROI sum (NormReference=1).
% If geometry/polarization are attenuation factors, CorrectionMode="divide".
% If they are already inverse correction multipliers, use "multiply" instead.
report = pilatus_stream_analysis(files,outDir,window, ...
    Workers=2, ...                    % benchmark 1,2,4 on your storage
    Reader="tiff", ...
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
    SaveQPhi=true, ...                % chunked DeltaI(q,phi)
    SaveDiff2D=false, ...             % avoid hundreds of GB of detector differences
    SaveRawSum=true);

%% Main outputs
S = load(fullfile(outDir,"summary.mat"),"summary");
R = S.summary;

% 1) Per-frame normalization monitor trace
figure;
plot(R.normalization.framePosition,R.normalization.value);
xlabel('Source frame position'); ylabel('Normalization ROI sum'); grid on;
title('Per-frame normalization trace');

% 2) Average normalized 2-D difference over the complete run
figure;
imagesc(R.avgDiff2D); axis image; colorbar;
title('Mean normalized I_{even}-I_{odd}');
colormap(jet(256))
colorbar
% 3) Sum of all normalized source images (NOT differences), useful for mask work
figure;
imagesc(log10(max(R.sumAll2D,eps))); axis image; colorbar;
title('log_{10} sum of all normalized source images');

% 4) Delta I(q,t/window-index): q x output-window matrix
figure;
imagesc(1:size(R.deltaIq,2),R.q,R.deltaIq);
axis xy; colorbar;
xlabel('Output window index'); ylabel('q'); title('\DeltaI(q)');
colormap(jet(256))
caxis([-1e-6 1e-6])
% 5) I(q) of the all-image sum
figure;
plot(R.q,R.sumAllIq); grid on;
xlabel('q'); ylabel('I(q)'); title('I(q) of all-image normalized sum');

% 6) I(q,phi) of the all-image sum
figure;
imagesc(R.phi,R.q,R.sumAllQPhi); axis xy; colorbar;
xlabel('\phi'); ylabel('q'); title('I(q,\phi) of all-image normalized sum');

% Read one saved per-window DeltaI(q,phi), if needed:
[qphi1,iq1,sourceFiles1] = read_pilatus_qphi(outDir,1);  


%%
% I    : 2D intensity image
% qmap : same size as I, radial q value at each pixel

dq = 0.01;   % q-bin width

qedges = min(qmap(:)):dq:max(qmap(:));
q = (qedges(1:end-1) + qedges(2:end))/2;

bin = discretize(qmap(:), qedges);

valid = ~isnan(bin) & isfinite(I(:));

% Mean intensity in each radial q bin
Iq = accumarray(bin(valid), I(valid), ...
    [numel(q), 1], @mean, NaN);

% Number of pixels contributing to each bin
Npix = accumarray(bin(valid), 1, ...
    [numel(q), 1], @sum, 0);

plot(q, Iq);
xlabel('q');
ylabel('I(q)');
