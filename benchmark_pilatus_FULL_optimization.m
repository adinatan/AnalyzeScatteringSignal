
============================================================
PILATUS FULL OPTIMIZATION BENCHMARK
============================================================
Files found : 14520
Input size  : 0.144 TB
Frame IDs   : 1 ... 14520
Readers     : tiff, imread
Workers     : [1 2 3 4 6 8]
Windows     : [1 3 10]
q bins      : [500 1000 1500]
GPU modes   : frame, batch
Detector     : 1679 x 1475
Norm ROI     : 188900 good pixels
GPU          : Quadro RTX 4000, 8.6 GB
Quick block  : 120 TIFFs (~1.19 GB)
Final block  : 600 TIFFs (~5.94 GB)

============================================================
STAGE 0: RADIAL REDUCTION KERNEL
============================================================
accumarray: 1.0833 s for 20 reductions
sparse    : 0.6005 s for 20 reductions
Selected CPU radial kernel: SPARSE

============================================================
STAGE 1: READER x PROCESS-WORKER SEARCH
Representative science setting: W=3, Nq=1000
============================================================

CPU reader=TIFF, workers=1
Parallel pool using the 'Processes' profile is shutting down.
10.79 s | 11.12 TIFF/s | 110.2 MB/s | projected 21.8 min

CPU reader=IMREAD, workers=1
17.32 s | 6.93 TIFF/s | 68.6 MB/s | projected 34.9 min

CPU reader=TIFF, workers=2
Starting parallel pool (parpool) using the 'Processes' profile ...
Connected to parallel pool with 2 workers.
10.62 s | 11.30 TIFF/s | 112.0 MB/s | projected 21.4 min

CPU reader=IMREAD, workers=2
10.71 s | 11.21 TIFF/s | 111.1 MB/s | projected 21.6 min

CPU reader=TIFF, workers=3
Parallel pool using the 'Processes' profile is shutting down.
Starting parallel pool (parpool) using the 'Processes' profile ...
Connected to parallel pool with 3 workers.
8.15 s | 14.73 TIFF/s | 146.0 MB/s | projected 16.4 min

CPU reader=IMREAD, workers=3
8.27 s | 14.52 TIFF/s | 143.8 MB/s | projected 16.7 min

CPU reader=TIFF, workers=4
Parallel pool using the 'Processes' profile is shutting down.
Starting parallel pool (parpool) using the 'Processes' profile ...
Connected to parallel pool with 4 workers.
6.83 s | 17.57 TIFF/s | 174.1 MB/s | projected 13.8 min

CPU reader=IMREAD, workers=4
6.91 s | 17.37 TIFF/s | 172.1 MB/s | projected 13.9 min

CPU reader=TIFF, workers=6
Parallel pool using the 'Processes' profile is shutting down.
Starting parallel pool (parpool) using the 'Processes' profile ...
SKIPPED: Minimum number of workers requested (6) must be less than or equal to the "Processes" NumWorkers value of 4 workers. Request fewer workers or increase the NumWorkers value by modifying the NumWorkers property for the "Processes" profile in the Cluster Profile Manager (up to a maximum of 512 for a Processes profile).

CPU reader=IMREAD, workers=6
Starting parallel pool (parpool) using the 'Processes' profile ...
SKIPPED: Minimum number of workers requested (6) must be less than or equal to the "Processes" NumWorkers value of 4 workers. Request fewer workers or increase the NumWorkers value by modifying the NumWorkers property for the "Processes" profile in the Cluster Profile Manager (up to a maximum of 512 for a Processes profile).

CPU reader=TIFF, workers=8
Starting parallel pool (parpool) using the 'Processes' profile ...
SKIPPED: Minimum number of workers requested (8) must be less than or equal to the "Processes" NumWorkers value of 4 workers. Request fewer workers or increase the NumWorkers value by modifying the NumWorkers property for the "Processes" profile in the Cluster Profile Manager (up to a maximum of 512 for a Processes profile).

CPU reader=IMREAD, workers=8
Starting parallel pool (parpool) using the 'Processes' profile ...
SKIPPED: Minimum number of workers requested (8) must be less than or equal to the "Processes" NumWorkers value of 4 workers. Request fewer workers or increase the NumWorkers value by modifying the NumWorkers property for the "Processes" profile in the Cluster Profile Manager (up to a maximum of 512 for a Processes profile).

STAGE 1 RANKING
     Reader     Workers    Seconds    TIFF_per_s    Input_MB_per_s    Projected_full_min
    ________    _______    _______    __________    ______________    __________________

    "tiff"         4       6.8284       17.574          174.12              13.771      
    "imread"       4       6.9088       17.369          172.09              13.933      
    "tiff"         3       8.1461       14.731          145.95              16.428      
    "imread"       3       8.2666       14.516          143.82              16.671      
    "tiff"         2       10.618       11.301          111.97              21.414      
    "imread"       2       10.705       11.209          111.06              21.589      
    "tiff"         1       10.788       11.124          110.21              21.755      
    "imread"       1       17.325       6.9266          68.627              34.938      

Carrying the top 3 I/O configurations into Stage 2.

============================================================
STAGE 2: WINDOW x q-RESOLUTION SEARCH
============================================================
Starting parallel pool (parpool) using the 'Processes' profile ...
Connected to parallel pool with 4 workers.

CPU TIFF, workers=4, W=1, Nq=500
7.65 s | 15.68 TIFF/s | 155.4 MB/s | projected 15.4 min

CPU TIFF, workers=4, W=1, Nq=1000
7.53 s | 15.93 TIFF/s | 157.8 MB/s | projected 15.2 min

CPU TIFF, workers=4, W=1, Nq=1500
7.82 s | 15.34 TIFF/s | 152.0 MB/s | projected 15.8 min

CPU TIFF, workers=4, W=3, Nq=500
6.77 s | 17.72 TIFF/s | 175.6 MB/s | projected 13.7 min

CPU TIFF, workers=4, W=3, Nq=1000
6.56 s | 18.30 TIFF/s | 181.3 MB/s | projected 13.2 min

CPU TIFF, workers=4, W=3, Nq=1500
6.66 s | 18.02 TIFF/s | 178.6 MB/s | projected 13.4 min

CPU TIFF, workers=4, W=10, Nq=500
7.44 s | 16.13 TIFF/s | 159.9 MB/s | projected 15.0 min

CPU TIFF, workers=4, W=10, Nq=1000
7.48 s | 16.04 TIFF/s | 158.9 MB/s | projected 15.1 min

CPU TIFF, workers=4, W=10, Nq=1500
7.42 s | 16.18 TIFF/s | 160.3 MB/s | projected 15.0 min

CPU IMREAD, workers=4, W=1, Nq=500
8.01 s | 14.99 TIFF/s | 148.5 MB/s | projected 16.1 min

CPU IMREAD, workers=4, W=1, Nq=1000
7.74 s | 15.50 TIFF/s | 153.6 MB/s | projected 15.6 min

CPU IMREAD, workers=4, W=1, Nq=1500
7.85 s | 15.28 TIFF/s | 151.4 MB/s | projected 15.8 min

CPU IMREAD, workers=4, W=3, Nq=500
6.94 s | 17.28 TIFF/s | 171.2 MB/s | projected 14.0 min

CPU IMREAD, workers=4, W=3, Nq=1000
6.65 s | 18.06 TIFF/s | 178.9 MB/s | projected 13.4 min

CPU IMREAD, workers=4, W=3, Nq=1500
6.68 s | 17.95 TIFF/s | 177.9 MB/s | projected 13.5 min

CPU IMREAD, workers=4, W=10, Nq=500
7.64 s | 15.71 TIFF/s | 155.6 MB/s | projected 15.4 min

CPU IMREAD, workers=4, W=10, Nq=1000
7.40 s | 16.22 TIFF/s | 160.7 MB/s | projected 14.9 min

CPU IMREAD, workers=4, W=10, Nq=1500
7.48 s | 16.04 TIFF/s | 158.9 MB/s | projected 15.1 min
Parallel pool using the 'Processes' profile is shutting down.
Starting parallel pool (parpool) using the 'Processes' profile ...
Connected to parallel pool with 3 workers.

CPU TIFF, workers=3, W=1, Nq=500
9.05 s | 13.26 TIFF/s | 131.4 MB/s | projected 18.2 min

CPU TIFF, workers=3, W=1, Nq=1000
8.95 s | 13.41 TIFF/s | 132.9 MB/s | projected 18.0 min

CPU TIFF, workers=3, W=1, Nq=1500
8.99 s | 13.35 TIFF/s | 132.3 MB/s | projected 18.1 min

CPU TIFF, workers=3, W=3, Nq=500
8.00 s | 15.01 TIFF/s | 148.7 MB/s | projected 16.1 min

CPU TIFF, workers=3, W=3, Nq=1000
8.20 s | 14.64 TIFF/s | 145.1 MB/s | projected 16.5 min

CPU TIFF, workers=3, W=3, Nq=1500
7.99 s | 15.03 TIFF/s | 148.9 MB/s | projected 16.1 min

CPU TIFF, workers=3, W=10, Nq=500
7.36 s | 16.31 TIFF/s | 161.6 MB/s | projected 14.8 min

CPU TIFF, workers=3, W=10, Nq=1000
7.35 s | 16.33 TIFF/s | 161.8 MB/s | projected 14.8 min

CPU TIFF, workers=3, W=10, Nq=1500
7.29 s | 16.46 TIFF/s | 163.1 MB/s | projected 14.7 min

STAGE 2 CPU RANKING
     Reader     Workers    Window     Nq     Seconds    TIFF_per_s    Input_MB_per_s    Projected_full_min
    ________    _______    ______    ____    _______    __________    ______________    __________________

    "tiff"         4          3      1000    6.5591       18.295          181.27              13.227      
    "imread"       4          3      1000    6.6452       18.058          178.91              13.401      
    "tiff"         4          3      1500    6.6579       18.024          178.58              13.427      
    "imread"       4          3      1500    6.6835       17.955          177.89              13.478      
    "tiff"         4          3       500    6.7714       17.722          175.58              13.656      
    "imread"       4          3       500    6.9437       17.282          171.22              14.003      
    "tiff"         3         10      1500    7.2903        16.46          163.09              14.702      
    "tiff"         3         10      1000    7.3492       16.328          161.78              14.821      
    "tiff"         3         10       500    7.3559       16.313          161.63              14.834      
    "imread"       4         10      1000    7.3976       16.221          160.72              14.919      
    "tiff"         4         10      1500    7.4175       16.178          160.29              14.959      
    "tiff"         4         10       500    7.4376       16.134          159.85              14.999      
    "tiff"         4         10      1000    7.4825       16.037          158.89               15.09      
    "imread"       4         10      1500    7.4828       16.037          158.89               15.09      
    "tiff"         4          1      1000    7.5333       15.929          157.82              15.192      
    "imread"       4         10       500    7.6392       15.708          155.64              15.406      
    "tiff"         4          1       500    7.6522       15.682          155.37              15.432      
    "imread"       4          1      1000    7.7409       15.502          153.59              15.611      
    "tiff"         4          1      1500    7.8213       15.343          152.01              15.773      
    "imread"       4          1      1500    7.8545       15.278          151.37               15.84      
    "tiff"         3          3      1500    7.9855       15.027          148.89              16.104      
    "tiff"         3          3       500    7.9969       15.006          148.67              16.127      
    "imread"       4          1       500    8.0061       14.989           148.5              16.146      
    "tiff"         3          3      1000    8.1953       14.643          145.08              16.527      
    "tiff"         3          1      1000    8.9461       13.414           132.9              18.041      
    "tiff"         3          1      1500    8.9873       13.352          132.29              18.124      
    "tiff"         3          1       500    9.0483       13.262           131.4              18.247      


============================================================
STAGE 3: GPU SEARCH
============================================================
Parallel pool using the 'Processes' profile is shutting down.

GPU reader=TIFF, transfer=FRAME, W=1, Nq=500
21.26 s | 5.64 TIFF/s | 55.9 MB/s | projected 42.9 min

GPU reader=TIFF, transfer=FRAME, W=1, Nq=1000
19.49 s | 6.16 TIFF/s | 61.0 MB/s | projected 39.3 min

GPU reader=TIFF, transfer=FRAME, W=1, Nq=1500
19.99 s | 6.00 TIFF/s | 59.5 MB/s | projected 40.3 min

GPU reader=TIFF, transfer=FRAME, W=3, Nq=500
17.94 s | 6.69 TIFF/s | 66.3 MB/s | projected 36.2 min

GPU reader=TIFF, transfer=FRAME, W=3, Nq=1000
17.91 s | 6.70 TIFF/s | 66.4 MB/s | projected 36.1 min

GPU reader=TIFF, transfer=FRAME, W=3, Nq=1500
17.95 s | 6.69 TIFF/s | 66.2 MB/s | projected 36.2 min

GPU reader=TIFF, transfer=FRAME, W=10, Nq=500
17.36 s | 6.91 TIFF/s | 68.5 MB/s | projected 35.0 min

GPU reader=TIFF, transfer=FRAME, W=10, Nq=1000
17.28 s | 6.94 TIFF/s | 68.8 MB/s | projected 34.9 min

GPU reader=TIFF, transfer=FRAME, W=10, Nq=1500
17.16 s | 6.99 TIFF/s | 69.3 MB/s | projected 34.6 min

GPU reader=TIFF, transfer=BATCH, W=1, Nq=500
20.32 s | 5.90 TIFF/s | 58.5 MB/s | projected 41.0 min

GPU reader=TIFF, transfer=BATCH, W=1, Nq=1000
20.15 s | 5.96 TIFF/s | 59.0 MB/s | projected 40.6 min

GPU reader=TIFF, transfer=BATCH, W=1, Nq=1500
20.31 s | 5.91 TIFF/s | 58.6 MB/s | projected 41.0 min

GPU reader=TIFF, transfer=BATCH, W=3, Nq=500
19.19 s | 6.25 TIFF/s | 62.0 MB/s | projected 38.7 min

GPU reader=TIFF, transfer=BATCH, W=3, Nq=1000
18.84 s | 6.37 TIFF/s | 63.1 MB/s | projected 38.0 min

GPU reader=TIFF, transfer=BATCH, W=3, Nq=1500
18.56 s | 6.47 TIFF/s | 64.1 MB/s | projected 37.4 min

GPU reader=TIFF, transfer=BATCH, W=10, Nq=500
18.78 s | 6.39 TIFF/s | 63.3 MB/s | projected 37.9 min

GPU reader=TIFF, transfer=BATCH, W=10, Nq=1000
20.48 s | 5.86 TIFF/s | 58.1 MB/s | projected 41.3 min

GPU reader=TIFF, transfer=BATCH, W=10, Nq=1500
18.57 s | 6.46 TIFF/s | 64.0 MB/s | projected 37.5 min

GPU reader=IMREAD, transfer=FRAME, W=1, Nq=500
21.10 s | 5.69 TIFF/s | 56.4 MB/s | projected 42.5 min

GPU reader=IMREAD, transfer=FRAME, W=1, Nq=1000
20.50 s | 5.85 TIFF/s | 58.0 MB/s | projected 41.3 min

GPU reader=IMREAD, transfer=FRAME, W=1, Nq=1500
20.71 s | 5.79 TIFF/s | 57.4 MB/s | projected 41.8 min

GPU reader=IMREAD, transfer=FRAME, W=3, Nq=500
18.53 s | 6.48 TIFF/s | 64.2 MB/s | projected 37.4 min

GPU reader=IMREAD, transfer=FRAME, W=3, Nq=1000
18.64 s | 6.44 TIFF/s | 63.8 MB/s | projected 37.6 min

GPU reader=IMREAD, transfer=FRAME, W=3, Nq=1500
18.47 s | 6.50 TIFF/s | 64.4 MB/s | projected 37.2 min

GPU reader=IMREAD, transfer=FRAME, W=10, Nq=500
17.68 s | 6.79 TIFF/s | 67.2 MB/s | projected 35.7 min

GPU reader=IMREAD, transfer=FRAME, W=10, Nq=1000
17.75 s | 6.76 TIFF/s | 67.0 MB/s | projected 35.8 min

GPU reader=IMREAD, transfer=FRAME, W=10, Nq=1500
17.69 s | 6.78 TIFF/s | 67.2 MB/s | projected 35.7 min

GPU reader=IMREAD, transfer=BATCH, W=1, Nq=500
21.52 s | 5.58 TIFF/s | 55.2 MB/s | projected 43.4 min

GPU reader=IMREAD, transfer=BATCH, W=1, Nq=1000
21.18 s | 5.66 TIFF/s | 56.1 MB/s | projected 42.7 min

GPU reader=IMREAD, transfer=BATCH, W=1, Nq=1500
21.55 s | 5.57 TIFF/s | 55.2 MB/s | projected 43.5 min

GPU reader=IMREAD, transfer=BATCH, W=3, Nq=500
20.17 s | 5.95 TIFF/s | 58.9 MB/s | projected 40.7 min

GPU reader=IMREAD, transfer=BATCH, W=3, Nq=1000
18.73 s | 6.41 TIFF/s | 63.5 MB/s | projected 37.8 min

GPU reader=IMREAD, transfer=BATCH, W=3, Nq=1500
18.66 s | 6.43 TIFF/s | 63.7 MB/s | projected 37.6 min

GPU reader=IMREAD, transfer=BATCH, W=10, Nq=500
18.73 s | 6.41 TIFF/s | 63.5 MB/s | projected 37.8 min

GPU reader=IMREAD, transfer=BATCH, W=10, Nq=1000
18.65 s | 6.43 TIFF/s | 63.7 MB/s | projected 37.6 min

GPU reader=IMREAD, transfer=BATCH, W=10, Nq=1500
18.89 s | 6.35 TIFF/s | 62.9 MB/s | projected 38.1 min

STAGE 3 GPU RANKING
     Reader     GPUTransfer    Window     Nq     Seconds    TIFF_per_s    Input_MB_per_s    Projected_full_min
    ________    ___________    ______    ____    _______    __________    ______________    __________________

    "tiff"        "frame"        10      1500     17.16       6.9929          69.284              34.606      
    "tiff"        "frame"        10      1000    17.285       6.9426          68.785              34.857      
    "tiff"        "frame"        10       500    17.355       6.9143          68.505                  35      
    "imread"      "frame"        10       500    17.682       6.7865          67.239              35.659      
    "imread"      "frame"        10      1500    17.693       6.7823          67.197              35.681      
    "imread"      "frame"        10      1000    17.747       6.7617          66.993               35.79      
    "tiff"        "frame"         3      1000    17.907       6.7012          66.394              36.113      
    "tiff"        "frame"         3       500    17.937       6.6902          66.285              36.172      
    "tiff"        "frame"         3      1500     17.95       6.6851          66.234                36.2      
    "imread"      "frame"         3      1500    18.467        6.498           64.38              37.242      
    "imread"      "frame"         3       500    18.532       6.4752          64.155              37.373      
    "tiff"        "batch"         3      1500    18.558       6.4661          64.065              37.426      
    "tiff"        "batch"        10      1500    18.574       6.4607          64.011              37.457      
    "imread"      "frame"         3      1000    18.645       6.4361          63.768                37.6      
    "imread"      "batch"        10      1000    18.652       6.4337          63.743              37.615      
    "imread"      "batch"         3      1500    18.661       6.4305          63.712              37.633      
    "imread"      "batch"        10       500    18.729       6.4071           63.48               37.77      
    "imread"      "batch"         3      1000     18.73       6.4068          63.477              37.772      
    "tiff"        "batch"        10       500    18.778       6.3904          63.315              37.869      
    "tiff"        "batch"         3      1000    18.836       6.3708           63.12              37.986      
    "imread"      "batch"        10      1500    18.889       6.3529          62.944              38.093      
    "tiff"        "batch"         3       500    19.187       6.2543          61.967              38.693      
    "tiff"        "frame"         1      1000    19.493       6.1559          60.991              39.312      
    "tiff"        "frame"         1      1500    19.986       6.0041          59.487              40.306      
    "tiff"        "batch"         1      1000    20.147       5.9563          59.013              40.629      
    "imread"      "batch"         3       500    20.175        5.948          58.932              40.686      
    "tiff"        "batch"         1      1500    20.306       5.9095           58.55              40.951      
    "tiff"        "batch"         1       500    20.324       5.9043          58.498              40.987      
    "tiff"        "batch"        10      1000     20.48       5.8594          58.053              41.301      
    "imread"      "frame"         1      1000    20.497       5.8544          58.004              41.336      
    "imread"      "frame"         1      1500    20.714       5.7933          57.399              41.772      
    "imread"      "frame"         1       500    21.097       5.6881          56.357              42.545      
    "imread"      "batch"         1      1000    21.184       5.6647          56.124              42.721      
    "tiff"        "frame"         1       500    21.259       5.6447          55.926              42.872      
    "imread"      "batch"         1       500    21.524       5.5752          55.238              43.406      
    "imread"      "batch"         1      1500    21.554       5.5675          55.162              43.466      


============================================================
STAGE 4: LARGE-BLOCK FINAL VALIDATION
============================================================

FINAL CPU: TIFF workers=4 W=3 Nq=1000
Starting parallel pool (parpool) using the 'Processes' profile ...
Connected to parallel pool with 4 workers.
30.44 s | 19.71 TIFF/s | 195.3 MB/s | projected 12.3 min

FINAL CPU: IMREAD workers=4 W=3 Nq=1000
31.36 s | 19.13 TIFF/s | 189.6 MB/s | projected 12.6 min
Parallel pool using the 'Processes' profile is shutting down.

FINAL GPU: TIFF FRAME W=10 Nq=1500
85.98 s | 6.98 TIFF/s | 69.1 MB/s | projected 34.7 min

FINAL GPU: TIFF FRAME W=10 Nq=1000
87.44 s | 6.86 TIFF/s | 68.0 MB/s | projected 35.3 min

FINAL VALIDATION RANKING
    Platform     Reader     Workers    Window     Nq       RadialMethod      GPUTransfer    Seconds    TIFF_per_s    Input_MB_per_s    Projected_full_min
    ________    ________    _______    ______    ____    ________________    ___________    _______    __________    ______________    __________________

     "CPU"      "tiff"         4          3      1000    "sparse"              ""           30.443       19.709          195.27              12.279      
     "CPU"      "imread"       4          3      1000    "sparse"              ""            31.36       19.132          189.56              12.649      
     "GPU"      "tiff"         0         10      1500    "gpu_accumarray"      "frame"      85.978       6.9785          69.142              34.678      
     "GPU"      "tiff"         0         10      1000    "gpu_accumarray"      "frame"       87.44       6.8619          67.986              35.267      


============================================================
WINNING PRODUCTION CONFIGURATION
============================================================
    Platform    Reader    Workers    Window     Nq     RadialMethod    GPUTransfer    TIFF_per_s    Input_MB_per_s    Projected_full_min
    ________    ______    _______    ______    ____    ____________    ___________    __________    ______________    __________________

     "CPU"      "tiff"       4         3       1000      "sparse"          ""           19.709          195.27              12.279      


============================================================
STAGE 5: X: DIRECT vs LOCAL C: STAGING
============================================================
Starting parallel pool (parpool) using the 'Processes' profile ...
Connected to parallel pool with 4 workers.
Copying 300 TIFFs to local staging folder...
    DirectX_TIFF_per_s    LocalC_TIFF_per_s    Copy_seconds    Copy_GB_per_s    Projected_DirectX_min    Projected_Copy_min    Projected_LocalProcess_min    Projected_CopyPlusProcess_min
    __________________    _________________    ____________    _____________    _____________________    __________________    __________________________    _____________________________

          18.788               38.799             8.4853          0.35029              12.881                  6.8448                    6.2373                         13.082            

DIRECT X: PROCESSING WINS when copy time is included.
There is no reason to stage the entire dataset locally for speed.

============================================================
BENCHMARK COMPLETE
============================================================
Full results saved to:
C:\Users\11idbuser\Downloads\Natan\pilatus_fast_diff_matlab\analysis\full_optimization_benchmark\FULL_BENCHMARK_RESULTS.mat

For production use the FINAL validation winner, not merely the
fastest 120-file quick candidate.
>> 
