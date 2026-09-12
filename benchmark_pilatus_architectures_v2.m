
============================================================
PILATUS ARCHITECTURE BENCHMARK v2
============================================================
Files found : 14520
Input size  : 0.144 TB
Frame IDs   : 1 ... 14520
Window      : 3 (6 TIFFs/output)
CPU workers : 2
Detector     : 1679 x 1475
Norm ROI     : 188900 good pixels
Timed/test   : 180 TIFFs (~1.78 GB)
Warm-up/test : 12 TIFFs
GPU          : Quadro RTX 4000, 8.6 GB

------------------------------------------------------------
TEST 1/3: CPU 1000q x 360phi + qphi(t)
------------------------------------------------------------
Warm-up...
Timing 180 TIFFs...
Elapsed           : 21.28 s
Throughput        : 8.46 TIFF/s
Input throughput  : 83.8 MB/s
Benchmark output  : 97.2 MB
Projected full run: 28.6 min

------------------------------------------------------------
TEST 2/3: CPU  500q x  90phi + qphi(t)
------------------------------------------------------------
Warm-up...
Timing 180 TIFFs...
Elapsed           : 20.15 s
Throughput        : 8.93 TIFF/s
Input throughput  : 88.5 MB/s
Benchmark output  : 55.5 MB
Projected full run: 27.1 min

------------------------------------------------------------
TEST 3/3: GPU  500q x  90phi, radial(t) only
------------------------------------------------------------
Warm-up...
Timing 180 TIFFs...
Elapsed           : 30.55 s
Throughput        : 5.89 TIFF/s
Input throughput  : 58.4 MB/s
Benchmark output  : 50.1 MB
Projected full run: 41.1 min


============================================================
FINAL RESULT
============================================================
    Test               Configuration                 Nq     Nphi     GPU     SaveQPhiTime    Seconds    TIFF_per_s    Input_MB_per_s    Output_MB    Projected_full_min
    ____    ____________________________________    ____    ____    _____    ____________    _______    __________    ______________    _________    __________________

     2      "CPU  500q x  90phi + qphi(t)"           500     90     false       true         20.151       8.9324            88.5         55.541            27.092      
     1      "CPU 1000q x 360phi + qphi(t)"          1000    360     false       true         21.283       8.4576          83.796         97.191            28.613      
     3      "GPU  500q x  90phi, radial(t) only"     500     90     true        false        30.545       5.8929          58.385         50.141            41.066      


*** FASTEST: TEST 2 - CPU  500q x  90phi + qphi(t) ***
8.93 TIFF/s, projected 27.1 min for all 14520 TIFFs.
Fastest is 5.3% quicker than second place.

How to interpret the tests:
  TEST 1 vs TEST 2: effect of reducing q/phi resolution.
  TEST 2 vs TEST 3: combined gain from GPU processing and not storing
                     the full time-resolved q/phi cube.
  If all Input_MB_per_s values are nearly the same, X: I/O is probably
  the dominant bottleneck and further compute optimization will help less.

Results saved to:
C:\Users\11idbuser\Downloads\Natan\pilatus_fast_diff_matlab\analysis\benchmark_fast\benchmark_result.mat
============================================================
>> 
