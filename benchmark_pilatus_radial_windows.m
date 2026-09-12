
============================================================
PILATUS RADIAL-ONLY WINDOW BENCHMARK
============================================================
Files found : 14520
Input size  : 0.144 TB
Frame IDs   : 1 ... 14520
Workers     : 2
q bins      : 1000
Windows     : [1 3 10]
Detector     : 1679 x 1475
Norm ROI     : 188900 good pixels
q pixels     : 2272563
Timed/test   : 180 TIFFs (~1.78 GB)

------------------------------------------------------------
WINDOW = 1
------------------------------------------------------------
Benchmark outputs : 90 DeltaI(q) curves
Full-run outputs  : 7260 DeltaI(q) curves
Elapsed           : 17.78 s
Throughput        : 10.12 TIFF/s
Input throughput  : 100.3 MB/s
Output file       : 20.2 MB
Projected full run: 23.9 min

------------------------------------------------------------
WINDOW = 3
------------------------------------------------------------
Benchmark outputs : 30 DeltaI(q) curves
Full-run outputs  : 2420 DeltaI(q) curves
Elapsed           : 15.22 s
Throughput        : 11.83 TIFF/s
Input throughput  : 117.2 MB/s
Output file       : 19.9 MB
Projected full run: 20.5 min

------------------------------------------------------------
WINDOW = 10
------------------------------------------------------------
Benchmark outputs : 9 DeltaI(q) curves
Full-run outputs  : 726 DeltaI(q) curves
Elapsed           : 15.40 s
Throughput        : 11.69 TIFF/s
Input throughput  : 115.8 MB/s
Output file       : 19.9 MB
Projected full run: 20.7 min


============================================================
FINAL WINDOW BENCHMARK
============================================================
    Window    FullRun_DeltaIq_Curves    Seconds    TIFF_per_s    Input_MB_per_s    Benchmark_Output_MB    Projected_full_min
    ______    ______________________    _______    __________    ______________    ___________________    __________________

       3               2420             15.219       11.827          117.18              19.941                 20.462      
      10                726               15.4       11.688           115.8              19.857                 20.705      
       1               7260             17.779       10.124          100.31              20.181                 23.903      


Interpretation:
Fastest measured window: W=3
Speed advantage over second place: 1.2%
The timing difference is <5%. Treat the windows as computationally
equivalent and choose W from the science/SNR/time-resolution tradeoff.

Remember:
  The final average 2-D difference is essentially independent of W
  (apart from an incomplete tail / masking details).
  W mainly changes how many pair-differences are averaged into each
  DeltaI(q) time point:

    W=1   -> 7260 DeltaI(q) curves over the full 14520-file run
    W=3   -> 2420 DeltaI(q) curves over the full 14520-file run
    W=10  -> 726 DeltaI(q) curves over the full 14520-file run

Results saved to:
C:\Users\11idbuser\Downloads\Natan\pilatus_fast_diff_matlab\analysis\benchmark_radial_windows\benchmark_radial_windows_result.mat
============================================================
>> 
