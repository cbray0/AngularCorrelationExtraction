# AngularCorrelationExtraction
A flexible GRSISort macro for extracting angular correlations from data.

This repo contains a variety of useful functions for automating simulation execution.

### Dependencies:

GRSISort (Commit 6423d91 or newer).

### Usage instructions:

After you have converted the g4out.root file using NTuple, launch GRSISort. Load the commands using `.L angCorr.cpp+` (or the full path if it is in another directory).

Then you can use the angCorr() command to calculate angular correlations. See the [full documentation](https://cbray0.github.io/AngularCorrelationExtraction/html/angCorr_8cpp.html) for specific syntax, but in general you can either provide a set of energies to find all angular correlations between, or tell the command to find the n largest peaks and calculate angular correlations between those.

Based on initial testing, the macro is capable of extracting all angular correlations between 80 peaks (3160 total angular correlations) in under an hour (single threaded test performed on an Intel(R) Xeon(R) CPU E5-2680 v3 @ 2.50GHz).

### TODO:

* Further testing and QA

* Determine if multithreading is worthwhile and possible without massive memory usage

Go to [Github pages](https://cbray0.github.io/AngularCorrelationExtraction/html/angCorr_8cpp.html) for full documentation.
