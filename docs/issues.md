# Known issues

## Warning about damaged file contents upon opening of the Excel report

If opening with Excel, a warning window to recover damaged contents might pop up.
No issue with Libre Office. Related issues found across GitHub:

* https://github.com/jmcnamara/XlsxWriter/issues/1109
* https://github.com/jmcnamara/XlsxWriter/issues/739

## Random failures of CNVkit

We have experienced several failures in the pipeline when CNVkit is run and we have not found any clear explanation yet for why this is happening.
Usually, merely restarting the workflow and let it run with the same setup is the solution to a successful execution.

## PureCN is not working

This is a silent error: it looks like the related jobs have completed successfully, however the output is empty.

## Bamsnap runs out of time

Bamsnap is tricky to get to work and computationally expensive. 
If the samples contain hundreds or thousands of variants that must be screenshot in IGV,
the time resources for Bamsnap must be increased.
If the samples were sequenced at very high coverage, the time also need to be increased.

Currently, the time is tuned to fit the needs of samples sequenced on NovaSeqX at ca. 1200x 
and we expect at most a few dozen of variants to screenshoot in each sample.

We use computing 8-core nodes of at least 128 GB (8 x 16 GB and 2.4 GHz). 