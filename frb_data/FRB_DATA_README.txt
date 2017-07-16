Evan Keane
15/12/2014

This directory contains data from 7 archival FRBs 

FRB010125
FRB010621
FRB010724
FRB110220
FRB110626
FRB110703
FRB120127

For these FRBs there are the following:

scamp_files   	For the first 3 FRBs only there are SCAMP dat files (not presto dat files!) for all 13 beams for the pointing, as produced by sc_td when reading these data from the original DLT magnetic tapes the survey was recorded on, or on the dd'd images from those same tapes.

fil_files	Full sigproc filterbank files for all 13 beams for the pointing

ar_files	2-second psrchive PSRFITS files centred around the FRB time, with full time and frequency resolution.

ascii_files	Ascii headers (the output from psredit run on the .ar files). For the first 3 FRBs there are 2-second frequency-time-amplitude data cubes at raw resolution (0.125, 0.250 and 1.000 ms, respectively), and 4-ms resolution. Ascii versions of the BPSR FRBs are not included, but a simple csh script entitled ascii_cube.csh, which is included, can be run to create them.


