# Debuncher
code for the estimation of the efficiency of a debuncher and production of root file from FASTER acquisition file

usage:

./debuncher <filename> <traping time [ms]> <extraction time [ms]>

the code estimates the number of ions registrated by the MCP

the efficientcy test:
the ions are injected in the debuncher trapped for the certain time and then ejected in a short time registrated by MCP and the FASTER acquisition. The code will estimate the number ions per bunch and the average number of ions perbunch

The debunching>

in this case for the trapping time, we should put injection time ~10us and for the extraction time we should put the real trapping time as in this case the ions are extracted during the trapping 

Prerequisites:
ROOT data analysis framework (https://root.cern.ch/)
FASTER analysis frame work (http://faster.in2p3.fr/index.php/download/category/2-software#)
