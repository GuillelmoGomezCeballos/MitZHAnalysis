#!/bin/sh

root -l -q -b MitZHAnalysis/macros/74x/zhAnalysis.C+\($1,0\);
root -l -q -b MitZHAnalysis/macros/74x/zhAnalysis.C+\($1,1\);
root -l -q -b MitZHAnalysis/macros/74x/zhAnalysis.C+\($1+1,0\);
root -l -q -b MitZHAnalysis/macros/74x/zhAnalysis.C+\($1+1,1\);
