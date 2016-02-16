#!/bin/sh

if [ $1 -lt 100 ];
then 

root -l -q -b MitZHAnalysis/macros/74x/zhAnalysis.C+\($1+2,0\);
root -l -q -b MitZHAnalysis/macros/74x/zhAnalysis.C+\($1+2,1\);
root -l -q -b MitZHAnalysis/macros/74x/zhAnalysis.C+\($1+1,0\);
root -l -q -b MitZHAnalysis/macros/74x/zhAnalysis.C+\($1+1,1\);
root -l -q -b MitZHAnalysis/macros/74x/zhAnalysis.C+\($1,0\);
root -l -q -b MitZHAnalysis/macros/74x/zhAnalysis.C+\($1,1\);

else

root -l -q -b MitZHAnalysis/macros/74x/zhAnalysis.C+\($1,0\);
root -l -q -b MitZHAnalysis/macros/74x/zhAnalysis.C+\($1,1\);

fi
