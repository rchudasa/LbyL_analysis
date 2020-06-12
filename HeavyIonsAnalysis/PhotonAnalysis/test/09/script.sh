#!/bin/sh
numOfLines=`wc -l outlist | awk '{print $1}'`
#echo $numOfLines
bulkSize=39
numOfFiles=`echo $numOfLines/$bulkSize | bc `
echo "Total number of files to be generated : " $numOfFiles
numOfFiles=`expr $numOfFiles - 1`
for i in $(seq 0 1 $numOfFiles)
do
   linefrom=`expr $i \\* $bulkSize`
   linefrom=`expr $linefrom + 1`
   ito=`expr $i + 1`
   lineto=`expr $ito \\* $bulkSize`
   #echo $linefrom"  :  "$lineto
   len=${#ito}
   #echo $len
   if [ "$len" -eq 1 ]
   then
      ito="0"$ito
   fi
   outfile="filelist_"$ito
   #echo $outfile
   #cat outlist | awk -v awklinefrom="$linefrom" -v awklineto="$lineto" 'NR==awklinefrom,NR==awklineto {print $5}' > $outfile
   cat outlist | awk -v awklinefrom="$linefrom" -v awklineto="$lineto" 'NR==awklinefrom,NR==awklineto' > $outfile
   echo "Generated file : "$outfile
done
