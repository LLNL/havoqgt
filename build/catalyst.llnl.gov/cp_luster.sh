#!/bin/bash
FILES=/l/ssd/*
OUTFILE="/p/lscratchf/mrdalek/"$1"/"
mkdir -p $OUTFILE
for f in $FILES
do
while :
do 
  fname=$(basename $f)
  path=$(dirname $f)
  echo "command: cp $f $OUTFILE$fname" 
  cp $f $OUTFILE$fname
  iv1=$(md5sum $f | awk '{print $1;}')
  iv2=$(md5sum $OUTFILE$fname | awk '{print $1;}')
  if [ "$iv1" = "$iv2" ]; then
     break
  else
    echo "did not copy file correctly( $iv1 != $iv2 ), trying again.."
  fi;
done
done
