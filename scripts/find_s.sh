#!/bin/bash
# Find the redshift of a given halo catalog

path=$1
loc_path=`pwd`
tmp=$loc_path'/tmp/sahf.tmp'
tmp_out=$loc_path'/tmp/out_s.tmp'
cd $path

if [ -f "$tmp_out" ]
then
	rm $tmp_out
fi

for f_ahf in `ls *AHF_halos`
do
	echo $f_ahf > $tmp
	this_s=`grep -o _[0-9][0-9][0-9] $tmp | sed 's/_//g ' `
	echo $this_s >> $tmp_out
done

rm $tmp
echo $tmp_out
