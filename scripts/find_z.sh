#!/bin/bash
# Find the redshift of a given halo catalog

path=$1
loc_path=`pwd`
tmp=$loc_path'/tmp/ahf.tmp'
tmp_out=$loc_path'/tmp/out_z.tmp'
cd $path

if [ -f "$tmp_out" ]
then
	rm $tmp_out
fi

for f_ahf in `ls *AHF_halos`
do
	echo $f_ahf > $tmp
	this_z=`grep -o z[0-9].[0-9][0-9][0-9] $tmp | sed 's/z//g ' `
	
	if [ "$this_z" == "" ]
	then
		this_z=`grep -o z[0-9][0-9].[0-9][0-9][0-9] $tmp | sed 's/z//g ' `
	fi

	echo $this_z >> $tmp_out
done

echo $tmp_out
rm $tmp
