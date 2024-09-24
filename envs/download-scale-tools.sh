#!/bin/sh
srcDir=$(dirname -- "$0")
dir=${1:-"$(readlink -f "${srcDir}/../bin")"}
mkdir -p $dir
echo "Downloading to $dir"
wget https://s3.us-east-2.amazonaws.com/scale.pub/scaleTools/cicd/bcParser/master/240918-g21ce2a3/bc_parser -P $dir && chmod 755 $dir/bc_parser
wget https://s3.us-east-2.amazonaws.com/scale.pub/scaleTools/cicd/scDedup/master/230908-ga1c0044/sc_dedup -P $dir && chmod 755 $dir/sc_dedup
