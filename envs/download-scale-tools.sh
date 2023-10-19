#!/bin/sh
srcDir=$(dirname -- "$0")
dir=${1:-"$(readlink -f "${srcDir}/../bin")"}
mkdir -p $dir
echo "Downloading to $dir"
curl https://s3.us-east-2.amazonaws.com/scale.pub/scaleTools/cicd/bcParser/master/230220-gf9e5e90/bc_parser -o $dir/bc_parser && chmod 755 $dir/bc_parser
curl https://s3.us-east-2.amazonaws.com/scale.pub/scaleTools/cicd/scDedup/master/230908-ga1c0044/sc_dedup -o $dir/sc_dedup && chmod 755 $dir/sc_dedup
