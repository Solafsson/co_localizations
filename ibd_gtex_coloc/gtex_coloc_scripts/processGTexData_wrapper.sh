#!/bin/bash

script=$1
function=$2
traitList=$3
gtex_path=$4
ibd_assoc_path=$5
working_dir=$6
jobFile=$7


${script} ${function} ${traitList} ${gtex_path} ${ibd_assoc_path} ${working_dir} ${jobFile} ${LSB_JOBINDEX}

exit $?


