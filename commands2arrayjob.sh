#!/bin/bash
# $Id: commands2arrayjob.sh,v 1.3 2006/05/16 14:17:31 sbender Exp $

if [ -z "$1" ] ; then
    echo "Usage: $0 <commandfile> [qsub options]"
    exit 1;
fi

CMD_FILE=$1
shift

if ! echo $CMD_FILE | grep -q ^/ ; then
    CMD_FILE=$PWD/$CMD_FILE
fi

MAX=$(wc -l $CMD_FILE | awk '{print $1}')

qsub -t 1-$MAX $@ /net/shendure/vol1/home/mkircher/bin/generic_array_job.csh "$CMD_FILE"
