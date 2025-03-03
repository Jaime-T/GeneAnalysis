#!/usr/bin/env bash
 
DIRECTORY=$(cd `dirname $0` && pwd)
 
port=8002
 
cmd='screen -dmS depdash '$DIRECTORY'/depdash.sh '$port
 
echo $cmd
echo 'start depdash on http://127.0.0.1:'$port
 
eval $cmd