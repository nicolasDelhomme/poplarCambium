#!/bin/sh
set -ex

## create the dirs
## safe
proj='24_Bellini_cambium'
projName='BelliniCambium'

mkdir -p /mnt/picea/storage/projects/$proj/$projName

## w/o backup
mkdir -p /mnt/picea/storage/projects/$proj/$projName/raw

## link the raw data
cd /mnt/picea/storage/projects/$proj/$projName/raw

## directories without re-runs
find /mnt/picea/storage/projects/24_Bellini_cambium/RNAseq_pop_cambium_bellini/ -name "*.gz" -exec mv "{}" . \;


