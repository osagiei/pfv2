#!/bin/sh

CODE=$(dirname "$0")

mkdir -p $CODE/classes

javac -cp $CODE/src:$CODE/classes:$CODE/lib/commons-lang3-3.2.1.jar $CODE/src/bio/igm/entities/*.java -d $CODE/classes

javac -cp $CODE/src:$CODE/classes:$CODE/lib/commons-lang3-3.2.1.jar $CODE/src/bio/igm/*/*/*.java -d $CODE/classes

jar cfvm $CODE/PFv2.jar $CODE/manifest.mf -C $CODE/classes/ .

cp $CODE/PFv2.jar $CODE/scripts

export PATH=$PATH:$CODE/scripts
