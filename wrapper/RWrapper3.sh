#!/bin/bash

SCRIPT=$1
VAR1=$2
VAR2=$3
VAR3=$4

Rscript --vanilla "${SCRIPT}" "${VAR1}" "${VAR2}" "${VAR3}"
