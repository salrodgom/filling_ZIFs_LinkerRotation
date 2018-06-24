#!/bin/bash
nCPUs=26
for CIFFile in $(ls struc/*.cif) ; do
  structure=$(echo $CIFFile | sed 's/\.cif//g' | sed 's/struc\///g')
  file=$(echo $CIFFile | sed 's/struc\///g')
  cp struc/$file .
  bash main.sh $file $nCPUs > salida.${structure}.txt 
done
exit 0
