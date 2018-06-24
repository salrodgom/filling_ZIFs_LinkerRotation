#!/bin/bash
nCPUs=8
for CIFFile in $(ls -rS struc/*.cif) ; do
  structure=$(echo $CIFFile | sed 's/\.cif//g' | sed 's/struc\///g')
  file=$(echo $CIFFile | sed 's/struc\///g')
  #if [ ! -d $structure_* ] ; then
   cp struc/$file .
   bash main.sh $file $nCPUs > salida.${structure}.txt
   mv $file struc/dones
  #elif [ -f $file ] ; then
  # mv $file struc/dones
  #fi
done
exit 0
