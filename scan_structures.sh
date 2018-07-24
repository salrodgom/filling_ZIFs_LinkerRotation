#!/bin/bash
nCPUs=16
for CIFFile in $(ls -rS struc/*.cif) ; do
  structure=$(echo $CIFFile | sed 's/\.cif//g' | sed 's/struc\///g')
  file=$(echo $CIFFile | sed 's/struc\///g')
  main_folder=${structure}_
  if [ ! -d ${main_folder}* ] ; then
   echo "${main_folder} does not exits"
   echo "[in] analysis"
   mv struc/$file .
   bash main.sh $file $nCPUs > salida.${structure}.txt
   mv $file struc/dones
   echo "[out] analysis"
  else
   echo "${main_folder} exits"
   if [ -f $file ] ; then mv $file struc/dones ; fi
  fi
done
exit 0
