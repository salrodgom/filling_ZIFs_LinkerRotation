#!/usr/bin/gnuplot -persist
set xtics nomirror
set xtics nomirror
set ytics nomirror
set ytics nomirror
unset x2tics
set y2tics nomirror 
set xlabel "Minimisation Step / -" 
set ylabel "(E - E_0) / N_Z_n / eV" 
set y2label "Grad. Norm." 
unset logscale
set logscale y 10
set logscale y2 10
set locale "en_GB.UTF-8"
set fit errorvariables nocovariancevariables
!e=$(grep "l p" logs/minimization.txt | sort -gk6 | awk '{print $6}' | head -n1) ; echo "e = $e" > gp_tmp
!n=$(grep "Zn" 00*topol.cif | wc -l | awk '{print $1}') ; echo "n = $n" >> gp_tmp
load "gp_tmp"
plot 'logs/minimization.txt' u 1:(($4-e)/n) w l t sprintf("E_0/N_Z_n = (%.3f) ",e/n) , '' u 1:5 axes x1y2 w l t 'Grad. Norm.'
