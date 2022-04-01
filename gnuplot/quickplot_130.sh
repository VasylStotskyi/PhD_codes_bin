#!/bin/bash
grep " MD| Cell lengths " tmp.*.out | grep ang |awk '{print NR/2"\t"$6"\t"$7"\t"$8}' > cell.dat
grep "Pressure" tmp.*.*|grep bar |awk '{print NR/2"\t"$5"\t"$6}' > pressure.dat
gnuplot < gnuplot_input_130.in
