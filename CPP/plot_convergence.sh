#!/bin/bash

# Fichier d'entrée
input_file="./Output/convergence_errors.dat"

# Fichier de sortie
output_file="./Plot/convergence_plot.png"

# Script Gnuplot pour générer le graphique
gnuplot -e "
    set logscale xy;
    set term pngcairo size 1280,720 enhanced font 'Verdana,10';
    set output '${output_file}';
    set title 'Convergence of Numerical Scheme with Order 1/2';
    set xlabel 'N';
    set ylabel 'Error';
    set grid;
    plot '${input_file}' using 1:2 with linespoints title 'Empirical error', \
         '${input_file}' using 1:3 with lines title 'Expected O(1/N)';
"

echo "Graphique de convergence généré : ${output_file}"