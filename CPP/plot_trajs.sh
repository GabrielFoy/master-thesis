#!/bin/bash

# Nombre de trajectoires à afficher
num_trajectories=5

# Nombre de simulations à traiter
num_simulations=3

# Boucle pour chaque simulation
for sim in $(seq 1 $num_simulations); do
    # Déterminer la valeur de la ligne en fonction de la simulation
    if [ "$sim" -eq 3 ]; then
        line_value=0.1
    else
        line_value=1.67
    fi

    # Boucle pour chaque trajectoire
    for i in $(seq 1 $num_trajectories); do
        gnuplot -e "
        set term pngcairo size 1280,720 enhanced font 'Verdana,10';
        set output './Plot/simulation_${sim}_traj_${i}.png';
        set title 'Sample ${i} - Simulation ${sim}';
        set xlabel 'Time';
        set ylabel 'Values';
        set key left top;
        plot './Output/simulation_${sim}_traj_${i}.dat' using 1:2 with lines title 'Y Trajectory ${i}', \
             './Output/simulation_${sim}_traj_${i}.dat' using 1:3 with lines title 'Z Trajectory ${i}', \
             './Output/simulation_${sim}_traj_${i}.dat' using 1:4 with lines title 'V Trajectory ${i}', \
             './Output/simulation_${sim}_traj_${i}.dat' using 1:5 with lines title 'S Trajectory ${i}', \
             $line_value with lines title 'E[Y_inf] = $line_value' lc rgb 'black' lw 2;"
    done
done