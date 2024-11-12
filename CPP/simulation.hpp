#ifndef SIMULATION_HPP
#define SIMULATION_HPP

// Importation des librairies utiles
#include <cmath>
#include <gsl/gsl_integration.h>
#include <Eigen/Dense>
#include <random>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <format>
#include <thread>
#include <future>

// Structure pour stocker les différentes simulations
struct SimulationResult
{
    Eigen::MatrixXd Y;
    Eigen::MatrixXd Z;
    Eigen::MatrixXd V;
    Eigen::MatrixXd S;
};

// Déclaration des fonctions
double b_ij(int i, int j, double alpha);

double integrand(double v, void *params);

double b_func(double y, double mu, double lam);

double sigma(double y, double a, double b, double c);

double deterministic_term(int l, double xi, int k, double mu, double lamb, double H, double dt);

Eigen::MatrixXd build_G(int N, 
                        double dt, 
                        double H, 
                        const Eigen::MatrixXd &L_B_int,
                        const Eigen::VectorXd &V, 
                        const Eigen::VectorXd &U,
                        const Eigen::MatrixXd &K, 
                        std::mt19937 &rng);

SimulationResult simulate_model_mp(double s0, 
                                   double eta, 
                                   double mu, 
                                   double lamb, 
                                   double a, 
                                   double b, 
                                   double c,
                                   double H, 
                                   double T, 
                                   double xi0, 
                                   int N, 
                                   int size_sample, 
                                   double rho,
                                   std::mt19937 &rng);

#endif
