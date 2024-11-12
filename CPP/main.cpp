/*
@author : Gabriel FOY.

L'objectif de ce fichier c++ est de lancer les simulations du modèle implémenté selon les trois contextes pris en exemples dans le papier de recherche.
*/

// Importation des librairies utiles
#include "simulation.hpp"
#include <fstream>


// Fonction pour enregistrer les données dans un fichier
void save_data(const std::string &filename, 
               const std::vector<double> &times,
               const Eigen::MatrixXd &Y, 
               const Eigen::MatrixXd &Z,
               const Eigen::MatrixXd &V, 
               const Eigen::MatrixXd &S, 
               int sample_index)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Erreur d'ouverture du fichier " << filename << std::endl;
        return;
    }

    for (std::size_t n = 0; n < times.size(); ++n)
    {
        file << times[n] << " " << Y(n, sample_index) << " " << Z(n, sample_index) << " "
             << V(n, sample_index) << " " << S(n, sample_index) << "\n";
    }

    file.close();
}

int main()
{
    // Paramètres de la simulation
    double s0 = 1;       // Prix initial de l'actif
    double T = 5;        // Maturité / Temps d'arrêt de la simulation
    int N = 500;        // N+1 points de discrétisation de l'intervalle [0, T]
    int size_sample = 5; // Nombre de trajectoires à générer

    double H = 0.1; // Paramètre de Hurst

    // Paramètres de la fonction sigma
    double a = 0.384;
    double b = 0.095;
    double c = 0.0025;

    // Taux de corrélation des browniens
    double rho = -0.7;

    // Discrétisation du temps
    std::vector<double> times(N + 1);
    for (int i = 0; i <= N; ++i)
    {
        times[i] = i * T / N;
    }

    // Simulation 1
    double eta = 0.1;
    double mu = 2;
    double lamb = 1.2;
    double xi0 = 0;

    std::mt19937 rng;
    SimulationResult result1 = simulate_model_mp(s0, eta, mu, lamb, a, b, c, H, T, xi0, N, size_sample, rho, rng);

    // Sauvegarde des données pour la première simulation
    for (int i = 0; i < size_sample; ++i)
    {
        save_data("./Output/simulation_1_traj_" + std::to_string(i + 1) + ".dat", times, result1.Y, result1.Z, result1.V, result1.S, i);
    }
    std::cout << "Fin de la première simulation." << std::endl << std::endl;

    // Simulation 2
    xi0 = mu / lamb;
    SimulationResult result2 = simulate_model_mp(s0, eta, mu, lamb, a, b, c, H, T, xi0, N, size_sample, rho, rng);

    // Sauvegarde des données pour la deuxième simulation
    for (int i = 0; i < size_sample; ++i)
    {
        save_data("./Output/simulation_2_traj_" + std::to_string(i + 1) + ".dat", times, result2.Y, result2.Z, result2.V, result2.S, i);
    }
    std::cout << "Fin de la deuxième simulation." << std::endl << std::endl;

    // Simulation 3
    eta = 0.01;
    lamb = 20;
    xi0 = mu / lamb;
    SimulationResult result3 = simulate_model_mp(s0, eta, mu, lamb, a, b, c, H, T, xi0, N, size_sample, rho, rng);

    // Sauvegarde des données pour la troisième simulation
    for (int i = 0; i < size_sample; ++i)
    {
        save_data("./Output/simulation_3_traj_" + std::to_string(i + 1) + ".dat", times, result3.Y, result3.Z, result3.V, result3.S, i);
    }
    std::cout << "Fin de la troisième simulation." << std::endl;

    return 0;
}
