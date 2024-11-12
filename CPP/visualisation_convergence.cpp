/*
@author : Gabriel FOY.

L'objectif de ce fichier c++ est d'implémenter les fonctions nécessaires à la production du graphique permettant de visualiser l'ordre de convergence
du schéma numérique du processus de Volterra pris en exemple (i.e le processus "éléphant").
*/
// Importation des librairies utiles
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include <fstream>

// Définition des fonctions utiles
// Fonction de simulation de l'EDS avec schéma d'Euler-Maruyama
Eigen::VectorXd simulate_eds(int N, double T)
{
    double dt = T / N;
    Eigen::VectorXd X = Eigen::VectorXd::Zero(N + 1);
    std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> dist(0.0, std::sqrt(dt));

    for (int k = 1; k <= N; ++k)
    {
        X[k] = X[k - 1] + std::sqrt(dt)*dist(rng);
    }

    return X;
}

// Fonction pour calculer les erreurs de convergence
std::vector<double> compute_errors(const std::vector<int> &Ns, double T, int num_simulations = 5)
{
    std::vector<std::vector<Eigen::VectorXd>> liste_trajs(num_simulations);

    // Simuler pour chaque N et chaque simulation
    for (int sim = 0; sim < num_simulations; ++sim)
    {
        for (int N : Ns)
        {
            liste_trajs[sim].push_back(simulate_eds(N, T));
        }
    }

    std::vector<double> errors(Ns.size() - 1, 0.0);

    // Calculer les erreurs pour chaque N
    for (int sim = 0; sim < num_simulations; ++sim)
    {
        for (size_t i = 0; i < Ns.size() - 1; ++i)
        {
            // Extraire un élément sur deux de liste_trajs[sim][i+1]
            Eigen::VectorXd half_segment(Ns[i] + 1);
            for (int j = 0; j < Ns[i]; ++j)
            {
                half_segment[j] = liste_trajs[sim][i + 1][2 * j];
            }

            Eigen::VectorXd diff = half_segment - liste_trajs[sim][i];
            double max_error = diff.cwiseAbs2().maxCoeff();
            errors[i] += max_error;
        }
    }

    // Moyenne des erreurs sur toutes les simulations
    for (double &error : errors)
    {
        error /= num_simulations;
    }

    return errors;
}

// Fonction pour sauvegarder les erreurs dans un fichier .dat
void save_errors(const std::string &filename, const std::vector<int> &Ns, const std::vector<double> &errors)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
        return;
    }

    for (size_t i = 0; i < errors.size(); ++i)
    {
        file << Ns[i] << " " << errors[i] << " " << (errors[0] * (static_cast<double>(Ns[0]) / Ns[i])) << "\n";
    }

    file.close();
}

// Execution des caluls / Point d'entrée de l'executable
int main()
{
    double T = 1.0;
    std::vector<int> rs = {7, 8, 9, 10};
    std::vector<int> Ns;

    // Générer les valeurs de Ns comme des puissances de 2
    for (int r : rs)
    {
        Ns.push_back(1 << r); // Opérateur de déclage binaire : très optimisé pour obtenir des puissances de 2
    }

    // Calculer les erreurs de convergence
    std::vector<double> errors = compute_errors(Ns, T, 1000);

    // Sauvegarder les erreurs dans un fichier .dat
    save_errors("./Output/convergence_errors.dat", Ns, errors);

    std::cout << "Les erreurs de convergence ont été sauvegardées dans convergence_errors.dat" << std::endl;

    return 0;
}
