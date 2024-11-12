/*
@author : Gabriel FOY.

L'objectif de ce fichier c++ est d'implémenter toutes les fonctions nécessaires à la simulation du modèle pris en exemple dans le papier de recherche.
*/

// Importation des librairies utiles
#include "simulation.hpp"
#include <chrono>

// Définition des fonctions utiles à la simulation des différents processus
/**
 * @brief Fonction intégrée utilisée pour calculer les coefficients de la matrice B_int.
 *
 * @param v La variable d'intégration.
 * @param params Les paramètres i, j et alpha de la fonction à intégrer.
 *
 * @return Valeur de l'intégrale pour les paramètres donnés.
 */
double integrand(double v, void *params)
{
    double *p = (double *)params;
    int i = (int)p[0];
    int j = (int)p[1];
    double alpha = p[2];

    return std::pow(i + v, alpha) * std::pow(j + v, alpha);
}

/**
 * @brief Calcule le coefficient (B_ij) de la matrice B_int.
 *
 * @param i L'indice de ligne pour le coefficient de la matrice : 0 <= i <= N-1.
 * @param j L'indice de colonne pour le coefficient de la matrice : 0 <= j <= N-1.
 * @param alpha Paramètre de forme qui influence la partie fractionnaire du noyau de convolution.
 *
 * @return Le coefficient B_ij de la matrice B_int.
 */
double b_ij(int i, int j, double alpha)
{
    if (i == j)
    {
        return (1.0 / (2 * alpha + 1)) * (std::pow(i + 1, 2 * alpha + 1) - std::pow(i, 2 * alpha + 1));
    }
    else
    {
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);

        double result, error;
        double params[3] = {double(i), double(j), alpha};

        gsl_function F;
        F.function = &integrand;
        F.params = &params;

        int status = gsl_integration_qag(&F, 0, 1, 0, 1e-7, 1000, 6, workspace, &result, &error);

        if (status)
        { // Erreur lors de l'intégration
            std::string formatted_str = std::format("Erreur lors de l'intégration de la fonction : i={}, j={}, alpha={}", i, j, alpha);
            throw std::runtime_error(formatted_str);
        }

        gsl_integration_workspace_free(workspace);

        return result;
    }
}

/**
 * @brief Calcule le terme de drift pour le processus éléphant.
 *
 * @param y La valeur actuelle du processus "poisson rouge" pour lequel le drift est calculé.
 * @param mu Le taux de rendement moyen attendu.
 * @param lam Taux de retour à la moyenne pour le processus.
 *
 * @return La valeur du drift pour le processus éléphant.
 */
double b_func(double y, double mu, double lamb)
{
    return mu - lamb * y;
}

/**
 * @brief Permet de calculer le terme déterministe dans la somme pour simuler le processus éléphant au temps k+1.
 *
 * @param l Indice de la somme.
 * @param xi La valeur du processus "poisson rouge" à l'indice l : 0 <= l <= k.
 * @param k Indice représentant la fin de l'intervalle temporel pour le calcul (i.e indice max de la somme).
 * @param mu Le taux de rendement moyen attendu.
 * @param lamb Taux de retour à la moyenne.
 * @param H Paramètre de Hurst, influençant le degré de mémoire du processus.
 * @param dt dt = T/N, le pas de temps.
 *
 * @return La contribution déterministe du terme pour les indices l et k.
 */
double deterministic_term(int l, double xi, int k, double mu, double lamb, double H, double dt)
{
    return b_func(xi, mu, lamb) * (std::pow(k - l + 1, H + 1.0 / 2) - std::pow(k - l, H + 1.0 / 2)) * std::pow(dt, H + 1.0 / 2) * (1.0 / (H + 1.0 / 2));
}

/**
 * @brief Permet de calculer le terme de dérive pour le processus éléphant.
 *
 * @param y La valeur actuelle du processus "poisson rouge".
 * @param a Coefficient influençant la dépendance de la volatilité par rapport à l'écart entre y et b.
 * @param b Valeur autour de laquelle la volatilité est ajustée.
 * @param c Terme constant assurant que la volatilité ne tombe jamais en dessous d'un certain seuil.
 *
 * @return La valeur du terme de dérive pour le processus éléphant selon la valeur du processus "poisson rouge"
 */
double sigma(double y, double a, double b, double c)
{
    return std::sqrt(a * std::pow(y - b, 2) + c);
}

/**
 * @brief Construit la matrice G utilisée pour simuler les processus stochastiques.
 *
 * @param N Nombre de points de discrétisation dans le temps.
 * @param dt Pas de temps, T/N où T est le temps total de simulation.
 * @param H Paramètre de Hurst dans le processus.
 * @param L_B_int Matrice obtenue après la décomposition de Cholesky de la matrice B_int.
 * @param V Vecteur V utilisé dans la simulation.
 * @param U Vecteur U utilisé dans la simulation.
 * @param K Matrice utilisée dans la décomposition de Cholesky de la matrice de covariance globale.
 * @param rng Générateur de nombres pseudo aléatoires.
 *
 * @return Eigen::MatrixXd
 */
Eigen::MatrixXd build_G(int N, double dt, double H, const Eigen::MatrixXd &L_B_int,
                        const Eigen::VectorXd &V, const Eigen::VectorXd &U,
                        const Eigen::MatrixXd &K, std::mt19937 &rng)
{
    Eigen::MatrixXd G_traj = Eigen::MatrixXd::Zero(N + 2, N);

    std::normal_distribution<> dist(0.0, 1.0); // Loi normale centrée réduite

    for (int l = 1; l <= N; ++l)
    {
        // Décomposition de Cholesky de la sous matrice B_l
        Eigen::MatrixXd L_B_l = std::pow(dt, H) * L_B_int.block(0, 0, N - l + 1, N - l + 1);

        // Sous vecteur V_l
        Eigen::VectorXd V_l = V.head(N - l + 1);

        // Sous vecteur U_l
        Eigen::VectorXd U_l = U.head(N - l + 1);

        // Sous matrice s_l
        Eigen::MatrixXd s_l(N - l + 1, 2);
        s_l.col(0) = V_l;
        s_l.col(1) = U_l;

        // Décomposition de Cholesky de la matrice de covariance globale Sigma_l
        Eigen::MatrixXd delta_l = L_B_l.lu().solve(s_l);

        Eigen::MatrixXd Temp = K - delta_l.transpose() * delta_l;

        // On s'assure que la matrice Temp est semi définie positive
        Temp += (1e-9) * Eigen::MatrixXd::Identity(2, 2);

        Eigen::MatrixXd M;
        try
        {
            M = Temp.llt().matrixL();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Erreur lors du calcul de la décomposition de Cholesky de la matrice K - delta_l^T * delta_l: "
                      << e.what() << std::endl;
            throw std::runtime_error("Cholesky decomposition failed");
        }

        // Initialisation de la matrice L_l
        Eigen::MatrixXd L_l = Eigen::MatrixXd::Zero(N - l + 3, N - l + 3);
        L_l.block(0, 0, N - l + 1, N - l + 1) = L_B_l;

        // Ajout de delta_l dans L_l
        L_l.block(N - l + 1, 0, 2, N - l + 1) = delta_l.transpose();

        // Ajout de M dans L_l
        L_l.block(N - l + 1, N - l + 1, 2, 2) = M;

        // Accroissements browniens
        Eigen::VectorXd w_acc = Eigen::VectorXd::NullaryExpr(N - l + 3, [&]()
                                                             { return dist(rng); });

        // l-ième colonne de G simulée
        G_traj.block(l - 1, l - 1, N + 2 - l + 1, 1) = L_l * w_acc;
    }

    std::cout << "Fin de la construction d'une matrice G" << std::endl;

    return G_traj;
}

// Fonction principale implémentant la simulation des différents processus
/**
 * @brief
 * Simule les trajectoires d'un actif financier selon le modèle présenté dans l'article.
 * Ce modèle est particulier dans le sens où la volatilité de l'actif possède un effet de mémoire.
 *
 * @param s0 Prix initial de l'actif financier.
 * @param eta Facteur d'échelle pour lisser la volatilité des processus "poisson rouge" et "éléphant".
 * @param mu  Taux de rendement moyen attendu de l'actif.
 * @param lamb Taux de retour à la moyenne pour le modèle de "poisson rouge".
 * @param a Coefficient de l'échelle de la variance dans le modèle de l'actif.
 * @param b Point moyen autour duquel la variance fluctue.
 * @param c Terme constant ajouté à la variance. Seuil minimal / Variance "planchée"
 * @param H Paramètre de Hurst, contrôlant le degré de mémoire du processus "éléphant".
 * @param T Durée totale de la période de simulation, en années.
 * @param xi0 Condition initiale pour le processus "éléphant".
 * @param N N+1 points de temps pour un pas de T/N.
 * @param size_sample Nombre de trajectoires indépendantes à simuler.
 * @param rho Taux de corrélation entre les mouvements browniens guidant les processus modélisant l'actif financier et la volatilité.
 * @param rng Générateur de nombres pseudo aléatoires.
 *
 * @return SimulationResult
 */
SimulationResult simulate_model_mp(double s0, double eta, double mu, double lamb, double a, double b, double c,
                                   double H, double T, double xi0, int N, int size_sample, double rho,
                                   std::mt19937 &rng)
{
    double dt = T / N;

    // -------- Construction de la matrice B pour la simulation des processus --------
    std::cout << "Construction de la matrice B_int" << std::endl;
    Eigen::MatrixXd B_int = Eigen::MatrixXd::Zero(N, N);

    for (int i = 0; i < N; ++i)
    {
        for (int j = i; j < N; ++j)
        {
            B_int(i, j) = b_ij(i, j, H - 0.5);
        }
    }

    B_int += B_int.transpose().eval() - Eigen::MatrixXd(B_int.diagonal().asDiagonal());

    // -------- Décomposition de Cholesky de la matrice B --------
    // On s'assure que la matrice B_int est semi définie positive en gommant les erreurs numériques
    B_int += 1e-9 * Eigen::MatrixXd::Identity(N, N);

    std::cout << "Décomposition de Cholesky de la matrice B_int" << std::endl;
    const auto start_b{std::chrono::steady_clock::now()};

    Eigen::MatrixXd L_B_int;
    try
    {
        L_B_int = B_int.llt().matrixL();
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erreur lors du calcul de la décomposition de Cholesky de la matrice B_int: "
                  << e.what() << std::endl;
        throw std::runtime_error("Cholesky decomposition failed");
    }

    const auto end_b{std::chrono::steady_clock::now()};
    std::cout << "Fin de la décomposition de Cholesky de la matrice B_int. Temps : " << std::chrono::duration<double>{end_b - start_b} << std::endl;

    Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(N, 0, N - 1);
    Eigen::VectorXd V = std::pow(T / N, H + 0.5) * (1.0 / (H + 0.5)) * ((t.array() + 1).pow(H + 0.5) - t.array().pow(H + 0.5));
    Eigen::VectorXd U = rho * V;

    Eigen::MatrixXd K(2, 2);
    K << T / N, rho * T / N,
        rho * T / N, T / N;

    // -------- Construction des matrices G : une matrice par trajectoire --------
    std::cout << "Construction des matrices G" << std::endl;
    const auto start_g{std::chrono::steady_clock::now()};

    std::vector<Eigen::MatrixXd> list_G(size_sample);

    auto build_task = [&](int index)
    {
        // Le code ici va appeler la fonction build_G pour chaque échantillon
        list_G[index] = build_G(N, dt, H, L_B_int, V, U, K, rng);
    };

    std::vector<std::future<void>> futures;
    for (int i = 0; i < size_sample; ++i)
    {
        futures.push_back(std::async(std::launch::async, build_task, i));
    }

    for (auto &f : futures)
    {
        f.get();
    }
    const auto end_g{std::chrono::steady_clock::now()};
    std::cout << "Fin de la construction des matrices G. Temps : " << std::chrono::duration<double>{end_g - start_g} << std::endl
              << "Construction des trajectoires des processus." << std::endl;

    // -------- Simulation des processus --------
    Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(N + 1, size_sample);
    Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(N + 1, size_sample);
    Z.row(0).setConstant(xi0);
    Eigen::MatrixXd V_matrix = Eigen::MatrixXd::Zero(N + 1, size_sample);
    V_matrix.row(0).setConstant(a * std::pow(xi0 - b, 2) + c);
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(N + 1, size_sample);
    S.row(0).setConstant(s0);

    for (int n = 1; n <= N; ++n)
    {
        for (int num_traj = 0; num_traj < size_sample; ++num_traj)
        {
            Eigen::MatrixXd &G_traj = list_G[num_traj];

            // Simulation du processus markovien "poisson rouge"
            Y(n, num_traj) = Y(n - 1, num_traj) + xi0 * (1 / std::tgamma(1.5 - H)) *
                            (std::pow(n * dt, 0.5 - H) - std::pow((n - 1) * dt, 0.5 - H)) +
                            dt * b_func(Y(n - 1, num_traj), mu, lamb) +
                            eta * sigma(Y(n - 1, num_traj), a, b, c) * G_traj(N, n - 1);

            // Simulation du processus "éléphant"
            Z(n, num_traj) = xi0;

            // Calcul de la somme pour le processus "éléphant"
            for (int l = 0; l < n; ++l)
            {
                // Terme déterministe
                double det_term = deterministic_term(l, Y(l, num_traj), n, mu, lamb, H, dt);
                
                // Terme stochastique
                double sto_term = eta * sigma(Y(l, num_traj), a, b, c) * G_traj(l, n - 1);

                Z(n, num_traj) += (det_term + sto_term);
            }

            Z(n, num_traj) *= 1 / std::tgamma(H + 0.5);

            // Simulation du processus de Variance
            V_matrix(n, num_traj) = a * std::pow(Z(n, num_traj) - b, 2) + c;

            // Simulation de l'actif financier
            S(n, num_traj) = S(n - 1, num_traj) + S(n - 1, num_traj) * std::sqrt(V_matrix(n - 1, num_traj)) * G_traj(N + 1, n - 1);
        }
    }

    return {Y, Z, V_matrix, S};
}