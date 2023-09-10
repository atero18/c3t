#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include <math.h>

/* Création de fonctions permettant de calculer
 * divers distances de liaison entre clusters
 * via un objet `AbstractSymMat` (R) contenant toutes
 * les distances (déjà calculées) inter-éléments.
*/

#include "distInter.h"

//' @description Calcule toutes les distances inter-clusters demandées
//' à partir de distances inter-éléments stockées dans un objet
//' `SymVMat` (R).
//' @name distanceInterSymVMat
//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector distanceInterSymVMat(const arma::vec& values,
                                         const arma::uvec& partition,
                                         const Rcpp::NumericMatrix& indexs,
                                         const unsigned int& dim,
                                         const bool& aDefautDiag,
                                         const double& defaultDiag,
                                         const Rcpp::String& comp)
{
  const unsigned int nbDistances = indexs.nrow();
  Rcpp::NumericVector distancesInter(nbDistances);
  for (unsigned int k = 0; k < nbDistances; k++)
  {

    const unsigned int cluster1 = indexs(k, 0);
    arma::uvec elementsCluster1 = find(partition == cluster1);

    const unsigned int cluster2 = indexs(k, 1);
    arma::uvec elementsCluster2 = find(partition == cluster2);

    distancesInter[k] = compSMSymVMat(values, dim, elementsCluster1,
                                      elementsCluster2, comp,
                                      aDefautDiag, defaultDiag);
  }

  return distancesInter;
}

//' @description Calcule toutes les distances inter-clusters demandées
//' à partir de distances inter-éléments stockées dans un objet
//' `SymVMat` (R).
//' @name distanceInterSymMMat
//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector distanceInterSymMMat(const arma::mat& values,
                                         const arma::uvec& partition,
                                         const Rcpp::NumericMatrix& indexs,
                                         const Rcpp::String& comp)
{
  const unsigned int nbDistances = indexs.nrow();
  Rcpp::NumericVector distancesInter(nbDistances);
  unsigned int cluster1, cluster2;

  for (unsigned int k = 0; k < nbDistances; k++)
  {
    cluster1 = indexs(k, 0);
    arma::uvec elementsCluster1 = find(partition == cluster1);

    cluster2 = indexs(k, 1);
    arma::uvec elementsCluster2 = find(partition == cluster2);

    distancesInter[k] = compSMSymMMat(values, elementsCluster1,
                                      elementsCluster2, comp);
  }

  return distancesInter;
}


//' @description Recherche de la position d'un élément d'une matrice symétrique
//' stockée sous forme d'un vecteur.
//'
//' @param i,j position en ligne et colonne (démarre à 0L) de l'élément.
//' (entiers compris entre 0 et `dim - 1`)
//' @param dim dimension de la matrice. (entier positif)
//' @param aDefautDiag `true` si la matrice a une valeur par défaut
//' sur sa diagonale (ex : 0.0 pour une matrice de distance), `false` sinon.
//' (booléen)
//' @noRd
unsigned int posElementSymVMat(unsigned int i,
                               unsigned int j,
                               const int& dim,
                               const bool& aDefautDiag)
{

  unsigned int pos = 0;

  if (i > j)
  {
    unsigned int temp = i;
    i = j;
    j = temp;
  }
  else if(aDefautDiag && i == j)
    return 0;


  // Avancée sur les lignes
  if (aDefautDiag)
    pos += i * dim - i * (i + 1) / 2;
  else
    pos += i * (dim + 1) - i * (i + 1) / 2;

  // Avancée sur les colonnes
  pos += aDefautDiag ? j - i - 1 : j - i;


  return pos;
}

const std::string ERRUNAVAILABLELINKAGE =
  "Asked linkage distance is unavailable";

//' @description Effectue le calcul d'une distance inter-clusters avec
//' une matrice de distances `SymVMat` (R).
//' @name compSMSymVMat
//' @noRd
float compSMSymVMat(const arma::vec& values,
                    const unsigned int& dim,
                    const arma::uvec& lignes,
                    const arma::uvec& colonnes,
                    const Rcpp::String& comp,
                    const bool& aDefautDiag,
                    const double& defaultDiag)
{
  float valeurElement = 0.0;

  if (comp == "min" || comp == "max")
  {
    float valeur;

    if (comp == "min")
      valeur = INFINITY;
    else
      valeur = -1.0;

    for(const auto& i: lignes)
    {

      for(const auto& j: colonnes)
      {
        if(aDefautDiag && i == j)
          valeurElement = defaultDiag;
        else
        {
          const unsigned int posElement = posElementSymVMat(i, j, dim, aDefautDiag);
          valeurElement = values(posElement);
        }

        if((comp == "min" && valeurElement < valeur) ||
           (comp == "max" && valeurElement > valeur))
        {
          valeur = valeurElement;
        }

      }
    }
    return valeur;

  }

  else if (comp == "mean")
  {

    float valeur = 0.0;

    for(const auto& i: lignes)
    {
      for(const auto& j: colonnes)
      {

        if(aDefautDiag && i == j)
          valeur = defaultDiag;
        else
        {
          const unsigned int posElement = posElementSymVMat(i, j, dim, aDefautDiag);
          valeurElement = values(posElement);
        }
        valeur += valeurElement;
      }
    }

    return valeur / (lignes.size() * colonnes.size());

  }

  else if (comp == "hausdorff")
  {

    // Maximum sur les lignes
    float maxLignes = -1.0;

    for(const auto& i: lignes)
    {
      float minLigne = INFINITY;

      for (const auto& j: colonnes)
      {
        if (aDefautDiag && i == j)
          valeurElement = defaultDiag;

        else
        {
          const unsigned int posElement = posElementSymVMat(i, j, dim, aDefautDiag);
          valeurElement = values(posElement);
        }

        if(valeurElement < minLigne)
          minLigne = valeurElement;
      }

      if(minLigne > maxLignes)
        maxLignes = minLigne;
    }

    // Maximum sur les colonnes
    float maxColonnes = -1.0;

    for(const auto& j: colonnes)
    {
      float minColonne = INFINITY;
      for(const auto& i: lignes)
      {
        if(aDefautDiag && i == j)
          valeurElement = defaultDiag;
        else
        {
          const unsigned int posElement = posElementSymVMat(i, j, dim, aDefautDiag);
          valeurElement = values(posElement);
        }

        if(valeurElement < minColonne)
          minColonne = valeurElement;

      }

      if(minColonne > maxColonnes)
        maxColonnes = minColonne;

    }

    return maxLignes > maxColonnes ? maxLignes : maxColonnes;
  }
  else
  {
    Rcpp::stop(ERRUNAVAILABLELINKAGE);
    return -1.;
  }
}


//' @description Effectue le calcul d'une distance inter-clusters avec
//' une matrice de distances `SymMMat` (R).
//' @name compSMSymMMat
//' @noRd
float compSMSymMMat(const arma::mat& values,
                    const arma::uvec& lignes,
                    const arma::uvec& colonnes,
                    const Rcpp::String& comp)
{
  if (comp == "min")
    return values.submat(lignes, colonnes).min();

  else if (comp == "max")
    return values.submat(lignes, colonnes).max();

  else if (comp == "mean")
  {
    float res = 0.;
    for (const auto& i : lignes)
    {
      for (const auto& j: colonnes)
        res += values(i,j);
    }
    return res / (lignes.size() * colonnes.size());
  }

  else if (comp == "hausdorff")
  {
    arma::mat sousMat = values.submat(lignes, colonnes);
    float max1 = arma::min(sousMat, 0).max();
    float max2 = arma::min(sousMat, 1).max();
    return max1 >= max2 ? max1 : max2;
  }
  else
  {
    Rcpp::stop(ERRUNAVAILABLELINKAGE);
    return -1.;
  }
}
