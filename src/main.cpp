#include "include/molecule.h"
#include <unordered_set>
#include <string>


int main(int argc, char** argv) {

    std::cout << "Hello World!" << std::endl;
    
    // Setup
    std::unordered_set<std::string> allowed_atoms = {"H", "C"};
    std::string filePath = argv[1];
    int p = std::atoi(argv[2]);
    int q = std::atoi(argv[3]);
    
    // Run the calculation
    Molecule molecule = moleculeFromTxt(filePath, allowed_atoms, p, q);
    arma::mat moleculeOverlapMatrix = molecule.overlapMatrix();
    arma::mat moleculeGammaMatrix = molecule.gammaMatrix();
    arma::mat hCoreMatrix = molecule.hCoreMatrix(
      moleculeGammaMatrix,
      moleculeOverlapMatrix
    );
    // std::tuple<arma::mat, arma::mat, 
    //            arma::mat, arma::mat,
    //            arma::mat, arma::mat, 
    //            arma::vec, arma::vec> result = molecule.directInversionIterativeSubspaceAlogrithm(
    //             moleculeGammaMatrix, 
    //             moleculeOverlapMatrix,
    //             hCoreMatrix,
    //             p,
    //             q,
    //             7,
    //             true,
    //             1.0
    // );
      std::tuple<arma::mat, arma::mat, 
               arma::mat, arma::mat,
               arma::mat, arma::mat, 
               arma::vec, arma::vec> result = molecule.selfConsistentFieldAlgorithm(
                moleculeGammaMatrix, 
                moleculeOverlapMatrix,
                hCoreMatrix,
                p,
                q,
                true
    );
    arma::mat molculeFockAlphaMatrix = std::get<0>(result);
    arma::mat moleculeFockBetaMatrix = std::get<1>(result);
    arma::mat moleculeDensityAlphaMatrix = std::get<2>(result);
    arma::mat moleculeDensityBetaMatrix = std::get<3>(result);
    arma::mat molecularOrbitalCoefficientsAlpha = std::get<4>(result);
    arma::mat molecularOrbitalCoefficientsBeta = std::get<5>(result);

    double moleculeNuclearRepulsionEnergy = molecule.nuclearRepulsionEnergy();
    std::cout << "Nuclear Repulsion Energy: " << moleculeNuclearRepulsionEnergy << " eV" << std::endl;

    double moleculeElectronicEnergy = molecule.electronicEnergy(
      molculeFockAlphaMatrix,
      moleculeFockBetaMatrix,
      hCoreMatrix,
      hCoreMatrix,
      moleculeDensityAlphaMatrix,
      moleculeDensityBetaMatrix
    );
    std::cout << "Electronic energy: " << moleculeElectronicEnergy << " eV" << std::endl;
    std::cout << "Testing totalEnergy(): " << molecule.totalEnergy() << std::endl;
    double moleculeTotalEnergy = molecule.totalEnergy(
      molculeFockAlphaMatrix,
      moleculeFockBetaMatrix,
      hCoreMatrix,
      hCoreMatrix,
      moleculeDensityAlphaMatrix,
      moleculeDensityBetaMatrix
    );
    std::cout << "Total energy: " << moleculeTotalEnergy << " eV" << std::endl;

    molecule.energyDerivative(moleculeDensityAlphaMatrix, moleculeDensityBetaMatrix).print("energyDerivative");

    molecule.coordinates.print("coordinates before opt.");

    std::string animationFilePath = filePath; 
    animationFilePath.replace(animationFilePath.length() - 3, animationFilePath.length(), "xyz");
    molecule.steepestDescentGeometryOptimizer(0.75, 1e-2, true, animationFilePath).print("coordinates after opt.");

    return 0;
}