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
    Molecule molecule = moleculeFromTxt(filePath, allowed_atoms);
    arma::mat moleculeOverlapMatrix = molecule.overlapMatrix();
    arma::mat moleculeGammaMatrix = molecule.gammaMatrix();
    arma::mat hCoreMatrix = molecule.hCoreMatrix(
      moleculeGammaMatrix,
      moleculeOverlapMatrix
    );
    std::tuple<arma::mat, arma::mat, 
               arma::mat, arma::mat, 
               arma::mat, arma::mat, 
               arma::vec, arma::vec> result = molecule.directInversionIterativeSubspaceAlogrithm(
                moleculeGammaMatrix, 
                moleculeOverlapMatrix,
                hCoreMatrix,
                p,
                q,
                7,
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

    double moleculeTotalEnergy = molecule.totalEnergy(
      molculeFockAlphaMatrix,
      moleculeFockBetaMatrix,
      hCoreMatrix,
      hCoreMatrix,
      moleculeDensityAlphaMatrix,
      moleculeDensityBetaMatrix
    );
    std::cout << "Total energy: " << moleculeTotalEnergy << " eV" << std::endl;

    std::cout << "Homework 5" << std::endl;

    arma::mat moleculeOverlapMatrixPositionDerivative = molecule.overlapMatrixPositionDerivative();
    // moleculeOverlapMatrixPositionDerivative.print("OV_RA");

    arma::mat moleculeGammaMatrixPositionDerivative = molecule.gammaMatrixPositionDerivative();
    // moleculeGammaMatrixPositionDerivative.print("gamma_RA");

    arma::mat moleculeNuclearRepulsionDerivative = molecule.nuclearRepulsionEnergyPositionDerivative();
    // moleculeNuclearRepulsionDerivative.print("gradient (Nuclear part)");

    arma::mat moleculeXCoefficientMatrix = molecule.xCoefficientMatrix(moleculeDensityAlphaMatrix, moleculeDensityBetaMatrix);
    // moleculeXCoefficientMatrix.print("X Coefficient Matrix");

    arma::mat moleculeYCoefficientMatrix = molecule.yCoefficientMatrix(moleculeDensityAlphaMatrix, moleculeDensityBetaMatrix);
    // moleculeYCoefficientMatrix.print("Y Coefficient Matrix");

    arma::mat moleculeEnergyDerivative = molecule.energyDerivative(
      moleculeOverlapMatrixPositionDerivative,
      moleculeGammaMatrixPositionDerivative,
      moleculeXCoefficientMatrix,
      moleculeYCoefficientMatrix,
      moleculeNuclearRepulsionDerivative
    );
    // moleculeEnergyDerivative.print("gradient");

    return 0;
}