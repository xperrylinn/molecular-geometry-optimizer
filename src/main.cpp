#include "include/molecule.h"
#include <unordered_set>
#include <string>


int main(int argc, char** argv) {
    
    // Setup
    if (argc != 5) {
      throw std::invalid_argument("Invalid number of arguments. Expected four arguments: input txt data file, number of alpha electrons, number of beta electrons, convergence algorithm"); 
    }

    std::unordered_set<std::string> allowed_atoms = {"H", "C"};
    std::string filePath = argv[1];
    int p = std::atoi(argv[2]);
    int q = std::atoi(argv[3]);
    std::string scfAlgo = argv[4];

    std::tuple<arma::mat, arma::mat, arma::mat, arma::mat, arma::mat, arma::mat, arma::vec, arma::vec, int> result;

    std::string animationFilePath = filePath; 
    animationFilePath.replace(animationFilePath.length() - 3, animationFilePath.length(), "xyz");

    std::cout << "Hello World! Running geometry optimization on input file " << filePath << std::endl;

    // Run the convergence
    Molecule molecule = moleculeFromTxt(filePath, allowed_atoms, p, q);
    arma::mat moleculeOverlapMatrix = molecule.overlapMatrix();
    arma::mat moleculeGammaMatrix = molecule.gammaMatrix();
    arma::mat moleculeHCoreMatrix = molecule.hCoreMatrix(
      moleculeGammaMatrix,
      moleculeOverlapMatrix
    );
    if (scfAlgo == "fixedPointIteration") {
      result = molecule.scfFixedPointIteration(
                moleculeGammaMatrix, 
                moleculeOverlapMatrix,
                moleculeHCoreMatrix,
                p,
                q,
                false
      );
    } else if (scfAlgo == "DIIS") {
      result = molecule.scfDIIS(
        moleculeGammaMatrix, 
        moleculeOverlapMatrix,
        moleculeHCoreMatrix,
        p,
        q,
        7,
        false
      );
    } else {
      throw std::invalid_argument("Received an invalid convergence algorithm option. Expected one of the following {\"DIIS\", \"fixedPointIteration\"}");
    }

    arma::mat molculeFockAlphaMatrix = std::get<0>(result);
    arma::mat moleculeFockBetaMatrix = std::get<1>(result);
    arma::mat moleculeDensityAlphaMatrix = std::get<2>(result);
    arma::mat moleculeDensityBetaMatrix = std::get<3>(result);
    arma::mat molecularOrbitalCoefficientsAlpha = std::get<4>(result);
    arma::mat molecularOrbitalCoefficientsBeta = std::get<5>(result);
    int scfIterationCount = std::get<8>(result);

    // Calculate energies
    double moleculeNuclearRepulsionEnergy = molecule.nuclearRepulsionEnergy();
    std::cout << "Nuclear Repulsion Energy: " << moleculeNuclearRepulsionEnergy << " eV" << std::endl;

    double moleculeElectronicEnergy = molecule.electronicEnergy(
      molculeFockAlphaMatrix,
      moleculeFockBetaMatrix,
      moleculeHCoreMatrix,
      moleculeHCoreMatrix,
      moleculeDensityAlphaMatrix,
      moleculeDensityBetaMatrix
    );
    std::cout << "Electronic energy: " << moleculeElectronicEnergy << " eV" << std::endl;

    double moleculeTotalEnergy = molecule.totalEnergy(
      molculeFockAlphaMatrix,
      moleculeFockBetaMatrix,
      moleculeHCoreMatrix,
      moleculeHCoreMatrix,
      moleculeDensityAlphaMatrix,
      moleculeDensityBetaMatrix
    );
    std::cout << "Total energy: " << moleculeTotalEnergy << " eV" << std::endl;

    std::cout << "SCF iteration count: " << scfIterationCount << std::endl;

    // Run geometry optimization
    molecule.steepestDescentGeometryOptimizer(
      0.75, 
      1e-2,
      scfAlgo,
      true, 
      animationFilePath
      ).print("coordinates after opt.");

    return 0;
}