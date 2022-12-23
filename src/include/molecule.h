#include <armadillo>
#include <unordered_set>


class Molecule {
    public:
      /* 
      Data members 
      */
      // Data structures
      int numberOfAtoms;
      arma::mat coordinates;
      arma::ivec atomicNumbers;
      arma::imat basisFunctions;
      int numberAtomicBasisFunctions;
      int numberElectronPairs;
      int numberAlphaElectrons, numberBetaElectrons;
      std::map<int, std::map<int, arma::vec>> contractionCoefficientMap;
      std::map<int, std::map<int, arma::vec>> contractionExponentMap;
      std::map<int, std::map<int, double>> ionizationEnergyAtomicNumberAngularMomentumMap;
      std::map<int, std::map<int, double>> ionizationAndElectronAffinityEnergyTermMap;
      std::map<int, int> atomicNumberBondingParameterMap;
      std::map<int, int> valenceAtomicNumberMap;
      std::map<int, std::string> atomicNumberSymbolMap;

      // Constants
      double electronVoltsToAtomicUnitsConversionFactor = 27.211;
      double atomicNumberUnitsToAngstromConversionFactor = 0.529177249;

      /* 
      Constructors 
      */
     // Data structure arguments
      Molecule(arma::mat coords, arma::ivec atomicNums, int p, int q);
      
      /* 
      Member functions 
      */
      arma::imat constructBasisFunctions();
      int numAtomicBasisFunctions();
      arma::mat overlapMatrix();
      arma::mat extendedHuckleHamiltonianMatrix(arma::mat overlapMatrix);
      arma::mat orthogonalizationTransfromation(arma::mat S);
      arma::mat hamiltonianOrthogonalizedBasis(arma::mat S_neg_sqrt, arma::mat hamiltonianMatrix);      
      arma::mat molecularOrbitalCoefficientMatrix(arma::mat S_neg_sqrt, arma::mat hamiltonianMatrix);
      arma::mat molecularOrbitalOverlapMatrix(
        arma::mat molecularOrbitalCoefficientMatrix, 
        arma::mat hamiltonianMatrixH2, 
        arma::mat overlapMatrix
      );
      arma::mat gammaMatrix();
      arma::mat fockOperatorMatrix(
        arma::mat gammaMatrix, 
        arma::mat densityMatrix, 
        arma::mat overlapMatrix,
        arma::vec pTot
      );
      arma::mat hCoreMatrix(
        arma::mat gammaMatrix,
        arma::mat overlapMatrix
      );
      arma::mat densityMatrix(
        arma::mat molecularOrbitalCoeffMatrix, 
        int numberOfElectrons
      );
      std::pair<arma::mat, arma::vec> huckleTheoryEigenValueProblem(
        arma::mat molecularOrbitalCoefficientMatrix, 
        arma::mat hamiltonianMatrix, 
        arma::mat overlapMatrix
      );
      double molecularEnergy(arma::mat molecularOrbitalCoefficientMatrix, arma::mat hamiltonianMatrix, arma::mat overlapMatrix);
      double bondingEnergy(double molecularEnergy);
      double gammaTwoCenterTwoElectronRepulsionIntegral(
        arma::vec contractionCoefficientsA,
        arma::vec contractionCoefficientsB,
        arma::vec contractionExponentsA,
        arma::vec contractionExponentsB,
        arma::vec Ra,
        arma::vec Rb,
        arma::ivec La,
        arma::ivec Lb
      );
      double fockOperatorMatrixElementOffDiagonal(
        double overlapMatrixElem, 
        double densityMatrixElem, 
        double gammaAB, 
        double betaA,
        double betaB
      );
      double fockOperatorMatrixElementDiagonal(
        int atomIndex,
        double ionizationAndElectronAffinityTerm,
        double pAlphaMuMu,
        arma::mat gammaMatrix,
        arma::vec pTot
      );
      double hCoreMatrixElementOffDiagonal(
        double overlapMatrixElem, 
        double densityMatrixElem,
        double betaA,
        double betaB
      );
      double hCoreMatrixElementDiagonal(
        int atomIndex,
        double ionizationAndElectronAffinityTerm,
        arma::mat gammaMatrix
      );
      double boysSigma(double alpha, double alphaPrime);
      double boysU(double boysSigma);
      double boysVSquared(double boysSigmaA, double boysSigmaB);
      double boysT(arma::vec Ra, arma::vec Rb, double boysVSquared);
      double atomToAtomEuclideanDistance(arma::vec Ra, arma::vec Rb);
      double boysIntegral(
        arma::vec Ra,
        arma::vec Rb,
        double boysUA, 
        double boysUB, 
        double boysT,
        double boysVSquared
      );
      arma::vec boysIntegralPositionDerivative(
        arma::vec Ra,
        arma::vec Rb,
        double boysUA, 
        double boysUB, 
        double boysT,
        double boysVSquared
      );
      std::tuple<
      arma::mat, arma::mat, 
      arma::mat, arma::mat, 
      arma::mat, arma::mat, 
      arma::vec, arma::vec> selfConsistentFieldAlgorithm(
        arma::mat gammaMatrix, 
        arma::mat overlapMatrix,
        arma::mat hCoreMatrix,
        int pElectrons,
        int qElectrons,
        bool consoleLog
      );
      std::tuple<
      arma::mat, arma::mat, 
      arma::mat, arma::mat, 
      arma::mat, arma::mat, 
      arma::vec, arma::vec> directInversionIterativeSubspaceAlogrithm(
        arma::mat gammaMatrix, 
        arma::mat overlapMatrix,
        arma::mat hCoreMatrix,
        int pElectrons,
        int qElectrons,
        int numPrevIters,
        bool consoleLog,
        double lambda
      );
      arma::mat generateLagrangeMultiplierMatrix(std::deque<arma::mat> errorMatrices);
      double electronicEnergy(
        arma::mat fockAlphaMatrix,
        arma::mat fockBetaMatrix,
        arma::mat hCoreAlphaMatrix,
        arma::mat hCoreBetaMatrix,
        arma::mat densityBetaMatrix,
        arma::mat densityAlphaMatrix
      );
      double nuclearRepulsionEnergy();
      arma::mat nuclearRepulsionEnergyPositionDerivative();
      double totalEnergy(
        arma::mat fockAlphaMatrix,
        arma::mat fockBetaMatrix,
        arma::mat hCoreAlphaMatrix,
        arma::mat hCoreBetaMatrix,
        arma::mat densityBetaMatrix,
        arma::mat densityAlphaMatrix
      );
      double totalEnergy();
      arma::vec densityPerAtom(arma::mat densityMatrixAlpha, arma::mat densityMatrixBeta);
      arma::mat overlapMatrixPositionDerivative();
      arma::vec gammaTwoCenterTwoElectronRepulsionIntegralPositionDerivative(
        arma::vec contractionCoefficientsA,
        arma::vec contractionCoefficientsB,
        arma::vec contractionExponentsA,
        arma::vec contractionExponentsB,
        arma::vec Ra,
        arma::vec Rb,
        arma::ivec La,
        arma::ivec Lb
      );
      arma::mat gammaMatrixPositionDerivative();
      arma::mat xCoefficientMatrix(
        arma::mat densityMatrixAlpha,
        arma::mat densityMatrixBeta
      );
      arma::mat yCoefficientMatrix(
        arma::mat densityMatrixAlpha,
        arma::mat densityMatrixBeta
      );
      arma::mat energyDerivative(
        arma::mat overlapMatrixPositionDerivative,
        arma::mat gammaMatrixPositionDerivative,
        arma::mat xCoefficientMatrix,
        arma::mat yCoefficientMatrix,
        arma::mat nuclearRepulsionEnergyPositionDerivative
      );
      arma::mat energyDerivative(arma::mat densityAlphaMatrix, arma::mat densityBetaMatrix);
      arma::mat geometryOptimizer(std::string optimizer);
      arma::mat steepestDescentGeometryOptimizer(double stepSize, double tolerance, bool logging, std::string animationPath);
      void xyzCoordinatesToStream(std::ofstream& ofs, std::string commentLine);
};

Molecule moleculeFromTxt(std::string rel_file_path, std::unordered_set<std::string> allowed_symbols, int p, int q);
