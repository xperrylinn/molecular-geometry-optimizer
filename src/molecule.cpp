#include "include/molecule.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <string>
#include <cstring>
#include <iterator>
#include <unordered_set>
#include <deque>
#include <map>
#include <stdexcept>
#include <math.h>
#include <armadillo>
#include <limits>
#include <cstdlib>
#include <tuple>


Molecule moleculeFromTxt(std::string rel_file_path, std::unordered_set<std::string> allowed_symbols, int p, int q) {
  
  // Attempt to read in file
  std::ifstream txtFile(rel_file_path);

  // If sucessfully open, parse into matrix
  if (txtFile.is_open()) {
    std::vector<int> disallowedSymbols;
    int atomCountActual = 0, atomCountRef, molecularCharge, atomicNumber;
    double x, y, z;
    std::string line;

    // Read line by line and convert to a row vector
    txtFile >> atomCountRef;
    txtFile >> molecularCharge;
    arma::mat coordinates(0, 3);
    arma::ivec atomicNumbers(atomCountRef);
    std::getline(txtFile, line);
    int i = 0;  // Index iterator for atomicNumbers
    for (;std::getline(txtFile, line);) {
      arma::rowvec rowVec(3);
      // Parse atomic number and coordinates from line
      std::istringstream ss(line);
      ss >> atomicNumber;
      ss >> x;
      ss >> y;
      ss >> z;
      rowVec[0] = x;
      rowVec[1] = y;
      rowVec[2] = z;
      // Load data into data structures
      atomicNumbers[i] = atomicNumber;
      coordinates = arma::join_cols(coordinates, rowVec);
      // Increment index and counter
      i += 1;
      atomCountActual += 1;

    }

    // Validate inputs, collect any errors, and throw at user
    bool invalid_atom_symbol_error = disallowedSymbols.size() > 0;
    bool invalid_atom_count_error = atomCountActual != atomCountRef;
    if (invalid_atom_count_error || invalid_atom_symbol_error)  {
      std::string error_string = "The following invalid arguments were encountered:\n";
      if (invalid_atom_symbol_error) {
        error_string += "- disallowed atoms were encountered: ";
        for (auto & symbol : disallowedSymbols) {
          error_string += std::to_string(symbol) + " ";
        }
        error_string += "\n";
      }
      if (invalid_atom_count_error) {
        error_string += ("- the txt file header stated " + std::to_string(atomCountRef) + " atoms, but " + std::to_string(atomCountActual) + " were found.\n");
      }
      throw std::invalid_argument(error_string);
    }

    // Instantiate Molecule and return
    Molecule molecule(coordinates, atomicNumbers, p, q);
    return molecule;
  } else {
    // Molecule molecule;
    throw std::invalid_argument("Unable to open file.");
  }
};

void starsAndBarsMatrix(int numStars, int numBars, std::vector<std::vector<int>>& combos, std::vector<int> prev) {
    if (numStars == 0) {
        for (int i = 0; i <= numBars; i += 1) {
            prev.push_back(0);
        }
        combos.push_back(prev);
    } else if (numBars == 0) {
        prev.push_back(numStars);
        combos.push_back(prev);
    } else {
        for (int i = 0; i <= numStars; i += 1) {
            std::vector<int> curr(prev);
            curr.push_back(i);
            starsAndBarsMatrix(numStars - i, numBars - 1, combos, curr);
        }
    }
}

arma::imat starsAndBarsMatrix(int numStars, int numBars) {
    
    // Compute each combo into vector of vectors, then transfer to imat
    std::vector<std::vector<int>> combos;
    std::vector<int> prev;

    // Run the recusrive alogrithm to generate combos
    starsAndBarsMatrix(numStars, numBars, combos, prev);

    // Create matrix for return
    arma::imat iMat(0, 3);
    for (auto & combo: combos) {
        arma::irowvec iRowVec =  arma::conv_to<arma::irowvec>::from(combo);
        iMat = arma::join_cols(iMat, iRowVec);
    }
    return iMat;
}

Molecule::Molecule(arma::mat coords, arma::ivec atomicNums, int p, int q) {
  numberOfAtoms = atomicNums.size();
  coordinates = coords;
  atomicNumbers = atomicNums;
  numberAlphaElectrons = p;
  numberBetaElectrons = q;
  /* 
  Populate maps 
  Contraction Coefficients and Constants are sourced from sto-3g.1.json, https://www.basissetexchange.org/
  */
  // Hyrogen
  arma::vec hydrogen1sContrCoeffs = {0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00};
  arma::vec hydrogen1sContrExps = {0.3425250914E+01, 0.6239137298E+00, 0.1688554040E+00};
  std::map<int, arma::vec> hydrogenContrCoeffMap = {
    {0, hydrogen1sContrCoeffs},
  };
  std::map<int, arma::vec> hydrogenContrExpMap = {
    {0, hydrogen1sContrExps},
  };
  // Carbon
  arma::vec carbon2sContrCoeffs = {-0.9996722919E-01, 0.3995128261E+00, 0.7001154689E+00};
  arma::vec carbon2pContrCoeffs = {0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00};
  arma::vec carbon2sContrExps = {0.2941249355E+01, 0.6834830964E+00, 0.2222899159E+00};
  arma::vec carbon2pContrExps = {0.2941249355E+01, 0.6834830964E+00, 0.2222899159E+00};
  std::map<int, arma::vec> carbonContrCoeffMap = {
    {0, carbon2sContrCoeffs},
    {1, carbon2pContrCoeffs},
  };
  std::map<int, arma::vec> carbonContrExpMap = {
    {0, carbon2sContrExps},
    {1, carbon2pContrExps},
  };

  // Aggregate map
  std::map<int, std::map<int, arma::vec>> contrCoeffMap = {
    {1, hydrogenContrCoeffMap},
    {6, carbonContrCoeffMap},
  };
  std::map<int, std::map<int, arma::vec>> contrExpMap = {
    {1, hydrogenContrExpMap},
    {6, carbonContrExpMap},
  };
  contractionCoefficientMap = contrCoeffMap;
  contractionExponentMap = contrExpMap;

  // Calculate number of basis functions
  numberAtomicBasisFunctions = numAtomicBasisFunctions();
  numberElectronPairs = numberAtomicBasisFunctions / 2;

  // Construct basis functions
  basisFunctions = constructBasisFunctions();

  // Create map of ionization energy angular momentum atomic number map
  std::map<int, double> hydrogenAnularMomentumIonizationMap = {{0, -13.6}};
  std::map<int, double> carbonAnularMomentumIonizationMap = {{0, -21.4}, {1, -11.4}};
  ionizationEnergyAtomicNumberAngularMomentumMap = {
    {1, hydrogenAnularMomentumIonizationMap},
    {6, carbonAnularMomentumIonizationMap},
  };

  // Crate ionization energy and electron affinity energy map
  std::map<int, double> hydrogenAnularIonizationAndElectronAffinityTermMap = {{0, 7.176}};
  std::map<int, double> carbonAnularIonizationAndElectronAffinityTermMap = {{0, 14.051}, {1, 5.572}};
  std::map<int, double> nitrogenAnularIonizationAndElectronAffinityTermMap = {{0, 19.316}, {1, 7.275}};
  std::map<int, double> oxygenAnularIonizationAndElectronAffinityTermMap = {{0, 25.390}, {1, 9.111}};
  std::map<int, double> flourineAnularIonizationAndElectronAffinityTermMap = {{0, 32.272}, {1, 11.080}};
  ionizationAndElectronAffinityEnergyTermMap = {
    {1, hydrogenAnularIonizationAndElectronAffinityTermMap},
    {6, carbonAnularIonizationAndElectronAffinityTermMap},
    {7, nitrogenAnularIonizationAndElectronAffinityTermMap},
    {8, oxygenAnularIonizationAndElectronAffinityTermMap},
    {9, flourineAnularIonizationAndElectronAffinityTermMap},
  };

  // Create atomic bonding parameters map
  atomicNumberBondingParameterMap = {
    {1, -9}, 
    {6, -21}, 
    {7, -25}, 
    {8, -31},
    {9, -39},
  };

  // Create map of valence atmoic number map
  valenceAtomicNumberMap = {
    {1, 1}, 
    {6, 4}, 
    {7, 5}, 
    {8, 6},
    {9, 7},
    };
}

arma::imat Molecule::constructBasisFunctions() {
  arma::imat basisFunctions(0, 6);  // Number of basis functions x 4 (atom index, atomic number, x, y, z)
  for (int i = 0; i < atomicNumbers.size(); i += 1) {
    int atomicNumber = atomicNumbers[i];
    for (std::pair<int, arma::vec> pair : contractionCoefficientMap[atomicNumber]) {
      int L = pair.first;
      arma::imat tempAngularMomentumCombos = starsAndBarsMatrix(L, 2);
      arma::imat tempAtomIndex(arma::size(tempAngularMomentumCombos, 0), 1, arma::fill::value(i));
      arma::imat tempAtomicNums(arma::size(tempAngularMomentumCombos, 0), 1, arma::fill::value(atomicNumber));
      arma::imat tempAngularMomentums(arma::size(tempAngularMomentumCombos, 0), 1, arma::fill::value(L));
      arma::imat tempBaisFunctions = arma::join_rows(tempAtomIndex, tempAtomicNums);
      tempBaisFunctions = arma::join_rows(tempBaisFunctions, tempAngularMomentums);
      tempBaisFunctions = arma::join_rows(tempBaisFunctions, tempAngularMomentumCombos);
      basisFunctions = arma::join_cols(tempBaisFunctions, basisFunctions);
    }
  }
  return basisFunctions;
}

int Molecule::numAtomicBasisFunctions() {
  std::map<int, int> atomicNumberBasisFunctionMap = {
    {1, 1},
    {6, 4},
  };
  int result = 0;
  for (int i = 0; i < atomicNumbers.size(); i += 1) {
    int atomicNumber = atomicNumbers[i];
    result += atomicNumberBasisFunctionMap[atomicNumber];
  }
  return result;
}

bool assertEquals(double val1, double val2, double tolerance) {
    double delta = abs(val1 - val2);
    return delta < tolerance;
}

arma::mat atomMatrixFromTxt(std::string rel_file_path, std::unordered_set<std::string> allowed_symbols) {
  
  // Attempt to read in file
  std::ifstream txtFile(rel_file_path);

  // If sucessfully open, parse into matrix
  if (txtFile.is_open()) {
    std::vector<int> disallowedSymbols;
    int atomCountActual = 0, atomCountRef, molecularCharge, atomicNumber;
    double x, y, z;
    std::string line;

    // Read line by line and convert to a row vector
    txtFile >> atomCountRef;
    txtFile >> molecularCharge;
    arma::mat result(0, 4);
    std::getline(txtFile, line);
    for (; std::getline(txtFile, line); ) {
        arma::rowvec rowVec(line);
        result = arma::join_cols(result, rowVec);
        atomCountActual += 1;
    }

    // Validate inputs, collect any errors, and throw at user
    bool invalid_atom_symbol_error = disallowedSymbols.size() > 0;
    bool invalid_atom_count_error = atomCountActual != atomCountRef;
    if (invalid_atom_count_error || invalid_atom_symbol_error)  {
      std::string error_string = "The following invalid arguments were encountered:\n";
      if (invalid_atom_symbol_error) {
        error_string += "- disallowed atoms were encountered: ";
        for (auto & symbol : disallowedSymbols) {
          error_string += std::to_string(symbol) + " ";
        }
        error_string += "\n";
      }
      if (invalid_atom_count_error) {
        error_string += ("- the txt file header stated " + std::to_string(atomCountRef) + " atoms, but " + std::to_string(atomCountActual) + " were found.\n");
      }
          throw std::invalid_argument(error_string);
    }
    return result;
  } else {
    arma::mat result;
    throw std::invalid_argument("Unable to open file.");
  }
};

int numBasisFunctions(arma::mat m) {
  std::map<int, int> atomicNumberBasisFunctionMap = {
    {1, 1},
    {6, 4},
  };
  int result = 0;
  for (int i = 0; i < arma::size(m, 0); i += 1) {
    int atomicNumber = m.at(i, 0);
    result += atomicNumberBasisFunctionMap[atomicNumber];
  }
  return result;
}

int numBasisFunctions(arma::ivec atomicNumbers) {
  std::map<int, int> atomicNumberBasisFunctionMap = {
    {1, 1},
    {6, 4},
  };
  int result = 0;
  for (int i = 0; i < atomicNumbers.size(); i += 1) {
    int atomicNumber = atomicNumbers[i];
    result += atomicNumberBasisFunctionMap[atomicNumber];
  }
  return result;
}

int numElectrons(arma::mat m) {
  int numElectrons = numBasisFunctions(m);
  return numElectrons;
}

int numElectronPairs(arma::mat m) {
  double intPart;
  if (std::modf(numElectrons(m), &intPart) != 0.0) {
    throw std::invalid_argument("Not integer number of electron pairs");
  }
  int numElectronsPairs = numBasisFunctions(m) / 2;
  return numElectronsPairs;
}

int doubleFactorial(int value) {
    if (value < -1) {
        throw std::invalid_argument("Double factorial of values less than -1 are not defined.");
        return 0;
    } else if (value <= 1) {
        return 1;
    } else if (value == 2) {
        return 2;
    } else {
        return value * doubleFactorial(value - 2);
    }
}

int binomialCoefficient(int n, int k) {
    if (k == 0 || k == n) {
        return 1;
    } else {
        return binomialCoefficient(n - 1, k - 1) + binomialCoefficient(n - 1, k);
    }
}

arma::vec centerProduct(arma::vec Ra, arma::vec Rb, double alpha, double beta) {
    arma::vec result = alpha / (alpha + beta) * Ra + beta / (alpha + beta) * Rb; 
    return result;
}

int numStarsAndBarsCombos(int numStars, int numBars) {
    return binomialCoefficient(numBars + numStars, numBars);
}

double gaussianStandardIntegral(double Xa, double Xb, double Xp, int la,  int lb, double alpha, double beta) {
    double exponentialPrefactor = exp(-1 * (alpha * beta * pow(Xa - Xb, 2)) / (alpha + beta));
    double sqrtFactor = sqrt(M_PI / (alpha + beta));
    double doubleSummmation = 0.0;
    for (int i = 0; i <= la; i += 1) {
        double innerSum = 0.0;
        for (int j = 0; j <= lb; j += 1) {
            if ((i + j) % 2 == 0) {
                innerSum += \
                  binomialCoefficient(la, i) * \
                  binomialCoefficient(lb, j) * \
                  (
                    doubleFactorial(i + j - 1) * \
                    pow(Xp - Xa, la - i) * \
                    pow(Xp - Xb, lb - j)) / \
                    pow(2 * (alpha + beta), double(i + j) / 2.0
                  ); 
            }
        }
        doubleSummmation += innerSum;
    }
    double indefiniteIntegral = exponentialPrefactor * sqrtFactor * doubleSummmation;
    return indefiniteIntegral;
}

double gaussianStandardIntegral(arma::vec Rp, arma::vec Ra, arma::vec Rb, arma::ivec La, arma::ivec Lb, double alpha, double beta) {

    double integral = 1.0;
    integral *= gaussianStandardIntegral(Ra[0], Rb[0], Rp[0], La[0],  Lb[0], alpha, beta);
    integral *= gaussianStandardIntegral(Ra[1], Rb[1], Rp[1], La[1],  Lb[1], alpha, beta);
    integral *= gaussianStandardIntegral(Ra[2], Rb[2], Rp[2], La[2],  Lb[2], alpha, beta);

    return integral;
}

arma::vec gaussianStandardIntegralPositionDerivative(arma::vec Rp, arma::vec Ra, arma::vec Rb, arma::ivec La, arma::ivec Lb, double alpha, double beta) {

    arma::vec result(3);

    double integralDerivativeX = 1.0;
    integralDerivativeX *= (
      -1 * La[0] * gaussianStandardIntegral(Ra[0], Rb[0], Rp[0], La[0] - 1,  Lb[0], alpha, beta) +
      2 * alpha * gaussianStandardIntegral(Ra[0], Rb[0], Rp[0], La[0] + 1,  Lb[0], alpha, beta)
    );
    integralDerivativeX *= gaussianStandardIntegral(Ra[1], Rb[1], Rp[1], La[1],  Lb[1], alpha, beta);
    integralDerivativeX *= gaussianStandardIntegral(Ra[2], Rb[2], Rp[2], La[2],  Lb[2], alpha, beta);

    double integralDerivativeY = 1.0;
    integralDerivativeY *= gaussianStandardIntegral(Ra[0], Rb[0], Rp[0], La[0],  Lb[0], alpha, beta);
    integralDerivativeY *= (
      -1 * La[1] * gaussianStandardIntegral(Ra[1], Rb[1], Rp[1], La[1] - 1,  Lb[1], alpha, beta) +
      2 * alpha * gaussianStandardIntegral(Ra[1], Rb[1], Rp[1], La[1] + 1,  Lb[1], alpha, beta)
    );
    integralDerivativeY *= gaussianStandardIntegral(Ra[2], Rb[2], Rp[2], La[2],  Lb[2], alpha, beta);

    double integralDerivativeZ = 1.0;
    integralDerivativeZ *= gaussianStandardIntegral(Ra[0], Rb[0], Rp[0], La[0],  Lb[0], alpha, beta);
    integralDerivativeZ *= gaussianStandardIntegral(Ra[1], Rb[1], Rp[1], La[1],  Lb[1], alpha, beta);
    integralDerivativeZ *= (
      -1 * La[2] * gaussianStandardIntegral(Ra[2], Rb[2], Rp[2], La[2] - 1,  Lb[2], alpha, beta) +
      2 * alpha * gaussianStandardIntegral(Ra[2], Rb[2], Rp[2], La[2] + 1,  Lb[2], alpha, beta)
    );

    result[0] = integralDerivativeX;
    result[1] = integralDerivativeY;
    result[2] = integralDerivativeZ;

    return result;
}

arma::mat overlapMatrix(arma::vec Ra, arma::vec Rb, int La, int Lb, double alpha, double beta) {
    
    // int numShellFxnsA = shellDimensionsMap[La], numShellFxnsB = shellDimensionsMap[Lb];
    arma::vec Rp = centerProduct(Ra, Rb, alpha, beta);
    int numShellFxnsA = numStarsAndBarsCombos(La, 2), numShellFxnsB = numStarsAndBarsCombos(Lb, 2);

    // Get the possible combinations of angular momentum
    arma::imat lACombosMat = starsAndBarsMatrix(La, 2);
    arma::imat lBCombosMat = starsAndBarsMatrix(Lb, 2);

    // Create a matrix for populating with overlap matrix elements
    arma::mat m(numShellFxnsA, numShellFxnsB);
    
    for (int i = 0; i < numShellFxnsA; i += 1) {
        int la = lACombosMat.at(i, 0), ma = lACombosMat.at(i, 1), na = lACombosMat.at(i, 2);
        for (int j = 0; j < numShellFxnsB; j += 1) {
            int lb = lBCombosMat.at(j, 0), mb = lBCombosMat.at(j, 1), nb = lBCombosMat.at(j, 2);
            double product = 1.0;
            product *= gaussianStandardIntegral(Ra[0], Rb[0], Rp[0], la,  lb, alpha, beta);
            product *= gaussianStandardIntegral(Ra[1], Rb[1], Rp[1], ma,  mb, alpha, beta);
            product *= gaussianStandardIntegral(Ra[2], Rb[2], Rp[2], na,  nb, alpha, beta);
            m.at(i, j) = product;
        }
    }
    m = arma::reverse(arma::reverse(m, 1));
    return m;
}

double primitiveGaussianNormConst(arma::vec R, arma::ivec L, double alpha) {
  arma::vec Rp = centerProduct(R, R, alpha, alpha);
  double integral = gaussianStandardIntegral(Rp, R, R, L, L, alpha, alpha);
  double normalizationConst = 1 / sqrt(integral);
  return normalizationConst;
}

double contractedGaussianOverlap(arma::vec contrCoeffsA, arma::vec contrExpsA, arma::vec contrCoeffsB, arma::vec contrExpsB, arma::vec Ra, arma::vec Rb, arma::ivec La, arma::ivec Lb) {
  double integral = 0.0;
  for (int i = 0; i < contrCoeffsA.size(); i += 1) {
    for (int j = 0; j < contrCoeffsB.size(); j += 1) {
      double contrCoeffA = contrCoeffsA[i], contrCoeffB = contrCoeffsB[j];
      double contrExpA = contrExpsA[i], contrExpB = contrExpsB[j];
      int la = La[0], ma = La[1], na = La[2], lb = Lb[0], mb = Lb[1], nb = Lb[2];
      arma::vec Rp = centerProduct(Ra, Rb, contrExpA, contrExpB);
      double contrNormConstA = primitiveGaussianNormConst(Ra, La, contrExpA);
      double contrNormConstB = primitiveGaussianNormConst(Rb, Lb, contrExpB);
      integral += contrCoeffA * contrCoeffB * contrNormConstA  * contrNormConstB * gaussianStandardIntegral(Rp, Ra, Rb, La, Lb, contrExpA, contrExpB);
    }
  }
  return integral;
}

arma::vec contractedGaussianOverlapPositionDerivative(arma::vec contrCoeffsA, arma::vec contrExpsA, arma::vec contrCoeffsB, arma::vec contrExpsB, arma::vec Ra, arma::vec Rb, arma::ivec La, arma::ivec Lb) {
  arma::vec result(3);
  for (int i = 0; i < contrCoeffsA.size(); i += 1) {
    for (int j = 0; j < contrCoeffsB.size(); j += 1) {
      double contrCoeffA = contrCoeffsA[i], contrCoeffB = contrCoeffsB[j];
      double contrExpA = contrExpsA[i], contrExpB = contrExpsB[j];
      int la = La[0], ma = La[1], na = La[2], lb = Lb[0], mb = Lb[1], nb = Lb[2];
      arma::vec Rp = centerProduct(Ra, Rb, contrExpA, contrExpB);
      double contrNormConstA = primitiveGaussianNormConst(Ra, La, contrExpA);
      double contrNormConstB = primitiveGaussianNormConst(Rb, Lb, contrExpB);
      result += contrCoeffA * contrCoeffB * contrNormConstA  * contrNormConstB * gaussianStandardIntegralPositionDerivative(Rp, Ra, Rb, La, Lb, contrExpA, contrExpB);
    }
  }
  return result;
}

arma::mat Molecule::fockOperatorMatrix(
  arma::mat gammaMatrix, 
  arma::mat densityMatrix, 
  arma::mat overlapMatrix,
  arma::vec pTot
) 
{
  arma::mat fockOperatorMatrix(numberAtomicBasisFunctions, numberAtomicBasisFunctions);
  for (int i = 0; i < numberAtomicBasisFunctions; i += 1) {
    arma::irowvec basisFunctionA = basisFunctions.row(i);
    int atomIndexA = basisFunctionA[0], atomicNumberA = basisFunctionA[1], totalLa = basisFunctionA[2];
    arma::vec Ra = {coordinates.at(atomIndexA, 0), coordinates.at(atomIndexA, 1), coordinates.at(atomIndexA, 2)};
    arma::ivec La = {basisFunctionA[3], basisFunctionA[4], basisFunctionA[5]};
    double pTotAA = pTot[atomIndexA];
    for (int j = 0; j < numberAtomicBasisFunctions; j += 1) {
      arma::irowvec basisFunctionB = basisFunctions.row(j);
      int atomIndexB = basisFunctionB[0], atomicNumberB = basisFunctionB[1], totalLb = basisFunctionB[2];
      arma::ivec Lb = {basisFunctionB[3], basisFunctionB[4], basisFunctionB[5]};
      arma::vec Rb = {coordinates.at(atomIndexB, 0), coordinates.at(atomIndexB, 1), coordinates.at(atomIndexB, 2)};
      double pTotBB = pTot[atomIndexB];
      double element;
      if (i != j) {
        double overlapMatrixElem = overlapMatrix.at(i, j);
        double densityMatrixElem = densityMatrix.at(i, j);
        double gammaAB = gammaMatrix.at(atomIndexA, atomIndexB);
        double betaA = atomicNumberBondingParameterMap[atomicNumberA];
        double betaB = atomicNumberBondingParameterMap[atomicNumberB];     
        element = fockOperatorMatrixElementOffDiagonal(
          overlapMatrixElem, 
          densityMatrixElem, 
          gammaAB, 
          betaA,
          betaB
        );
      } else {
        double ionizationAndElectronAffinityTerm = ionizationAndElectronAffinityEnergyTermMap[atomicNumberA][totalLb];
        double pAlphaMuMu = densityMatrix.at(i, j);
        element = fockOperatorMatrixElementDiagonal(
          atomIndexA,
          ionizationAndElectronAffinityTerm,
          pAlphaMuMu,
          gammaMatrix,
          pTot
        );
      }
      fockOperatorMatrix.at(i, j) = element;
    }
  }
  return fockOperatorMatrix;
}

arma::mat Molecule::hCoreMatrix(
  arma::mat gammaMatrix,
  arma::mat overlapMatrix
) 
{
  arma::mat zeroDensityMatrix(numberAtomicBasisFunctions, numberAtomicBasisFunctions, arma::fill::zeros);
  arma::vec zeroDensityPerAtom(numberOfAtoms, arma::fill::zeros);
  return fockOperatorMatrix(
    gammaMatrix, 
    zeroDensityMatrix, 
    overlapMatrix,
    zeroDensityPerAtom
  );
}

arma::mat Molecule::extendedHuckleHamiltonianMatrix(arma::mat overlapMatrix) {
  arma::mat result(numberAtomicBasisFunctions, numberAtomicBasisFunctions);
  // Compute diagonal elements
  double diagonalElements[numberAtomicBasisFunctions];
  for (int i = 0; i < numberAtomicBasisFunctions; i += 1) {
    int atomicNumber = basisFunctions.at(i, 1);
    int angularMomentum = basisFunctions.at(i, 2);
    diagonalElements[i] = ionizationEnergyAtomicNumberAngularMomentumMap[atomicNumber][angularMomentum];
  }

  // Compute off diagonal elements
  for (int i = 0; i < numberAtomicBasisFunctions; i += 1) {
    for (int j = 0; j < numberAtomicBasisFunctions; j += 1) {
      if (i == j) {
        result.at(i, j) = diagonalElements[i];
      } else {
        double averageIonizationOverlapWeighted = 1.75 / 2 * (diagonalElements[i] + diagonalElements[j]) * overlapMatrix.at(i, j);
        result.at(i, j) = averageIonizationOverlapWeighted;
      }
    }
  }
  return result;
}

arma::mat Molecule::orthogonalizationTransfromation(arma::mat S) {
  arma::vec S_eval;
  arma::mat S_evec;
  arma::eig_sym(S_eval, S_evec, S);
  arma::mat S_eval_mat = arma::diagmat(S_eval.transform([](double val) { return 1 / sqrt(val); }));
  arma::mat result = S_evec * S_eval_mat * S_evec.t();
  return result;
}

arma::mat Molecule::hamiltonianOrthogonalizedBasis(arma::mat S_neg_sqrt, arma::mat hamiltonianMatrix) {
  arma::mat result = S_neg_sqrt.t() * hamiltonianMatrix * S_neg_sqrt;
  return result;
}

arma::mat Molecule::molecularOrbitalCoefficientMatrix(arma::mat S_neg_sqrt, arma::mat hamiltonianMatrixH2) {
  arma::vec epsilon;
  arma::mat cPrime;
  arma::eig_sym(epsilon, cPrime, hamiltonianMatrixH2);
  arma::mat result = S_neg_sqrt * cPrime;
  return result; 
}

 std::pair<arma::mat, arma::vec>  Molecule::huckleTheoryEigenValueProblem(arma::mat molecularOrbitalCoefficientMatrix, arma::mat hamiltonianMatrix, arma::mat overlapMatrix) {
  arma::vec epsilon;
  arma::mat molecularOrbitalOverlapMat;
  arma::mat tempMat = arma::inv(overlapMatrix * molecularOrbitalCoefficientMatrix) * hamiltonianMatrix * molecularOrbitalCoefficientMatrix;
  arma::eig_sym(epsilon, molecularOrbitalOverlapMat, tempMat);
  std::pair<arma::mat, arma::vec> moOrbitalsAndEnergiesPair(molecularOrbitalOverlapMat, epsilon);
  return moOrbitalsAndEnergiesPair; 
 }

arma::mat Molecule::molecularOrbitalOverlapMatrix(arma::mat molecularOrbitalCoefficientMatrix, arma::mat hamiltonianMatrix, arma::mat overlapMatrix) {
  arma::mat molecularOritbalOverlapMat = huckleTheoryEigenValueProblem(molecularOrbitalCoefficientMatrix, hamiltonianMatrix, overlapMatrix).first;
  return molecularOritbalOverlapMat;
}

double Molecule::molecularEnergy(arma::mat molecularOrbitalCoefficientMatrix, arma::mat hamiltonianMatrix, arma::mat overlapMatrix) {
  arma::vec energyEigenValues = huckleTheoryEigenValueProblem(molecularOrbitalCoefficientMatrix, hamiltonianMatrix, overlapMatrix).second;
  double molecEnergy = 0.0;
  for (int i = numberElectronPairs, j = 0; i > 0;i -= 1, j += 1) {
    molecEnergy += (2 * energyEigenValues[j]);
  }
  return molecEnergy;
}

double Molecule::bondingEnergy(double molecularEnergy) {
  double isolatedAtomicEnergy = 0.0;
  for (int i = 0; i < arma::size(basisFunctions, 0); i += 1) {
    arma::irowvec basisFunction = basisFunctions.row(i);
    int atomicNumber = basisFunction[1], angularMomentum = basisFunction[2];
    isolatedAtomicEnergy += ionizationEnergyAtomicNumberAngularMomentumMap[atomicNumber][angularMomentum];
  }
  double bondingEnergy = molecularEnergy - isolatedAtomicEnergy;
  return bondingEnergy;
}

arma::mat Molecule::overlapMatrix() {
  arma::mat overlapMat(numberAtomicBasisFunctions, numberAtomicBasisFunctions);
  for (int i = 0; i < numberAtomicBasisFunctions; i += 1) {
    arma::irowvec basisFunctionA = basisFunctions.row(i);
    int atomIndexA = basisFunctionA[0], atomicNumberA = basisFunctionA[1], totalLa = basisFunctionA[2];
    arma::vec Ra = {coordinates.at(atomIndexA, 0), coordinates.at(atomIndexA, 1), coordinates.at(atomIndexA, 2)};
    arma::ivec La = {basisFunctionA[3], basisFunctionA[4], basisFunctionA[5]};
    arma::vec contrCoeffsA = contractionCoefficientMap[atomicNumberA][totalLa];
    arma::vec contrExpsA = contractionExponentMap[atomicNumberA][totalLa];
    for (int j = 0; j < numberAtomicBasisFunctions; j += 1) {
      arma::irowvec basisFunctionB = basisFunctions.row(j);
      int atomIndexB = basisFunctionB[0], atomicNumberB = basisFunctionB[1], totalLb = basisFunctionB[2];
      arma::ivec Lb = {basisFunctionB[3], basisFunctionB[4], basisFunctionB[5]};
      arma::vec Rb = {coordinates.at(atomIndexB, 0), coordinates.at(atomIndexB, 1), coordinates.at(atomIndexB, 2)};
      arma::vec contrCoeffsB = contractionCoefficientMap[atomicNumberB][totalLb];
      arma::vec contrExpsB = contractionExponentMap[atomicNumberB][totalLb];
      double overlap = contractedGaussianOverlap(contrCoeffsA, contrExpsA, contrCoeffsB, contrExpsB, Ra, Rb, La, Lb);
      overlapMat.at(i, j) = overlap;
    }
  }
  return overlapMat;
}

arma::mat Molecule::overlapMatrixPositionDerivative() {
  arma::mat overlapMat(3, numberAtomicBasisFunctions * numberAtomicBasisFunctions);
  for (int i = 0; i < numberAtomicBasisFunctions; i += 1) {
    arma::irowvec basisFunctionA = basisFunctions.row(i);
    int atomIndexA = basisFunctionA[0], atomicNumberA = basisFunctionA[1], totalLa = basisFunctionA[2];
    arma::vec Ra = {coordinates.at(atomIndexA, 0), coordinates.at(atomIndexA, 1), coordinates.at(atomIndexA, 2)};
    arma::ivec La = {basisFunctionA[3], basisFunctionA[4], basisFunctionA[5]};
    arma::vec contrCoeffsA = contractionCoefficientMap[atomicNumberA][totalLa];
    arma::vec contrExpsA = contractionExponentMap[atomicNumberA][totalLa];
    for (int j = 0; j < numberAtomicBasisFunctions; j += 1) {
      arma::irowvec basisFunctionB = basisFunctions.row(j);
      int atomIndexB = basisFunctionB[0], atomicNumberB = basisFunctionB[1], totalLb = basisFunctionB[2];
      arma::ivec Lb = {basisFunctionB[3], basisFunctionB[4], basisFunctionB[5]};
      arma::vec Rb = {coordinates.at(atomIndexB, 0), coordinates.at(atomIndexB, 1), coordinates.at(atomIndexB, 2)};
      arma::vec contrCoeffsB = contractionCoefficientMap[atomicNumberB][totalLb];
      arma::vec contrExpsB = contractionExponentMap[atomicNumberB][totalLb];
      arma::vec overlap = contractedGaussianOverlapPositionDerivative(contrCoeffsA, contrExpsA, contrCoeffsB, contrExpsB, Ra, Rb, La, Lb);
      for (int k = 0; k < 3; k += 1) {
        overlapMat.at(k, numberAtomicBasisFunctions * i + j) = overlap[k];
      }
    }
  }
  return overlapMat;
}

double Molecule::gammaTwoCenterTwoElectronRepulsionIntegral(
  arma::vec contractionCoefficientsA,
  arma::vec contractionCoefficientsB,
  arma::vec contractionExponentsA,
  arma::vec contractionExponentsB,
  arma::vec Ra,
  arma::vec Rb,
  arma::ivec La,
  arma::ivec Lb
)
{
  double gamma = 0.0;
  for (int k = 0; k < 3; k += 1) {
    double contractionCoefficientK = contractionCoefficientsA[k];
    double contractionExponentK = contractionExponentsA[k];
    double contractedGaussianNormConstA = primitiveGaussianNormConst(Ra, La, contractionExponentK);
    contractionCoefficientK *= contractedGaussianNormConstA;
    for (int kPrime = 0; kPrime < 3; kPrime += 1) {
      double contractionCoefficientKPrime = contractionCoefficientsA[kPrime];
      double contractionExponentKPrime = contractionExponentsA[kPrime];
      double contractedGaussianNormConstAPrime = primitiveGaussianNormConst(Ra, La, contractionExponentKPrime);
      contractionCoefficientKPrime *= contractedGaussianNormConstAPrime;
      double sigmaA = boysSigma(contractionExponentK, contractionExponentKPrime);
      double UA = boysU(sigmaA);
      for (int l = 0; l < 3; l += 1) {
        double contractionCoefficientL = contractionCoefficientsB[l];
        double contractionExponentL = contractionExponentsB[l];
        double contractedGaussianNormConstB = primitiveGaussianNormConst(Rb, Lb, contractionExponentL);
        contractionCoefficientL *= contractedGaussianNormConstB;
        for (int lPrime = 0; lPrime < 3; lPrime += 1) {
          // Contracted Gassian setup
          double contractionCoefficientLPrime = contractionCoefficientsB[lPrime];
          double contractionExponentLPrime = contractionExponentsB[lPrime];
          double contractedGaussianNormConstBPrime = primitiveGaussianNormConst(Rb, Lb, contractionExponentLPrime);
          contractionCoefficientLPrime *= contractedGaussianNormConstBPrime;
          double sigmaB = boysSigma(contractionExponentL, contractionExponentLPrime);
          double UB = boysU(sigmaB);
          double vSquared = boysVSquared(sigmaA, sigmaB);
          double T = boysT(Ra, Rb, vSquared);
          // Boys function setup
          double byzIntegral = boysIntegral(Ra, Rb, UA, UB, T, vSquared);
          gamma += contractionCoefficientK * \
          contractionCoefficientKPrime * \
          contractionCoefficientL * \
          contractionCoefficientLPrime * \
          byzIntegral;
          
          // std::cout << "d'k=" << contractionCoefficientK << " d'k'=" << contractionCoefficientKPrime << " d'l=" << contractionCoefficientL << " d'l'=" << contractionCoefficientLPrime << std::endl;
          // std::cout << "k=" << k << " k'=" << kPrime << " l=" << l << " l'=" << lPrime << " [o](0)=" << byzIntegral << std::endl;
        }
      }
    }
  }
  // Convert units from a.u. to eV
  gamma *= electronVoltsToAtomicUnitsConversionFactor;
  return gamma;
}

arma::vec Molecule::gammaTwoCenterTwoElectronRepulsionIntegralPositionDerivative(
  arma::vec contractionCoefficientsA,
  arma::vec contractionCoefficientsB,
  arma::vec contractionExponentsA,
  arma::vec contractionExponentsB,
  arma::vec Ra,
  arma::vec Rb,
  arma::ivec La,
  arma::ivec Lb
)
{
  arma::vec result(3, arma::fill::zeros);
  for (int k = 0; k < 3; k += 1) {
    double contractionCoefficientK = contractionCoefficientsA[k];
    double contractionExponentK = contractionExponentsA[k];
    double contractedGaussianNormConstA = primitiveGaussianNormConst(Ra, La, contractionExponentK);
    contractionCoefficientK *= contractedGaussianNormConstA;
    for (int kPrime = 0; kPrime < 3; kPrime += 1) {
      double contractionCoefficientKPrime = contractionCoefficientsA[kPrime];
      double contractionExponentKPrime = contractionExponentsA[kPrime];
      double contractedGaussianNormConstAPrime = primitiveGaussianNormConst(Ra, La, contractionExponentKPrime);
      contractionCoefficientKPrime *= contractedGaussianNormConstAPrime;
      double sigmaA = boysSigma(contractionExponentK, contractionExponentKPrime);
      double UA = boysU(sigmaA);
      for (int l = 0; l < 3; l += 1) {
        double contractionCoefficientL = contractionCoefficientsB[l];
        double contractionExponentL = contractionExponentsB[l];
        double contractedGaussianNormConstB = primitiveGaussianNormConst(Rb, Lb, contractionExponentL);
        contractionCoefficientL *= contractedGaussianNormConstB;
        for (int lPrime = 0; lPrime < 3; lPrime += 1) {
          // Contracted Gassian setup
          double contractionCoefficientLPrime = contractionCoefficientsB[lPrime];
          double contractionExponentLPrime = contractionExponentsB[lPrime];
          double contractedGaussianNormConstBPrime = primitiveGaussianNormConst(Rb, Lb, contractionExponentLPrime);
          contractionCoefficientLPrime *= contractedGaussianNormConstBPrime;
          double sigmaB = boysSigma(contractionExponentL, contractionExponentLPrime);
          double UB = boysU(sigmaB);
          double vSquared = boysVSquared(sigmaA, sigmaB);
          double T = boysT(Ra, Rb, vSquared);
          // Boys function setup
          arma::vec byzIntegral = boysIntegralPositionDerivative(Ra, Rb, UA, UB, T, vSquared);
          result += \
            contractionCoefficientK * \
            contractionCoefficientKPrime * \
            contractionCoefficientL * \
            contractionCoefficientLPrime * \
            byzIntegral;
          // std::cout << "d'k=" << contractionCoefficientK << " d'k'=" << contractionCoefficientKPrime << " d'l=" << contractionCoefficientL << " d'l'=" << contractionCoefficientLPrime << std::endl;
          // std::cout << "k=" << k << " k'=" << kPrime << " l=" << l << " l'=" << lPrime << " [o](0)=" << byzIntegral << std::endl;
        }
      }
    }
  }
  // Convert units from a.u. to eV
  result *= electronVoltsToAtomicUnitsConversionFactor;
  return result;
}

double Molecule::fockOperatorMatrixElementOffDiagonal(
  double overlapMatrixElem, 
  double densityMatrixElem, 
  double gammaAB, 
  double betaA,
  double betaB
)
{
  double fockOperatorMatricElement = 0.5 * (betaA + betaB) * overlapMatrixElem - densityMatrixElem * gammaAB;
  return fockOperatorMatricElement;
}

double Molecule::fockOperatorMatrixElementDiagonal(
  int atomIndex,
  double ionizationAndElectronAffinityTerm,
  double pAlphaMuMu,
  arma::mat gammaMatrix,
  arma::vec densityPerAtom
)
{
  double gammaAA = gammaMatrix.at(atomIndex, atomIndex);
  int atomicNumberA = atomicNumbers[atomIndex];
  int zA = valenceAtomicNumberMap[atomicNumberA];
  double pTotAA = densityPerAtom[atomIndex];
  double repulsionSummationTerm = 0.0;
  
  for (int j = 0; j < atomicNumbers.size(); j += 1) {
    if (j != atomIndex) {
      double gammaAB = gammaMatrix.at(atomIndex, j);
      double pTotBB = densityPerAtom[j];
      int atomicNumberB = atomicNumbers[j];
      int zB = valenceAtomicNumberMap[atomicNumberB];
      repulsionSummationTerm += (pTotBB - zB) * gammaAB;
    }
  }
  double permutationTerm = ((pTotAA - zA) - (pAlphaMuMu - 0.5)) * gammaAA;
  double fockOperatorMatrixElement = \
    -1 * ionizationAndElectronAffinityTerm + \
    permutationTerm + \
    repulsionSummationTerm;
  return fockOperatorMatrixElement;
}

double Molecule::hCoreMatrixElementOffDiagonal(
  double overlapMatrixElem, 
  double densityMatrixElem,
  double betaA,
  double betaB
)
{
  double gammaAB = 0.0;
  return fockOperatorMatrixElementOffDiagonal(
    overlapMatrixElem, 
    densityMatrixElem, 
    gammaAB, 
    betaA,
    betaB
  );
}

// double Molecule::hCoreMatrixElementDiagonal(
//   int atomIndex,
//   double ionizationAndElectronAffinityTerm,
//   arma::mat gammaMatrix
// )
// {
//   double pAlphaMuMu = 0.0;
//   double pTotAA = 0.0;
//   double pTotBB = 0.0;
  
//   return fockOperatorMatrixElementDiagonal(
//     int atomIndex,
//     double ionizationAndElectronAffinityTerm,
//     double pAlphaMuMu,
//     gammaMatrix,
//     zeroDensityPerAtom
//   );
// }

double Molecule::atomToAtomEuclideanDistance(arma::vec Ra, arma::vec Rb) {
  double euclideanDistance = sqrt(pow(Ra[0] - Rb[0], 2) + pow(Ra[1] - Rb[1], 2) + pow(Ra[2] - Rb[2], 2));
  return euclideanDistance;
}

double Molecule::boysSigma(double alpha, double alphaPrime) {
  double sigma = pow(alpha + alphaPrime, -1);
  return sigma;
}

double Molecule::boysU(double boysSigma) { 
  double U = pow(sqrt(M_PI * boysSigma), 3);
  // double U = pow(M_PI * boysSigma, 1.5);
  return U;
}

double Molecule::boysVSquared(double boysSigmaA, double boysSigmaB) {
  double vSquared = pow(boysSigmaA + boysSigmaB, -1);
  return vSquared;
}

double Molecule::boysT(arma::vec Ra, arma::vec Rb, double boysVSquared) {
  double atomDistance = atomToAtomEuclideanDistance(Ra, Rb);
  double T = boysVSquared * pow(atomDistance, 2);
  return T;
}

double Molecule::boysIntegral(
        arma::vec Ra,
        arma::vec Rb,
        double boysUA, 
        double boysUB, 
        double boysT,
        double boysVSquared
)
{
  double boysIntegral = 0.0;
  double atomicDistance = atomToAtomEuclideanDistance(Ra, Rb);
  if (boysT == 0.0) {
    boysIntegral = boysUA * boysUB * sqrt(2 * boysVSquared) * sqrt(2 / M_PI);
  } else {
    boysIntegral = boysUA * boysUB * erf(sqrt(boysT)) / atomicDistance;
  }
  return boysIntegral;
}

arma::vec Molecule::boysIntegralPositionDerivative(
        arma::vec Ra,
        arma::vec Rb,
        double boysUA, 
        double boysUB, 
        double boysT,
        double boysVSquared
)
{
  double distance = atomToAtomEuclideanDistance(Ra, Rb);
  arma::vec result(3, arma::fill::zeros);
  if (boysT != 0.0) {
    arma::vec temp = boysUA * boysUB / pow(distance, 2) * ((2 * sqrt(boysVSquared) / sqrt(M_PI)) * exp(-1 * boysT) - erf(sqrt(boysT)) / distance) * (Ra - Rb);
    for (int i = 0; i < 3; i += 1) {
      result[i] = temp[i];
    }
  }
  return result;
}

arma::mat Molecule::gammaMatrix() {
  arma::mat gammaMat(numberOfAtoms, numberOfAtoms);
  for (int i = 0; i < numberOfAtoms; i += 1) { // Todo: remove iteration over basis functions, as it will break for H2. Iterate over the atoms.
    arma::vec Ra = {coordinates.at(i, 0), coordinates.at(i, 1), coordinates.at(i, 2)};
    arma::ivec La(3, arma::fill::zeros);
    int atomicNumberA = atomicNumbers[i];
    int totalLa = 0;
    arma::vec contrCoeffsA = contractionCoefficientMap[atomicNumberA][totalLa];
    arma::vec contrExpsA = contractionExponentMap[atomicNumberA][totalLa];
    for (int j = 0; j < numberOfAtoms; j += 1) {
      arma::ivec Lb(3, arma::fill::zeros);
      int atomicNumberB = atomicNumbers[j];
      int totalLb = 0;
      arma::vec Rb = {coordinates.at(j, 0), coordinates.at(j, 1), coordinates.at(j, 2)};
      arma::vec contrCoeffsB = contractionCoefficientMap[atomicNumberB][totalLb];
      arma::vec contrExpsB = contractionExponentMap[atomicNumberB][totalLb];
      double gamma = gammaTwoCenterTwoElectronRepulsionIntegral(
        contrCoeffsA,
        contrCoeffsB,
        contrExpsA,
        contrExpsB,
        Ra,
        Rb,
        La,
        Lb
      );
      gammaMat.at(i, j) = gamma;
    }
  }
  return gammaMat;
}

arma::mat Molecule::gammaMatrixPositionDerivative() {
  arma::mat gammaMat(3, numberOfAtoms * numberOfAtoms);
  for (int i = 0; i < numberOfAtoms; i += 1) { // Todo: remove iteration over basis functions, as it will break for H2. Iterate over the atoms.
    arma::vec Ra = {coordinates.at(i, 0), coordinates.at(i, 1), coordinates.at(i, 2)};
    arma::ivec La(3, arma::fill::zeros);
    int atomicNumberA = atomicNumbers[i];
    int totalLa = 0;
    arma::vec contrCoeffsA = contractionCoefficientMap[atomicNumberA][totalLa];
    arma::vec contrExpsA = contractionExponentMap[atomicNumberA][totalLa];
    for (int j = 0; j < numberOfAtoms; j += 1) {
      arma::ivec Lb(3, arma::fill::zeros);
      int atomicNumberB = atomicNumbers[j];
      int totalLb = 0;
      arma::vec Rb = {coordinates.at(j, 0), coordinates.at(j, 1), coordinates.at(j, 2)};
      arma::vec contrCoeffsB = contractionCoefficientMap[atomicNumberB][totalLb];
      arma::vec contrExpsB = contractionExponentMap[atomicNumberB][totalLb];
      arma::vec gamma = gammaTwoCenterTwoElectronRepulsionIntegralPositionDerivative(
        contrCoeffsA,
        contrCoeffsB,
        contrExpsA,
        contrExpsB,
        Ra,
        Rb,
        La,
        Lb
      );
      for (int k = 0; k < 3; k += 1) {
        gammaMat.at(k, i * numberOfAtoms + j) = gamma[k];
      }
    }
  }
  return gammaMat;
}
        
std::tuple<arma::mat, arma::mat, arma::mat, arma::mat, arma::mat, arma::mat, arma::vec, arma::vec> Molecule::selfConsistentFieldAlgorithm(
  arma::mat gammaMatrix,
  arma::mat overlapMatrix,
  arma::mat hCoreMatrix,
  int pElectrons,
  int qElectrons,
  bool consoleLog=true
) 
{
  arma::mat fockAlpha;
  arma::mat fockBeta;
  arma::vec epsilonAlpha;
  arma::vec epsilonBeta;
  arma::mat molecularOrbitalCoefficientsAlpha;
  arma::mat molecularOrbitalCoefficientsBeta;
  double totalE = 0, totalEnergyOld = std::numeric_limits<double>::max();
  arma::vec pTot(numberOfAtoms, arma::fill::zeros);
  arma::mat densityMatrixAlphaOld;
  arma::mat densityMatrixBetaOld;
  arma::mat densityMatrixAlpha(numberAtomicBasisFunctions, numberAtomicBasisFunctions, arma::fill::zeros);
  arma::mat densityMatrixBeta(numberAtomicBasisFunctions, numberAtomicBasisFunctions, arma::fill::zeros);
  arma::mat (numberAtomicBasisFunctions, numberAtomicBasisFunctions, arma::fill::zeros);
  double tolerance = 1e-6;
  int iterationCount = 0;

  if (consoleLog) {
    std::cout << "gamma" << std::endl;
    gammaMatrix.print();
    std::cout << "Overlap" << std::endl;
    overlapMatrix.print();
    std::cout << "H_core" << std::endl;
    hCoreMatrix.print();
    std::cout << "p = " << pElectrons << " q = " << qElectrons << std::endl;
  }

  // while (maxDelta > tolerance) {
  while (sqrt(pow(totalE - totalEnergyOld, 2)) > tolerance && iterationCount < 300) {

    fockAlpha = fockOperatorMatrix(gammaMatrix, densityMatrixAlpha, overlapMatrix, pTot);
    fockBeta = fockOperatorMatrix(gammaMatrix, densityMatrixBeta, overlapMatrix, pTot);

    arma::eig_sym(epsilonAlpha, molecularOrbitalCoefficientsAlpha, fockAlpha);
    arma::eig_sym(epsilonBeta, molecularOrbitalCoefficientsBeta, fockBeta);
      
    densityMatrixAlphaOld = densityMatrixAlpha;
    densityMatrixBetaOld = densityMatrixBeta;

    if (pElectrons > 0) {
      densityMatrixAlpha = molecularOrbitalCoefficientsAlpha.cols(0, pElectrons - 1) * molecularOrbitalCoefficientsAlpha.cols(0, pElectrons - 1).t();
    }
    if (qElectrons > 0) {
      densityMatrixBeta = molecularOrbitalCoefficientsBeta.cols(0, qElectrons - 1) * molecularOrbitalCoefficientsBeta.cols(0, qElectrons - 1).t();
    }

    pTot = densityPerAtom(densityMatrixAlpha, densityMatrixBeta);

    if (consoleLog) {
      std::cout << "Iteration: " << iterationCount << std::endl;
    
      std::cout << "Fa" << std::endl;
      fockAlpha.print();
      std::cout << "Fb" << std::endl;
      fockBeta.print();

      std::cout << "Ca" << std::endl;
      molecularOrbitalCoefficientsAlpha.print();
      std::cout << "Cb" << std::endl;
      molecularOrbitalCoefficientsBeta.print();

      std::cout << "Ea" << std::endl;
      epsilonAlpha.print();
      std::cout << "Eb" << std::endl;
      epsilonBeta.print();

      std::cout << "Pa_old" << std::endl;
      densityMatrixAlphaOld.print();
      std::cout << "Pb_old" << std::endl;
      densityMatrixBetaOld.print();
    
      std::cout << "Pa_new" << std::endl;
      densityMatrixAlpha.print();
      std::cout << "Pb_new" << std::endl;
      densityMatrixBeta.print();

      std::cout << "P_t" << std::endl;
      pTot.print();
  
      std::cout << "Ga" << std::endl;
      (fockAlpha - hCoreMatrix).print();
      std::cout << "Gb" << std::endl;
      (fockBeta - hCoreMatrix).print();
    }

    totalEnergyOld = totalE;
    totalE = totalEnergy(
      fockAlpha,
      fockBeta,
      hCoreMatrix,
      hCoreMatrix,
      densityMatrixAlpha,
      densityMatrixBeta
    );

    iterationCount += 1;
  }

  std::tuple<
      arma::mat, arma::mat, 
      arma::mat, arma::mat, 
      arma::mat, arma::mat, 
      arma::vec, arma::vec> result = std::make_tuple(
        fockAlpha, fockBeta,
        densityMatrixAlpha, densityMatrixBeta,
        molecularOrbitalCoefficientsAlpha, molecularOrbitalCoefficientsBeta,
        epsilonAlpha, epsilonBeta
      );
  return result;
}

std::tuple<arma::mat, arma::mat, arma::mat, arma::mat, arma::mat, arma::mat, arma::vec, arma::vec> Molecule::directInversionIterativeSubspaceAlogrithm(
  arma::mat gammaMatrix, 
  arma::mat overlapMatrix,
  arma::mat hCoreMatrix,
  int pElectrons,
  int qElectrons,
  int numPrevIters,  
  bool consoleLog,
  double lambda
) {
  // Initialize data structures
  std::deque<arma::mat> kPrevMolOrbitalCoeffMatricesAlpha;
  std::deque<arma::mat> kPrevMolOrbitalCoeffMatricesBeta;
  std::deque<arma::mat> kPrevDensityMatricesAlpha;
  std::deque<arma::mat> kPrevDensityMatricesBeta;
  std::deque<arma::mat> kPrevFockMatricesAlpha;
  std::deque<arma::mat> kPrevFockMatricesBeta;
  std::deque<arma::mat> kPrevErrorMatricesAlpha;
  std::deque<arma::mat> kPrevErrorMatricesBeta;
  
  arma::vec epsilonAlpha;
  arma::vec epsilonBeta;
  arma::vec pTot(numberOfAtoms, arma::fill::zeros);

  arma::mat molecularOrbitalCoefficientsAlpha;
  arma::mat molecularOrbitalCoefficientsBeta;
  arma::mat densityMatrixAlpha(numberAtomicBasisFunctions, numberAtomicBasisFunctions, arma::fill::zeros);
  arma::mat densityMatrixBeta(numberAtomicBasisFunctions, numberAtomicBasisFunctions, arma::fill::zeros);
  arma::mat fockAlpha = fockOperatorMatrix(gammaMatrix, densityMatrixAlpha, overlapMatrix, pTot);
  arma::mat fockBeta = fockOperatorMatrix(gammaMatrix, densityMatrixBeta, overlapMatrix, pTot);

  arma::mat prevDensityMatrixAlpha(numberAtomicBasisFunctions, numberAtomicBasisFunctions, arma::fill::value(std::numeric_limits<double>::max()));

  double totalE = 0, totalEnergyOld = std::numeric_limits<double>::max();
  double tolerance = 1e-6;
  
  int iterationCount = 0;
  int queueSize = 0;

  // Initialize the queues
  kPrevFockMatricesAlpha.push_front(fockAlpha);
  kPrevFockMatricesBeta.push_front(fockBeta);
  kPrevDensityMatricesAlpha.push_front(densityMatrixAlpha);
  kPrevDensityMatricesBeta.push_front(densityMatrixBeta);
  kPrevErrorMatricesAlpha.push_front(fockAlpha * densityMatrixAlpha - densityMatrixAlpha * fockAlpha);
  kPrevErrorMatricesBeta.push_front(fockBeta * densityMatrixBeta - densityMatrixBeta * fockBeta);
  queueSize += 1;

  // Check for convergence
  while (sqrt(pow(totalE - totalEnergyOld, 2)) > tolerance) {
    
    // Create DIIS equations and solve for extrapolation coefficients
    arma::mat lagrangeMatixAlpha = generateLagrangeMultiplierMatrix(kPrevErrorMatricesAlpha);
    arma::mat lagrangeMatixBeta = generateLagrangeMultiplierMatrix(kPrevErrorMatricesBeta);
    arma::vec lagrangeContantsVectorAlpha(kPrevErrorMatricesBeta.size() + 1, arma::fill::zeros);
    lagrangeContantsVectorAlpha[kPrevErrorMatricesBeta.size()] = -1;
    arma::vec lagrangeContantsVectorBeta(kPrevErrorMatricesBeta.size() + 1, arma::fill::zeros);
    lagrangeContantsVectorBeta[kPrevErrorMatricesBeta.size()] = -1;
    arma::vec lagrangeSolutionAlpha = arma::solve(lagrangeMatixAlpha, lagrangeContantsVectorAlpha);
    arma::vec lagrangeSolutionBeta = arma::solve(lagrangeMatixBeta, lagrangeContantsVectorBeta);
    if (consoleLog) {
      lagrangeMatixAlpha.print("lagrangeMatixAlpha");
      lagrangeMatixBeta.print("lagrangeMatixBeta");
      lagrangeContantsVectorAlpha.print("lagrangeContantsVectorAlpha");
      lagrangeContantsVectorBeta.print("lagrangeContantsVectorBeta");
      lagrangeSolutionAlpha.print("lagrangeSolutionAlpha");
      lagrangeSolutionBeta.print("lagrangeSolutionBeta");
    }

    // Compute new extrapolated Matrices
    // Define extrapolated matrices to go to back of queue
    arma::mat extrapolatedErrorMatrixAlpha(numberAtomicBasisFunctions, numberAtomicBasisFunctions, arma::fill::zeros);
    arma::mat extrapolatedErrorMatrixBeta(numberAtomicBasisFunctions, numberAtomicBasisFunctions, arma::fill::zeros);
    arma::mat extrapolatedFockMatrixAlpha(numberAtomicBasisFunctions, numberAtomicBasisFunctions, arma::fill::zeros);
    arma::mat extrapolatedFockMatrixBeta(numberAtomicBasisFunctions, numberAtomicBasisFunctions, arma::fill::zeros);
    for (int i = 0; i < queueSize; i += 1) {
      arma::mat currErrorMatrixAlpha = kPrevErrorMatricesAlpha[i];
      arma::mat currErrorMatrixBeta = kPrevErrorMatricesBeta[i];
      arma::mat currFockMatrixAlpha = kPrevFockMatricesAlpha[i];
      arma::mat currFockMatrixBeta = kPrevFockMatricesBeta[i];
      extrapolatedErrorMatrixAlpha += lagrangeSolutionAlpha.at(i, 0) * currErrorMatrixAlpha;
      extrapolatedErrorMatrixBeta += lagrangeSolutionBeta.at(i, 0) * currErrorMatrixBeta;
      extrapolatedFockMatrixAlpha += lagrangeSolutionBeta.at(i, 0) * currFockMatrixAlpha;
      extrapolatedFockMatrixBeta += lagrangeSolutionBeta.at(i, 0) * currFockMatrixBeta;
    }

      // Logging information
    if (consoleLog) {
      std::cout << "extrapolated matrices - linear combos of k previous matrices" << std::endl;
      extrapolatedFockMatrixAlpha.print("extrapolatedFockMatrixAlpha");
      extrapolatedFockMatrixBeta.print("extrapolatedFockMatrixBeta");
      extrapolatedErrorMatrixAlpha.print("extrapolatedErrorMatrixAlpha");
      extrapolatedErrorMatrixBeta.print("extrapolatedErrorMatrixBeta");
    }

    // Solve the SCF euqations to create a molecular orbital coefficient matrix
    arma::mat newExtrapMolOrbCoefficientsAlpha;
    arma::mat newExtrapMolOrbCoefficientsBeta;
    arma::eig_sym(epsilonAlpha, newExtrapMolOrbCoefficientsAlpha, extrapolatedFockMatrixAlpha);
    arma::eig_sym(epsilonBeta, newExtrapMolOrbCoefficientsBeta, extrapolatedFockMatrixBeta);

    // Calculate a new density from extrapolated molecular orbital coefficient matrix
    arma::mat newExtrapolatedDensityMatrixAlpha;
    arma::mat newExtrapolatedDensityMatrixBeta;
    if (pElectrons > 0) {
      newExtrapolatedDensityMatrixAlpha = newExtrapMolOrbCoefficientsAlpha.cols(0, pElectrons - 1) * \
                            newExtrapMolOrbCoefficientsAlpha.cols(0, pElectrons - 1).t();
    }
    if (qElectrons > 0) {
      newExtrapolatedDensityMatrixBeta = newExtrapMolOrbCoefficientsBeta.cols(0, qElectrons - 1) * \
                            newExtrapMolOrbCoefficientsBeta.cols(0, qElectrons - 1).t();
    }    
    pTot = densityPerAtom(newExtrapolatedDensityMatrixAlpha, newExtrapolatedDensityMatrixBeta);

    // Calculate new fock matrices
    arma::mat newExtrapolatedFockMatrixAlpha = fockOperatorMatrix(gammaMatrix, newExtrapolatedDensityMatrixAlpha, overlapMatrix, pTot);
    arma::mat newExtrapolatedFockMatrixBeta = fockOperatorMatrix(gammaMatrix, newExtrapolatedDensityMatrixBeta, overlapMatrix, pTot);

    // Using newly calculated fock and density matrices, calculate new error matrices
    arma::mat newExtrapolatedErrorMatrixAlpha = newExtrapolatedFockMatrixAlpha * newExtrapolatedDensityMatrixAlpha - newExtrapolatedDensityMatrixAlpha * newExtrapolatedFockMatrixAlpha;
    arma::mat newExtrapolatedErrorMatrixBeta = newExtrapolatedFockMatrixBeta * newExtrapolatedDensityMatrixBeta - newExtrapolatedDensityMatrixBeta * newExtrapolatedFockMatrixBeta;

    // Logging information
    if (consoleLog) {
      std::cout << "new extrapolated matrices" << std::endl;
      newExtrapolatedDensityMatrixAlpha.print("newExtrapolatedDensityMatrixAlpha " + std::to_string(iterationCount));
      newExtrapolatedDensityMatrixBeta.print("newExtrapolatedDensityMatrixBeta " + std::to_string(iterationCount));
      newExtrapolatedFockMatrixAlpha.print("newExtrapolatedFockMatrixAlpha " + std::to_string(iterationCount));
      newExtrapolatedFockMatrixBeta.print("newExtrapolatedFockMatrixBeta " + std::to_string(iterationCount));
      newExtrapolatedErrorMatrixAlpha.print("newExtrapolatedErrorMatrixAlpha " + std::to_string(iterationCount));
      newExtrapolatedErrorMatrixBeta.print("newExtrapolatedErrorMatrixBeta " + std::to_string(iterationCount));
    }

    // Decide how to manage queue
    kPrevMolOrbitalCoeffMatricesAlpha.push_back(newExtrapMolOrbCoefficientsAlpha);
    kPrevMolOrbitalCoeffMatricesBeta.push_back(newExtrapMolOrbCoefficientsBeta);
    kPrevDensityMatricesAlpha.push_back(newExtrapolatedDensityMatrixAlpha);
    kPrevDensityMatricesBeta.push_back(newExtrapolatedDensityMatrixBeta);
    kPrevErrorMatricesAlpha.push_back(newExtrapolatedErrorMatrixAlpha);
    kPrevErrorMatricesBeta.push_back(newExtrapolatedErrorMatrixBeta);
    kPrevFockMatricesAlpha.push_back(newExtrapolatedFockMatrixAlpha);
    kPrevFockMatricesBeta.push_back(newExtrapolatedFockMatrixBeta);
    queueSize += 1;
    if (queueSize > numPrevIters) {
      kPrevDensityMatricesAlpha.pop_front();
      kPrevDensityMatricesBeta.pop_front();
      kPrevErrorMatricesAlpha.pop_front();
      kPrevErrorMatricesBeta.pop_front();
      kPrevFockMatricesAlpha.pop_front();
      kPrevFockMatricesBeta.pop_front();
      queueSize -= 1;
    }
    
    // Logging information
    if (consoleLog) {
      std::cout << "Density matrices queue" << std::endl;
      for (int i = 0; i < queueSize; i += 1) {
        kPrevDensityMatricesAlpha[i].print("kPrevDensityMatricesAlpha " + std::to_string(i));
      }
      for (int i = 0; i < queueSize; i += 1) {
        kPrevDensityMatricesBeta[i].print("kPrevDensityMatricesBeta " + std::to_string(i));
      }
      std::cout << "Fock matrices queue" << std::endl;
      for (int i = 0; i < queueSize; i += 1) {
        kPrevFockMatricesAlpha[i].print("kPrevFockMatricesAlpha " + std::to_string(i));
      }
      for (int i = 0; i < queueSize; i += 1) {
        kPrevFockMatricesBeta[i].print("kPrevFockMatricesBeta " + std::to_string(i));
      }
      std::cout << "Error matrices queue" << std::endl;
      for (int i = 0; i < queueSize; i += 1) {
        kPrevErrorMatricesAlpha[i].print("kPrevErrorMatricesAlpha " + std::to_string(i));
      }
      for (int i = 0; i < queueSize; i += 1) {
        kPrevErrorMatricesBeta[i].print("kPrevErrorMatricesBeta " + std::to_string(i));
      }
    }

    // Evaluate total energy
    totalEnergyOld = totalE;
    totalE = totalEnergy(
      kPrevFockMatricesAlpha.back(),
      kPrevFockMatricesBeta.back(),
      hCoreMatrix,
      hCoreMatrix,
      kPrevDensityMatricesAlpha.back(),
      kPrevDensityMatricesBeta.back()
    );    

    // Logging information
    if (consoleLog) {
      std::cout << "Iteration count: " << iterationCount << ", queueSize: " << queueSize - 1 << ", totaLEnergy" << totalE << ", totalEnergyOld: " << totalEnergyOld << std::endl;
    }

    iterationCount += 1;
  }

  std::tuple<
      arma::mat, arma::mat, 
      arma::mat, arma::mat, 
      arma::mat, arma::mat, 
      arma::vec, arma::vec> result = std::make_tuple(
        kPrevFockMatricesAlpha.back(), kPrevFockMatricesBeta.back(),
        kPrevDensityMatricesAlpha.back(), kPrevDensityMatricesBeta.back(),
        kPrevMolOrbitalCoeffMatricesAlpha.back(), kPrevMolOrbitalCoeffMatricesBeta.back(),
        epsilonAlpha, epsilonBeta
      );
  return result;
}

arma::mat Molecule::generateLagrangeMultiplierMatrix(std::deque<arma::mat> errorMatrices) {
  arma::mat lagrangeMultMatrix(errorMatrices.size() + 1, errorMatrices.size() + 1, arma::fill::zeros);
  for (int i = 0; i < errorMatrices.size(); i += 1) {
    arma::mat leftErrorMatrix = errorMatrices[i];
    for (int j = 0; j < errorMatrices.size(); j += 1) {
      arma::mat rightErrorMatrix = errorMatrices[j];
      double element = arma::dot(arma::vectorise(leftErrorMatrix), arma::vectorise(rightErrorMatrix));
      lagrangeMultMatrix.at(i, j) = element;
    }
  }
  for (int i = 0; i < errorMatrices.size(); i += 1) {
    lagrangeMultMatrix.at(i, errorMatrices.size()) = -1;
    lagrangeMultMatrix.at(errorMatrices.size(), i) = -1;
  }
  return lagrangeMultMatrix;
}

double Molecule::nuclearRepulsionEnergy() {
  double nuclearRepulsion = 0.0;
  for (int i = 0; i < atomicNumbers.size(); i += 1) {
    int atomicNumberA = atomicNumbers[i];
    int zA = valenceAtomicNumberMap[atomicNumberA];
    arma::vec Ra = {coordinates.at(i, 0), coordinates.at(i, 1), coordinates.at(i, 2)};
    for (int j = i + 1; j < atomicNumbers.size(); j += 1) {
      int atomicNumberB = atomicNumbers[j];
      int zB = valenceAtomicNumberMap[atomicNumberB];
      arma::vec Rb = {coordinates.at(j, 0), coordinates.at(j, 1), coordinates.at(j, 2)};
      double atomicDistance = atomToAtomEuclideanDistance(Ra, Rb);
      nuclearRepulsion += (zA * zB) / atomicDistance * electronVoltsToAtomicUnitsConversionFactor;
    }
  }
  return nuclearRepulsion;
}

arma::mat Molecule::nuclearRepulsionEnergyPositionDerivative() {
  arma::mat nuclearRepulsionEnergyPositionDerivative(3, 0);
  for (int i = 0; i < atomicNumbers.size(); i += 1) {
    int atomicNumberA = atomicNumbers[i];
    int zA = valenceAtomicNumberMap[atomicNumberA];
    arma::vec Ra = {coordinates.at(i, 0), coordinates.at(i, 1), coordinates.at(i, 2)};
    arma::mat ithNuclearRepulsionPositionDerivative(3, 1, arma::fill::zeros);
    for (int j = 0; j < atomicNumbers.size(); j += 1) {
      if (i == j) {
        continue;
      }
      int atomicNumberB = atomicNumbers[j];
      int zB = valenceAtomicNumberMap[atomicNumberB];
      arma::vec Rb = {coordinates.at(j, 0), coordinates.at(j, 1), coordinates.at(j, 2)};
      double atomicDistance = atomToAtomEuclideanDistance(Ra, Rb);
      ithNuclearRepulsionPositionDerivative += -1 * zA * zB * (Ra - Rb) * pow(atomicDistance, -3);
    }
    nuclearRepulsionEnergyPositionDerivative = arma::join_rows(nuclearRepulsionEnergyPositionDerivative, ithNuclearRepulsionPositionDerivative);
  }
  return nuclearRepulsionEnergyPositionDerivative * electronVoltsToAtomicUnitsConversionFactor;
}

double Molecule::electronicEnergy(
  arma::mat fockAlphaMatrix,
  arma::mat fockBetaMatrix,
  arma::mat hCoreAlphaMatrix,
  arma::mat hCoreBetaMatrix,
  arma::mat densityAlphaMatrix,
  arma::mat densityBetaMatrix
) 
{
  double firstTerm = 0.0;
  for (int i = 0; i < numberAtomicBasisFunctions ; i += 1) {
    for (int j = 0; j < numberAtomicBasisFunctions; j += 1) {
      firstTerm += 0.5 * densityAlphaMatrix.at(i, j) * (hCoreAlphaMatrix.at(i, j) + fockAlphaMatrix.at(i, j));
    }
  }
  double secondTerm = 0.0;
  for (int i = 0; i < numberAtomicBasisFunctions; i += 1) {
    for (int j = 0; j < numberAtomicBasisFunctions; j += 1) {
      secondTerm += 0.5 * densityBetaMatrix.at(i, j) * (hCoreBetaMatrix.at(i, j) + fockBetaMatrix.at(i, j));
    }
  }
  double energy = firstTerm + secondTerm;
  return energy;
}

double Molecule::totalEnergy(
  arma::mat fockAlphaMatrix,
  arma::mat fockBetaMatrix,
  arma::mat hCoreAlphaMatrix,
  arma::mat hCoreBetaMatrix,
  arma::mat densityBetaMatrix,
  arma::mat densityAlphaMatrix
) {
  double totalEnergy = electronicEnergy(
    fockAlphaMatrix, 
    fockBetaMatrix, 
    hCoreAlphaMatrix,
    hCoreBetaMatrix,
    densityBetaMatrix,
    densityAlphaMatrix
  ) + nuclearRepulsionEnergy();
  return totalEnergy;
}

arma::vec Molecule::densityPerAtom(arma::mat densityMatrixAlpha, arma::mat densityMatrixBeta) {
  arma::mat totalDensityMatrix = densityMatrixAlpha + densityMatrixBeta;
  arma::vec energyPerAtom(numberOfAtoms, arma::fill::zeros);
  for (int i = 0; i < numberAtomicBasisFunctions; i += 1) {
    arma::irowvec basisFunction = basisFunctions.row(i);
    int atomIndex = basisFunction[0];
    energyPerAtom[atomIndex] += totalDensityMatrix(i, i);
  }
  return energyPerAtom;
}

arma::mat Molecule::xCoefficientMatrix(arma::mat densityMatrixAlpha, arma::mat densityMatrixBeta) {
  arma::mat result = densityMatrixAlpha + densityMatrixBeta;
  for (int i = 0; i < numberAtomicBasisFunctions; i += 1) {
    arma::irowvec basisFunctionA = basisFunctions.row(i);
    int atomicNumberA = basisFunctionA[1];
    for (int j = 0; j < numberAtomicBasisFunctions; j += 1) {
      arma::irowvec basisFunctionB = basisFunctions.row(j);
      int atomicNumberB = basisFunctionB[1];
      int betaA = atomicNumberBondingParameterMap[atomicNumberA], betaB = atomicNumberBondingParameterMap[atomicNumberB];
      result.at(i, j) = (betaA + betaB) * result.at(i, j);
    }
  }
  return result;
}

arma::mat Molecule::yCoefficientMatrix(arma::mat densityMatrixAlpha, arma::mat densityMatrixBeta) {
  arma::mat pTot = densityPerAtom(densityMatrixAlpha, densityMatrixBeta);
  arma::mat result(numberOfAtoms, numberOfAtoms);
  for (int i = 0; i < numberOfAtoms; i += 1) {
    int atomicNumberA = atomicNumbers[i];
    int zA = valenceAtomicNumberMap[atomicNumberA];
    double pTotAA = pTot[i];
    for (int j = 0; j < numberOfAtoms; j += 1) {
      int atomicNumberB = atomicNumbers[j];
      int zB = valenceAtomicNumberMap[atomicNumberB];
      double pTotBB = pTot[j];
      double densityProductsDoubleSum = 0.0;
      for (int k = 0; k < numberAtomicBasisFunctions; k += 1) {
        arma::irowvec basisFunctionA = basisFunctions.row(k);
        int indexAtomA = basisFunctionA[0];
        if (i != indexAtomA) {
          continue;
        }
        for (int l = 0; l < numberAtomicBasisFunctions; l += 1) {
          arma::irowvec basisFunctionB = basisFunctions.row(l);
          int indexAtomB = basisFunctionB[0];
          if (j != indexAtomB) {
            continue;
          }
          densityProductsDoubleSum += (
            densityMatrixAlpha.at(k, l) * densityMatrixAlpha.at(k, l) + \
            densityMatrixBeta.at(k, l) * densityMatrixBeta.at(k, l)
          );
        }
      }
      result.at(i, j) = pTotAA * pTotBB - zB * pTotAA - zA * pTotBB - densityProductsDoubleSum;
    }
  }
  return result;
}

arma::mat Molecule::energyDerivative(
  arma::mat overlapMatrixPositionDerivative,
  arma::mat gammaMatrixPositionDerivative,
  arma::mat xCoefficientMatrix,
  arma::mat yCoefficientMatrix,
  arma::mat nuclearRepulsionEnergyPositionDerivative
) {
  // Calculate overlap term
  arma::mat overlapSummationTerm(3, numberOfAtoms, arma::fill::zeros);
  for (int i = 0; i < numberAtomicBasisFunctions; i += 1) {
    int atomIndex = basisFunctions.at(i, 0);
    for (int j = 0; j < numberAtomicBasisFunctions; j += 1) {
      double coefficientX = xCoefficientMatrix.at(i, j);
      arma::vec overlapDerivative = overlapMatrixPositionDerivative.col(i * numberAtomicBasisFunctions + j);
      overlapDerivative *= coefficientX;
      for (int k = 0; k < 3; k += 1) {
        overlapSummationTerm.at(k, atomIndex) += overlapDerivative[k];
      }
    }
  }

  // Calculate gamma term
  arma::mat gammaSummationTerm(3, numberOfAtoms, arma::fill::zeros);
  for (int i = 0; i < numberOfAtoms; i += 1) {
    for (int j = 0; j < numberOfAtoms; j += 1) {
      if (i == j) {
        continue;
      }
      double coefficientY = yCoefficientMatrix.at(i, j);
      arma::vec gammaDerivative = gammaMatrixPositionDerivative.col(i * numberOfAtoms + j);
      gammaDerivative *= coefficientY;
      for (int k = 0; k < 3; k += 1) {
        gammaSummationTerm.at(k, i) += gammaDerivative[k];
      }
    }
  }

  // Add overlap, gamma term,
  arma::mat result = overlapSummationTerm + gammaSummationTerm + nuclearRepulsionEnergyPositionDerivative;

  return result;
}

arma::mat Molecule::energyDerivative(arma::mat densityAlphaMatrix, arma::mat densityBetaMatrix) {
  arma::mat moleculeOverlapMatrixPositionDerivative = overlapMatrixPositionDerivative();
  arma::mat moleculeGammaMatrixPositionDerivative = gammaMatrixPositionDerivative();
  arma::mat moleculeNuclearRepulsionDerivative = nuclearRepulsionEnergyPositionDerivative();
  arma::mat moleculeXCoefficientMatrix = xCoefficientMatrix(densityAlphaMatrix, densityBetaMatrix);
  arma::mat moleculeYCoefficientMatrix = yCoefficientMatrix(densityAlphaMatrix, densityBetaMatrix);

  arma::mat energyGradient = energyDerivative(
    moleculeOverlapMatrixPositionDerivative,
    moleculeGammaMatrixPositionDerivative,
    moleculeXCoefficientMatrix,
    moleculeYCoefficientMatrix,
    moleculeNuclearRepulsionDerivative
  );

  return energyGradient;
}

double Molecule::totalEnergy() {
  // Compute necessary input matrices for DIIS algorithm
  arma::mat overlapMat = overlapMatrix();
  arma::mat gammaMat = gammaMatrix();
  arma::mat hCoreMat = hCoreMatrix(
    gammaMat,
    overlapMat
  );
  // Run DIIS algorithm
    // std::tuple<arma::mat, arma::mat, 
    //            arma::mat, arma::mat,
    //            arma::mat, arma::mat, 
    //            arma::vec, arma::vec> result = directInversionIterativeSubspaceAlogrithm(
    //             gammaMat, 
    //             overlapMat,
    //             hCoreMat,
    //             numberAlphaElectrons,
    //             numberBetaElectrons,
    //             7,
    //             false,
    //             1.0
    // );
    std::tuple<arma::mat, arma::mat, 
               arma::mat, arma::mat,
               arma::mat, arma::mat, 
               arma::vec, arma::vec> result = selfConsistentFieldAlgorithm(
                gammaMat, 
                overlapMat,
                hCoreMat,
                numberAlphaElectrons,
                numberBetaElectrons,
                false
    );
    arma::mat fockAlphaMatrix = std::get<0>(result);
    arma::mat fockBetaMatrix = std::get<1>(result);
    arma::mat densityAlphaMatrix = std::get<2>(result);
    arma::mat densityBetaMatrix = std::get<3>(result);
    arma::mat orbitalCoefficientsAlpha = std::get<4>(result);
    arma::mat orbitalCoefficientsBeta = std::get<5>(result);

  // Calculate total energy
  double totalE = totalEnergy(
                    fockAlphaMatrix,
                    fockBetaMatrix,
                    hCoreMat,
                    hCoreMat,
                    densityAlphaMatrix,
                    densityBetaMatrix
                  );
  return totalE;
}

arma::mat Molecule::steepestDescentGeometryOptimizer(double stepSize, double tolerance, bool logging=true) {
  // Setup
  double previousEnergy, currentEnergy;
  arma::mat previousCoordinates;
  arma::mat energyGradient;

  std::tuple<arma::mat, arma::mat, arma::mat, arma::mat, arma::mat, arma::mat, arma::vec, arma::vec> diisResult;

  arma::mat overlapMat;
  arma::mat gammaMat;
  arma::mat hCoreMat;
  arma::mat fockAlphaMatrix;
  arma::mat fockBetaMatrix;
  arma::mat densityAlphaMatrix;
  arma::mat densityBetaMatrix;
  arma::mat orbitalCoefficientsAlpha;
  arma::mat orbitalCoefficientsBeta;

  // Initialize optimization data
  previousEnergy = totalEnergy();
  overlapMat = overlapMatrix();
  gammaMat = gammaMatrix();
  hCoreMat = hCoreMatrix(
    gammaMat,
    overlapMat
  );
  diisResult = selfConsistentFieldAlgorithm(
    gammaMat, 
    overlapMat,
    hCoreMat,
    numberAlphaElectrons,
    numberBetaElectrons,
    false
  );
  fockAlphaMatrix = std::get<0>(diisResult);
  fockBetaMatrix = std::get<1>(diisResult);
  densityAlphaMatrix = std::get<2>(diisResult);
  densityBetaMatrix = std::get<3>(diisResult);
  orbitalCoefficientsAlpha = std::get<4>(diisResult);
  orbitalCoefficientsBeta = std::get<5>(diisResult);

  // Evaulate forces on initial cooridantes (CNDO/2 Energy Derivative)
  energyGradient = energyDerivative(densityAlphaMatrix, densityBetaMatrix).t();

  // While the the norm of the gradient is above some tolerance
  int iterationCount = 0;
  // while (arma::norm(energyGradient, "fro") > tolerance && iterationCount < 6) {
  while (arma::norm(energyGradient, "fro") > tolerance) {
    // Calculate new coordinates
    previousCoordinates = coordinates;
    coordinates = previousCoordinates - stepSize * energyGradient / arma::norm(energyGradient, 2);
    currentEnergy = totalEnergy();

    // Logging ingormation
    if (logging) {
      std::cout << "iterationCount: " << iterationCount << ", previousEnergy: " << previousEnergy << ", currentEnergy: " << currentEnergy << ", stepSize: " << stepSize << std::endl;
      previousCoordinates.print("previous coordinates:");
      energyGradient.print("energy gradient:");
      coordinates.print("new coordinates from applied forces:");
    } 
    if (currentEnergy < previousEnergy) { // Productive step, increase step size
      // Update energy to newly accepted lower energy
      previousEnergy = currentEnergy;

      // Calculate a new gradient
      overlapMat = overlapMatrix();
      gammaMat = gammaMatrix();
      hCoreMat = hCoreMatrix(
        gammaMat,
        overlapMat
      );
      diisResult = selfConsistentFieldAlgorithm(
        gammaMat, 
        overlapMat,
        hCoreMat,
        numberAlphaElectrons,
        numberBetaElectrons,
        false
      );
        densityAlphaMatrix = std::get<2>(diisResult);
        densityBetaMatrix = std::get<3>(diisResult);

        energyGradient = energyDerivative(densityAlphaMatrix, densityBetaMatrix).t();

        // Increase setpsize as step was productive
        stepSize *= 1.25;

    } else {  // Counter productive step, decrease step size
      // Restore coordinates to previous state as step was not productive in lower energy
      coordinates = previousCoordinates;
      
      // Descreate step size
      stepSize /= 2;
    }
    iterationCount += 1;
  }
  return coordinates;
}
