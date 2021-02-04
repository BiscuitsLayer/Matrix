#pragma once

//  SYSTEM
#include <random>
#include <fstream>

//  MATRIX
#include "../Matrix/Matrix.hpp"

const int DEFAULT_SIZE = 50;
const double UNIFORM_MIN = -0.5;
const double UNIFORM_MAX = 0.5;

class Generator {
    private:
        std::mt19937 generator_ {};
        std::uniform_real_distribution <> uniformDistribution_ { UNIFORM_MIN, UNIFORM_MAX };

        std::ofstream* OpenFile () {
            static int testIdx = 0;
            std::stringstream filenameStream {};
            filenameStream << "Generated/" << testIdx++;
            std::ofstream* file  = new std::ofstream { filenameStream.str () };
            return file;
        }

        Linear::Matrix <double> GenerateUpperTriangular (int size) {
            Linear::Matrix <double> ans { size };
            for (int i = 0; i < size; ++i) {
                for (int j = i; j < size; ++j) {
                    ans.At (i, j) = uniformDistribution_ (generator_);
                }
            }
            return ans;
        }

        bool UpperTriangularTest (int size = DEFAULT_SIZE) {
            Linear::Matrix <double> m1 = GenerateUpperTriangular (size);
            Linear::Matrix <double> m2 = GenerateUpperTriangular (size);
            m2.Transpose ();

            double correctAns = 1.0;
            for (int i = 0; i < size; ++i) {
                correctAns *= m1.At (i, i);
                correctAns *= m2.At (i, i);
            }

            m1 *= m2;
            double myAns = m1.Determinant (Linear::Determinant::Type::GAUSS);
            bool result = std::fabs (myAns - correctAns) < EPS;
            //if (!result) {
                std::ofstream* outfile = OpenFile ();
                *outfile << m1 << std::endl;
                *outfile << "Correct Answer = " << correctAns << ", My Answer = " << myAns << std::endl;
                *outfile << "Difference = " << std::fabs (myAns - correctAns) << std::endl;
                delete outfile;
            //}
            return result;
        }

        Linear::Matrix <double> GenerateGivenDeterminant (double determinant, int size) {
            Linear::Matrix <double> m1 = GenerateUpperTriangular (size);
            Linear::Matrix <double> m2 = GenerateUpperTriangular (size);
            m2.Transpose ();

            for (int i = 0; i < size; ++i) {
                m1.At (i, i) = m2.At (i, i) = 1;
            }

            int row1 = std::abs (static_cast <int> (uniformDistribution_ (generator_))) % size,
                row2 = std::abs (static_cast <int> (uniformDistribution_ (generator_))) % size;
            m1.At (row1, row1) = m2.At (row2, row2) = std::sqrt (determinant);
            return m1 * m2;
        }

        bool GivenDeterminantTest (double determinant, int size = DEFAULT_SIZE) {
            Linear::Matrix <double> m1 = GenerateGivenDeterminant (determinant, size);
            double myAns = m1.Determinant (Linear::Determinant::Type::GAUSS);
            bool result = std::fabs (myAns - determinant) < EPS;
            //if (!result) {
                std::ofstream* outfile = OpenFile ();
                *outfile << m1 << std::endl;
                *outfile << "Correct Answer = " << determinant << ", My Answer = " << myAns << std::endl;
                *outfile << "Difference = " << std::fabs (myAns - determinant) << std::endl;
                delete outfile;
            //}
            return result;
        }

        Linear::Matrix <double> GenerateSingular (int size) {
            Linear::Matrix <double> ans { size };
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    ans.At (i, j) = uniformDistribution_ (generator_);
                }
            }

            int row1 = std::abs (static_cast <int> (uniformDistribution_ (generator_))) % size,
                row2 = std::abs (static_cast <int> (uniformDistribution_ (generator_))) % size;
            if (row1 == row2) {
                row2 = (row2 + 1) % size;
            }
            for (int i = 0; i < size; ++i) {
                ans.At (row1, i) = ans.At (row2, i);
            }
            return ans;
        }

        bool SingularTest (int size = DEFAULT_SIZE) {
            Linear::Matrix <double> m1 = GenerateSingular (size);
            //std::cout << m1 << std::endl;
            double myAns = m1.Determinant (Linear::Determinant::Type::GAUSS);
            bool result = std::fabs (myAns) < EPS;
            //if (!result) {
                std::ofstream* outfile = OpenFile ();
                *outfile << m1 << std::endl;
                *outfile << "Correct Answer = " << 0 << ", My Answer = " << myAns << std::endl;
                *outfile << "Difference = " << std::fabs (myAns) << std::endl;
                delete outfile;
            //}
            return result;
        }

    public:
        void Execute () {
            std::cerr << "----------------------------------" << std::endl;
            std::cerr << "UPPER TRIANGULAR TESTS" << std::endl;
            std::cerr << "----------------------------------" << std::endl;
            for (int i = 0; i < 10; ++i) {
                std::cout << std::boolalpha << UpperTriangularTest () << std::endl;
            }
            std::cerr << "----------------------------------" << std::endl;
            std::cerr << "SINGULAR TESTS" << std::endl;
            std::cerr << "----------------------------------" << std::endl;
            for (int i = 0; i < 10; ++i) {
                std::cout << std::boolalpha << SingularTest () << std::endl;
            }
            std::cerr << "----------------------------------" << std::endl;
            std::cerr << "GIVEN DETERMINANT TESTS" << std::endl;
            std::cerr << "----------------------------------" << std::endl;
            for (int i = 0; i < 10; ++i) {
                std::cout << std::boolalpha << GivenDeterminantTest (std::fabs (uniformDistribution_ (generator_))) << std::endl;
            }
            std::cout << std::boolalpha << GivenDeterminantTest (42) << std::endl;
        }
};