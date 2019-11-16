#include <iostream>
#include <cmath>
#include <memory>
#include <vector>
#include "Eigen/Dense"
#include "Aletas.hpp"
#include "AletaTask.hpp"
#include "FunDer.hpp"
#include "GeratrizesTask.hpp"

void test_AletaTask();

int main()
{
	test_AletaTask();
	std::getchar();
}

void test_AletaTask()
{
	const int numberOfPoints = 9;
	std::unique_ptr<AletaTask> aleta;
	std::vector<std::shared_ptr<AletaTaskInput>> geratrizes;
	std::string filename;

	// ========== Altere aqui! ==========
#pragma region  Geratrizes
	geratrizes.push_back(std::make_shared<GeratrizA>());
	geratrizes.push_back(std::make_shared<GeratrizB>());
	geratrizes.push_back(std::make_shared<GeratrizC>());
	geratrizes.push_back(std::make_shared<GeratrizD>());
	geratrizes.push_back(std::make_shared<GeratrizE>());
	geratrizes.push_back(std::make_shared<GeratrizF>());
	geratrizes.push_back(std::make_shared<GeratrizG>());
	geratrizes.push_back(std::make_shared<GeratrizH>());
	geratrizes.push_back(std::make_shared<GeratrizI>());
	geratrizes.push_back(std::make_shared<GeratrizJ>());
	geratrizes.push_back(std::make_shared<GeratrizK>());
	geratrizes.push_back(std::make_shared<GeratrizL>());
	geratrizes.push_back(std::make_shared<GeratrizM>());
#pragma endregion

	// ==== não altere mais! ======= !

	// Para cada geratriz
	for (int i = 0; i < geratrizes.size(); i++) {
		aleta = std::make_unique<AletaTask>(numberOfPoints, geratrizes[i]);

		std::cout << "Geratriz " << aleta->ID() << std::endl;
		Eigen::VectorXd T = aleta->getT();

		std::cout << T << std::endl;
		//std::cout << aleta.getA() << std::endl;
		std::cout << "\n " <<   std::endl;

		filename = "geratriz_" + aleta->ID() + ".txt";
		aleta->writeTtoFile(filename);
		aleta.release();
	}
}
