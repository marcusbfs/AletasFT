#include <memory>
#include <Array>
#include <string>

#include "AletaTask.hpp"
#include "GeratrizesTask.hpp"

int main(int argc, char** argv)
{
	const int numberOfPoints = argc > 1 ? atoi(argv[1]) : 9;

	std::array<std::unique_ptr<AletaTaskInput>, 13> geratrizes;
	std::string filename;

	// ========== Altere aqui! ==========
	geratrizes[0]  = std::make_unique<GeratrizA>();
	geratrizes[1]  = std::make_unique<GeratrizB>();
	geratrizes[2]  = std::make_unique<GeratrizC>();
	geratrizes[3]  = std::make_unique<GeratrizD>();
	geratrizes[4]  = std::make_unique<GeratrizE>();
	geratrizes[5]  = std::make_unique<GeratrizF>();
	geratrizes[6]  = std::make_unique<GeratrizG>();
	geratrizes[7]  = std::make_unique<GeratrizH>();
	geratrizes[8]  = std::make_unique<GeratrizI>();
	geratrizes[9]  = std::make_unique<GeratrizJ>();
	geratrizes[10] = std::make_unique<GeratrizK>();
	geratrizes[11] = std::make_unique<GeratrizL>();
	geratrizes[12] = std::make_unique<GeratrizM>();

	 //==== nao altere mais! ======= !

	 //Para cada geratriz
	for (const auto& geratriz : geratrizes) {

		AletaTask aleta(numberOfPoints, geratriz.get());

		aleta.build();
		aleta.solve();

		filename = "geratriz_" + aleta.ID() + ".txt";

		aleta.writeDataToFile(filename);
	}
}
