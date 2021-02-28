#include <cstdlib>
#include <chrono>
#include "modelePartiel.h"

int main(int argc, char* argv[])
{
	IloEnv env;

	env.setDeleter(IloSafeDeleterMode);

	auto start = std::chrono::system_clock::now();

	ModelePartiel::relaxAndFix(env, argv[1], atoi(argv[2]), atoi(argv[3]));
	
	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	double duration(elapsed_seconds.count());

	std::cout << "time=" << duration << std::endl;

	return EXIT_SUCCESS;
}