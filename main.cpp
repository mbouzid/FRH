#include <cstdlib>
#include <chrono>
#include "utils.h"

int main(int argc, char* argv[])
{

	auto start = std::chrono::system_clock::now();

	utils::process(argv, argc);

	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	double duration(elapsed_seconds.count());

	std::cout << "time=" << duration << std::endl;

	return EXIT_SUCCESS;
}