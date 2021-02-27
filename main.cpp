#include <cstdlib>
#include "modelpulse.h"

int main(int argc, char* argv[])
{
	IloEnv env;
	ModelPulse* p = ModelPulse::load(env,argv[1]);

	p->relaxAndFix(env, 20, 10);

	delete p;
	return EXIT_SUCCESS;
}