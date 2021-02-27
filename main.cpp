#include <cstdlib>
#include "modelpulse.h"
#include "ilobenders.h"

int main(int argc, char* argv[])
{
	IloEnv env;
	ModelPulse* p = ModelPulse::load(env,argv[1]);

	try
	{
		p->relaxAndFix(env, 20, 5);
	}
	catch (const IloException& e)
	{
		std::cerr << e << std::endl;
	}
	delete p;

	//exportLP(env, "pulse.mod", argv[1], "export.LP");
	
	return EXIT_SUCCESS;
}