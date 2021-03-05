#pragma once


#include "types.h"
#include <ilopl/iloopl.h>
#include <vector>


namespace utils
{

	int process(char ** argv, int argc);

	SETUP parseArg(const char* setup);

	void usage(std::ostream& os, const char* progname);

	std::vector<IloInt> getSetups(const IloInt & j, const std::vector<IloInt> & orders, IloIntMap & s);

	IloInt getMinSetup(const std::vector<IloInt> & setups);
	IloInt getMedSetup(const std::vector<IloInt>& setups);
	IloInt getMaxSetup(const std::vector<IloInt>& setups);

	IloInt computeSetup(SETUP stype, const IloInt& j, const std::vector<IloInt>& orders, IloIntMap& s);
};

