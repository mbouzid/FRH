#include "utils.h"
#include <algorithm>
#include "modelePartiel.h"
#include "ModelePartielOnOff.h"

int utils::process(char** argv, int argc)
{

	IloEnv env;

	env.setDeleter(IloSafeDeleterMode);

	try
	{

		if (argc < 5)
		{
			usage(std::cerr, argv[0]);
			return EXIT_FAILURE;
		}


		const char* datfile(argv[1]);
		IloInt sigma(atoi(argv[2]));
		IloInt delta(atoi(argv[3]));
		SETUP stype(parseArg(argv[4]));

		ModelePartiel::relaxAndFix(env, datfile, sigma, delta, stype);
		//ModelePartielOnOff::relaxAndFix(env, datfile,sigma, delta,stype);
	}
	catch (const IloException& e)
	{
		std::cout << e << std::endl;

	}

}

SETUP utils::parseArg(const char* setup)
{
	if (strcmp(setup, "MIN") == 0)
	{
		std::cout << "MIN" << std::endl;
		return SETUP::MIN;
	}
	else
	{
		if (strcmp(setup, "MED") == 0)
		{
			std::cout << "MED" << std::endl;
			return SETUP::MED;
		}
		else
		{
			if (strcmp(setup, "MAX") == 0)
			{
				std::cout << "MAX" << std::endl;
				return SETUP::MAX;
			}
			else
			{
				return SETUP::INVALID;
			}
		}

	}
}

void utils::usage(std::ostream & os, const char* progname)
{
	os << progname << " <datfile> <sigma> <delta> <setup> " << std::endl;
	os << "<datfile>:" << "string | .dat filename" << std::endl;
	os << "<sigma>:" << "int | observation window length " << std::endl;
	os << "<delta>:" << "int | number of overlapping periods between two steps" << std::endl;
	os << "<setup>:"<< "string : MIN,MED,MAX | setup approximation strategy" << std::endl;
}

std::vector<IloInt> utils::getSetups(const IloInt& j, const std::vector<IloInt>& orders, IloIntMap& s)
{
	std::vector<IloInt> setups;
	for (IloInt i : orders)
	{
		if (i != j)
		{
			setups.push_back(s.getSub(i).get(j));
		}
	}


	return setups;
}

IloInt utils::getMinSetup(const std::vector<IloInt>& setups)
{
	return *std::min_element(setups.cbegin(), setups.cend());
}

IloInt utils::getMedSetup(const std::vector<IloInt>& setups)
{
	
	std::vector<IloInt> setups_sorted(setups);
	std::sort(setups_sorted.begin(), setups_sorted.end());
	
	size_t nb(setups.size());
	if (nb % 2 == 0)
	{
		IloNum res((setups_sorted.at(nb / 2) + setups_sorted.at((nb+1) / 2))/2 );
		
		return IloRound(res);
	}
	else
	{
		IloInt res(setups_sorted.at(nb / 2));
		return res;
	}
}

IloInt utils::getMaxSetup(const std::vector<IloInt>& setups)
{
	return *std::max_element(setups.cbegin(), setups.cend());
}

IloInt utils::computeSetup(SETUP stype, const IloInt& j, const std::vector<IloInt>& orders, IloIntMap& s)
{
	IloInt setup(0);
	std::vector<IloInt> setups(getSetups(j, orders, s));


	switch (stype)
	{
		case SETUP::MIN:
		{
			setup = getMinSetup(setups);
			break;
		}

		case SETUP::MED:
		{
			setup = getMedSetup(setups);
			break;
		}

		case SETUP::MAX:
		{
			setup = getMaxSetup(setups);
			break;
		}

		case SETUP::INVALID:
		{
			break;
		}

	}

	return setup;
}
