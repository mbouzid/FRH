#pragma once
#include <ilopl/iloopl.h>
#include <vector>
#include <map>

#include "types.h"

typedef IloArray<IloNumVarArray>  NumVarMatrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<IloIntArray> IntMatrix;
typedef std::vector<IloInt> Orders;



typedef std::map <IloInt, std::map <IloInt, IloNumVar> > VarMatrix;
typedef std::map <IloInt, IloNumVar> VarArray;



class ModelePartiel
{

	static const char* s_params;

	private:

		IloModel _model;
		IloOplRunConfiguration _dat;

		SETUP _setup;


		IloInt _a;
		IloInt _b;

		Orders _nonProcessed;
		IloInt _precOrder;
		// starting time of precOrder
		IloInt _tt;

		VarMatrix _x;
		VarMatrix _u;

		VarArray _alpha;
		VarArray _omega;

		VarArray _CES;




		static IloOplRunConfiguration loadRC(IloEnv& env, const char* datfile);

		void initVars(IloEnv& env);
		void initObj(IloEnv& env);

		void initConstraints(IloEnv& env);

	public:

		ModelePartiel
		(
			IloEnv& env,
			const IloOplRunConfiguration& dat,
			SETUP setup,
			const IloInt & a,
			const IloInt & b   ,
			const Orders & nonProcessed,
			const IloInt & precOrder,
			const IloInt & tt
		) :
			_model(env),
			_dat(dat),
			_setup(setup),
			_a(a),
			_b(b),
			_nonProcessed(nonProcessed),
			_precOrder(precOrder),
			_tt(tt),
			_x(),
			_u(),
			_alpha(),
			_omega(),
			_CES()

		{
			//std::cout << "initVars" << std::endl;
			initVars(env);
			//std::cout << "initObj" << std::endl;
			initObj(env);
			//std::cout << "initConstraints" << std::endl;
			initConstraints(env);
		}


		// get the number of released orders between a and b
		std::vector<IloInt> getNumberOfReleasedOrders(IloInt a, IloInt b, const std::vector <IloInt> & nonProc);

		std::vector<IloInt> getNumberOfDeadlineOrders(IloInt a, IloInt b, const std::vector <IloInt>& nonProc);

		static ModelePartiel* load(IloEnv& env, IloOplRunConfiguration& rc, SETUP setup, const IloInt& a, const IloInt& b, const Orders & nonProcessed, const IloInt & precOrder, const IloInt & tt);

		static void relaxAndFix(IloEnv & env,  const char * datfile, const IloInt& sigma, const IloInt & delta, SETUP setup);

		static void relaxAndFixLoop(IloOplRunConfiguration& rc, SETUP setup, const IloInt& k, const IloInt& a, const IloInt& b, const IloInt& delta, IloInt& prec, Orders& nonProc, IloInt& tt, IntMatrix & vals);


		void fix(IloEnv& env, const IntMatrix & vals, const IloInt & from, const IloInt & to);

		void get(IloCplex& cplx, IntMatrix & vals, const IloInt& from, const IloInt& to, std::vector<IloInt> & orders, IloInt& precOrder, IloInt & tt);

		static double computeObj(IloOplRunConfiguration& opl, const IntMatrix & vals, const Orders & sequence);

		void printSol(std::ostream & os, const IntMatrix & vals);

		~ModelePartiel()
		{
			IloOplModel opl(_dat.getOplModel());
			IloInt n(opl.getElement("n").asInt());

			for (IloInt i : _nonProcessed)
			{
				_x[i].clear();
				_u[i].clear();
			}
	
			_x[_precOrder].clear();
			_u[_precOrder].clear();

			_alpha.clear();
			_omega.clear();

		}
};

