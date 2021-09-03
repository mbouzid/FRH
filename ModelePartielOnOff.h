#pragma once


#include <ilopl/iloopl.h>
#include <vector>
#include <map>

#include "types.h"

typedef IloArray<IloInt> IntArray;
typedef IloArray<IloNumVarArray>  NumVarMatrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<IloIntArray> IntMatrix;
typedef std::vector<IloInt> Orders;



typedef std::map <IloInt, std::map <IloInt, IloNumVar> > VarMatrix;
typedef std::map <IloInt, IloNumVar> VarArray;



class ModelePartielOnOff
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
		// prod times of precOrder
		std::vector<IloInt> _tt;


		/* decision variables */
		// prod 
		VarMatrix _x;
		// setup
		VarMatrix _y;
		// seq
		VarMatrix _u;

		// acceptation
		VarArray _A;

		// completion times
		VarArray _C;
		// tardiness
		VarArray _T;

		// specific to fr heuristic
		VarArray _alpha;
		VarArray _omega;



		static IloOplRunConfiguration loadRC(IloEnv& env, const char* datfile);

		void initVars(IloEnv& env);
		void initObj(IloEnv& env);

		void initConstraints(IloEnv& env);

	public:

		ModelePartielOnOff
		(
			IloEnv& env,
			const IloOplRunConfiguration& dat,
			SETUP setup,
			const IloInt& a,
			const IloInt& b,
			const Orders& nonProcessed,
			const IloInt& precOrder,
			const std::vector<IloInt>& tt
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
			_y(),
			_u(),
			_A(),
			_C(),
			_T(),
			_alpha(),
			_omega()

		{
			std::cout << "initVars" << std::endl;
			initVars(env);
			std::cout << "initObj" << std::endl;
			initObj(env);
			std::cout << "initConstraints" << std::endl;
			initConstraints(env);
		}

		// factory
		static ModelePartielOnOff * load(IloEnv& env, IloOplRunConfiguration& rc, SETUP setup, const IloInt& a, const IloInt& b, const Orders& nonProcessed, const IloInt& precOrder, const std::vector<IloInt>& tt);

		static void relaxAndFix(IloEnv& env, const char* datfile, const IloInt& sigma, const IloInt& delta, SETUP setup);

		static void relaxAndFixLoop(IloOplRunConfiguration& rc, SETUP setup, const IloInt& k, const IloInt& a, const IloInt& b, const IloInt& delta, IloInt& precOrder, Orders& nonProc, std::vector< IloInt>& tt, IntMatrix& vals, IntMatrix & yvals);


		void fix(IloEnv& env, const IntMatrix& xvals, const IntMatrix & yvals, const IloInt& from, const IloInt& to);

		void get(IloCplex& cplx, IntMatrix& vals, IntMatrix & yvals, const IloInt& from, const IloInt& to, std::vector<IloInt>& orders, IloInt& precOrder, std::vector<IloInt>& tt);

		static double computeObj(IloOplRunConfiguration& rc, const IntMatrix& xvals, const IntMatrix & yvals, const IntArray & A, const IntArray & Tardiness, const Orders& sequence);

		void printSol(std::ostream& os, const IntMatrix& vals);
		void printCplx(std::ostream& os, IloCplex & cplx);

		~ModelePartielOnOff()
		{
			IloOplModel opl(_dat.getOplModel());
			IloInt n(opl.getElement("n").asInt());

			for (IloInt i : _nonProcessed)
			{
				_x[i].clear();
				_y[i].clear();
				_u[i].clear();
			}

			_x[_precOrder].clear();
			_y[_precOrder].clear();
			_u[_precOrder].clear();

			_alpha.clear();
			_omega.clear();

		}
};


