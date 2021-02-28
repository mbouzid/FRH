#pragma once
#include <ilopl/iloopl.h>


typedef IloArray<IloNumVarArray>  NumVarMatrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<IloIntArray> IntMatrix;
#include <vector>

class ModelePartiel
{

	static const char* s_params;

	private:

		IloModel _model;
		IloOplRunConfiguration _dat;

		IloInt _a;
		IloInt _b;


		NumVarMatrix _x;
		NumVarMatrix _u;

		IloNumVarArray _alpha;
		IloNumVarArray _omega;


		static IloOplRunConfiguration loadRC(IloEnv& env, const char* datfile);

		void initVars(IloEnv& env);
		void initObj(IloEnv& env);

		void initConstraints(IloEnv& env);

	public:

		ModelePartiel
		(
			IloEnv& env,
			const IloOplRunConfiguration& dat,
			const IloInt & a,
			const IloInt & b
		) :
			_model(env),
			_dat(dat),
			_a(a),
			_b(b),
			_x(env),
			_u(env),
			_alpha(env),
			_omega(env)

		{
			initVars(env);
			initObj(env);
			initConstraints(env);
		}


		static ModelePartiel* load(IloEnv& env, IloOplRunConfiguration& rc, const IloInt& a, const IloInt& b);

		static void relaxAndFix(IloEnv & env,  const char * datfile, const IloInt& sigma, const IloInt & delta);

		void fix(IloEnv& env, const IntMatrix & vals, const IloInt & from, const IloInt & to);

		void get(IloCplex& cplx, IntMatrix & vals, const IloInt& from, const IloInt& to, std::vector<IloInt> & orders);

		~ModelePartiel()
		{
			IloOplModel opl(_dat.getOplModel());
			IloInt n(opl.getElement("n").asInt());

			_x.end();
			_u.end();
		

			
			_alpha.end();
			
			_omega.end();
			

		}
};

