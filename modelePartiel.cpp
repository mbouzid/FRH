#include "modelePartiel.h"


#include <algorithm>

const char* ModelePartiel::s_params = "params.mod";

IloOplRunConfiguration ModelePartiel::loadRC(IloEnv& env, const char* datfile)
{
	IloOplErrorHandler errHdlr(env, std::cerr);
	IloOplSettings settings(env, errHdlr);
	IloOplModelSource src(env, s_params);

	IloOplModelDefinition def(src, settings);
	IloOplDataSource dat(env, datfile);

	IloOplDataElements elts(def, dat);
	IloOplRunConfiguration rc(def, elts);
	rc.getOplModel().generate();

	return rc;
}


ModelePartiel* ModelePartiel::load(IloEnv& env, IloOplRunConfiguration& rc, const IloInt& a, const IloInt& b, const Orders& nonProcessed, const IloInt & precOrder)
{

	IloInt T(rc.getOplModel().getElement("T").asInt());
	IloInt n(rc.getOplModel().getElement("n").asInt());


	return new ModelePartiel(env, rc, a,b,nonProcessed, precOrder);

}

void ModelePartiel::relaxAndFix(IloEnv& env, const char* datfile, const IloInt & sigma, const IloInt & delta)
{
	IloOplRunConfiguration rc(loadRC(env, datfile));

	IloInt T(rc.getOplModel().getElement("T").asInt());
	IloInt n(rc.getOplModel().getElement("n").asInt());


	
	IloInt b(sigma - 1);
	IloInt a(0);
	IloInt k(0);

	Orders nonProc;
	IloInt precOrder(0);

	IntMatrix vals(env, n + 1);
	for (IloInt i(0); i <= n; ++i)
	{
		vals[i] = IloIntArray(env, T + 1);
		for (IloInt t(0); t <= T; ++t)
		{
			vals[i][t] = 0;
		}

		if (i != 0)
		{
			nonProc.push_back(i);
		}
	}
	

	

	double timeBudg(3600);
	
	while (b < T)
	{
		std::cout << "non processed orders: " << std::endl;
		for (IloInt i : nonProc)
		{
			std::cout << i << " ";
		}
		std::cout << std::endl;

		std::cout << "precOrder=" << precOrder << std::endl;

		IloEnv env1;
		std::cout << "b=" << b << std::endl;
		ModelePartiel* subProblem = load(env1, rc, 0, b,nonProc,precOrder);


		if (k > 0)
		{
			IloInt to(std::max(IloInt(0), a - 1));
			std::cout << "fix from " << 0 << " to " << to << std::endl;
			subProblem->fix(env1,vals, 0, to);
		}

	
		IloCplex cplx(env1);
		cplx.extract(subProblem->_model.getImpl());

		//cplx.setOut(env.getNullStream());
		cplx.setParam(IloCplex::Param::Threads, 4);
		cplx.setParam(IloCplex::Param::TimeLimit, 60);

		/*std::ostringstream oss;
		oss << "subModel_" << a << "_" << b << "_pulse" << ".LP";
		std::string fname(oss.str());
		cplx.exportModel(fname.c_str());*/

		if (cplx.solve())
		{

			std::cerr << cplx.getStatus() << std::endl;

			IloInt to(b - delta - 1);
			std::cout << "get from " << a << " to " << to << std::endl;
			subProblem->get(cplx,vals, a, to,nonProc,precOrder);

			for (IloInt i : nonProc)
			{
				std::cout << "omega#" << i << "=" << cplx.getValue(subProblem->_omega[i]) << std::endl;
			}

			for (IloInt i : nonProc)
			{
				std::cout << "alpha#" << i << "=" << cplx.getValue(subProblem->_alpha[i]) << std::endl;
			}	 
			
			std::cout << "partial obj=" << cplx.getObjValue() << std::endl;

		
		}
		else
		{
			std::cerr << cplx.getStatus() << std::endl;


			//printConflict(env, cplx, subProblem);

			break;
		}

		a = b - delta;
		b = b + sigma - delta;
		++k;

		if (b > T)
		{
			b = T;
		}

		//subProblem->_model.end();
		delete subProblem;
		env1.end();


	}

	for (IloInt i(0); i <= n; ++i)
	{
		for (IloInt t(0); t <= T; ++t)
		{
			std::cout << vals[i][t];
		}
		std::cout << std::endl;
	}

}

void ModelePartiel::fix(IloEnv& env, const IntMatrix& vals, const IloInt& from, const IloInt& to)
{
	IloOplModel opl(_dat.getOplModel());
	IloInt n(opl.getElement("n").asInt());


	Orders E(_nonProcessed);
	E.push_back(_precOrder);

	for (IloInt i : E)
	{
		for (IloInt t(from); t <= to; ++t)
		{
			if (vals[i][t] == IloInt(1))
				std::cout << "x#" << i << "#" << t << std::endl;
			_model.add(_x[i][t] == vals[i][t]);
		}
	}


}

void ModelePartiel::get(IloCplex& cplx, IntMatrix& vals, const IloInt& from, const IloInt& to, std::vector <IloInt> & orders, IloInt& precOrder)
{
	IloOplModel opl(_dat.getOplModel());
	IloInt n(opl.getElement("n").asInt());

	Orders E(orders);
	if (_precOrder == 0)
	{
		E.push_back(0);
	}

	IloInt prevOrder(_precOrder);

	for (IloInt t(from); t <= to; ++t)
	{
		for (IloInt i : E)
		{
			IloNumVar x(_x[i][t]);

			vals[i][t] = IloRound(cplx.getValue(_x[i][t]));
			
			if (vals[i][t] == 1 and i != 0)
			{
				orders.erase(std::find(orders.begin(), orders.end(), i));
				precOrder = i;
			}

			if (vals[i][t] and i == 0)
			{
				precOrder = i;
			}
		}
	}

	if (prevOrder != precOrder and prevOrder != 0)
	{
	    if (std::find(orders.begin(), orders.end(), prevOrder) != orders.end())
			orders.erase(std::find(orders.begin(), orders.end(), prevOrder));
	}


}


void ModelePartiel::initVars(IloEnv& env)
{
	IloOplModel opl(_dat.getOplModel());
	IloInt n(opl.getElement("n").asInt());
	IloInt N(_nonProcessed.size() + 1);

	IloInt T(opl.getElement("T").asInt());
	IloIntMap db(opl.getElement("db").asIntMap());
	IloIntMap r(opl.getElement("r").asIntMap());
	IloIntMap p(opl.getElement("p").asIntMap());
	IloIntMap s(opl.getElement("s").asIntMap());



	// E = non processed + previous order
	Orders E(_nonProcessed);
	E.push_back(_precOrder);

	for (IloInt i : E)
	{
		_x.emplace(i, std::map<IloInt, IloNumVar>());
		_u.emplace(i, std::map<IloInt, IloNumVar>());

		for (IloInt t(0); t <= _b; ++t)
		{
			std::ostringstream oss;
			oss << "x#" << i << "#" << t;
			std::string repr(oss.str());

			_x[i].emplace(t,IloNumVar(env, 0, 1, ILOBOOL, repr.c_str()));
		}

		for (IloInt j : E)
		{
			std::ostringstream oss;
			oss << "u#" << i << "#" << j;
			std::string repr(oss.str());
			_u[i].emplace(j,IloNumVar(env, IloFalse, IloTrue, ILOBOOL, repr.c_str()));
		}
	}

	

	// add alpha variables
	for (IloInt j : E)
	{
		std::ostringstream oss;
		oss << "alpha#" << j;
		std::string repr(oss.str());
		_alpha.emplace(j,IloNumVar(env, IloFalse, IloTrue, ILOBOOL, repr.c_str()));

		std::ostringstream oss1;
		oss1 << "omega#" << j;
		std::string repr1(oss1.str());
		_omega.emplace(j,IloNumVar(env, IloFalse, IloTrue, ILOBOOL, repr1.c_str()));
	}


}


void ModelePartiel::initObj(IloEnv& env)
{
	IloOplModel opl(_dat.getOplModel());

	IloInt n(opl.getElement("n").asInt());
	IloInt T(opl.getElement("T").asInt());
	IloIntMap db(opl.getElement("db").asIntMap());
	IloIntMap r(opl.getElement("r").asIntMap());
	IloIntMap p(opl.getElement("p").asIntMap());
	IloNumMap f(opl.getElement("f").asNumMap());
	IloIntMap s(opl.getElement("s").asIntMap());
	IloNumMap c(opl.getElement("c").asNumMap());
	IloIntMap e(opl.getElement("e").asIntMap());


	Orders E(_nonProcessed);
	E.push_back(_precOrder);

	IloExpr obj(env);
	for (IloInt j : _nonProcessed)
	{
		IloExpr profit(env);
		for (IloInt t(r.get(j)); t <= std::min(_b,db.get(j) - p.get(j) + 1); ++t)
		{
			IloNum energyCostProcessing(0.0);
			for (IloInt tp(0); tp < p.get(j); ++tp)
			{
				energyCostProcessing += c.getSub(j).get(t + tp);

			}

			IloExpr	energyCostSetup(env);
			for (IloInt i : E)
			{
				if (i != j)
				{
					IloNum cost(0.0);
					for (IloInt tp(1); tp <= s.getSub(i).get(j); ++tp)
					{
						if (t - tp >= r.get(j))
						{
							cost += c.getSub(j).get(t - tp);
						}

					}

					energyCostSetup += _u[i][j] * cost;
				}
			}


			profit += _x[j][t] * (f.getSub(j).get(t + p.get(j) - 1) - energyCostProcessing - energyCostSetup);

		}

		IloNum range(((db.get(j) - r.get(j)) / (T+p.get(j))));

		profit +=  
			
			//range
			0.8
			* _alpha[j] * /*e.get(j)*/ f.getSub(j).get(std::min(db.get(j) - p.get(j)+1, _b )) ;


		obj += profit;

	}


	_model.add(IloMaximize(env, obj));

}


void ModelePartiel::initConstraints(IloEnv& env)
{
	IloOplModel opl(_dat.getOplModel());

	IloInt n(opl.getElement("n").asInt());
	IloInt T(opl.getElement("T").asInt());
	IloIntMap r(opl.getElement("r").asIntMap());
	IloIntMap s(opl.getElement("s").asIntMap());
	IloIntMap p(opl.getElement("p").asIntMap());
	IloIntMap db(opl.getElement("db").asIntMap());
	IloIntMap d(opl.getElement("d").asIntMap());

	Orders E(_nonProcessed);
	E.push_back(_precOrder);
								   

	// capacity	 


	for (IloInt t(0); t <= _b; ++t)
	{
		IloExpr capacity(env);
		for (IloInt j : E)
		{
			capacity += _x[j][t];
		}
		_model.add(capacity <= 1);
	}


	// capacity 1
	for (IloInt j : E)
	{
		IloExpr sum(env);
		for (IloInt t(r.get(j)); t <= std::min(_b,db.get(j) - p.get(j) + 1); ++t)
		{
			sum += _x[j][t];
		}
		_model.add(sum <= 1);
	}



	// def1	& def2
	for (IloInt j : E)
	{
		IloExpr sum1(env);
		for (IloInt t(0); t < std::min(_b,r.get(j)); ++t)
		{
			sum1 += _x[j][t];
		}
		_model.add(sum1 == 0);

		IloExpr sum2(env);
		for (IloInt t(db.get(j) - p.get(j) + 2); t <= _b; ++t)
		{
			sum1 += _x[j][t];
		}
		_model.add(sum1 == 0);

	}

	if (_precOrder == 0)
	{	 
		_model.add(_x[0][0] == 1);

		IloExpr sumZ(env);

		for (IloInt t(1); t <= _b; ++t)
		{
			sumZ += _x[0][t];
		}

		_model.add(sumZ == 0);
	}


	// prec
	for (IloInt i : E)
	{
		for (IloInt j : _nonProcessed)
		{
			if (i != j)
			{
				IloExpr sumA(env);
				IloInt ub1(std::min(_b,db.get(i) - p.get(i) + 1));
				for (IloInt t(r.get(i)); t <= ub1; ++t)
				{
					sumA += (t)*_x[i][t];
				}
				sumA += ((s.getSub(i).get(j) + p.get(i)) * _u[i][j]) - (db.get(i)) * (1 - _u[i][j])
					- (db.get(i)) * (1 - _omega[j]) - db.get(i) * (1 - _omega[i]) + db.get(i) * (_omega[i] - _omega[j])

					;

				/*if (i == 0)
					sumA += _x[0][0]; */

				IloInt lb(r.get(j) + s.getSub(i).get(j) + 1);


				IloInt ub(std::min(db.get(j) - p.get(j) + 1, _b));

				IloExpr sumB(env);
				for (IloInt t(lb); t <= ub; ++t)
				{
					sumB += t * _x[j][t];
				}

				IloConstraint prec(sumA <= sumB);
				std::ostringstream oss;
				oss << "prec#" << i << "#" << j;
				std::string name(oss.str());

				prec.setName(name.c_str());

				_model.add(prec);



			}
		}
	}


	for (IloInt i : E)
	{
		IloExpr sum1(env);
		for (IloInt j : _nonProcessed)
		{
			if (i != j)
			{
				sum1 += _u[i][j];
			}
		}

		IloExpr sum2(env);
		for (IloInt t(r.get(i)); t <= std::min(_b,db.get(i) - p.get(i) + 1); ++t)
		{
			sum2 += _x[i][t];
		}

		_model.add(sum1 <= _omega[i]);
	}


	// pred
	for (IloInt i : _nonProcessed)
	{
		IloExpr sum1(env);
		for (IloInt j : E)
		{
			if (i != j)
			{
				sum1 += _u[j][i];
			}
		}


		_model.add(sum1 == _omega[i]);
	}

	// approx

	IloExpr sum(env);
	for (IloInt j : E)
	{
		

		// max
		std::vector<IloInt> setups;
		for (IloInt i(0); i <= n; ++i)
		{
			if (i!=j)
				setups.push_back(s.getSub(i).get(j));
		}

		IloInt setup(*std::max_element(setups.begin(),setups.end()));

		

		sum += (p.get(j) + setup )* _alpha[j];
	}

	IloConstraint approx(sum <= T - _b);
	approx.setName("approx");
	_model.add(approx);

	std::cout << "T-b=" << T - _b << std::endl;


	// define omega
	for (IloInt j : E)
	{
		IloExpr subsum2(env);

		for (IloInt t(0); t <= _b; ++t)
		{
			subsum2 += _x[j][t];
		}

		_model.add(_omega[j] >= subsum2);
	}

	// link alpha omega
	for (IloInt j : E)
	{
		_model.add(_alpha[j] <= 1 - _omega[j]);
		_model.add(_omega[j] <= 1 - _alpha[j]);
	}

	_model.add(_omega[_precOrder] == 1);

	// si b + p_j >= db_j - p_j => omega[j] = 0

	/*for (IloInt j(1); j <= n; ++j)
	{
		if (_b + p.get(j) >= db.get(j) - p.get(j) + 1)
		{
			_model.add(_omega[j] == 0);
			_model.add(_alpha[j] == 0);
		}
	}*/
	
	for (IloInt j : E)
	{
		if (_b >= db.get(j))
		{
			_model.add(_alpha[j] == 0);
		}
	}

}
