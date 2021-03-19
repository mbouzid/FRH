#include "modelePartiel.h"

#include "utils.h"
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


ModelePartiel* ModelePartiel::load(IloEnv& env, IloOplRunConfiguration& rc, SETUP setup, const IloInt& a, const IloInt& b, const Orders& nonProcessed, const IloInt & precOrder, const IloInt & tt)
{
	IloInt T(rc.getOplModel().getElement("T").asInt());
	IloInt n(rc.getOplModel().getElement("n").asInt());


	return new ModelePartiel(env, rc,setup, a,b,nonProcessed, precOrder,tt);

}

void ModelePartiel::relaxAndFix(IloEnv& env, const char* datfile, const IloInt & sigma, const IloInt & delta, SETUP setup)
{
	IloOplRunConfiguration rc(loadRC(env, datfile));

	IloInt T(rc.getOplModel().getElement("T").asInt());
	IloInt n(rc.getOplModel().getElement("n").asInt());

 
	//IloNumMap CES(rc.getOplModel().getElement("C").asNumMap());

	//std::cout << CES.getSub(IloInt(0)).getSub(IloInt(2)).get(IloInt(30)) << std::endl;
	
	IloInt b(sigma - 1);
	IloInt a(0);
	IloInt k(0);

	Orders nonProc;

	IloInt precOrder(0);
	IloInt tt(0);

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
		
		relaxAndFixLoop(rc, setup, k, a, b, delta, precOrder, nonProc, tt, vals);

		a = b - delta;
		b = b + sigma - delta;
		++k;

		if (b > T)
		{
			b = T;
		}
	}


	std::cout << "last loop" << std::endl;
	relaxAndFixLoop(rc, setup, k, a, b, delta, precOrder, nonProc, tt, vals);
	


	/*for (IloInt i(0); i <= n; ++i)
	{
		for (IloInt t(0); t <= T; ++t)
		{
			std::cout << vals[i][t];
		}
		std::cout << std::endl;
	}*/

	Orders seq;
	std::map<IloInt, IloInt> u;
	// get sequence dirty
	IloInt previous(0);
	for (IloInt t(0); t <= T; ++t)
	{
		for (IloInt i(1); i <= n; ++i)
		{
			if (vals[i][t])
			{
				u.emplace(previous, i);
				previous = i;
				seq.push_back(i);
			}
		}
	}

	IloOplModel opl(rc.getOplModel());
	IloIntMap p(opl.getElement("p").asIntMap());
	IloIntMap s(opl.getElement("s").asIntMap());
	for (IloInt i(0); i <= n; ++i)
	{
		for (IloInt t(0); t <= T; ++t)
		{
			if (vals[i][t] == 1)
			{
				if (u.find(i) != u.end())
				{
					IloInt prev(u.at(i));
					//std::cout << "ST[" << i << "]==" << t - s.getSub(prev).get(i) -1 << ";" << std::endl;
					std::cout << "C[" << i << "]==" << t + p.get(i) -1 << ";" << std::endl;
				}
			}
		}

	}



	for (const std::pair<IloInt,IloInt> & p : u)
	{
		std::cout << "u[" << p.first << "][" << p.second << "]==" << 1 << ";" << std::endl;
	}
	std::cout << std::endl;

	std::cout << "obj=" << computeObj(rc,vals, seq) << std::endl;

}

void ModelePartiel::relaxAndFixLoop(IloOplRunConfiguration& rc, SETUP setup, const IloInt& k, const IloInt& a, const IloInt& b, const IloInt& delta, IloInt & precOrder, Orders & nonProc, IloInt & tt, IntMatrix& vals)
{
	std::cout << "non processed orders: " << std::endl;
	for (IloInt i : nonProc)
	{
		std::cout << i << " ";
	}
	std::cout << std::endl;

	std::cout << "precOrder=" << precOrder << std::endl;

	IloEnv env1;
	std::cout << "a=" << a << "b=" << b << std::endl;
	ModelePartiel* subProblem = load(env1, rc, setup, a, b, nonProc, precOrder, tt);


	if (k > 0)
	{
		IloInt to(std::max(IloInt(0), a - 1));

		std::cout << "fix from " << tt << " to " << to << std::endl;
		subProblem->fix(env1, vals, tt, to);
	}


	IloCplex cplx(env1);
	cplx.extract(subProblem->_model.getImpl());

	cplx.setOut(env1.getNullStream());
	cplx.setParam(IloCplex::Param::Threads, 4);
	cplx.setParam(IloCplex::Param::TimeLimit, 15);

	/*std::ostringstream oss;
	oss << "subModel_" << a << "_" << b << "_pulse" << ".LP";
	std::string fname(oss.str());
	cplx.exportModel(fname.c_str());*/

	if (cplx.solve())
	{

		std::cout << cplx.getStatus() << std::endl;

		IloInt to(b - delta - 1);
		std::cout << "get from " << a << " to " << to << std::endl;
		subProblem->get(cplx, vals, a, to, nonProc, precOrder, tt);

		for (IloInt i : nonProc)
		{
			std::cout << "omega#" << i << "=" << cplx.getValue(subProblem->_omega[i]) << std::endl;
		}

		for (IloInt i : nonProc)
		{
			std::cout << "alpha#" << i << "=" << cplx.getValue(subProblem->_alpha[i]) << std::endl;
		}

		std::cout << "partial obj=" << cplx.getObjValue() << std::endl;
		std::cout << "tt=" << tt << std::endl;

	}
	else
	{
		std::cerr << cplx.getStatus() << std::endl;


		//printConflict(env, cplx, subProblem);

		return;
	}


	delete subProblem;
	cplx.end();
	env1.end();

}


void ModelePartiel::fix(IloEnv& env, const IntMatrix& vals, const IloInt& from, const IloInt& to)
{
	IloOplModel opl(_dat.getOplModel());
	IloInt n(opl.getElement("n").asInt());


	/*for (IloInt i : _nonProcessed)
	{
		for (IloInt t(from); t <= to; ++t)
		{
			if (vals[i][t] == IloInt(1))
				std::cout << "x#" << i << "#" << t << std::endl;
			_model.add(_x[i][t] == vals[i][t]);
		}
	}*/

	_model.add(_x[_precOrder][_tt] == 1);

}

void ModelePartiel::get(IloCplex& cplx, IntMatrix& vals, const IloInt& from, const IloInt& to, std::vector <IloInt> & orders, IloInt& precOrder, IloInt & tt)
{
	IloOplModel opl(_dat.getOplModel());
	IloInt n(opl.getElement("n").asInt());

	Orders E(orders);

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
				tt = t;
			}

			if (vals[i][t] and i == 0)
			{
				precOrder = i;
				tt = 0;
			}
		}
	}

	if (precOrder == 0)
	{
		vals[0][0] = 1;
	}

	/*if (prevOrder != precOrder and prevOrder != 0)
	{
	    if (std::find(orders.begin(), orders.end(), prevOrder) != orders.end())
			orders.erase(std::find(orders.begin(), orders.end(), prevOrder));
	}*/


}

double ModelePartiel::computeObj(IloOplRunConfiguration& rc, const IntMatrix& vals, const Orders& sequence)
{

	IloOplModel opl(rc.getOplModel());

	IloInt n(opl.getElement("n").asInt());
	IloInt T(opl.getElement("T").asInt());
	IloIntMap db(opl.getElement("db").asIntMap());
	IloIntMap r(opl.getElement("r").asIntMap());
	IloIntMap p(opl.getElement("p").asIntMap());
	IloNumMap f(opl.getElement("f").asNumMap());
	IloIntMap s(opl.getElement("s").asIntMap());
	IloNumMap c(opl.getElement("c").asNumMap());
	IloIntMap e(opl.getElement("e").asIntMap());


	IloNum obj(0.0);
	IloInt prev(0);

	for (IloInt j : sequence)
	{
		IloNum profit(0.0);
		for (IloInt t(r.get(j)); t <=  db.get(j) - p.get(j) + 1; ++t)
		{
			IloNum energyCostProcessing(0.0);
			for (IloInt tp(0); tp < p.get(j); ++tp)
			{
				energyCostProcessing += c.getSub(j).get(t + tp);

			}

			IloNum	energyCostSetup(0.0);
			for (IloInt tp(1); tp <= s.getSub(prev).get(j); ++tp)
			{
				if (t - tp >= r.get(j))
				{
					energyCostSetup += c.getSub(j).get(t - tp);
				}
			}

			profit += vals[j][t] * (f.getSub(j).get(t + p.get(j) - 1) - energyCostProcessing - energyCostSetup);
		}

		prev = j;

		obj += profit;

	}



	return obj;
}

void ModelePartiel::printSol(std::ostream& os, const IntMatrix& vals)
{



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



	for (IloInt i : _nonProcessed)
	{
		_x.emplace(i, std::map<IloInt, IloNumVar>());
		_u.emplace(i, std::map<IloInt, IloNumVar>());

		IloInt lb(std::max(_a, _tt));
		for (IloInt t(lb); t <= _b; ++t)
		{
			std::ostringstream oss;
			oss << "x#" << i << "#" << t;
			std::string repr(oss.str());

			_x[i].emplace(t,IloNumVar(env, 0, 1, ILOBOOL, repr.c_str()));
		}

		for (IloInt j : _nonProcessed)
		{
			std::ostringstream oss;
			oss << "u#" << i << "#" << j;
			std::string repr(oss.str());
			_u[i].emplace(j,IloNumVar(env, IloFalse, IloTrue, ILOBOOL, repr.c_str()));
		}
	}


	Orders E(_nonProcessed);
	E.push_back(_precOrder);
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

		std::ostringstream oss2;
		oss2 << "CES#" << j;
		std::string repr2(oss2.str());
		_CES.emplace(j, IloNumVar(env, 0, IloInfinity, ILOFLOAT, repr1.c_str()));

	}


	// prec Order
	_x.emplace(_precOrder, std::map<IloInt, IloNumVar>());
	std::ostringstream ossPrec;
	ossPrec << "x#" << _precOrder << "#" << _tt;
	std::string reprPrec(ossPrec.str());

	// x[_precOrder][_tt]
	_x[_precOrder].emplace(_tt, IloNumVar(env, 1, 1, ILOBOOL, reprPrec.c_str()));

	// u[_precOrder][j]
	for (IloInt j : _nonProcessed)
	{
		std::ostringstream oss;
		oss << "u#" << _precOrder << "#" << j;
		std::string repr(oss.str());
		_u[_precOrder].emplace(j, IloNumVar(env, IloFalse, IloTrue, ILOBOOL, repr.c_str()));
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
	IloNumMap w(opl.getElement("w").asNumMap());
	IloIntMap d(opl.getElement("d").asIntMap());


	Orders E(_nonProcessed);
	E.push_back(_precOrder);

	// compute pavg
	IloNum psum(0.0);
	for (IloInt i : E)
	{
		psum += p.get(i);
	}
	IloNum pavg(psum / (IloNum)E.size());

	
	// compute savg
	IloNum ssum(0.0);
	for (IloInt i : E)
	{
		for (IloInt j : _nonProcessed)
		{
			if (i != j)
			{
				ssum += s.getSub(i).get(j);
			}
		}
	}

	IloNum savg(ssum /(IloNum) (_nonProcessed.size()*E.size()));

	



	IloExpr obj(env);
	for (IloInt j : _nonProcessed)
	{
		IloInt lb(std::max(_tt, _a));
		IloExpr profit(env);
		for (IloInt t(std::max(lb,r.get(j))); t <= std::min(_b,db.get(j) - p.get(j) + 1); ++t)
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


			profit += _x[j][t] * (f.getSub(j).get(t + p.get(j) - 1) - energyCostProcessing 
				//- energyCostSetup
				);

		}

		IloNum rangeTito2021(((db.get(j) - r.get(j)) / (T+p.get(j))));

	/*	IloInt setup(utils::computeSetup(_setup, j, _nonProcessed, s));
		IloNum coeff((w.get(j) / p.get(j)) * std::exp(-((std::max(db.get(j) - p.get(j) - _b,  _a) / (pavg * 0.9))))
			* std::exp(-((setup) / (0.9*savg))));
	 */
		profit +=
			
			/*_alpha[j] * (1/coeff);
			*/ 
			
			/*0.8
			* _alpha[j] * f.getSub(j).get(std::min(db.get(j) - p.get(j)+1, _b )) ;	 */

					
			(0.8) * _alpha[j] * f.getSub(j).get(std::min(db.get(j) - p.get(j) + 1, _b));
		

		obj += profit - _CES[j];

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


	long int nbConstraints(0);

	// capacity	 

	//std::cout << "capacity" << std::endl;
	for (IloInt t(std::max(_a, _tt)); t <= _b; ++t)
	{
		IloExpr capacity(env);
		for (IloInt j : _nonProcessed)
		{
			capacity += _x[j][t];

		}
		_model.add(capacity <= 1);
		++nbConstraints;
	}



	//std::cout << "capacity1" << std::endl;
	// capacity 1
	for (IloInt j : _nonProcessed)
	{
		IloExpr sum(env);
		IloInt lb(std::max(_a, _tt));
		for (IloInt t(std::max(lb, r.get(j))); t <= std::min(_b, db.get(j) - p.get(j) + 1); ++t)
		{
			sum += _x[j][t];
		}
		_model.add(sum <= 1);
		++nbConstraints;
	}




	//std::cout << "def1 & def2" << std::endl;
	// def1	& def2
	for (IloInt j : _nonProcessed)
	{
		IloExpr sum1(env);
		IloInt lb(std::max(_a, _tt));
		for (IloInt t(lb); t < std::min(_b, r.get(j)); ++t)
		{
			sum1 += _x[j][t];
		}
		_model.add(sum1 == 0);
		++nbConstraints;

		IloExpr sum2(env);
		for (IloInt t(std::max(db.get(j) - p.get(j) + 2, lb)); t <= _b; ++t)
		{
			sum1 += _x[j][t];
		}
		_model.add(sum1 == 0);
		++nbConstraints;

	}

	IloInt LB(std::max(_a, _tt));
	_model.add(_x[_precOrder][_tt] == 1);
	++nbConstraints;


	for (IloInt i : _nonProcessed)
	{
		IloExpr sumI(env);
		for (IloInt t(LB); t <= std::min(_b, _tt + p.get(_precOrder) + s.getSub(_precOrder).get(i) - 1); ++t)
		{
			sumI += _x[i][t];
		}
		_model.add(sumI == 0);
		++nbConstraints;
	}


	//}


	// prec
	for (IloInt i : E)
	{
		for (IloInt j : _nonProcessed)
		{
			if (i != j)
			{
				IloExpr sumA(env);
				IloInt lb1(std::max(LB, r.get(i)));
				IloInt ub1(std::min(_b, db.get(i) - p.get(i) + 1));
				if (i == _precOrder)
				{
					sumA += _tt * _x[_precOrder][_tt];
				}
				else
				{
					for (IloInt t(lb1); t <= ub1; ++t)
					{
						sumA += (t)*_x[i][t];
					}
				}
				sumA += ((s.getSub(i).get(j) + p.get(i)) * _u[i][j]) - (db.get(i)) * (1 - _u[i][j])
					- (db.get(i)) * (1 - _omega[j]) - db.get(i) * (1 - _omega[i]) + db.get(i) * (_omega[i] - _omega[j])

					;

				/*if (i == 0)
					sumA += _x[0][0]; */

				IloInt lb(std::max(LB, r.get(j) + s.getSub(i).get(j) + 1));


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
				++nbConstraints;


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


		_model.add(sum1 <= _omega[i]);
		++nbConstraints;
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
		++nbConstraints;
	}

	// approx

	IloExpr sum(env);
	for (IloInt j : _nonProcessed)
	{

		IloInt setup(utils::computeSetup(_setup,j,E,s));



		sum += (p.get(j) + setup) * _alpha[j];
	}

	IloConstraint approx(sum <= T - _b);
	approx.setName("approx");
	_model.add(approx);
	++nbConstraints;

	std::cout << "T-b=" << T - _b << std::endl;


	// define omega
	for (IloInt j : _nonProcessed)
	{
		IloExpr subsum2(env);

		for (IloInt t(LB); t <= _b; ++t)
		{
			subsum2 += _x[j][t];
		}

		_model.add(_omega[j] >= subsum2);
		++nbConstraints;
	}

	// link alpha omega
	for (IloInt j : _nonProcessed)
	{
		_model.add(_alpha[j] <= 1 - _omega[j]);
		++nbConstraints;
		_model.add(_omega[j] <= 1 - _alpha[j]);
		++nbConstraints;
	}

	_model.add(_alpha[_precOrder] == 0);
	++nbConstraints;
	_model.add(_omega[_precOrder] == 1);
	++nbConstraints;

	// si b + p_j >= db_j - p_j => omega[j] = 0

	/*for (IloInt j(1); j <= n; ++j)
	{
		if (_b + p.get(j) >= db.get(j) - p.get(j) + 1)
		{
			_model.add(_omega[j] == 0);
			_model.add(_alpha[j] == 0);
		}
	}*/

	for (IloInt j : _nonProcessed)
	{
		if (_b >= db.get(j))
		{
			_model.add(_alpha[j] == 0);
			++nbConstraints;
		}
	}

	std::map<IloInt, std::map<IloInt, std::map<IloInt, IloNum>>> CES;

	IloNumMap c(opl.getElement("c").asNumMap());
	for (IloInt i : E)
	{
		CES.emplace(i, std::map<IloInt, std::map<IloInt, IloNum>>());
		for (IloInt j : _nonProcessed)
		{
			if (i != j)
			{
				CES[i].emplace(j, std::map<IloInt, IloNum>());
				for (IloInt t(LB); t <= _b; ++t)
				{
					if ((t >= r.get(j) + s.getSub(i).get(j) +1) and (t <= db.get(j) - p.get(j) + 1))
					{
						IloNum sumCout(0.0);
						for (IloInt tp(t - s.getSub(i).get(j)); tp <= t - 1; ++tp)
						{

							sumCout += c.getSub(j).get(tp);
						}
						CES[i][j].emplace(t, sumCout);
					}
				}
			}
		}
	}

	for (IloInt j : _nonProcessed)
	{
		for (IloInt t(LB); t <= _b; ++t)
		{
			IloExpr sumCost(env);
			for (IloInt i : E)
			{
				sumCost += _x[j][t] * CES[i][j][t] - CES[i][j][t] * (1 - _u[i][j])
					;

			}
			//sumCost -= 100 * (1 - _omega[j]);
			_model.add(_CES[j] >= sumCost);
			++nbConstraints;
		}
	}

	std::cout << "NB CONSTRAINTS = " << nbConstraints << std::endl;

	for (IloInt i : E)
	{
		for (IloInt j : _nonProcessed)
		{
			if (i != j)
			{
				CES[i][j].clear();
			}

		}
		CES[i].clear();
	}

	CES.clear();

}