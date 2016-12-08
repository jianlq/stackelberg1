#ifndef LB_H
#define LB_H
#include"CGraph.h"
#include <ilcplex/ilocplex.h>

// 4+x^2
void LBdictor(CGraph *G,vector<demand> & req,int ornum){
	G->clearOcc();
	IloEnv env;
	IloModel model(env);
	IloCplex EEsolver(model);

	int num = req.size();	
	IloArray<IloIntVarArray> x(env, num); 
	for(int d = 0; d < num; d++)
		x[d] = IloIntVarArray(env, G->m, 0, 1); 	

	//优化目标
	IloNumVar z(env, 0, 1);	
	model.add(IloMinimize(env, z));

	// 对每个点进行流量守恒约束  
	for(int d = 0; d < num; d++){
		for(int i = 0; i < G->n; i++){    // n为顶点数
			IloExpr constraint(env);
			for(unsigned int k = 0; k < G->adjL[i].size(); k++) // 出度边
				constraint += x[d][G->adjL[i][k]->id];
			for(unsigned int k = 0; k < G->adjRL[i].size(); k++) // 入度边
				constraint -= x[d][G->adjRL[i][k]->id];

			if(i == req[d].org)
				model.add(constraint == 1);
			else if(i == req[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}
	}

	for(int i = 0; i < G->m; i++){
		IloExpr load(env);
		for(int d = 0; d <  num; d++)
			load += req[d].flow*x[d][i];
		model.add( load <= G->Link[i]->capacity );  
		model.add( load <= z*G->Link[i]->capacity);	
	}	

	EEsolver.setOut(env.getNullStream());
	G->mlu = INF;
	if(EEsolver.solve()){
		G->mlu = EEsolver.getObjValue(); //mlu
		
		for(int i=0;i<G->m;i++){  
			double loadc = 0;
			for(int d=0;d < num;d++)
				loadc += EEsolver.getValue(x[d][i])*req[d].flow;
			G->Link[i]->latency = linearCal(loadc,G->Link[i]->capacity); //拟合
		}
		
		double del = 0;
		for(int d = 0; d <ornum; d++){  
			for(int i=0;i<G->m;i++)
				del +=  EEsolver.getValue(x[d][i])*G->Link[i]->latency;
		}
		G->delay = del;
		cout << "LB : "<<G->mlu<<"\t"<<del<<endl;
	}
	else{
		cout << "LB unfeasible"<<endl;
	}
	for(int i = 0; i < req.size(); i++)
		x[i].end();
	x.end();
	env.end();
}

void ORdictor(CGraph *G,vector<demand>& background,vector<demand>& overlay){    
	G->clearOcc();
	IloEnv env;
	IloModel model(env);
	IloCplex ORsolver(model);
	
	// background : underlying traffic,other overlay traffic
	int backnum = background.size();
	IloArray<IloNumVarArray> y(env, backnum); 
	for(int d = 0; d < backnum; d++)
		y[d] = IloNumVarArray(env, G->m, 0, 1); 

	//specific overlay
	int ornum = overlay.size();
	IloArray<IloNumVarArray> x(env, ornum); 
	for(int d = 0; d < ornum; d++)
		x[d] = IloNumVarArray(env, G->m, 0, 1);

	IloNumVarArray cost(env,G->m,0.0 ,INF); 
	IloExpr res(env);
	for(int i = 0;i < G->m; i++)
	{
		double c = G->Link[i]->capacity;
		IloExpr load(env);
		IloIntVar one(env,0,1);
		IloExpr load1(env);
		
		for(int d = 0; d <  ornum; d++){
			load1 += overlay[d].flow*x[d][i];
			load += overlay[d].flow*x[d][i];
		}
		model.add(one >= load1/G->Link[i]->capacity);

		for(int d = 0; d <  backnum; d++)
			load += background[d].flow*y[d][i];
		
		model.add(cost[i] >= 0 );		
		model.add(cost[i] >= (load-(1-one)*INF) );
		model.add(cost[i] >= ( 3*load-2*c/3-(1-one)*INF ) );
		model.add(cost[i] >= ( 10*load-16*c/3 -(1-one)*INF ) );
		model.add(cost[i] >= ( 70*load-178*c/3 -(1-one)*INF ) );
		res += cost[i];
	}
	model.add(IloMinimize(env, res));

	// overlay流量守恒约束
	for(int d = 0; d < ornum; d++){
		for(int i = 0; i < G->n; i++){  
			IloExpr constraint(env);
			for(unsigned int k = 0; k < G->adjL[i].size(); k++) // 出度边
				constraint += x[d][G->adjL[i][k]->id]; /////// y
			for(unsigned int k = 0; k < G->adjRL[i].size(); k++) // 入度边
				constraint -= x[d][G->adjRL[i][k]->id];
			// 出 - 入
			if(i == overlay[d].org)
				model.add(constraint == 1);
			else if(i == overlay[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}
	}

	//background 流守恒
	for(int d = 0; d < backnum; d++){
		for(int i = 0; i < G->n; i++){    // n为顶点数
			IloExpr constraint(env);
			for(unsigned int k = 0; k < G->adjL[i].size(); k++) // 出度边
				constraint += y[d][G->adjL[i][k]->id];
			for(unsigned int k = 0; k < G->adjRL[i].size(); k++) // 入度边
				constraint -= y[d][G->adjRL[i][k]->id];
			// 出 - 入
			if(i == background[d].org)
				model.add(constraint == 1);
			else if(i == background[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}
	}

	//带宽约束
	for(int i = 0; i < G->m; i++){
		IloExpr load(env);
		for(int d = 0; d <  ornum; d++)
			load += overlay[d].flow*x[d][i];
		for(int d = 0; d <  backnum; d++)
			load += background[d].flow*y[d][i];
		model.add( load <= G->Link[i]->capacity );  
	}

	ORsolver.setOut(env.getNullStream());
	G->delay = INF;

	if(ORsolver.solve()){
	
		// mlu
		double util = 0;
		for(int i = 0; i < G->m; i++){  
			double load = 0;
			for(int d = 0; d < ornum; d++)
				load += ORsolver.getValue(x[d][i])*overlay[d].flow;
			for(int d = 0; d < backnum; d++)
				load += ORsolver.getValue(y[d][i])*background[d].flow;
			util = max(util,load/G->Link[i]->capacity);
			G->Link[i]->use = load;
		}
		G->mlu = util;

		// delay
		double del = 0;
		for(int i = 0; i < G->m; i++){  
			double cur = linearCal(G->Link[i]->use,G->Link[i]->capacity);
			for(int d = 0; d < ornum; d++)
				del += ORsolver.getValue(x[d][i])*cur;
		}
		G->delay = del;
		cout << "OR\t利用率 "<<util<< "  \t延时 "<<del<<endl;
	}
	else{
		env.out() << "OR unfeasible" << endl;
	}

	for(int i = 0; i < overlay.size(); i++)
		x[i].end();
	x.end();
	for(int i = 0; i < background.size(); i++)
		y[i].end();
	y.end();
	env.end();
}

#endif