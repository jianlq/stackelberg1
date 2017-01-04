#include"LB.h"
#include"evolutionbit.h"
#include"evolution.h"

vector<demand> eqOR;
vector<demand> eqTE;
vector<demand> eqbase;

 bool selfishrouting = false;
 int ORNUM ;
void SelfishRouting(CGraph *G, double &maxLinkUtil, double &sumDelay){
	selfishrouting = true;
	maxLinkUtil = sumDelay = 0;
	G->clearOcc();
	int NUMEVENT = eqTE.size();
	for(int i = 0; i < NUMEVENT; i++){
		double can = G->dijkstraLB(i,eqTE[i].org, eqTE[i].des, eqTE[i].flow);
		if(can + SMALL >= INF){ //无解
			maxLinkUtil = INF;
			sumDelay = INF;
			return;
		}
	}

	for(int i = 0; i < G->m; i++){
		maxLinkUtil = max(maxLinkUtil, (double)G->Link[i]->use);
	}

	for(int d = 0; d < eqOR.size(); d++){
		for(unsigned int ij = 0;ij < G->reqPathID[d].size(); ij++)
			sumDelay +=  eqOR[d].flow*linearCal( G->Link[G->reqPathID[d][ij]]->use,G->Link[G->reqPathID[d][ij]]->capacity);
	}
	selfishrouting = false;
}

void CentralizedOptimization(CGraph *G, double &maxLinkUtil, double &sumDelay){
	G->clearOcc();
	ilpSolver(G,eqTE,ORNUM);
	maxLinkUtil = G->mlu;
	sumDelay = G->delay;
	if(G->mlu + SMALL > INF){
		maxLinkUtil = INF;
		sumDelay = INF;
	}
}

void Stackelberg(CGraph *G, CGraph *GOR,double &maxLinkUtil, double &sumDelay){
	evoluPopu popubit(50,G->m,G,GOR,&eqTE);
	evoluDiv ret = popubit.evolution();
	maxLinkUtil = INF;
	sumDelay = INF;
	if(ret.ability + SMALL < INF){
		maxLinkUtil = ret.ability;
		sumDelay = ret.delay;
	}
}

int TestCase(CGraph *G, CGraph *GOR,double &alg1, double &alg1w,double &alg2, double &alg2w,double &alg3, double &alg3w){	

	printf("\n*****************************************\n");
	double mlu, sd;

	SelfishRouting(G, mlu, sd);  
	double dbase = sd;
	printf("SelfishRouting:\t\t%f\t%f\n", mlu, sd);
	double gremlu = mlu;
	
	CentralizedOptimization(G, mlu, sd);  
	double ubase = mlu;
	printf("Centralized Optimization:\t%f\t%f\n", mlu, sd);
	
	alg1 = gremlu/ubase;
	alg1w = 1;

	////ilp
	alg2 = mlu/ubase;
	alg2w = sd/dbase;

	Stackelberg(G, GOR,mlu, sd);
	printf("Stackelberg:\t\t%f\t%f\n", mlu, sd);
	alg3 = mlu/ubase;
	alg3w = sd/dbase;

	if(ubase + SMALL > INF || mlu >= INF)
		return 0;
	else
		return 1;
}

int main(){
	srand((unsigned)time(NULL));
	int Time = 3;
	int CASEnum= 200;	

	for(int i = 0;i < Time;i++){
		FILE *res = fopen("outputFile//result.csv", "a");		
		fprintf(res,"%d case average\n",CASEnum);
		fprintf(res,",SelfishRouting,,,Traffic Engineering,,,Stackelberg\n");
		fclose(res);
		for(int NUMDEMAND = 5; NUMDEMAND <= 200; NUMDEMAND += 5){
			int sucCase = 0 ;
			double   alg1 = 0, alg2 = 0, alg3 = 0;
			double   alg1w = 0, alg2w = 0, alg3w = 0;

			for(int casenum = 0; casenum < CASEnum; casenum++){

				genGraph(18,98,"inputFile//graph.txt");
				genGraphOR(18,9,45,"inputFile//graphOR.txt");
				CGraph *G = new CGraph("inputFile//graph.txt");
				CGraph *GOR = new CGraph("inputFile//graphOR.txt");
			
				
				eqbase.clear();
				int BGNUM;
				if( NUMDEMAND > 100){
					ORNUM = NUMDEMAND*(2.0/5.0);
					BGNUM = NUMDEMAND - ORNUM;
					
				}
				else{
					ORNUM = rand()%NUMDEMAND + 1;
					BGNUM = NUMDEMAND - ORNUM;
					
				}
				
				for(int i = 0; i < BGNUM; i++){
					int s , t;
					do{
						s = rand()%G->n;
						t = rand()%G->n;
					}while( s == t || G->canNotReach(s,t));
					eqbase.push_back(demand(s, t, rand()%(MAXDEMAND)+1));
				}

			
				eqOR.clear();
				for(int i = 0; i < ORNUM; i++){
					int n = GOR->ver.size();
					int s1, t1;
					int org , dst;
					do{
						s1 = rand()% n;
						t1 = rand()% n;
						org = GOR->ver[s1];
						dst = GOR->ver[t1];
					}while( s1 == t1 || G->canNotReach(org, dst) || GOR->canNotReach( org, dst) );	
					eqOR.push_back(demand(org, dst, rand()%(MAXDEMAND)+1));
				}
				
				eqTE.clear();
				for(unsigned int i=0;i<eqOR.size();i++)
					eqTE.push_back(eqOR[i]);
				for(unsigned int i=0;i<eqbase.size();i++)
					eqTE.push_back(eqbase[i]);
				
				double a, b, c, d, e, f;
				int success = TestCase(G,GOR, a, b, c, d, e, f);
				sucCase += success;
			
				if(success){
					alg1 += a;
					alg1w += b;
					alg2 += c;
					alg2w += d;
					alg3 += e;
					alg3w += f;
				}
				delete G;
				delete GOR;
				
			} // end of CASENum for

			res = fopen("outputFile//result.csv", "a");		
			fprintf(res, "%d,%lf,%lf,,%lf,%lf,,%lf,%lf \n",NUMDEMAND, alg1/sucCase, alg1w/sucCase, alg2/sucCase, alg2w/sucCase, alg3/sucCase,alg3w/sucCase); 
			fclose(res);
		
		}
	}
	system("pause");
	return 0;
}