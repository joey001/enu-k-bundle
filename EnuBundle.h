#pragma once
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <ctime>
#include "Defines.h"
#include "RandList.h"
#include "LinearHeap.h"
using namespace std;

#ifdef DBGMOD
class DBGNode {
public:
	unsigned long long nodeid;
	ui vtx;
	//ui orvtx;
	vector<int> P;
	vector<int> cand;
	vector<int> excl;
	DBGNode *take;
	DBGNode *notake;
	DBGNode(ui _vtx, unsigned long long _nodeid,
		RandList &_P,
		RandList &_Cand,
		RandList &_Excl) :vtx(_vtx), nodeid(_nodeid) {
		for (ui i = 0; i < _P.getSize(); i++) {
			P.push_back(_P.get(i));
		}
		for (ui i = 0; i < _Cand.getSize(); i++) {
			cand.push_back(_Cand.get(i));
		}
		for (ui i = 0; i < _Excl.getSize(); i++) {
			excl.push_back(_Excl.get(i));
		}
	}
	void showNode() {
		printf("ID: %llu, vtx: %u\n", nodeid, vtx);
		printf("P:");
		for (auto u : P) {
			printf("%u ", u);
		}
		printf("\nCand:");
		for (auto u : cand) {
			printf("%u ", u);
		}
		printf("\nExcl:");
		for (auto u : excl) {
			printf("%u ", u);
		}
	}
	void addTake(DBGNode *_take) {
		take = _take;
	}
	void addNoTake(DBGNode *_notake) {
		notake = _notake;
	}
};
#endif // DBGMOD

class Solution {
private:
	set<ui > elems;
public:
	Solution(vector<ui> unorderedSol) {
		elems.insert(unorderedSol.begin(), unorderedSol.end());
	}
	bool operator<(const Solution &s) const {
		set<ui>::iterator mit = elems.begin();
		set<ui>::iterator sit = s.elems.begin();
		while (mit != elems.end() && sit != s.elems.end()) {
			if (*mit < *sit) 
				return true;
			else if(*mit > *sit) {
				return false;
			}
			else {
				mit++, sit++;
			}
		}
		if (mit == elems.end() && sit != s.elems.end()) 
			return true;
		return false;
	}
	void printSol() {
		for (auto u : elems)
			printf("%u ", u);
		printf("\n");
	}
};

class EnuBundle
{
	//orignial 
	ui n, m;
	ui* pstart;
	ui* edges;
	ui* reverse;
	//map<ui, ui> idmp; //idmp[orin] maps orignal id to new consectutive id
	ui* dseq;
	ui *dpos;
	ui* core;
	ui k;
	ui lb;

	//induced subgraph
	ui bvtx;
	ui *bmark;
	ui* bID; //bID[x] is the original vertex id of x in block(v)
	ui* nID;// nID[u] is the new id of u in G
	ui* bstart;
	ui* bedges;
	ui bm;
	ui bn; // block vertex number
	MBitSet* badc;
	MBitSet* binv;

	RandList P;		//current solution
	RandList Cand;
	RandList Excl;

	ui* neiInP; // neiInP[u] is the number of neighbors in P	
	ui* neiInG; //number of neibors in the graph induced by G[P\cup Cand]	
	ListLinearHeap *degHeap;

	//void reduceCand(int v);	//Reduce a vertex from candidate set
	//void backToCand(int v); //Move a vertex back to candidate set

	void removeFrCand(ui u);

	void addToCand(ui u);

	int canMoveToP(ui u);

	//void addToP(ui u);

	//void removeFrP(ui u);

	void CandToP(ui u);

	void PToCand(ui u);


	int isGlobalMaximal();

	void branch();

	void recurSearch(ui start);

	void showSolution();


	ui checkSolution();

	void stopAsSolution();

	int cntplex;
	set<Solution> allsols;

	unsigned long long nnodes;
	clock_t startclk;
	clock_t sortclk;
	clock_t enumclk;
	
public:
	int readRawDIM10Text(const char * filepath);
	int readRawSNAPText(const char* filePath);
	int readBinaryGraph(const char * filepath);
	int writeBinaryGraph(const char * filepath);

	int degeneracyOrder(ui * seq, ui * core, ui * pos);
	int writeBlockToBin(char * filepath);
	int buildBlock(int v);
	
	void enumPlex(ui _k, ui _lb);

	ui checkMaximal(vector<ui>& S, ui * degS);

	void enumBruteforce(vector<ui>& CurS, vector<ui>& CandS, vector<ui>& VisitS, ui * degCur);

	EnuBundle();
	~EnuBundle();
};

