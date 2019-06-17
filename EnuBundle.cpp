#include <fstream>
#include "EnuBundle.h"
#include "Utility.h"
using namespace std;

#ifdef DBGMOD
//#define BFSEARCH
//#define SHOWBLK
//#define SHOWSOL
//#define	SHOWRECUR
#endif

int EnuBundle::writeBinaryGraph(const char* filepath) {
	vector<ui> nodes;

	FILE *f = Utility::open_file(filepath, "wb");
	ui tt = sizeof(ui);
	fwrite(&tt, sizeof(ui), 1, f); //length of ui
	fwrite(&n, sizeof(ui), 1, f);
	fwrite(&m, sizeof(ui), 1, f);
	ui *degree = new ui[n];
	for (ui i = 0; i < n; i++)
		degree[i] = pstart[i + 1] - pstart[i];
	fwrite(degree, sizeof(ui), n, f);
	fwrite(edges, sizeof(ui), m, f);
	fclose(f);
	return 0;
}


/**
TODO: reverse[] is not supported.
*/
int EnuBundle::readBinaryGraph(const char* filepath) {
	FILE *f = Utility::open_file(filepath, "rb");
	ui tt;
	fread(&tt, sizeof(ui), 1, f);
	if (tt != sizeof(ui)) {
		printf("sizeof unsigned int is different: file %u, machine %lu\n", tt, sizeof(ui));
	}
	fread(&n, sizeof(ui), 1, f);
	fread(&m, sizeof(ui), 1, f);
	printf("n=%u, m=%u (undirected)\n", n, m);
	ui* degree = new ui[n];
	fread(degree, sizeof(ui), n, f);

	if (pstart != nullptr) delete[] pstart;
	pstart = new ui[n + 1];
	if (edges != nullptr) delete[] edges;
	edges = new ui[m];

	//if (reverse != nullptr) delete[] reverse;
	//reverse = new ui[m];

	pstart[0] = 0;
	for (ui i = 0; i < n; i++) {
		if (degree[i] > 0) fread(edges + pstart[i], sizeof(ui), degree[i], f);
		pstart[i + 1] = pstart[i] + degree[i];
	}
	fclose(f);
	delete[] degree;
	return 0;
}

/**
TODO: reverse[] is not supported.
*/
int EnuBundle::readRawDIM10Text(const char* filepath) {
	ifstream infile;
	char buf[1024];
	vector<pair<ui, ui> > epairs;
	vector<ui> nodes;
	//FILE *f = Utility::open_file(filepath, "r");
	infile.open(filepath, ios::in);
	if (!infile.is_open()) {
		fprintf(stderr, "can not find file %s\n", filepath);
		exit(1);
	}

	infile.getline(buf, 1024);
	while (buf[0] == '%') infile.getline(buf, 1024);

	stringstream ss(buf);
	int fmt;
	ss >> n >> m >> fmt;
	m *= 2;
	pstart = new ui[n + 1];
	edges = new ui[m];
	reverse = new ui[m];
	ui j = 0;
	for (ui u = 0; u < n; u++) {
		pstart[u] = j;
		infile.getline(buf, 1024);
		stringstream ss(buf);
		int nei;
		while (ss >> nei) {
			edges[j] = nei - 1;
			reverse[j] = u;
			j++;
		}
		sort(edges + pstart[u], edges + j);
	}
	pstart[n] = j;
	assert(j == m);
}


int EnuBundle::readRawSNAPText(const char* filepath) {
	ifstream infile;
	char buf[1024];
	vector<pair<ui, ui> > epairs;
	vector<ui> nodes;
	//FILE *f = Utility::open_file(filepath, "r");
	infile.open(filepath, ios::in);
	if (!infile.is_open()) {
		fprintf(stderr, "can not find file %s\n", filepath);
		exit(1);
	}
	int max_id = 0;
	int from, to;
	while (infile.getline(buf, 1024)) {
		char *p = buf;
		while (*p == ' ' && *p != '\0') p++;
		if (*p == '#' || p == '\0') continue;
		stringstream ss(buf);
		ss >> from >> to;
		epairs.push_back(make_pair(from, to));
		epairs.push_back(make_pair(to, from));
		nodes.push_back(from);
		nodes.push_back(to);
	}
	infile.close();

	//去除重点
	sort(nodes.begin(), nodes.end());
	nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());

	//去除重边
	sort(epairs.begin(), epairs.end());
	epairs.erase(unique(epairs.begin(), epairs.end()), epairs.end());

	//去除孤立点。对剩下的点重新编号。 in p2p-Gnutella04, 10452 is a isolate vertex
	ui contn = 1;
	map<ui, ui> idmp;
	for (ui i = 0; i < nodes.size(); i++) {
		idmp[nodes[i]] = i;
		if (nodes[i] != i) {
			contn = 0;
		}
	}
	if (contn == 0) printf("Node ids are not preserved! \n");

	n = nodes.size();
	m = epairs.size();
	printf("n = %s, (undirected) m = %s\n",
		Utility::integer_to_string(n).c_str(),
		Utility::integer_to_string(m / 2).c_str());

	pstart = new ui[n + 1];
	edges = new ui[m];
	reverse = new ui[m];
	ui j = 0;
	for (ui i = 0; i < n; i++) {
		pstart[i] = j;
		while (j < m && epairs[j].first == nodes[i]) {
			edges[j] = idmp[epairs[j].second];
			reverse[j] = i;
			++j;
		}
	}
	pstart[n] = j;
}

int EnuBundle::degeneracyOrder(ui *seq, ui *core, ui* pos) {
	ui *id_s = seq, *degree = core;
	for (ui i = 0; i < n; i++) {
		id_s[i] = i;
		degree[i] = pstart[i + 1] - pstart[i];
	}

	ui max_core = 0;
	ListLinearHeap *linear_heap = new ListLinearHeap(n, n - 1);
	linear_heap->init(n, n - 1, id_s, degree);
	memset(core, 0, sizeof(ui)*n);
	for (ui i = 0; i < n; i++) {
		ui u, key;
		linear_heap->pop_min(u, key);
		if (key > max_core)
			max_core = key;
		seq[i] = u;
		core[u] = max_core;
		pos[u] = i;

		for (ui j = pstart[u]; j < pstart[u + 1]; j++)
			if (core[edges[j]] == 0) {
				linear_heap->decrement(edges[j]);
			}
	}
	delete linear_heap;

	return max_core;
}

#ifdef DBGMOD
int EnuBundle::writeBlockToBin(char* filepath) {
	FILE *f = Utility::open_file(filepath, "wb");
	ui tt = sizeof(ui);
	fwrite(&tt, sizeof(ui), 1, f); //length of ui
	fwrite(&bn, sizeof(ui), 1, f);
	fwrite(&bm, sizeof(ui), 1, f);
	ui *degree = new ui[bn];
	for (ui i = 0; i < bn; i++)
		degree[i] = bstart[i + 1] - bstart[i];
	fwrite(degree, sizeof(ui), bn, f);
	fwrite(edges, sizeof(ui), bm, f);
	fclose(f);
	return 0;
}
#else
int EnuBundle::writeBlockToBin(char* filepath) {}
#endif

int EnuBundle::buildBlock(ui v) {
	//mark all the inducing vertices
	memset(bmark, (ui)0, sizeof(ui)* n);
	bmark[v] = 1;
	bvtx = v;
	bn = 1;
	vector<ui> nei1s;
	vector<ui> blkset;
	for (ui j = pstart[v]; j < pstart[v + 1]; j++) {
		ui nei = edges[j];
		if (!bmark[nei] && dpos[nei] > dpos[v]) {
			bmark[nei] = 1;	//N(v)
			blkset.push_back(nei);
			nei1s.push_back(nei);
			bn++;
		}
	}	
	for (auto u : nei1s) {
		for (ui k = pstart[u]; k < pstart[u + 1]; k++) {
			ui nei2 = edges[k];
			if (!bmark[nei2] && dpos[nei2] > dpos[v]) {
				bmark[nei2] = 2; //N2(v)
				blkset.push_back(nei2);
				bn++;
			}
		}
	}
	//size pruning
	for (ui u : blkset) {
		if (bmark[u]) {
			//set<ui> vnei(edges + pstart[v], edges + pstart[v + 1]);
			//caclculate intersections
			ui cnt = 0;
			ui p1 = pstart[v], p2 = pstart[u];
			while (p1 != pstart[v + 1] && p2 != pstart[u + 1]) {
				if (edges[p1] < edges[p2]) {
					p1++;
				}
				else if (edges[p1] > edges[p2]) {
					p2++;
				}
				else {//equal
					if (bmark[edges[p1]]) {
						cnt++;
					}
					p1++, p2++;
				}
			}
			if (bmark[u] == 1) {
				if (cnt + 2 * k < lb) {
					bmark[u] = 0;
					//canprune = 1;
					bn--;
				}
			}else{
				if (cnt + 2 * k - 2 < lb) {
					bmark[u] = 0;
					//canprune = 1;
					bn--;
				}
			}
			
		}
	}
		//renumber the vertices, bid[x] is the new number of vertex x
	if (bn < lb) {
		//delete[] bmark;
		return 0;
	}
	ui cnt = 0;	
	memset(bID, n, sizeof(ui)*n);
	memset(nID, n, sizeof(ui)*n);
	for (ui j = dpos[v]; j < n; j++) {
		if (bmark[dseq[j]]) {
			bID[cnt] = dseq[j];
			nID[dseq[j]] = cnt;
			cnt++;
		}
	}
	assert(cnt == bn);
	//initialize the subgraph	
	/*if (bstart != nullptr) delete[] bstart;
	if (bedges != nullptr) delete[] bedges;	
	bstart = new ui[bn+1];
	bedges = new ui[bn*bn];
	*/

	bm = 0;
	for (ui j = 0; j < bn; j++) {
		int oruid = bID[j];
		badc[j].init(bn);
		//binv[j].init(bn);
		bstart[j] = bm;
		for (ui k = pstart[oruid]; k < pstart[oruid + 1]; k++) {
			if (bmark[edges[k]]) {
				bedges[bm++] = nID[edges[k]];
				badc[j].set(nID[edges[k]]);
			}
		}
		//binv[j].flip();
		//sort(bedges + bstart[j], bedges + bm); //we can optimize here by sorting the vertex
	}
	bstart[bn] = bm;
	return 1;

}

void EnuBundle::CandToP(ui u) {
	assert(Cand.contains(u));
	P.add(u);
	Cand.remove(u);

	//update neiInP
	for (ui i =  bstart[u]; i < bstart[u + 1]; i++) {
		ui nei = bedges[i];
		neiInP[nei]++;
	}	

}

void EnuBundle::PToCand(ui u) {
	assert(P.contains(u));
	P.remove(u);
	Cand.add(u);

	//update neiInP	
	for (ui i = bstart[u]; i < bstart[u + 1]; i++) {
		ui nei = bedges[i];
		neiInP[nei]--;		
	}
}

void EnuBundle::removeFrCand(ui u) {
	assert(Cand.contains(u));
	Cand.remove(u);
	for (ui i = bstart[u]; i < bstart[u + 1]; i++) {
		ui nei = bedges[i];
		neiInG[nei]--;
	}
}

void EnuBundle::addToCand(ui u) {
	assert(!Cand.contains(u)); 
	Cand.add(u);
	for (ui i = bstart[u]; i < bstart[u + 1]; i++) {
		ui nei = bedges[i];
		neiInG[nei]++;
	}
}

int EnuBundle::canMoveToP(ui u) {
	//assert(Cand.contains(u));
	if (neiInP[u] + k < P.getSize() + 1) {
		return 0;
	}
	for (ui i = 0; i < P.getSize(); i++) {
		ui v = P.get(i);
		assert(neiInP[v] + k >= P.getSize());
		if (neiInP[v] + k == P.getSize() && (!badc[v].test(u))) {
			return 0;
		}
	}
	return 1;	
}

void EnuBundle::updateSatSet() {
	satInP.clear();
	for (ui i = 0; i < P.getSize(); i++) {
		ui u = P.get(i);
		if (neiInP[u] + k == P.getSize())
			satInP.insert(u);
	}
}


int EnuBundle::isGlobalMaximal() {
	vector<ui> candidates;
	ui init = 0;
	updateSatSet();
	if (!satInP.empty()) {
		for (auto u : satInP) {
			ui oru = bID[u];
			if (!init) {
				candidates.insert(candidates.begin(), edges + pstart[oru], edges + pstart[oru + 1]);
				init = 1;
			}
			else {
				//vector<ui> tmp;
				set_intersection(edges + pstart[oru], edges + pstart[oru + 1], candidates.begin(), candidates.end(), std::back_inserter(candidates));
				//candidates = tmp;
			}
		}
	}
	if (!satInP.empty()) {
		for (auto oru : candidates) {
			if (dpos[oru] < dpos[bvtx] && canGloablaAdd(oru)) {
				return 0;
			}
		}
		return 1;
	}
	else {
		return isGlobalMaximal2();
	}
}
int EnuBundle::canGloablaAdd(ui oru) {
	if (core[oru] + k < P.getSize() + 1)
		return 0; //can not add
	int scon = 0;
	int pcon = 0;
	for (ui j = pstart[oru]; j < pstart[oru + 1]; j++) {
		ui nei = edges[j];
		if (bmark[nei]) {
			ui idInBlk = nID[nei];
			if (P.contains(idInBlk)) {
				pcon++;
				if (satInP.find(idInBlk) != satInP.end())
					scon++;
			}
		}
	}
	if (scon == satInP.size() && pcon + k >= P.getSize() + 1)
		return 1; //can add
	return 0;	//can not add
}
//TODO
int EnuBundle::isGlobalMaximal2() {
	updateSatSet();
	for (ui i = 0; i < dpos[bvtx]; i++) {		
		ui u = dseq[i];		
		if (canGloablaAdd(u))
			return 0;
	}
	//printf("SAT vs P %u: %u\n", sat.size(), P.getSize());
	return 1;
}

void EnuBundle::recurSearch(ui start) {
	//reduce
	assert(Cand.contains(start));
	//assert(canMoveToP(start));
	vector<ui> rcand;
	vector<ui> rexcl;

	if (Cand.getSize() + P.getSize() < lb) {
		return;
	}	
	CandToP(start);
	//update candidate
	for (ui i = 0; i < Cand.getSize(); i++) {
		ui u = Cand.get(i);
		if (!canMoveToP(u)) {			
			rcand.push_back(u);
		}
	}

	for (auto u : rcand) {
		removeFrCand(u);
	}
	//update excl
	for (ui i = 0; i < Excl.getSize(); i++) {
		ui u = Excl.get(i);
		if (!canMoveToP(u))
			rexcl.push_back(u);
	}
	for (auto u: rexcl) {
		Excl.remove(u);
	}
	branch(); // branch search	
	
	//recover
	PToCand(start);

	for (auto u : rcand) {
		addToCand(u);
	}
	for (auto u : rexcl) {
		Excl.add(u);
	}
}


void EnuBundle::checkSolution() {
	ui ismax = isGlobalMaximal2();
	if (ismax) {
		cntplex++;
#ifdef SHOWSOL
		printf("Sol:");
		for (ui i = 0; i < P.getSize(); i++) {
			printf("%u ", bID[P.get(i)]);
		}
		printf("\n");
		//ui rt = dbgCheckSolution();
		//if (!rt)
		//	printf("Wrong solution pause\n");
#endif
	}
}

ui EnuBundle::dbgCheckSolution() {
	//check k-plex
	if (P.getSize() < lb) {
		printf("ERROR:Samller than lower bound\n");
		exit(-1);
	}
	//get sat
	set<ui> sat;
	for (ui i = 0; i < P.getSize(); i++) {
		ui u = P.get(i);
		ui deg = 0;
		for (ui j = 0; j < P.getSize(); j++) {
			if (P.get(j) != u && badc[u].test(P.get(j))) //trust block adjacency matrix
				deg++;
		}
		if (deg == P.getSize() - k)
			sat.insert(u);
		else if (deg < P.getSize() - k) {
			printf("ERROR:Not a k-plex\n");
			exit(0); //
		}
	}
	//check maximality
	for (ui u = 0; u < n; u++) {
		if (bmark[u] && P.contains(nID[u])) continue;// in subgraph
		if (core[u] + k >= P.getSize() + 1) {
			int scon = 0;
			int pcon = 0;
			for (ui j = pstart[u]; j < pstart[u + 1]; j++) {
				ui nei = edges[j];
				if (bmark[nei]) {
					ui idInBlk = nID[nei];
					if (P.contains(idInBlk))
						pcon++;
					if (sat.find(idInBlk) != sat.end())
						scon++;
				}
			}
			if (scon == sat.size() && pcon + k >= P.getSize() + 1) {
				printf("ERROR: Not maximal\n");
			}
				
		}
	}
	//check uniqueness
	vector<ui> vecsol;
	for (ui i = 0; i < P.getSize(); i++) {
		vecsol.push_back(bID[P.get(i)]);
	}
	Solution sol(vecsol);
	if (allsols.find(sol) != allsols.end()) {
		printf("ERROR: Not unique\n");
	}
	else {
		allsols.insert(sol);
	}
}

/**
Stop as the whole graph is a solution.
*/
void EnuBundle::stopAsSolution() { 
	vector<ui> rcand;
	vector<ui> rexcl;
	nnodes++;
	for (ui i = 0; i < Cand.getSize(); i++)
		rcand.push_back(Cand.get(i));

	for (auto u : rcand) {
		CandToP(u);
	}
	for (ui i = 0; i < Excl.getSize(); i++) {
		ui u = Excl.get(i);
		if (!canMoveToP(u)) {
			rexcl.push_back(u);
		}
	}
#ifdef SHOWRECUR
	printf("-------Nodes: %u (LEAF)-------------\n", nnodes);
	printf("P:");
	for (ui i = 0; i < P.getSize(); i++)
		printf("%u ", P.get(i));
	printf("\nCand:");
	for (ui i = 0; i < Cand.getSize(); i++)
		printf("%u ", Cand.get(i));
	printf("\nExcl:");
	for (ui i = 0; i < Excl.getSize(); i++)
		printf("%u ", Excl.get(i));
	printf("\n");
#endif
	//Check maximal, Excl is not changed
	if (rexcl.size() == Excl.getSize()) { // Excel will be empty
		checkSolution();
	}
	for (auto u : rcand) {
		PToCand(u);
	}
}

 void EnuBundle::branch() {
	nnodes++;
#ifdef SHOWRECUR
	printf("-------Nodes: %u-------------\n", nnodes);
	printf("P:");
	for (ui i = 0; i < P.getSize(); i++)
		printf("%u ", P.get(i));
	printf("\nCand:");
	for (ui i = 0; i < Cand.getSize(); i++)
		printf("%u ", Cand.get(i));
	printf("\nExcl:");
	for (ui i = 0; i < Excl.getSize(); i++)
		printf("%u ", Excl.get(i));
	printf("\n");
#endif
	if (Cand.empty() && Excl.empty()) {
		//maximal;
		//P.printList();
		if (P.getSize() >= lb) {
			checkSolution();
			return;
		}
	}
	if (Cand.empty()) return;

	if (Cand.getSize() + P.getSize() < lb) {
		return;
	}	

	//find mindeg 
	ui minu = bn;
	ui branchu = bn;
	for (ui u = 0; u < bn; u++) {			
		if (Cand.contains(u) || P.contains(u)) {
			if(minu == bn || neiInG[u] < neiInG[minu])
				minu = u;
		}
		if (Cand.contains(u)) {
			if (branchu == bn || neiInG[u] < neiInG[branchu])
				branchu = u;
		}

	}

#ifdef SHOWRECUR
	printf("Mindeg vertex %u ,degG %u, degP %u\n", minu, neiInP[minu], neiInG[minu]);
#endif
	assert(minu != bn && branchu != bn);
	//The whole graph is a k-plex
	if (neiInG[minu] >= Cand.getSize() + P.getSize() - k) { 			
		stopAsSolution();
	}
	else if (Cand.getSize() + P.getSize() <= lb) {
#ifdef SHOWRECUR
		printf("No better solution as G=lb, G is not kplex \n");
#endif
		return; // the whole graph is not a k-plex and it is not possible to find a better one 
	}
	else {
#ifdef SHOWRECUR
		printf("Recur %u\n", branchu);
#endif
		recurSearch(branchu);	

		removeFrCand(branchu);
		Excl.add(branchu);
		branch();

		Excl.remove(branchu);
		addToCand(branchu);		
	}
	
}
 
void EnuBundle::enumPlex(ui _k, ui _lb)
{	
	startclk = clock();
	k = _k;
	lb = _lb;
	cntplex = 0;
	
	Excl.init(n);
	P.init(n);
	Cand.init(n);
	neiInG = new ui[n];
	neiInP = new ui[n];
	
	dseq = new ui[n];
	core = new ui[n];
	dpos = new ui[n];	//vtxmp[v] is the position of v in degeneracy order 
	ui maxcore = degeneracyOrder(dseq, core, dpos); // fast degeneracy order
	sortclk = clock();
		//build subgraph		
	nID = new ui[n];
	bmark = new ui[n];
	bID = new ui[n];
	bstart = new ui[n+1];
	bedges = new ui[m];
	
	badc = new MBitSet[min(maxcore*maxcore, n)];
	//binv = new MBitSet[min(maxcore*maxcore, n)];

	for (ui i = 0; i < n - 1; i++) {		
		ui v = dseq[i];
		if (core[v] +k >= lb ) {
			int built = buildBlock(v);
			if (!built) {
				printf("Vertex %u [%u] discard \n", i, v);
				continue;
			}
			printf("Block %u[%u]: %u %u\n", i, v, bn, bm);
#ifdef SHOWBLK	
			for (ui i = 0; i < bn; i++) {

				printf("%u(%u):", i, bID[i]);
				for (ui j = bstart[i]; j < bstart[i + 1]; j++) {
					printf("%u(%u) ", bedges[j], bID[bedges[j]]);
				}
				printf("\n");
			}
#endif

#ifdef BFSEARCH
			//Enum jazz.bin induced by vtx 178 is long, (when i = 168)
			vector<ui> CurS;
			vector<ui> CandS;
			vector<ui> VisitS;
			ui *degS = new ui[bn];
			CurS.push_back(0);
			memset(degS, 0, sizeof(ui)*bn);
			for (ui j = bstart[0]; j < bstart[1]; j++) {
				degS[bedges[j]] = 1;
			}
			for (ui u = bn - 1; u >= 1 ; u--) {
				CandS.push_back(u);
			}
			nnodes = 0;

			//if (i == 7)
				//printf("pause\n");
			enumBruteforce(CurS, CandS, VisitS, degS);
#else
			//build cand;
			P.clear();
			//P.add(0);
			Cand.clear();

			memset(neiInP, 0, sizeof(int) * bn);
			//for (ui j = bstart[0]; j < bstart[1]; j++) {
			//	neiInP[bedges[j]] = 1;
			//}

			for (ui u = 0; u < bn; u++) {
				Cand.add(u);
			}
			
			for (ui j = 0; j < bn; j++) {
				neiInG[j] = bstart[j + 1] - bstart[j];
			}
			Excl.clear();
			nnodes = 0;
			recurSearch(0);
			//branch();
#endif// BFSEARCH
			printf("Node number: %llu\n", nnodes);
		}
	}
	
	enumclk = clock();
	printf("Number of %u-cplex larger than %u:  %u\n", k, lb, cntplex);
	printf("Total search time %.2f\n", Utility::elapse_seconds(startclk, enumclk));
	printf("Sort time %.2f\n", Utility::elapse_seconds(startclk, sortclk));
	//printf("Totoal nodes %u \n", nnodes);
}
#ifdef BFSEARCH
ui EnuBundle::checkMaximal(vector<ui> &S, ui* degS) {
	set<ui> sat;
	for (auto u: S) {
		if (degS[u] + k == S.size())
			sat.insert(u);
	}
	for (ui i = 0; i < dpos[bvtx]; i++) {
		ui u = dseq[i];
		int scon = 0;
		int pcon = 0;
		for (ui j = pstart[u]; j < pstart[u + 1]; j++) {
			ui nei = edges[j];
			//TODO:WRONG
			ui nidInBlk = dpos[nei] - dpos[bvtx];
			if (dpos[nei] >= dpos[bvtx] && binary_search(S.begin(), S.end(), nidInBlk)) //nidBilk is in S
				pcon++;
			if (sat.find(nidInBlk) != sat.end())
				scon++;
		}
		if (scon == sat.size() && pcon + k - 1 >= S.size())
			return 0;
	}
	return 1;

}
void EnuBundle::enumBruteforce(vector<ui> &CurS, vector<ui> &CandS, vector<ui> &VisitS, ui* degCur) {
	nnodes++;
#ifdef SHOWRECUR
	printf("-------Nodes: %u-------------\n", nnodes);
	printf("P:");
	for (auto u : CurS)
		printf("%u ", u);
	printf("\nCand:");
	for (auto u : CandS)
		printf("%u ", u);
	printf("\nExcl:");
	for (auto u : VisitS)
		printf("%u ", u);
	printf("\n");
#endif
	if (CandS.empty() && VisitS.empty()) {
		if (checkMaximal(CurS, degCur)) {
			cntplex++;
#ifdef SHOWSOL
			printf("SOL:");
			for (auto u : CurS)
				printf("%u ", bID[u]);
			printf("\n");
#endif
		}
	}
	
	while(!CandS.empty()) {
		ui c = CandS.back();
		CandS.pop_back();

		CurS.push_back(c);

		for (ui i = bstart[c]; i < bstart[c+1]; i++) {
			ui u = bedges[i];
			degCur[u]++;
		}
		//Update sat
		vector<ui> sat;
		if (degCur[c] +k == CurS.size()) {
			sat.push_back(c);
		}
		for (auto v: CurS) {	
			if (v!= c && degCur[v] + k == CurS.size())
				sat.push_back(v);
		}
		

		vector<ui> newCand;
		vector<ui> newVisit;		
		for (auto u : CandS) {			
			if (degCur[u] + k< CurS.size()) {				
				continue;
			}
			ui keep = 1;
			for (auto s : sat) {
				if (!badc[u].test(s)) {
					keep = 0;
					break;
				}
			}
			if (keep)
				newCand.push_back(u);
		}
		for (auto u : VisitS) {
			if (degCur[u] + k < CurS.size())
				continue;
			ui keep = 1;
			for (auto s : sat) {
				if (!badc[u].test(s)) {
					keep = 0;
					break;
				}
			}
			if (keep)
				newVisit.push_back(u);
		}

		//prune
		if (newCand.size() + CurS.size() < lb) {
#ifdef SHOWRECUR
			printf("Prune by size < lb\n");
#endif // DEBUG
			//return;
		}
		else {

			enumBruteforce(CurS, newCand, newVisit, degCur);
		}
		//remove c from P
		CurS.pop_back();
		for (ui i = bstart[c]; i < bstart[c + 1]; i++) {
			ui u = bedges[i];
			degCur[u]--;
		}

		VisitS.push_back(c);
	}

}
#endif

EnuBundle::EnuBundle()
{
	n = 0;
	m = 0; 
	pstart = nullptr;
	edges = nullptr;
	reverse = nullptr;
	bstart = nullptr;
	bID = nullptr;
	bedges = nullptr;	
	badc = nullptr;
	binv = nullptr;	
}


EnuBundle::~EnuBundle()
{
	if (pstart != nullptr)
		delete[] pstart;
	if (edges != nullptr)
		delete[] edges;
	if (reverse != nullptr)
		delete[] reverse;
}

