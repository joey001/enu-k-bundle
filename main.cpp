#include "EnuBundle.h"
#define FILELEN 1024
//#define TRANSFER
//#include<vld.h>
void showUsage() {
	printf("enubundle [-f filename] [-k k] [-t maxsecond]\n");
}
int main(int argc, char** argv) {
	//p2p-Gnutella04,wiki-vote.txt
	//char filepath[1024] = "D:\\Home\\benchmarks\\splex\\10th_dimacs\\jazz.graph";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\10th_dimacs\\jazz.bin";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\snap\\soc-Slashdot0902.bin";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\snap\\wiki-Vote.bin";
	char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\snap\\email-EuAll.bin";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\snap\\amazon0505.bin";
	//char filepath[1024] = "graph1.bin";
	//char filepath[1024] = "graph2.bin";
	ui k = 2;
	ui lb = 10;
	ui maxsec = 100;
	for (int i = 1; i < argc; i += 2) {
		if (argv[i][0] != '-' || argv[i][2] != 0) {
			showUsage();
			exit(0);
		}
		else if (argv[i][1] == 'f') {
			strncpy_s(filepath, argv[i + 1], FILELEN);
		}
		else if (argv[i][1] == 'k') {
			k = atoi(argv[i + 1]);
		}
		else if (argv[i][1] == 't') {
			maxsec = atoi(argv[i + 1]);
		}
	}

	EnuBundle enbundle;
	enbundle.readBinaryGraph(filepath);
	enbundle.enumPlex(k,lb,maxsec);
	//_CrtDumpMemoryLeaks();
	return 0;
}

