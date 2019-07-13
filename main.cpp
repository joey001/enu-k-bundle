#include "EnuBundle.h"
#include "XGetopt.h"
#define FILELEN 1024
//#define TRANSFER
//#include<vld.h>
void showUsage() {
	fprintf(stderr, "enplex [-f filename] [--small] [-k k] [-q lb] [-t maxsecond] \n");
}
int main(int argc, char** argv) {
	//p2p-Gnutella04,wiki-vote.txt
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\10th_dimacs\\jazz.bin";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\snap\\soc-Slashdot0902.bin";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\snap\\wiki-Vote.bin";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\snap\\email-EuAll.bin";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\snap\\amazon0505.bin";
	//char filepath[1024] = "D:\\Home\\benchmarks\\splex\\10th_dimacs\\celegans_metabolic.bin";
	char filepath[1024] = "D:\\Home\\benchmarks\\splex\\snap\\as-caida20071105.bin";
	//char filepath[1024] = "graph1.bin";
	//char filepath[1024] = "graph2.bin";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\2nd_dimacs\\brock200_1.bin";

	ui k = 2;
	ui lb = 1;
	uli maxsec = 600;
	ui decompose = 0;
	ui isquiete = 0;
	while (true){
		int option_index = 0;
		static struct option long_options[] =
		{
			/* Set the trim flag. */
			{"graphfile",	required_argument, 0, 'f'},
			{"k",				required_argument, 0, 'k'},
			{"maxsecond",		required_argument, 0, 't'},
			{"lowerbound",	required_argument, 0, 'l'},
			{"decompose",	no_argument, 0, 'd'},
			{"quiete", no_argument, 0, 'q'},
			{0, 0, 0, 0}
		};
		int c = getopt_long(argc, argv, "f:k:t:l:dq",
			long_options, &option_index);
		if (c == -1)
			break;
		switch (c)
		{
		case 'f':
			strncpy(filepath, optarg, FILELEN);
			break;
		case 'k':
			k = atoi(optarg);
			break;
		case 't':
			maxsec = atoi(optarg);
			break;
		case 'l':
			lb = atoi(optarg);
			break;
		case 'd':
			decompose = 1;
			break;
		case 'q':
			isquiete = 1;
			break;
		default:
			showUsage();
			exit(-1);

		}
	}	
	if (decompose && lb < 2*k-2) {
		fprintf(stderr, "lb is at least 2k-2 in decompose mode\n");
		exit(-1);
	}
	EnuBundle enbundle;
	enbundle.readBinaryGraph(filepath);
	enbundle.enumPlex(k,lb,maxsec, decompose,isquiete);
	//_CrtDumpMemoryLeaks();
	return 0;
}

