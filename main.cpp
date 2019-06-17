#include "EnuBundle.h"
#define FILELEN 1024
//#define TRANSFER


int fileSuffixPos(char* filepath) {
	int j = strlen(filepath)-1;
	while (filepath[j] != '.')
		j--;
	return j + 1;
}

#ifdef TRANSFER
int transferFile(char* txtFilePath, char* binFilePath){
	EnuBundle enbdl;

	int pos = fileSuffixPos(txtFilePath);
	if (strcmp(txtFilePath + pos, "bin") == 0) {

		printf("Bin file\n");
		return 0;
	}
	
	if (strcmp(txtFilePath + pos, "graph") == 0) {//dimacs 10
		enbdl.readRawDIM10Text(txtFilePath);
	}
	else if (strcmp(txtFilePath + pos, "txt") == 0) { //snap
		enbdl.readRawSNAPText(txtFilePath);
	}
	else {
		printf("File not support\n");
		return 0;
	}
	strncpy_s(binFilePath,FILELEN, txtFilePath, strlen(txtFilePath)+1);
	binFilePath[pos++] = 'b';
	binFilePath[pos++] = 'i';
	binFilePath[pos++] = 'n';
	binFilePath[pos++] = '\0';
	enbdl.writeBinaryGraph(binFilePath);
	return 1;
}
int main(int argc, char** argv) {
	char filepath[FILELEN]="\0";
	char binFile[FILELEN];
	if (argc < 2) {
		printf("text2bin textfile\n");
	}
	else {
		strncpy_s(filepath, argv[1], FILELEN);
		transferFile(filepath, binFile);
		printf("File %s -> %s\n", filepath, binFile);
	}
#else

void showUsage() {
	printf("enubundle [-f filename] [-k k] [-t maxsecond]\n");
}
int main(int argc, char** argv) {
	//p2p-Gnutella04,wiki-vote.txt
	//char filepath[1024] = "D:\\Home\\benchmarks\\splex\\10th_dimacs\\jazz.graph";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\10th_dimacs\\jazz.bin";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\snap\\amazon0505.bin";
	char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\snap\\wiki-Vote.bin";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\snap\\email-EuAll.bin";
	//char filepath[1024] = "graph1.bin";
	//char filepath[1024] = "graph2.bin";
	ui k = 0;
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
	enbundle.enumPlex(2,20);
#endif	

	return 0;
}

