/*
 * global.h
 *
 *  Created on: May 22, 2012
 *      Author: zhuxiao
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_


// =================== global file names and paths ==================
char outputPathStr[256];				// default './'
char subjectsFile[256];
char mergedSegFile[256];

char inputQueryFileInit[256];
char inputBlastnFile[256];
char inputQueryFile[256];
char queryMatchInfoFile[256];

char newQueryFile[256];

int32_t singleCellFlag;					// default is 0. 0 -- non single-cell data; 1 -- single-cell data
int32_t minQueryLenThres;				// default 100
int32_t shortQueryLenThres;				// default 200
double matchPercentThres;				// default 0.95
int32_t threadNum;						// default is the number of CPU cores
int32_t indelSizeThres;					// default 5 base pairs

readFile_t *readFileList;

char configFile[256];

// =================== global variables ========================
pthread_t *threadArr;
threadPara_t *threadParaArr;

double matchPercentFactor;
int minTotalMatchLenThres;
int varyEndLenThres;
int endIgnoreLen;
int minDisjunctDistanceThres;
int minAlignedSegLenThres;

// row -- query, column -- reference
int8_t baseFlagArr[7][7] =
			{ // Ref.: A     C     G     T     N     .     -
					 { 1,    2,    3,    4,   18,   18,   27 }, // A
					 { 5,    1,    6,    7,   19,   19,   28 }, // C
					 { 8,    9,    1,   10,   20,   20,   29 }, // G
					 { 11,  12,   13,    1,   21,   21,   30 }, // T
					 { 14,  15,   16,   17,    1,    1,   31 }, // N
					 { 14,  15,   16,   17,    1,    1,   31 }, // .
					 { 22,  23,   24,   25,    26,  26,   -1 }  // -
			};
int32_t colNumBaseFlagArr = 7;

//================================================
uint64_t hashTableSizeReadseq;

readBlock_t *pReadBlockTmp;
read_t *pReadTmpDoing;

readseqBlock_t *pReadseqBlockTmp;
uint64_t *pReadseqTmpDoing;

readseqHashItemBlock_t *pReadseqHashItemBlockTmp;
readseqHashItem_t *pReadseqHashItemTmpDoing;


#endif /* GLOBAL_H_ */
