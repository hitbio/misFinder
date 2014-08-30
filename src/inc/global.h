/*
 * global.h
 *
 *  Created on: May 22, 2012
 *      Author: zhuxiao
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_


// =================== global file names and paths ==================
char outputPathStr[256];				// can be user-defined, default './'
char mergeSubjectsFile[256];			// must be user-defined
char mergedSegFile[256];

char inputQueryFileInit[256];			// must be user-defined
char inputBlastnFile[256];				// must be user-defined
char inputQueryFile[256];

char parseResultFile[256];
char queryMatchInfoFile[256];
char queryMatchInfoFileNew[256];
char queryStatisticsFile[256];			// "queryStatistics"
char perfectQueryFile[256];
char matchedQueryFile[256];
char disjunctQueryFile[256];
char unmatchedQueryFile[256];
char linkedQueryMatchInfoFile[256];
char sortedQueryFile[256];
char refDeletionFile[256];

char errorsFile[256];
char svFile[256];
char misUncertainFile[256];
char gapFile[256];

char newQueryFile[256];

int32_t minQueryLenThres;				// can be user-defined, default 100
int32_t shortQueryLenThres;				// can be user-defined, default 200
double matchPercentThres;				// can be user-defined, default 0.95
int64_t threadNum;						// can be user-defined, default is the number of CPU cores
int32_t indelSizeThres;					// can be user-defined, default 5 bp

char *readFilesInput[256];
int32_t readFileNum;
int32_t readsFileFormatType;
int32_t pairedMode;

// =================== global variables ========================
pthread_t *threadArr;
threadPara_t *threadParaArr;


char **segFileArr;
int64_t segFileNum;

FILE *fpParseResult, *fpPerfectQuery, *fpMatchedQuery, *fpDisjunctQuery, *fpUnmatchedQuery, *fpQueryMatch;

tmpQuery_t *tmpQueryArr;
int64_t itemNumTmpQueryArr;

queryMatchInfo_t *queryMatchInfoSet;


double matchPercentFactor;
int minTotalMatchLenThres;
int varyEndLenThres;
int endIgnoreLen;
int minDisjunctDistanceThres;
int minAlignedSegLenThres;

FILE *fpStatistics;
queryLenStatistic_t *lenStatisticArr, *lenStatisitcBufArr;
int64_t itemNumLenStatisticArr, maxItemNumLenStatisticArr;

segLinkSet_t *segLinkSet;

metrics_t *queryMetrics;

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
int colNumBaseFlagArr = 7;
//int8_t baseFlagArr[7][7];


//================================================
readSet_t *readSet;
int32_t reserveHashItemBlocksFlag;

uint64_t hashTableSizeReadseq;

readBlock_t *pReadBlockTmp;
read_t *pReadTmpDoing;

readseqBlock_t *pReadseqBlockTmp;
uint64_t *pReadseqTmpDoing;

readseqHashItemBlock_t *pReadseqHashItemBlockTmp;
readseqHashItem_t *pReadseqHashItemTmpDoing;

//==================== queryIndex ===============
queryIndex_t *queryIndex;

double SP_ratio_Thres, SMinus_ratio_Thres, SPlus_ratio_Thres, insertSize, standDev;
double discorRatio_Thres;


#endif /* GLOBAL_H_ */
