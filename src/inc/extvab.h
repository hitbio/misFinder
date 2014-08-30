/*
 * extvab.h
 *
 *  Created on: May 22, 2012
 *      Author: zhuxiao
 */

#ifndef EXTVAB_H_
#define EXTVAB_H_ 1

#include "constants.h"


// =================== global file names and paths ==================
extern char outputPathStr[256];					// can be user-defined, default './'
extern char mergeSubjectsFile[256];					// must be user-defined
extern char mergedSegFile[256];

extern char inputQueryFileInit[256];				// must be user-defined
extern char inputBlastnFile[256];				// must be user-defined
extern char inputQueryFile[256];

extern char parseResultFile[256];
extern char queryMatchInfoFile[256];
extern char queryMatchInfoFileNew[256];
extern char queryStatisticsFile[256];			// "queryStatistics"
extern char perfectQueryFile[256];
extern char matchedQueryFile[256];
extern char disjunctQueryFile[256];
extern char unmatchedQueryFile[256];
extern char linkedQueryMatchInfoFile[256];
extern char sortedQueryFile[256];
extern char refDeletionFile[256];

extern char errorsFile[256];
extern char svFile[256];
extern char misUncertainFile[256];
extern char gapFile[256];

extern char newQueryFile[256];

extern int32_t minQueryLenThres;				// can be user-defined, default 100
extern int32_t shortQueryLenThres;				// can be user-defined, default 200
extern double matchPercentThres;				// can be user-defined, default 0.95
extern int64_t threadNum;						// can be user-defined, default is the number of CPU cores
extern int32_t indelSizeThres;					// can be user-defined, default 5 bp

extern char *readFilesInput[256];
extern int32_t readFileNum;
extern int32_t readsFileFormatType;
extern int32_t pairedMode;

// =================== global variables ========================
extern pthread_t *threadArr;
extern threadPara_t *threadParaArr;


extern char **segFileArr;
extern int64_t segFileNum;

extern FILE *fpParseResult, *fpPerfectQuery, *fpMatchedQuery, *fpDisjunctQuery, *fpUnmatchedQuery, *fpQueryMatch;

extern tmpQuery_t *tmpQueryArr;
extern int64_t itemNumTmpQueryArr;

extern queryMatchInfo_t *queryMatchInfoSet;

extern double matchPercentFactor;
extern int minTotalMatchLenThres;
extern int varyEndLenThres;
extern int endIgnoreLen;
extern int minDisjunctDistanceThres;
extern int minAlignedSegLenThres;

extern FILE *fpStatistics;
extern queryLenStatistic_t *lenStatisticArr, *lenStatisitcBufArr;
extern int64_t itemNumLenStatisticArr, maxItemNumLenStatisticArr;

extern segLinkSet_t *segLinkSet;

// ====================================== New
extern metrics_t *queryMetrics;

// row -- query, column -- reference
extern int8_t baseFlagArr[7][7];
extern int colNumBaseFlagArr;


//================================================
extern readSet_t *readSet;
extern int32_t reserveHashItemBlocksFlag;

extern uint64_t hashTableSizeReadseq;

extern readBlock_t *pReadBlockTmp;
extern read_t *pReadTmpDoing;

extern readseqBlock_t *pReadseqBlockTmp;
extern uint64_t *pReadseqTmpDoing;

extern readseqHashItemBlock_t *pReadseqHashItemBlockTmp;
extern readseqHashItem_t *pReadseqHashItemTmpDoing;


extern queryIndex_t *queryIndex;

extern double SP_ratio_Thres, SMinus_ratio_Thres, SPlus_ratio_Thres, insertSize, standDev;
extern double discorRatio_Thres;


#endif /* EXTVAB_H_ */
