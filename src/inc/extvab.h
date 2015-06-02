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
extern char subjectsFile[256];					// must be user-defined
extern char mergedSegFile[256];

extern char inputQueryFileInit[256];				// must be user-defined
extern char inputBlastnFile[256];				// must be user-defined
extern char inputQueryFile[256];
extern char queryMatchInfoFile[256];

extern char newQueryFile[256];

extern int32_t minQueryLenThres;				// can be user-defined, default 100
extern int32_t shortQueryLenThres;				// can be user-defined, default 200
extern double matchPercentThres;				// can be user-defined, default 0.95
extern int32_t threadNum;						// can be user-defined, default is the number of CPU cores
extern int32_t indelSizeThres;					// can be user-defined, default 5 bp

extern readFile_t *readFileList;

extern char configFile[256];

// =================== global variables ========================
extern pthread_t *threadArr;
extern threadPara_t *threadParaArr;

extern double matchPercentFactor;
extern int minTotalMatchLenThres;
extern int varyEndLenThres;
extern int endIgnoreLen;
extern int minDisjunctDistanceThres;
extern int minAlignedSegLenThres;

// row -- query, column -- reference
extern int8_t baseFlagArr[7][7];
extern int32_t colNumBaseFlagArr;


//================================================
extern uint64_t hashTableSizeReadseq;

extern readBlock_t *pReadBlockTmp;
extern read_t *pReadTmpDoing;

extern readseqBlock_t *pReadseqBlockTmp;
extern uint64_t *pReadseqTmpDoing;

extern readseqHashItemBlock_t *pReadseqHashItemBlockTmp;
extern readseqHashItem_t *pReadseqHashItemTmpDoing;

extern double SP_ratio_Thres, SMinus_ratio_Thres, SPlus_ratio_Thres;
extern double discorRatio_Thres;


#endif /* EXTVAB_H_ */
