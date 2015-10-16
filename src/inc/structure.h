/*
 * structure.h
 *
 *  Created on: May 23, 2012
 *      Author: zhuxiao
 */

#ifndef STRUCTURE_H_
#define STRUCTURE_H_ 1


typedef struct ctrlLockNode
{
	int32_t readerCnt, writerCnt;
	pthread_mutex_t accessReaderCnt, accessWriterCnt;
	pthread_mutex_t writeLock, readerLock, outerLock;
}ctrlLock_t;

typedef struct
{
	int32_t threadID, validThreadNum;

	// parameters for blastn alignment
	char queryFileName[256];
	char subjectFileName[256];
	char blastnFileName[256];
	int64_t initQuerySubSum, querySubSum;

	// ====== global control variables ======
	int8_t successFlag, validFlag;

	// ====== global variables for mapping ======
	int8_t uniqueMapOpFlag;

	struct queryMatchInfoNode *queryMatchInfoSet;
	struct readSetNode *readSet;
	struct queryIndexNode *queryIndex;

	// control locks
	ctrlLock_t *lockArray;
	int32_t lockArraySize;

	// parameters for reads mapping
	int64_t startReadID, endReadID, itemNumRead;

	// parameters for mis-assembly identification, SV validation, error correction
	//int32_t startQueryID, endQueryID, itemNumQuery;
}threadPara_t;

typedef struct
{
	char *subjectSeq;
	char *subjectTitle;				// subject title
	uint32_t subjectTitleLen;
	uint32_t subjectID;				// starts from 1
	uint32_t subjectLen: 30;		// subject length
	uint32_t circularFlag: 2;		// circularFlag: YES (default) or NO
}subject_t;

typedef struct
{
	char *queryTitle;			// query title
	int32_t queryLen;
	int32_t subjectID;			// +++++++++++++++++++++++++++
	int32_t matchItemNum;
	int32_t firstRow;			// new added
}tmpQuery_t;

typedef struct
{
	int32_t startSubPos;
	int32_t endSubPos;
	int32_t startQueryPos;
	int32_t endQueryPos;
	int32_t strand;
	int32_t matchLen;
	int32_t totalMatchLen;
	int32_t gapNum;
	double matchPercent;
}matchItem_t, validSegment_t;

typedef struct querySubjectNode
{
	int32_t subjectID;
	int32_t matchItemNum;
	int32_t firstRow;
	int32_t validSegmentNum: 26;
	int32_t matchKind: 3;
	int32_t circularFlag: 3;
	validSegment_t *validSegArray;
}querySubject_t;

// alignment information for inter-chromosomes
typedef struct globalValidSegNode
{
	int32_t subjectID;
	int32_t startSubPos;
	int32_t endSubPos;
	int32_t startQueryPos;
	int32_t endQueryPos;
	int32_t strand;
	int32_t matchLen;
	int32_t totalMatchLen;
	int32_t gapNum;
	double matchPercent;
}globalValidSeg_t;

// queryReadNode
typedef struct queryReadNode
{
	uint32_t readID;				// read ID
	uint16_t seqlen;
	uint16_t startReadPos;
	uint16_t orientation;
	uint16_t alignSize;
	uint32_t queryPos;			// the query position
}queryRead_t;

typedef struct queryReadSetNode
{
	queryRead_t *queryReadArray;
	int32_t queryReadNum, setID;	// setID: starts from 1
}queryReadSet_t;

// new querySeq
typedef struct querySeqNode
{
	char *querySeq;
	int32_t queryLen;
	int32_t startQueryPos;
	int32_t endQueryPos;
	struct querySeqNode *next;
}querySeq_t;

// query node
typedef struct
{
	int32_t queryID;					// query ID
	int32_t queryLen;
	char *querySeq;
	char *queryTitle;					// query title
	querySubject_t *querySubArray;		// +++++++++++++++++++++++++++
	int32_t querySubjectNum;
	int32_t bestMatchRow;				// the row of the best match querySubject item array
	int32_t globalMatchKind : 10;
	int32_t circularFlag: 4;			// circular flag
	int32_t misassFlag: 6;
	int32_t queryTitleLen: 12;

	globalValidSeg_t *globalValidSegArray;
	int32_t globalValidSegNum;

	querySeq_t *newQuerySeqList, *tailNewQuerySeqList;
	int32_t newQuerySeqNum;

	queryReadSet_t *queryReadSetArray;
	int32_t queryReadSetNum;

	uint64_t *covFlagArray;

	struct misInfoNode *misInfoList, *tailMisInfo;
	int32_t misInfoItemNum;

	int32_t subRegSize;
	double meanRegCov, sdevRegCov;
}query_t;

typedef struct queryMatchInfoNode
{
	query_t *queryArray;
	subject_t *subjectArray;
	matchItem_t *matchItemArray;
	int64_t itemNumQueryArray, itemNumSubjectArray, itemNumMatchItemArray;
	int32_t potentMisassNum, trueMisassNum, SVNum;
	int64_t maxQueryID, secQueryID;
}queryMatchInfo_t;

typedef struct segmentLinkNode
{
	int32_t arrRow;		// the row of matchItemArr
	int32_t addedOrder; // starts from 1
	int16_t validFlag, usedFlag;
	int32_t previous;
	int32_t next;
}segmentLink_t;

typedef struct segLinkSetNode
{
	segmentLink_t *linkArray;
	int32_t headRow, tailRow;
	int32_t itemNum, maxArraySize;
}segLinkSet_t;

typedef struct lenStatisticNode
{
	int32_t queryID;
	int32_t queryLen;
}queryLenStatistic_t;

//==================================== New
typedef struct lengthMetricsNode
{
	int64_t totalLen;
	int32_t totalNum;
	int32_t maxSize;
	int32_t N50;
	int32_t medianSize;
	double meanSize;
	double lengthRatio;					// lengthRatio = totalLen / totalRefBaseNum;
	double coveredRatio;				// coveredRatio = coveredBaseNum / totalRefBaseNum;
	int64_t totalRefBaseNum;
	int64_t coveredBaseNum;
	double GCRatio;
}lengthMetrics_t;

typedef struct accuracyMetricsNode
{
	int32_t totalNum;
	int64_t totalLen;
	double meanSize;
	double lengthRatio;					// lengthRatio = totalLen / totalRefBaseNum;
}accuracyMetrics_t;

typedef struct baseBitArrayNode
{
	uint8_t *subjectBitArray;
	//uint8_t *subjectBitArrayDeletion;
	int32_t totalRefBaseNum;
	int32_t coveredBaseNum;
	double coveredRatio;
}baseBitArray_t;

typedef struct metricsNode
{
	lengthMetrics_t lengthMetrics;
	accuracyMetrics_t accuracyMetrics[4];
	baseBitArray_t *baseBitArrays;
	int32_t subjectNum;
}metrics_t;


//======================================================
//read
typedef struct readNode
{
	uint32_t rowReadseqInBlock;		// point to the row of readseqBlock.readseqArr
	uint16_t readseqBlockID;
	uint16_t nBaseNum;		// unknown base count: it will be changed to 'C' automatically
	uint16_t seqlen;
	uint16_t validFlag: 2;
	uint16_t successMapFlag: 2;
	uint16_t uniqueMapFlag: 2;
}read_t;

// read block
typedef struct readBlockNode
{
	uint32_t blockID;
	uint32_t itemNum;
	read_t *readArr;
}readBlock_t;

// readseq block
typedef struct readseqBlockNode
{
	uint32_t blockID;
	uint32_t rowsNum;
	uint64_t *readseqArr;
}readseqBlock_t;

//readseq hash node
typedef struct readseqHashItemNode
{
	uint32_t rowReadseqInBlock;		// point to the row of readseqBlock.readseqArr
	uint16_t readseqBlockID;
	uint16_t seqlen;
	uint32_t nextHashItemBlockID;
	uint32_t nextItemRowHashItemBlock;
}readseqHashItem_t;

// readseq hash block
typedef struct readseqHashItemBlockNode
{
	uint32_t blockID;
	uint32_t itemNum;
	readseqHashItem_t *readseqHashItemArr;
}readseqHashItemBlock_t;

//readseq hash node
typedef struct readseqHashBucketNode
{
	uint32_t hashItemBlockID;
	uint32_t itemRowHashItemBlock;
}readseqHashBucket_t;

// readMatchInfo
typedef struct readMatchInfoNode
{
	uint32_t queryID;
	uint32_t queryPos;			// the query position
	uint32_t readID;
	uint16_t seqlen;
	uint16_t alignSize;
	uint16_t startReadPos;
	uint16_t readOrientation;
}readMatchInfo_t;

// read block
typedef struct readMatchInfoBlockNode
{
	uint32_t blockID;
	uint32_t itemNum;
	readMatchInfo_t *readMatchInfoArr;
}readMatchInfoBlock_t;

//readSet
typedef struct readSetNode
{
	// read blocks
	readBlock_t *readBlockArr;		// point to kmerSeqBlock array
	int16_t blocksNumRead;
	int16_t maxBlocksNumRead;
	int16_t bytesPerRead;
	int16_t maxReadLen;
	int64_t totalItemNumRead;
	int64_t totalValidItemNumRead;
	int32_t maxItemNumPerReadBlock;

	// readseq blocks
	readseqBlock_t *readseqBlockArr;		// point to kmerSeqBlock array
	int32_t blocksNumReadseq;
	int16_t maxBlocksNumReadseq;
	int16_t bytesPerEntryReadseq;
	int64_t totalItemNumReadseq;
	int32_t maxEntryNumReadseqBlock;

	// readseq hash table
	readseqHashBucket_t *readseqHashtable;
	int32_t hashTableSizeReadseq;

	// readseq hash item blocks
	readseqHashItemBlock_t *readseqHashItemBlockArr;		// point to readseqHashBlock array
	int32_t blocksNumReadseqHashItem;
	int16_t maxBlocksNumReadseqHashItem;
	int16_t bytesPerReadseqHashItem;
	int64_t totalItemNumReadseqHashItem;
	int32_t maxItemNumPerReadseqHashItemBlock;

	// read blocks
	readMatchInfoBlock_t *readMatchInfoBlockArr;		// point to kmerSeqBlock array
	int16_t blocksNumReadMatchInfo;
	int16_t bytesPerReadMatchInfo;
	int64_t totalItemNumReadMatchInfo;
	int64_t totalValidItemNumReadMatchInfo;
	int32_t maxItemNumPerReadMatchInfoBlock;

	int8_t setID;								// starts from 1
	int8_t readType, pairedMode;
	double insertSize, standDev;

}readSet_t;

// readFile
typedef struct readFileNode
{
	char **readFiles;
	int8_t readFileNum;
	int8_t readsFileFormatType;
	int8_t readType;
	int8_t pairedMode;
	//int8_t shortReadsFlag;
	struct readFileNode *next;
}readFile_t;

//readSet array
typedef struct readSetArrNode
{
	readSet_t *readSetArray;
	int32_t readSetNum;
}readSetArr_t;

typedef struct readBufNode
{
	char *seq;
	char *qual;
	int16_t len, nBaseNum, validFlag;
	uint64_t *pReadseqInt, hashcode;
	readseqHashItem_t *pReadseqHashItem;
}readBuf_t;

//====================== structures for queryIndex =========
//kmer hash node
typedef struct kmerHashBucketNode
{
	uint32_t kmerBlockID;
	uint32_t itemRowKmerBlock;
}kmerHashBucket_t;

//queryPosNode
typedef struct queryPosNode
{
	int32_t queryID;
	int32_t queryPos;
}queryPos_t;

//queryKmer
typedef struct queryKmerNode
{
	queryPos_t *ppos;
	uint32_t multiplicity;
	uint32_t arraysize;
	uint32_t kmerseqBlockID;
	uint32_t itemRowKmerseqBlock;
	uint32_t nextKmerBlockID;
	uint32_t nextItemRowKmerBlock;
}queryKmer_t;

// queryKmer block
typedef struct queryKmerBlockNode
{
	uint32_t blockID;
	uint32_t itemNum;
	queryKmer_t *kmerArr;
}queryKmerBlock_t;

// kmerseq block
typedef struct kmerseqBlockNode
{
	uint32_t blockID;
	uint32_t itemNum;
	uint64_t *kmerSeqArr;
}kmerseqBlock_t;

// queryPos block
typedef struct queryPosBlockNode
{
	uint32_t blockID;
	uint32_t itemNum;			// unused
	queryPos_t *queryPosArr;
}queryPosBlock_t;

// queryIndex
typedef struct queryIndexNode
{
	// k-mer blocks
	kmerHashBucket_t *kmerHashtable;
	int32_t hashTableSize;
	int32_t kmerSize;

	queryKmerBlock_t *kmerBlockArr;				// point to queryKmerBlock array
	int32_t blocksNumKmer;
	int16_t maxBlocksNumKmer;
	int16_t bytesPerKmer;
	int64_t totalItemNumKmer;
	int32_t maxItemNumPerKmerBlock;

	// kmerseq blocks
	kmerseqBlock_t *kmerSeqBlockArr;		// point to kmerSeqBlock array
	int32_t blocksNumKmerSeq;
	int16_t maxBlocksNumKmerSeq;
	int8_t bytesPerKmerseq;
	int8_t entriesPerKmer;
	int64_t totalItemNumKmerSeq;
	int64_t maxItemNumPerKmerSeqBlock;
	uint64_t lastEntryMask;
	int32_t lastEntryBaseNum;		// the base number of the last entry of its sequence array

	// ridpos blocks
	queryPosBlock_t *queryPosBlockArr;		// point to queryPosBlock array
	int16_t blocksNumQuerypos;
	int16_t maxBlocksNumQuerypos;
	int16_t bytesPerQuerypos;
	int64_t totalItemNumQuerypos;
	int64_t maxItemNumPerQueryposBlock;

	uint64_t *pKmerSeqAdding;
}queryIndex_t;


typedef struct alignMatchItemNode
{
	int32_t queryID;
	int32_t queryPos;
	int16_t mismatchNum;
	int16_t pairRow;
	int8_t orientation;
	int8_t validFlag;
	int16_t headSkipNum, tailSkipNum;
	int16_t startReadPos;  // start with 1
	int16_t alignSize;
	int16_t fragSize;
}alignMatchItem_t;


typedef struct baseCovNode
{
	int32_t baseNumArray[6];	// [0-4]: A,C,G,T,N; [5]: total
}baseCov_t;


typedef struct queryMarginNode
{
	int32_t leftMargin, rightMargin, misassFlag, queryEndFlag;
	int32_t disagreeNum, zeroCovNum, discorNum, highCovRegNum, lowCovRegNum, disagreeRegSize;
	double discorRatio, multiReadsRatio, ratioMateNumOtherQuery;
}queryMargin_t;


typedef struct ratioRegionNode
{
	int32_t disagreeNum, zeroCovNum;
	int32_t startQPosLHalf, endQPosLHalf, startQPosRHalf, endQPosRHalf, midQPos;

	int32_t pairedNum, singleNum, singleNumLeftHalf, singleNumRightHalf, discorNum;
	int32_t multiMapReadsNum, totalReadsNum, mateNumOtherQuery;
	double discorRatio, multiReadsRatio, ratioMateNumOtherQuery;
}ratioRegion_t;


typedef struct queryIndelNode
{
	int32_t leftSegRow, rightSegRow;
	int16_t queryIndelKind, misassFlag;
	int32_t leftMargin, rightMargin, subjectID, leftMarginSubject, rightMarginSubject;
	int32_t pairNumLeft, discorNumLeft, pairNumRight, discorNumRight;
	int32_t disagreeNum, zeroCovNum, difSubject, difQuery, highCovRegNum, lowCovRegNum, disagreeRegSize;
	double averFragSizeLeft, averFragSizeRight, difFragSizeLeft, difFragSizeRight, discorRatioLeft, discorRatioRight;
	int32_t disagreeNumInRatioRegion, zeroCovNumInRatioRegion, discorRegNumInRatioRegion, multiReadsRegNumInRatioRegion;
}queryIndel_t;


typedef struct misassSeqNode
{
	int32_t startQueryPos, endQueryPos, misassKind;	// misasskind: ERR_MISJOIN, ERR_INSERT, ERR_DEL, COR_SV, ERR_UNCERTAIN
	int32_t subjectID, startSubjectPos, endSubjectPos;
	struct misassSeqNode *next;
}misassSeq_t;


typedef struct misInfoNode
{
	int32_t leftSegRow, rightSegRow;
	int8_t misType, misassFlag, gapFlag, innerFlag;
	queryMargin_t *queryMargin;
	queryIndel_t *queryIndel;
	misassSeq_t *misassSeqList, *tailMisassSeq;
	int32_t misassSeqNum;
	struct misInfoNode *next;
}misInfo_t;


#endif /* STRUCTURE_H_ */
