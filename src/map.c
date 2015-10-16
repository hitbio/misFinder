/*
 * map.c
 *
 *  Created on: May 30, 2014
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Map the reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mapReads(queryMatchInfo_t *queryMatchInfoSet, readSetArr_t *readSetArray, queryIndex_t *queryIndex, int32_t threadNum)
{
	struct timeval tpstart, tpend;
	double timeused_align;
	gettimeofday(&tpstart, NULL);

	int32_t i;

	for(i=0; i<readSetArray->readSetNum; i++)
	{
		if(mapReadsSingleReadset(queryMatchInfoSet, readSetArray->readSetArray+i, queryIndex, threadNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot align reads to queries, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// calculate aligning time
	gettimeofday(&tpend, NULL);
	timeused_align = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Align reads used time: %.2f seconds.\n", timeused_align);

	return SUCCESSFUL;
}

/**
 * Map the reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mapReadsSingleReadset(queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet, queryIndex_t *queryIndex, int32_t threadNum)
{
	// initialize read blocks
	if(initReadMatchInfoBlockInReadset(readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize read blocks, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize the coverage flag array
	if(initCovFlagArray(queryMatchInfoSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize coverage flag array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize the data for threads
	if(initThreadParasMap(&threadArr, &threadParaArr, threadNum, queryMatchInfoSet, readSet, queryIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize map threads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// map reads uniquely
	if(mapReadsOp(threadArr, threadParaArr, threadNum, YES)==FAILED)
	{
		printf("line=%d, In %s(), cannot map reads uniquely, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// map reads non-uniquely
	if(mapReadsOp(threadArr, threadParaArr, threadNum, NO)==FAILED)
	{
		printf("line=%d, In %s(), cannot map reads uniquely, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// free memory for thread parameters
	if(freeThreadParasMap(&threadArr, &threadParaArr, threadNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot free the thread parameters, error.\n", __LINE__, __func__);
		return FAILED;
	}

	// free coverage flag array
	freeCovFlagArray(queryMatchInfoSet);

	return SUCCESSFUL;
}

/**
 * Initialize the coverage flag array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initCovFlagArray(queryMatchInfo_t *queryMatchInfoSet)
{
	int32_t i, itemNum;
	query_t *queryArray;

	queryArray = queryMatchInfoSet->queryArray;
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		itemNum = (queryArray[i].queryLen - 1) / 64 + 1;
		queryArray[i].covFlagArray = (uint64_t*) calloc (itemNum, sizeof(uint64_t));
		if(queryArray[i].covFlagArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Free the coverage flag array.
 */
void freeCovFlagArray(queryMatchInfo_t *queryMatchInfoSet)
{
	int32_t i;
	query_t *queryArray;

	queryArray = queryMatchInfoSet->queryArray;
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		free(queryArray[i].covFlagArray);
		queryArray[i].covFlagArray = NULL;
	}
}

/**
 * Initialize the thread parameters for reads mapping.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short initThreadParasMap(pthread_t **threadArray, threadPara_t **threadParaArray, int32_t threadNum, queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet, queryIndex_t *queryIndex)
{
	int32_t i, subItemNum, totalReadsNum, total;

	*threadArray = (pthread_t*) calloc (threadNum, sizeof(pthread_t));
	if((*threadArray)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error.\n", __LINE__, __func__);
		return FAILED;
	}

	*threadParaArray = (threadPara_t*) calloc (threadNum, sizeof(threadPara_t));
	if((*threadParaArray)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error.\n", __LINE__, __func__);
		return FAILED;
	}

	totalReadsNum = readSet->totalItemNumRead;
	subItemNum = (totalReadsNum - 1) / threadNum + 1;
	if(subItemNum%2==1)
		subItemNum ++;

	total = 0;
	for(i=0; i<threadNum; i++)
	{
		if(total<totalReadsNum)
		{
			(*threadParaArray)[i].startReadID = total + 1;
			if(total+subItemNum<=totalReadsNum)
				(*threadParaArray)[i].endReadID = total + subItemNum;
			else
				(*threadParaArray)[i].endReadID = totalReadsNum;
			(*threadParaArray)[i].itemNumRead = (*threadParaArray)[i].endReadID - (*threadParaArray)[i].startReadID + 1;

			(*threadParaArray)[i].queryMatchInfoSet = queryMatchInfoSet;
			(*threadParaArray)[i].readSet = readSet;
			(*threadParaArray)[i].queryIndex = queryIndex;

			(*threadParaArray)[i].validFlag = YES;

			total += subItemNum;
		}else
			(*threadParaArray)[i].validFlag = NO;

		(*threadParaArray)[i].successFlag = NO;
	}

	return SUCCESSFUL;
}

/**
 * Free the thread parameters for reads mapping.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short freeThreadParasMap(pthread_t **threadArray, threadPara_t **threadParaArray, int32_t threadNum)
{
	int32_t i;

	for(i=0; i<threadNum; i++)
	{
		if((*threadParaArray)[i].validFlag==YES)
		{
			(*threadParaArray)[i].queryMatchInfoSet = NULL;
			(*threadParaArray)[i].readSet = NULL;
			(*threadParaArray)[i].queryIndex = NULL;

			(*threadParaArray)[i].lockArray = NULL;
			(*threadParaArray)[i].lockArraySize = 0;
		}
	}

	free(*threadArray);
	*threadArray = NULL;

	free(*threadParaArray);
	*threadParaArray = NULL;

	return SUCCESSFUL;
}

/**
 * Start the reads map.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mapReadsOp(pthread_t *threadArray, threadPara_t *threadParaArray, int32_t threadNum, int32_t uniqueMapOpFlag)
{
	readSet_t *readSet;

	readSet = threadParaArray[0].readSet;

	// initialize mutex
	if(initMutexMemMap(threadParaArray, threadNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize memory thread mutex for reads mapping, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// start multiple threads
	if(createThreadsMap(threadArray, threadParaArray, threadNum, uniqueMapOpFlag)==FAILED)
	{
		printf("line=%d, In %s(), cannot create threads, error.\n", __LINE__, __func__);
		return FAILED;
	}

	if(waitThreads(threadArray, threadParaArray, threadNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot wait threads, error.\n", __LINE__, __func__);
		return FAILED;
	}

	// free mutex
	if(freeMutexMemMap(threadParaArray, threadNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot free the memory for thread mutex for reads mapping, error.\n", __LINE__, __func__);
		return FAILED;
	}

	// calculate the total aligned reads
	if(uniqueMapOpFlag==NO)
	{
		if(computeTotalAlignedReadsNum(readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the total aligned reads count, error.\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Initialize the mutex memory for threads for reads mapping.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short initMutexMemMap(threadPara_t *threadParaArray, int32_t threadNum)
{
	ctrlLock_t *lockArray;
	int32_t i, lockArraySize;

	// initialize mutex
	lockArraySize = threadParaArray[0].queryMatchInfoSet->itemNumQueryArray;
	lockArray = (ctrlLock_t*) calloc (lockArraySize, sizeof(ctrlLock_t));
	if(lockArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<lockArraySize; i++)
	{
		lockArray[i].readerCnt = 0;
		lockArray[i].writerCnt = 0;

		pthread_mutex_init(&lockArray[i].accessReaderCnt, NULL);
		pthread_mutex_init(&lockArray[i].accessWriterCnt, NULL);
		pthread_mutex_init(&lockArray[i].writeLock, NULL);
		pthread_mutex_init(&lockArray[i].readerLock, NULL);
		pthread_mutex_init(&lockArray[i].outerLock, NULL);
	}

	for(i=0; i<threadNum; i++)
	{
		threadParaArray[i].lockArray = lockArray;
		threadParaArray[i].lockArraySize = lockArraySize;
	}

	return SUCCESSFUL;
}

/**
 * Free the mutex memory for threads for reads mapping.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short freeMutexMemMap(threadPara_t *threadParaArray, int32_t threadNum)
{
	ctrlLock_t *lockArray;
	int32_t i, lockArraySize;

	lockArray = threadParaArray[0].lockArray;
	lockArraySize = threadParaArray[0].lockArraySize;

	for(i=0; i<lockArraySize; i++)
	{
		lockArray[i].readerCnt = 0;
		lockArray[i].writerCnt = 0;

		pthread_mutex_init(&lockArray[i].accessReaderCnt, NULL);
		pthread_mutex_init(&lockArray[i].accessWriterCnt, NULL);
		pthread_mutex_init(&lockArray[i].writeLock, NULL);
		pthread_mutex_init(&lockArray[i].readerLock, NULL);
		pthread_mutex_init(&lockArray[i].outerLock, NULL);
	}

	for(i=0; i<threadNum; i++)
	{
		threadParaArray[0].lockArray = NULL;
		threadParaArray[0].lockArraySize = 0;
	}
	free(lockArray);

	return SUCCESSFUL;
}

/**
 * Compute total aligned reads count.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short computeTotalAlignedReadsNum(readSet_t *readSet)
{
	int32_t i, j;
	readBlock_t *readBlock;
	readMatchInfoBlock_t *readMatchInfoBlock;
	readMatchInfo_t *pReadMatchInfo;

	readSet->totalValidItemNumReadMatchInfo = 0;
	for(i=0; i<readSet->blocksNumRead; i++)
	{
		readBlock = readSet->readBlockArr + i;
		readMatchInfoBlock = readSet->readMatchInfoBlockArr + i;
		for(j=0; j<readBlock->itemNum; j++)
		{
			pReadMatchInfo = readMatchInfoBlock->readMatchInfoArr + j;
			if(pReadMatchInfo->queryID>0)
				readSet->totalValidItemNumReadMatchInfo ++;
		}
	}

	return SUCCESSFUL;
}

/**
 * Create the threads for reads mapping.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short createThreadsMap(pthread_t *threadArray, threadPara_t *threadParaArray, int32_t threadNum, int32_t uniqueMapOpFlag)
{
	int32_t i, ret, validNum;

	if(uniqueMapOpFlag==YES)
	{
		validNum = 0;
		for(i=0; i<threadNum; i++)
			if(threadParaArray[i].validFlag==YES)
				validNum ++;
		printf("The reads alignment will be performed using %d threads.\n", validNum);
	}

	if(uniqueMapOpFlag==YES)
		printf("Aligning reads uniquely ...\n");
	else
		printf("Aligning reads of multiple occurrences ...\n");

	for(i=0; i<threadNum; i++)
	{
		if(threadParaArray[i].validFlag==YES)
		{
			threadParaArray[i].uniqueMapOpFlag = uniqueMapOpFlag;
			ret = pthread_create(threadArray+i, NULL, (void  *) mapReadsOpSingleThread, threadParaArray+i);
			if(ret!=0)
			{
				printf("line=%d, In %s(), cannot create threads for reads mapping, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Start the reads map using single thread.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
void mapReadsOpSingleThread(threadPara_t *threadPara)
{
	int32_t i, j, k, seqLen, maxArraySize, *matchItemNum, matchItemNums[2], autoMismatchThreshold;
	readBlock_t *pReadBlockArray;
	read_t *pRead, *pReadArray;
	uint64_t rid, tmpRid, *readSeqInt, *readSeqIntRev;
	alignMatchItem_t *matchResultArray, *matchResultArrays[2], *matchResultArrayBuf1, *matchResultArrayBuf2;
	char readseqTmp[MAX_READ_LEN_IN_BUF+1];
	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo;
	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;

	queryMatchInfo_t *queryMatchInfoSet;
	readSet_t *readSet;
	queryIndex_t *queryIndex;
	double insertSize, standDev;
	int32_t uniqueMapOpFlag, matedNum, uniqueMapFlagArray[2];

	int64_t startReadID, endReadID;
	int32_t startReadBlockID, endReadBlockID, startRowInReadBlock, endRowInReadBlock;

	// get parameters from threadPara
	queryMatchInfoSet = threadPara->queryMatchInfoSet;
	readSet = threadPara->readSet;
	queryIndex = threadPara->queryIndex;
	uniqueMapOpFlag = threadPara->uniqueMapOpFlag;
	insertSize = readSet->insertSize;
	standDev = readSet->standDev;
	startReadID = threadPara->startReadID;
	endReadID = threadPara->endReadID;

	// the match result array
	if(getMaxArraySizeFromQueryIndex(&maxArraySize, queryIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the maximal array size from query index, error!\n", __LINE__, __func__);
		return;
	}
	if(maxArraySize<=0)
	{
		printf("line=%d, In %s(), invalid array size %d, error!\n", __LINE__, __func__, maxArraySize);
		return;
	}

	matchResultArrays[0] = (alignMatchItem_t*) calloc (4*queryIndex->kmerSize*maxArraySize, sizeof(alignMatchItem_t));
	matchResultArrays[1] = (alignMatchItem_t*) calloc (4*queryIndex->kmerSize*maxArraySize, sizeof(alignMatchItem_t));
	matchResultArrayBuf1 = (alignMatchItem_t*) calloc (4*queryIndex->kmerSize*maxArraySize, sizeof(alignMatchItem_t));
	matchResultArrayBuf2 = (alignMatchItem_t*) calloc (4*queryIndex->kmerSize*maxArraySize, sizeof(alignMatchItem_t));
	if(matchResultArrays[0]==NULL || matchResultArrays[1]==NULL || matchResultArrayBuf1==NULL || matchResultArrayBuf2==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return;
	}


	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

	// map reads
	startReadBlockID = (startReadID - 1) / readSet->maxItemNumPerReadBlock + 1;
	endReadBlockID = (endReadID - 1) / readSet->maxItemNumPerReadBlock + 1;

	rid = startReadID;
	for(i=startReadBlockID-1; i<=endReadBlockID-1; i++)
	{
		pReadBlockArray = readSet->readBlockArr + i;
		pReadArray = pReadBlockArray->readArr;

		startRowInReadBlock = (rid - 1) % readSet->maxItemNumPerReadBlock;
		if(i==endReadBlockID-1)
			endRowInReadBlock = (endReadID - 1) % readSet->maxItemNumPerReadBlock;
		else
			endRowInReadBlock = pReadBlockArray->itemNum - 1;

		if(endRowInReadBlock>pReadBlockArray->itemNum-1)
		{
			printf("line=%d, In %s(), endRowInReadBlock=%d, lastRow=%d, cannot allocate memory, error!\n", __LINE__, __func__, endRowInReadBlock, pReadBlockArray->itemNum-1);
			return;
		}

		for(j=startRowInReadBlock; j<=endRowInReadBlock; j++, rid++)
		{
			pRead = pReadArray + j;
			seqLen = pRead->seqlen;
			autoMismatchThreshold = seqLen * 0.05;

			//if(rid==31722926)
			//{
			//	printf("rid=%ld, seqLen=%d\n", rid, seqLen);
			//}

			if(rid%2==1)
			{
				matchResultArray = matchResultArrays[0];
				matchItemNum = matchItemNums;
			}else
			{
				matchResultArray = matchResultArrays[1];
				matchItemNum = matchItemNums + 1;
			}
			*matchItemNum = 0;

			if(pRead->validFlag==YES && (uniqueMapOpFlag==YES || (uniqueMapOpFlag==NO && pRead->successMapFlag==NO)))
			{
				// map single read
				if(mapSingleReadToQueries(rid, pRead, readSet, matchResultArray, matchResultArrayBuf1, matchResultArrayBuf2, matchItemNum, queryIndex, queryMatchInfoSet->queryArray, autoMismatchThreshold)==FAILED)
				{
					printf("line=%d, In %s(), cannot map the read %lu, error!\n", __LINE__, __func__, rid);
					return;
				}
			}

			if(rid%2==0)
			{
				if(getMatedNumAlignResult(&matedNum, matchResultArrays[0], matchResultArrays[1], matchItemNums[0], matchItemNums[1], (pRead-1)->seqlen, pRead->seqlen, insertSize, standDev)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the match information of paired-ends, error!\n", __LINE__, __func__);
					return;
				}

				if(matedNum==0 && (pRead-1)->successMapFlag==NO && pRead->successMapFlag==NO)
				{
					// get more align items from unmated paired ends
					if(getMoreAlignItemsFromUnmatedPE(&matedNum, rid-1, rid, matchResultArrays[0], matchResultArrays[1], matchItemNums, matchItemNums+1, 3, queryMatchInfoSet, readSet)==FAILED)
					{
						printf("line=%d, In %s(), cannot get more align items from unmated paired ends, error!\n", __LINE__, __func__);
						return;
					}
				}

				// compute unique map flag
				if(computeUniqueMapFlags(uniqueMapFlagArray, matedNum, matchResultArrays[0], matchResultArrays[1], matchItemNums[0], matchItemNums[1])==FAILED)
				{
					printf("line=%d, In %s(), cannot compute unique map flags for ends, error!\n", __LINE__, __func__);
					return;
				}

				// select the mated match information
				if(matchItemNums[0]>0 && matchItemNums[1]>0 && (matchItemNums[0]>=2 || matchItemNums[1]>=2))
				{
					if(selectMatedMatchPosPE(matchResultArrays[0], matchResultArrays[1], matchItemNums, matchItemNums+1, matedNum, (pRead-1)->seqlen, pRead->seqlen, insertSize, standDev, uniqueMapOpFlag, queryMatchInfoSet, threadPara)==FAILED)
					{
						printf("line=%d, In %s(), cannot adjust the match information of paired-ends, error!\n", __LINE__, __func__);
						return;
					}
				}else if((matchItemNums[0]>=2 && matchItemNums[1]==0) || (matchItemNums[0]==0 && matchItemNums[1]>=2))
				{
					// select best alignment for single read
					if(selectBestAlignInfoSingleRead(matchResultArrays[0], matchItemNums, uniqueMapOpFlag, queryMatchInfoSet, threadPara)==FAILED)
					{
						printf("line=%d, In %s(), cannot select the best alignment information for single read, error!\n", __LINE__, __func__);
						return;
					}

					// select best alignment for single read
					if(selectBestAlignInfoSingleRead(matchResultArrays[1], matchItemNums+1, uniqueMapOpFlag, queryMatchInfoSet, threadPara)==FAILED)
					{
						printf("line=%d, In %s(), cannot select the best alignment information for single read, error!\n", __LINE__, __func__);
						return;
					}
				}

				// save match information
				tmpRid = rid - 1;
				for(k=0; k<2; k++, tmpRid++)
				{
					if(k==0)
					{
						pRead --;
						matchResultArray = matchResultArrays[0];
					}else
					{
						pRead ++;
						matchResultArray = matchResultArrays[1];
					}

					if(matchItemNums[k]==1)
					{
						// save the match result
						readMatchInfoBlockID = (tmpRid - 1) / maxItemNumPerReadMatchInfoBlock;
						rowNumInReadMatchInfoBlock = (tmpRid - 1) % maxItemNumPerReadMatchInfoBlock;
						pReadMatchInfo = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;
						readMatchInfoBlockArr[readMatchInfoBlockID].itemNum ++;

						pReadMatchInfo->queryID = matchResultArray[0].queryID;
						pReadMatchInfo->queryPos = matchResultArray[0].queryPos;
						pReadMatchInfo->alignSize = matchResultArray[0].alignSize;
						pReadMatchInfo->startReadPos = matchResultArray[0].startReadPos;
						pReadMatchInfo->readID = tmpRid;
						pReadMatchInfo->seqlen = pRead->seqlen;
						pReadMatchInfo->readOrientation = matchResultArray[0].orientation;

						pRead->successMapFlag = YES;
						pRead->uniqueMapFlag = uniqueMapFlagArray[k];

						// update coverage flag
						if(updateQueryCovFlag(queryMatchInfoSet, matchResultArray, threadPara)==FAILED)
						{
							printf("line=%d, In %s(), cannot update query coverage flag, error!\n", __LINE__, __func__);
							return;
						}

/*
						// output the match result to file
						if(pReadMatchInfo->queryID==52 || pReadMatchInfo->queryID==1046)
						{
							readSeqInt = readSet->readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
							seqLen = pRead->seqlen;

							if(matchResultArray[0].orientation==ORIENTATION_PLUS)
							{
								getReadBaseFromPosByInt(readseqTmp, readSeqInt, seqLen, 0, seqLen);
								printf("%ld\t%d\t%d\t+\t%d\t%s\n", rid, pReadMatchInfo->queryID, pReadMatchInfo->queryPos, seqLen, readseqTmp);
							}
							else
							{
								readSeqIntRev = (uint64_t*) calloc (((seqLen-1)>>5)+1, sizeof(uint64_t));
								if(readSeqIntRev==NULL)
								{
									printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
									return FAILED;
								}

								// generate the reverse complements
								if(getReverseReadseqInt(readSeqIntRev, readSeqInt, seqLen)==FAILED)
								{
									printf("line=%d, In %s(), cannot get the reverse complements of a read, error!\n", __LINE__, __func__);
									return FAILED;
								}

								getReadBaseFromPosByInt(readseqTmp, readSeqIntRev, seqLen, 0, seqLen);
								printf("%ld\t%d\t%d\t-\t%d\t%s\n", rid, pReadMatchInfo->queryID, pReadMatchInfo->queryPos, seqLen, readseqTmp);

								free(readSeqIntRev);
							}
						}
*/
					}
				}
			}
		}
	}

	free(matchResultArrays[0]);
	free(matchResultArrays[1]);
	free(matchResultArrayBuf1);
	free(matchResultArrayBuf2);

	threadPara->successFlag = SUCCESSFUL;
}

/**
 * Get the maximal array size for the match result array for reads alignment to queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxArraySizeFromQueryIndex(int32_t *maxArraySize, queryIndex_t *queryIndex)
{
	int32_t i, j, maxValue;
	queryKmerBlock_t *pKmerBlock;
	queryKmer_t *kmer;

	maxValue = 0;
	for(i=0; i<queryIndex->blocksNumKmer; i++)
	{
		pKmerBlock = queryIndex->kmerBlockArr + i;
		for(j=0; j<pKmerBlock->itemNum; j++)
		{
			kmer = pKmerBlock->kmerArr + j;
			if(maxValue<kmer->arraysize)
				maxValue = kmer->arraysize;
		}
	}

	*maxArraySize = maxValue;

	return SUCCESSFUL;
}

/**
 * Estimate the insert size.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short estimateInsertSize(queryMatchInfo_t *queryMatchInfoSet, readSetArr_t *readSetArray, queryIndex_t *queryIndex)
{
	int32_t i;

	for(i=0; i<readSetArray->readSetNum; i++)
	{
		if(readSetArray->readSetArray[i].pairedMode==1 || readSetArray->readSetArray[i].pairedMode==2)
		{
			if(estimateInsertSizeSingleReadset(queryMatchInfoSet, readSetArray->readSetArray+i, queryIndex)==FAILED)
			{
				printf("line=%d, In %s(), cannot estimate the insert size, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Estimate the insert size.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short estimateInsertSizeSingleReadset(queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet, queryIndex_t *queryIndex)
{
	int32_t i, j, k, seqLen, maxArraySize, *matchItemNum, matchItemNums[2];
	readBlock_t *pReadBlockArray;
	read_t *pRead, *pReadArray;
	uint64_t rid;
	alignMatchItem_t *matchResultArray, *matchResultArrays[2], *matchResultArrayBuf1, *matchResultArrayBuf2;
	double insertSizeTmp, standDevTmp, fragSizeSum, fragSize, fragDif, fragDifSum;
	int32_t pairNum, pairNumTmp, pairNumMax, pairNumEnoughFlag, uniqueMapOpFlag;

	// the match result array
	if(getMaxArraySizeFromQueryIndex(&maxArraySize, queryIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the maximal array size from query index, error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(maxArraySize<=0)
	{
		printf("line=%d, In %s(), invalid array size %d, error!\n", __LINE__, __func__, maxArraySize);
		return FAILED;
	}

	matchResultArrays[0] = (alignMatchItem_t*) calloc (4*queryIndex->kmerSize*maxArraySize, sizeof(alignMatchItem_t));
	matchResultArrays[1] = (alignMatchItem_t*) calloc (4*queryIndex->kmerSize*maxArraySize, sizeof(alignMatchItem_t));
	matchResultArrayBuf1 = (alignMatchItem_t*) calloc (4*queryIndex->kmerSize*maxArraySize, sizeof(alignMatchItem_t));
	matchResultArrayBuf2 = (alignMatchItem_t*) calloc (4*queryIndex->kmerSize*maxArraySize, sizeof(alignMatchItem_t));
	if(matchResultArrays[0]==NULL || matchResultArrays[1]==NULL || matchResultArrayBuf1==NULL || matchResultArrayBuf2==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// estimate the insert size and standard deviation

	uniqueMapOpFlag = YES;
	pairNumMax = 2000;

	for(k=0; k<2; k++)
	{
		fragSizeSum = 0;
		pairNum = 0;

		pairNumTmp = 0;
		fragDifSum = 0;
		pairNumEnoughFlag = NO;

		rid = 0;
		for(i=0; i<readSet->blocksNumRead; i++)
		{
			pReadBlockArray = readSet->readBlockArr + i;
			pReadArray = pReadBlockArray->readArr;

			for(j=0; j<pReadBlockArray->itemNum; j++)
			{
				rid ++;
				pRead = pReadArray + j;
				seqLen = pRead->seqlen;

				if(rid%2==1)
				{
					matchResultArray = matchResultArrays[0];
					matchItemNum = matchItemNums;
				}else
				{
					matchResultArray = matchResultArrays[1];
					matchItemNum = matchItemNums + 1;
				}
				*matchItemNum = 0;

				if(pRead->validFlag==YES)
				{
					// map single read
					if(mapSingleReadToQueries(rid, pRead, readSet, matchResultArray, matchResultArrayBuf1, matchResultArrayBuf2, matchItemNum, queryIndex, queryMatchInfoSet->queryArray, 5)==FAILED)
					{
						printf("line=%d, In %s(), cannot map the read %lu, error!\n", __LINE__, __func__, rid);
						return FAILED;
					}
				}

				if(rid%2==0)
				{
					if(matchItemNums[0]==1 && matchItemNums[1]==1)
					{
						if(matchResultArrays[0][0].queryID==matchResultArrays[1][0].queryID && matchResultArrays[0][0].orientation!=matchResultArrays[1][0].orientation && matchResultArrays[0][0].mismatchNum==0 && matchResultArrays[1][0].mismatchNum==0)
						{
							fragSize = -1;
							if(matchResultArrays[0][0].orientation==ORIENTATION_PLUS && matchResultArrays[0][0].queryPos<=matchResultArrays[1][0].queryPos)
								fragSize = matchResultArrays[1][0].queryPos + pRead->seqlen - matchResultArrays[0][0].queryPos;
							else if(matchResultArrays[0][0].orientation==ORIENTATION_MINUS && matchResultArrays[0][0].queryPos>=matchResultArrays[1][0].queryPos)
								fragSize = matchResultArrays[0][0].queryPos + (pRead-1)->seqlen - matchResultArrays[1][0].queryPos;

							if(fragSize>=0)
							{
								if(k==0)
								{
									fragSizeSum += fragSize;
									pairNumTmp ++;
								}else
								{
									if(fragSize<2*insertSizeTmp)
									{
										fragDif = (fragSize - insertSizeTmp) * (fragSize - insertSizeTmp);
										fragDifSum += fragDif;
										pairNumTmp ++;
									}
								}

								if(pairNumTmp>=pairNumMax)
								{
									pairNumEnoughFlag = YES;
									break;
								}
							}
						}
					}
				}
			}

			if(pairNumEnoughFlag==YES)
				break;
		}

		if(k==0)
		{
			if(pairNumTmp>0)
				insertSizeTmp = fragSizeSum / pairNumTmp;
			else
			{
				printf("line=%d, In %s(), no valid paired reads for computing insert size of library, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else
		{
			if(pairNumTmp>0)
				standDevTmp = sqrt((fragDifSum) / pairNumTmp);
			else
			{
				printf("line=%d, In %s(), no valid paired reads for computing insert size of library, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	free(matchResultArrays[0]);
	free(matchResultArrays[1]);
	free(matchResultArrayBuf1);
	free(matchResultArrayBuf2);

	readSet->insertSize = insertSizeTmp;
	readSet->standDev = standDevTmp;

	printf("insertSize=%.4f, SDev=%.4f\n", readSet->insertSize, readSet->standDev);

	return SUCCESSFUL;
}

/**
 * Perfectly map single read to queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mapSingleReadToQueries(uint64_t rid, read_t *pRead, readSet_t *readSet, alignMatchItem_t *matchResultArray, alignMatchItem_t *matchResultArrayBuf1, alignMatchItem_t *matchResultArrayBuf2, int32_t *matchItemNum, queryIndex_t *queryIndex, query_t *queryArray, int32_t mismatchNumThreshold)
{
	int32_t tmpItemNum;

	*matchItemNum = 0;

/*
	// perfectly map single read
	if(mapSingleReadToQueriesPerfect(rid, pRead, readSet, matchResultArray, matchResultArrayBuf1, matchResultArrayBuf2, &tmpItemNum, queryIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot map the read %lu, error!\n", __LINE__, __func__, rid);
		return FAILED;
	}
	(*matchItemNum) += tmpItemNum;
*/

	//if((*matchItemNum)<=0)
	{
		// map single read with some mismatches
		if(mapSingleReadToQueriesWithMismatch(rid, pRead, readSet, matchResultArray+(*matchItemNum), matchResultArrayBuf1, &tmpItemNum, queryIndex, queryArray, mismatchNumThreshold)==FAILED)
		{
			printf("line=%d, In %s(), cannot map the read %lu, error!\n", __LINE__, __func__, rid);
			return FAILED;
		}
		(*matchItemNum) += tmpItemNum;
	}

	return SUCCESSFUL;
}

/**
 * Perfectly map single read to queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mapSingleReadToQueriesPerfect(int64_t rid, read_t *pRead, readSet_t *readSet, alignMatchItem_t *matchResultArray, alignMatchItem_t *matchResultArrayBuf1, alignMatchItem_t *matchResultArrayBuf2, int32_t *matchItemNum, queryIndex_t *queryIndex)
{
	int32_t i, seqLen, entriesNumRead, tmpItemNum;
	uint64_t *readSeqInt, *readSeqIntRev;

	readSeqInt = readSet->readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
	seqLen = pRead->seqlen;

	entriesNumRead = ((seqLen - 1) >> 5) + 1;

	*matchItemNum = 0;
	if(seqLen>=queryIndex->kmerSize)
	{
		if(getMatchedQueryPosPerfect(matchResultArray, matchResultArrayBuf1, matchResultArrayBuf2, &tmpItemNum, readSeqInt, seqLen, queryIndex)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the match query position, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(tmpItemNum>0)
		{
			for(i=0; i<tmpItemNum; i++)
				matchResultArray[i].orientation = ORIENTATION_PLUS;
			*matchItemNum = tmpItemNum;
		}

		// the reverse complements
		readSeqIntRev = (uint64_t*) calloc (entriesNumRead, sizeof(uint64_t));
		if(readSeqIntRev==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// generate the reverse complements
		if(getReverseReadseqInt(readSeqIntRev, readSeqInt, seqLen)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the reverse complements of a read, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// align the read
		if(getMatchedQueryPosPerfect(matchResultArray+(*matchItemNum), matchResultArrayBuf1, matchResultArrayBuf2, &tmpItemNum, readSeqIntRev, seqLen, queryIndex)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the match query position, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(tmpItemNum>0)
		{
			for(i=0; i<tmpItemNum; i++)
				matchResultArray[i+(*matchItemNum)].orientation = ORIENTATION_MINUS;
			*matchItemNum += tmpItemNum;
		}

		free(readSeqIntRev);
	}

	return SUCCESSFUL;
}

/**
 * Get the matched query position of a read sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMatchedQueryPosPerfect(alignMatchItem_t *matchResultArray, alignMatchItem_t *matchResultArrayBuf1, alignMatchItem_t *matchResultArrayBuf2, int32_t *matchItemNum, uint64_t *readSeqInt, int32_t seqLen, queryIndex_t *queryIndex)
{
	int32_t i, j, k, entriesNumRead, baseNumLastEntryRead, matchItemNum1, matchItemNum2, newMatchItemNum1, validItemNum, mapTimes, startBasePos;
	uint64_t hashcode, kmerSeqInt[queryIndex->entriesPerKmer];
	queryKmer_t *kmer;
	int32_t matchDistance, matchFlag, startRow, endRow, shiftBaseNum;

	entriesNumRead = ((seqLen - 1) >> 5) + 1;
	baseNumLastEntryRead = ((seqLen - 1) % 32) + 1;

	mapTimes = (seqLen - 1) / queryIndex->kmerSize + 1;
	matchItemNum1 = 0;
	newMatchItemNum1 = 0;
	if(seqLen>=queryIndex->kmerSize)
	{
		for(shiftBaseNum=0; shiftBaseNum<queryIndex->kmerSize; shiftBaseNum++)
		{
			startBasePos = shiftBaseNum;

			// generate the kmer integer sequence
			if(generateKmerSeqIntFromReadset(kmerSeqInt, readSeqInt, startBasePos, queryIndex->kmerSize, entriesNumRead, baseNumLastEntryRead)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// get the kmer
			hashcode = kmerhashInt(kmerSeqInt, queryIndex->entriesPerKmer, queryIndex->lastEntryBaseNum, queryIndex->hashTableSize);
			kmer = getQueryKmerByHash(hashcode, kmerSeqInt, queryIndex);

			// record the match information
			if(kmer==NULL)
			{
				matchItemNum1 = 0;
			}else
			{
				matchItemNum1 = kmer->arraysize;
				for(j=0; j<matchItemNum1; j++)
				{
					matchResultArrayBuf1[j].queryID = kmer->ppos[j].queryID;
					matchResultArrayBuf1[j].queryPos = kmer->ppos[j].queryPos;
					matchResultArrayBuf1[j].mismatchNum = 0;
					matchResultArrayBuf1[j].validFlag = YES;
				}

				// process the other bases
				validItemNum = matchItemNum1;
				for(i=1; i<mapTimes && validItemNum>0; i++)
				{
					if(i<mapTimes-1)
					{
						startBasePos += queryIndex->kmerSize;
						matchDistance = startBasePos;
					}else
					{
						startBasePos = seqLen - queryIndex->kmerSize;
						matchDistance = startBasePos;
					}

					// generate the kmer integer sequence
					if(generateKmerSeqIntFromReadset(kmerSeqInt, readSeqInt, startBasePos, queryIndex->kmerSize, entriesNumRead, baseNumLastEntryRead)==FAILED)
					{
						printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
						return FAILED;
					}

					// get the kmer
					hashcode = kmerhashInt(kmerSeqInt, queryIndex->entriesPerKmer, queryIndex->lastEntryBaseNum, queryIndex->hashTableSize);
					kmer = getQueryKmerByHash(hashcode, kmerSeqInt, queryIndex);

					// record the match information
					if(kmer==NULL)
					{
						matchItemNum2 = 0;
					}else
					{
						matchItemNum2 = kmer->arraysize;
						for(j=0; j<matchItemNum2; j++)
						{
							matchResultArrayBuf2[j].queryID = kmer->ppos[j].queryID;
							matchResultArrayBuf2[j].queryPos = kmer->ppos[j].queryPos;
							matchResultArrayBuf2[j].mismatchNum = 0;
						}
					}

					// get the intersection
					if(matchItemNum2==0)
					{
						matchItemNum1 = 0;
						validItemNum = 0;
						break;
					}else if(matchItemNum2>0)
					{
						for(j=0; j<matchItemNum1; j++)
						{
							if(matchResultArrayBuf1[j].validFlag==YES)
							{
								matchFlag = NO;

								// get the startRow and endRow of Buf2
								if(matchItemNum2>20)
								{
									if(getRowRangeMatchArray(&startRow, &endRow, matchResultArrayBuf2, matchItemNum2, matchResultArrayBuf1[j].queryID)==FAILED)
									{
										printf("line=%d, In %s(), cannot get the array range, error!\n", __LINE__, __func__);
										return FAILED;
									}
								}
								else
								{
									startRow = 0;
									endRow = matchItemNum2 - 1;
								}

								if(startRow>=0 && endRow>=0)
								{
									for(k=startRow; k<=endRow; k++)
									{
										if(matchResultArrayBuf1[j].queryID==matchResultArrayBuf2[k].queryID && matchResultArrayBuf1[j].queryPos+matchDistance==matchResultArrayBuf2[k].queryPos)
										{
											matchFlag = YES;
											break;
										}
									}
								}

								if(matchFlag==NO)
								{
									matchResultArrayBuf1[j].validFlag = NO;
									validItemNum --;
								}
							}
						}
					}
				} // for(i=1; i<mapTimes; i++)

				// copy the valid items
				if(validItemNum>0)
				{
					newMatchItemNum1 = 0;
					for(j=0; j<matchItemNum1; j++)
					{
						if(matchResultArrayBuf1[j].validFlag==YES)
						{
							if(memcpy(matchResultArray+newMatchItemNum1, matchResultArrayBuf1+j, sizeof(alignMatchItem_t))==NULL)
							{
								printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
								return FAILED;
							}
							newMatchItemNum1 ++;
						}
					}

					if(newMatchItemNum1!=validItemNum)
					{
						printf("newMatchItemNum1=%d, validItemNum=%d, error!\n", newMatchItemNum1, validItemNum);
					}
				}
			}
		}
	}

	*matchItemNum = newMatchItemNum1;

	return SUCCESSFUL;
}

/**
 * Map single read to queries with some mismatches.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mapSingleReadToQueriesWithMismatch(int64_t rid, read_t *pRead, readSet_t *readSet, alignMatchItem_t *matchResultArray, alignMatchItem_t *matchResultArrayBuf, int32_t *matchItemNum, queryIndex_t *queryIndex, query_t *queryArray, int32_t mismatchNumThreshold)
{
	int32_t seqLen, entriesNumRead, tmpItemNum;
	uint64_t *readSeqInt, *readSeqIntRev;

	readSeqInt = readSet->readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
	seqLen = pRead->seqlen;

	entriesNumRead = ((seqLen - 1) >> 5) + 1;

	*matchItemNum = 0;
	if(seqLen>=queryIndex->kmerSize)
	{
		if(getMatchedQueryPosWithMismatch(matchResultArray, &tmpItemNum, readSeqInt, seqLen, ORIENTATION_PLUS, queryIndex, queryArray, mismatchNumThreshold)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the match query position, error!\n", __LINE__, __func__);
			return FAILED;
		}
		*matchItemNum = tmpItemNum;

		// the reverse complements
		readSeqIntRev = (uint64_t*) calloc (entriesNumRead, sizeof(uint64_t));
		if(readSeqIntRev==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// generate the reverse complements
		if(getReverseReadseqInt(readSeqIntRev, readSeqInt, seqLen)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the reverse complements of a read, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// align the read
		if(getMatchedQueryPosWithMismatch(matchResultArray+(*matchItemNum), &tmpItemNum, readSeqIntRev, seqLen, ORIENTATION_MINUS, queryIndex, queryArray, mismatchNumThreshold)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the match query position, error!\n", __LINE__, __func__);
			return FAILED;
		}
		*matchItemNum += tmpItemNum;

		free(readSeqIntRev);

		// sort the results according to mismatch
		if((*matchItemNum)>=2)
		{
			if(sortMatchResults(matchResultArray, matchResultArrayBuf, *matchItemNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot sort the match results, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the matched query position of a read sequence with some mismatches.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMatchedQueryPosWithMismatch(alignMatchItem_t *matchResultArray, int32_t *matchItemNum, uint64_t *readSeqInt, int32_t seqLen, int32_t orientation, const queryIndex_t *queryIndex, query_t *queryArray, int32_t mismatchNumThreshold)
{
	int32_t i, j, entriesNumRead, baseNumLastEntryRead, startBasePos, processedFlag;
	uint64_t hashcode, kmerSeqInt[queryIndex->entriesPerKmer];
	queryKmer_t *kmer;
	char *querySeq, readSeq[MAX_READ_LEN_IN_BUF+1];
	int32_t queryID, queryLen, queryPos, maxStartBasePos, remainBaseNum, misNum, startCmpRpos, startCmpQpos;
	int32_t startCmpRposLeft, startCmpQposLeft, startCmpRposRight, startCmpQposRight, remainBaseNumLeft, remainBaseNumRight;
	int32_t startAlignReadPos, endAlignReadPos, startAlignQueryPos, queryEndIgnoreLen, shiftBaseNum;
	int32_t headSkipBaseNumRead, tailSkipBaseNumRead;

	entriesNumRead = ((seqLen - 1) >> 5) + 1;
	baseNumLastEntryRead = ((seqLen - 1) % 32) + 1;
	queryEndIgnoreLen = 0.15 * seqLen;

	*matchItemNum = 0;
	if(seqLen>=queryIndex->kmerSize)
	{
		// get the read base sequence
		getReadBaseFromPosByInt(readSeq, readSeqInt, seqLen, 0, seqLen);

		for(shiftBaseNum=0; shiftBaseNum<queryIndex->kmerSize; shiftBaseNum++)
		{
			maxStartBasePos = seqLen - queryIndex->kmerSize;
			for(startBasePos = shiftBaseNum; startBasePos<=maxStartBasePos; startBasePos += queryIndex->kmerSize)
			{
				// generate the kmer integer sequence
				if(generateKmerSeqIntFromReadset(kmerSeqInt, readSeqInt, startBasePos, queryIndex->kmerSize, entriesNumRead, baseNumLastEntryRead)==FAILED)
				{
					printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// get the kmer
				hashcode = kmerhashInt(kmerSeqInt, queryIndex->entriesPerKmer, queryIndex->lastEntryBaseNum, queryIndex->hashTableSize);
				kmer = getQueryKmerByHash(hashcode, kmerSeqInt, queryIndex);

				// record the match information
				if(kmer)
				{
					for(i=0; i<kmer->arraysize; i++)
					{
						queryID = kmer->ppos[i].queryID;
						queryPos = kmer->ppos[i].queryPos;

						querySeq = queryArray[queryID-1].querySeq;
						queryLen = queryArray[queryID-1].queryLen;

						if(queryPos-startBasePos>=1)
						{
							startAlignReadPos = 1;
							startAlignQueryPos = queryPos - startBasePos;
						}else
						{
							startAlignReadPos = startBasePos + 2 - queryPos;
							startAlignQueryPos = 1;
						}

						remainBaseNum = seqLen - startBasePos;
						if(queryPos+remainBaseNum-1<=queryLen)
							endAlignReadPos = seqLen;
						else
							endAlignReadPos = queryLen - (queryPos - startBasePos - 1);

						if((endAlignReadPos>=seqLen-queryEndIgnoreLen+1) && (startAlignReadPos<=queryEndIgnoreLen))
						{
							// check whether it is processed already
							if(isProcessedMatchItem(&processedFlag, queryID, startAlignQueryPos, matchResultArray, *matchItemNum)==FAILED)
							{
								printf("line=%d, In %s(), cannot get the processed flag of alignment item, error!\n", __LINE__, __func__);
								return FAILED;
							}

							if(processedFlag==NO)
							{
								// mismatch on the right side
								misNum = 0;

								startCmpRposRight = startBasePos + queryIndex->kmerSize;
								startCmpQposRight = queryPos + queryIndex->kmerSize - 1;
								if(queryLen-startCmpQposRight<seqLen-startCmpRposRight)
									remainBaseNumRight = queryLen - startCmpQposRight;
								else
									remainBaseNumRight = seqLen - startCmpRposRight;

								// compute the mismatches on the right side
								for(j=0; j<remainBaseNumRight; j++)
								{
									if(readSeq[startCmpRposRight+j]!=querySeq[startCmpQposRight+j])
									{
										misNum ++;
										if(misNum>mismatchNumThreshold)
											break;
									}
								}

								// compute the skip mismatched bases at the head of the read
								tailSkipBaseNumRead = 0;
								for(j=remainBaseNumRight-1; j>=0; j--)
								{
									if(readSeq[startCmpRposRight+j]!=querySeq[startCmpQposRight+j])
										tailSkipBaseNumRead ++;
									else
										break;
								}

								// mismatch on the left side
								startCmpRposLeft = startAlignReadPos - 1;
								if(startAlignReadPos==1)
								{
									startCmpQposLeft = queryPos - startBasePos - 1;
									remainBaseNumLeft = startBasePos;
								}else
								{
									startCmpQposLeft = 0;
									remainBaseNumLeft = startBasePos + 1 - startAlignReadPos;
								}

								// compute the mismatches on the left side
								for(j=0; j<remainBaseNumLeft; j++)
								{
									if(readSeq[startCmpRposLeft+j]!=querySeq[startCmpQposLeft+j])
									{
										misNum ++;
										if(misNum>mismatchNumThreshold)
											break;
									}
								}

								// compute the skip mismatched bases at the head of the read
								headSkipBaseNumRead = 0;
								for(j=0; j<remainBaseNumLeft; j++)
								{
									if(readSeq[startCmpRposLeft+j]!=querySeq[startCmpQposLeft+j])
										headSkipBaseNumRead ++;
									else
										break;
								}

								// save the match result to matchResultArray
								if(misNum<=mismatchNumThreshold)
								{ // add the match information
									matchResultArray[*matchItemNum].queryID = queryID;
									matchResultArray[*matchItemNum].queryPos = startAlignQueryPos;
									matchResultArray[*matchItemNum].mismatchNum = misNum;
									matchResultArray[*matchItemNum].startReadPos = startAlignReadPos;
									matchResultArray[*matchItemNum].alignSize = endAlignReadPos - startAlignReadPos + 1;
									matchResultArray[*matchItemNum].orientation = orientation;
									matchResultArray[*matchItemNum].headSkipNum = headSkipBaseNumRead;
									matchResultArray[*matchItemNum].tailSkipNum = tailSkipBaseNumRead;
									matchResultArray[*matchItemNum].validFlag = YES;
									(*matchItemNum) ++;
								}
							}
						}
					}
				}
			}
		}
	}

	// adjust the startReadPos according to headSkipNum and tailSkipNum
	for(i=0; i<*matchItemNum; i++)
	{
		if(matchResultArray[i].headSkipNum>0)
		{
			matchResultArray[i].startReadPos += matchResultArray[i].headSkipNum;
			matchResultArray[i].queryPos += matchResultArray[i].headSkipNum;
			matchResultArray[i].mismatchNum -= matchResultArray[i].headSkipNum;
			matchResultArray[i].alignSize -= matchResultArray[i].headSkipNum;
			matchResultArray[i].headSkipNum = -1;
		}

		if(matchResultArray[i].tailSkipNum>0)
		{
			matchResultArray[i].mismatchNum -= matchResultArray[i].tailSkipNum;
			matchResultArray[i].alignSize -= matchResultArray[i].tailSkipNum;
			matchResultArray[i].tailSkipNum = -1;
		}
	}

	return SUCCESSFUL;
}

/**
 * check whether the alignment item is processed already.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short isProcessedMatchItem(int32_t *processedFlag, int32_t queryID, int32_t startAlignQueryPos, alignMatchItem_t *matchResultArray, int32_t matchItemNum)
{
	int32_t i;

	*processedFlag = NO;
	for(i=0; i<matchItemNum; i++)
	{
		if(matchResultArray[i].queryID==queryID && matchResultArray[i].queryPos==startAlignQueryPos)
		{
			*processedFlag = YES;
			break;
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the mated pair count from paired-end reads alignment result.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMatedNumAlignResult(int32_t *matedNum, alignMatchItem_t *matchResultArray1, alignMatchItem_t *matchResultArray2, int32_t matchItemNum1, int32_t matchItemNum2, int32_t seqLen1, int32_t seqLen2, double insertSize, double standDev)
{
	int32_t i, j, queryID1, queryPos1, queryID2, queryPos2, orient1, orient2, matedRow1, matedRow2, pairRow;
	double fragSize, fragDif, fragDifExist;

	for(i=0; i<matchItemNum1; i++) matchResultArray1[i].pairRow = -1;
	for(i=0; i<matchItemNum2; i++) matchResultArray2[i].pairRow = -1;

	*matedNum = 0;
	matedRow1 = matedRow2 = -1;
	for(i=0; i<matchItemNum1; i++)
	{
		queryID1 = matchResultArray1[i].queryID;
		queryPos1 = matchResultArray1[i].queryPos;
		orient1 = matchResultArray1[i].orientation;

		for(j=0; j<matchItemNum2; j++)
		{
			queryID2 = matchResultArray2[j].queryID;
			queryPos2 = matchResultArray2[j].queryPos;
			orient2 = matchResultArray2[j].orientation;

			if(queryID1==queryID2 && orient1!=orient2)
			{
				fragSize = INT_MAX;
				if(orient1==ORIENTATION_PLUS)
					fragSize = queryPos2 + seqLen2 - queryPos1;
				else
					fragSize = queryPos1 + seqLen1 - queryPos2;
				fragDif = fragSize - insertSize;
				if(matchResultArray1[i].pairRow==-1)
				{
					if(fragDif>=-6*standDev && fragDif<=6*standDev)
					{
						matchResultArray1[i].pairRow = j;
						matchResultArray1[i].fragSize = fragSize;
						matedRow1 = i;

						matchResultArray2[j].pairRow = i;
						matchResultArray2[j].fragSize = fragSize;
						matedRow2 = j;

						(*matedNum) ++;
					}
				}else
				{
					pairRow = matchResultArray1[i].pairRow;
					if(matchResultArray2[j].mismatchNum<=matchResultArray2[pairRow].mismatchNum)
					{
						fragDifExist = matchResultArray1[i].fragSize - insertSize;
						if(fragDifExist<0)
							fragDifExist = -fragDifExist;
						if(fragDif<0)
							fragDif = -fragDif;
						if(fragDif<fragDifExist)
						{ // change the pair
							matchResultArray2[pairRow].pairRow = -1;
							matchResultArray2[pairRow].fragSize = 0;

							matchResultArray1[i].pairRow = j;
							matchResultArray1[i].fragSize = fragSize;

							matchResultArray2[j].pairRow = i;
							matchResultArray2[j].fragSize = fragSize;
							matedRow2 = j;
						}
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Compute unique map flag for paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeUniqueMapFlags(int32_t *uniqueMapFlagArray, int32_t matedNum, alignMatchItem_t *matchResultArray1, alignMatchItem_t *matchResultArray2, int32_t matchItemNum1, int32_t matchItemNum2)
{
	if(matedNum==1)
	{
		uniqueMapFlagArray[0] = YES;
		uniqueMapFlagArray[1] = YES;
	}else
	{
		if(matchItemNum1<=1)
			uniqueMapFlagArray[0] = YES;
		else
			uniqueMapFlagArray[0] = NO;
		if(matchItemNum2<=1)
			uniqueMapFlagArray[1] = YES;
		else
			uniqueMapFlagArray[1] = NO;
	}

	return SUCCESSFUL;
}

/**
 * Select the mated query position of paired reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short selectMatedMatchPosPE(alignMatchItem_t *matchResultArray1, alignMatchItem_t *matchResultArray2, int32_t *matchItemNum1, int32_t *matchItemNum2, int32_t matedNum, int32_t seqLen1, int32_t seqLen2, double insertSize, double standDev, int32_t uniqueMapOpFlag, queryMatchInfo_t *queryMatchInfoSet, threadPara_t *threadPara)
{
	struct pairMismatch
	{
		int32_t arrRow1, arrRow2, mismatchNum1, mismatchNum2, totalMismatchNum, fragSize, fragDif;
	};

	int32_t i, matedRow1, matedRow2, pairRow;
	double fragDif;
	int32_t *zeroCovNumArray, zeroNum, itemRow, maxValue, maxRow, maxNum;

	struct pairMismatch *pairMismatchArray;
	int32_t minValue, minRow, minNum, itemNum, arrRow1, arrRow2, minID;

	if(matedNum==1)
	{ // unique mated
		// get matedRow
		for(i=0; i<(*matchItemNum1); i++)
		{
			if(matchResultArray1[i].pairRow>=0)
			{
				matedRow1 = i;
				matedRow2 = matchResultArray1[i].pairRow;
				break;
			}
		}

		if((*matchItemNum1)>1)
		{
			if(matedRow1>=1)
			{
				if(memcpy(matchResultArray1, matchResultArray1+matedRow1, sizeof(alignMatchItem_t))==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				matchResultArray2[matchResultArray1[0].pairRow].pairRow = 0;
			}
			*matchItemNum1 = 1;
		}

		if((*matchItemNum2)>1)
		{
			if(matedRow2>=1)
			{
				if(memcpy(matchResultArray2, matchResultArray2+matedRow2, sizeof(alignMatchItem_t))==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				matchResultArray1[matchResultArray2[0].pairRow].pairRow = 0;
			}
			*matchItemNum2 = 1;
		}
	}else if(matedNum>=2)
	{ // multiple mated, then select the minimal mismatch pair

		pairMismatchArray = (struct pairMismatch *) calloc (matedNum, sizeof(struct pairMismatch));
		if(pairMismatchArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		itemNum = 0;
		for(i=0; i<(*matchItemNum1); i++)
		{
			pairRow = matchResultArray1[i].pairRow;
			if(pairRow>=0)
			{
				pairMismatchArray[itemNum].arrRow1 = i;
				pairMismatchArray[itemNum].arrRow2 = pairRow;
				pairMismatchArray[itemNum].mismatchNum1 = matchResultArray1[i].mismatchNum;
				pairMismatchArray[itemNum].mismatchNum2 = matchResultArray2[pairRow].mismatchNum;
				pairMismatchArray[itemNum].totalMismatchNum = pairMismatchArray[itemNum].mismatchNum1 + pairMismatchArray[itemNum].mismatchNum2;
				pairMismatchArray[itemNum].fragSize = matchResultArray1[i].fragSize;
				fragDif = matchResultArray1[i].fragSize - insertSize;
				if(fragDif<0)
					fragDif = -fragDif;
				pairMismatchArray[itemNum].fragDif = fragDif;
				itemNum ++;
			}
		}

		minValue = INT_MAX;
		minRow = -1;
		minNum = 0;
		for(i=0; i<matedNum; i++)
		{
			if(pairMismatchArray[i].totalMismatchNum<minValue)
			{
				minValue = pairMismatchArray[i].totalMismatchNum;
				minRow = i;
				minNum = 1;
			}else if(pairMismatchArray[i].totalMismatchNum==minValue)
			{
				minNum ++;
			}
		}

		if(minNum>1)
		{
			if(uniqueMapOpFlag==NO)
			{
				// select the item with maximum zero coverage
				zeroCovNumArray = (int32_t*) calloc (matedNum, sizeof(int32_t));
				if(zeroCovNumArray==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}

				for(i=0; i<matedNum; i++)
				{
					if(pairMismatchArray[i].totalMismatchNum==minValue)
					{
						itemRow = pairMismatchArray[i].arrRow1;
						if(getZeroCovBaseNumSingleRead(&zeroNum, matchResultArray1+itemRow, queryMatchInfoSet, threadPara)==FAILED)
						{
							printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
						zeroCovNumArray[i] += zeroNum;

						itemRow = pairMismatchArray[i].arrRow2;
						if(getZeroCovBaseNumSingleRead(&zeroNum, matchResultArray2+itemRow, queryMatchInfoSet, threadPara)==FAILED)
						{
							printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
						zeroCovNumArray[i] += zeroNum;
					}
				}

				maxValue = INT_MIN;
				maxRow = -1;
				maxNum = 0;
				for(i=0; i<matedNum; i++)
				{
					if(pairMismatchArray[i].totalMismatchNum==minValue)
					{
						if(maxValue<zeroCovNumArray[i])
						{
							maxValue = zeroCovNumArray[i];
							maxRow = i;
							maxNum = 1;
						}else if(maxValue==zeroCovNumArray[i])
						{
							maxNum ++;
						}
					}
				}

				if(maxNum==1)
				{
					minRow = maxRow;
				}else if(maxNum>1)
				{ // select the item randomly
					minRow = -1;
					minID = rand() % maxNum;
					itemNum = 0;
					for(i=0; i<matedNum; i++)
					{
						if(pairMismatchArray[i].totalMismatchNum==minValue && zeroCovNumArray[i]==maxValue)
						{
							if(itemNum==minID)
							{
								minRow = i;
								break;
							}
							itemNum ++;
						}
					}
				}else
				{
					printf("line=%d, In %s(), maxNum=%d, error!\n", __LINE__, __func__, maxNum);
					return FAILED;
				}

				free(zeroCovNumArray);
			}else
			{
				minRow = -1;
			}
		}

		if(minRow>=0)
		{
			arrRow1 = pairMismatchArray[minRow].arrRow1;
			arrRow2 = pairMismatchArray[minRow].arrRow2;
			if(arrRow1>0)
			{
				if(memcpy(matchResultArray1, matchResultArray1+arrRow1, sizeof(alignMatchItem_t))==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				matchResultArray2[arrRow2].pairRow = 0;
				arrRow1 = pairMismatchArray[minRow].arrRow1 = 0;
			}

			if(arrRow2>0)
			{
				if(memcpy(matchResultArray2, matchResultArray2+arrRow2, sizeof(alignMatchItem_t))==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				matchResultArray1[arrRow1].pairRow = 0;
				arrRow2 = pairMismatchArray[minRow].arrRow2 = 0;
			}

			*matchItemNum1 = 1;
			*matchItemNum2 = 1;
		}

		free(pairMismatchArray);

	}else if(matedNum==0)
	{
		// select best alignment for single read
		if(selectBestAlignInfoSingleRead(matchResultArray1, matchItemNum1, uniqueMapOpFlag, queryMatchInfoSet, threadPara)==FAILED)
		{
			printf("line=%d, In %s(), cannot select the best alignment information for single read, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// select best alignment for single read
		if(selectBestAlignInfoSingleRead(matchResultArray2, matchItemNum2, uniqueMapOpFlag, queryMatchInfoSet, threadPara)==FAILED)
		{
			printf("line=%d, In %s(), cannot select the best alignment information for single read, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the zero covered base count for single read.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getZeroCovBaseNumSingleRead(int32_t *zeroCovBaseNum, alignMatchItem_t *alignMatchItem, queryMatchInfo_t *queryMatchInfoSet, threadPara_t *threadPara)
{
	int32_t i, queryID, queryPos, queryLen, alignSize, arrRow, arrCol, baseCovFlag;
	uint64_t *covFlagArray;
	query_t *queryItem;
	ctrlLock_t *lockItem;

	queryID = alignMatchItem->queryID;
	queryPos = alignMatchItem->queryPos;
	queryItem = queryMatchInfoSet->queryArray + (queryID - 1);
	queryLen = queryItem->queryLen;
	covFlagArray = queryItem->covFlagArray;
	alignSize = alignMatchItem->alignSize;

	arrRow = (queryPos - 1) / 64;
	arrCol = (queryPos - 1) % 64;

	lockItem = threadPara->lockArray + queryID - 1;

	// increase the reader count
	pthread_mutex_lock(&lockItem->outerLock);
	{
		pthread_mutex_lock(&lockItem->readerLock);
		{
			pthread_mutex_lock(&lockItem->accessReaderCnt);
			{
				lockItem->readerCnt ++;
				if(lockItem->readerCnt == 1)
				{
					pthread_mutex_lock(&lockItem->writeLock);
				}
			}
			pthread_mutex_unlock(&lockItem->accessReaderCnt);
		}
		pthread_mutex_unlock(&lockItem->readerLock);
	}
	pthread_mutex_unlock(&lockItem->outerLock);

	*zeroCovBaseNum = 0;
	for(i=0; i<alignSize; i++)
	{
		baseCovFlag = (covFlagArray[arrRow] >> (63 - arrCol)) & 1;
		if(baseCovFlag==0)
			(*zeroCovBaseNum) ++;
		arrCol ++;
		if(arrCol==64)
		{
			arrCol = 0;
			arrRow ++;
		}

		if(queryPos+i>queryLen)
		{
			printf("line=%d, In %s(), basePos=%d, queryLen=%d, error!\n", __LINE__, __func__, queryPos+i, queryLen);
			return FAILED;
		}
	}

	// decrease the reader count
	pthread_mutex_lock(&lockItem->accessReaderCnt);
	{
		lockItem->readerCnt --;
		if(lockItem->readerCnt==0)
		{
			pthread_mutex_unlock(&lockItem->writeLock);
		}
	}
	pthread_mutex_unlock(&lockItem->accessReaderCnt);

	return SUCCESSFUL;
}

/**
 * Select best alignment for single read.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short selectBestAlignInfoSingleRead(alignMatchItem_t *matchResultArray, int32_t *matchItemNum, int32_t uniqueMapOpFlag, queryMatchInfo_t *queryMatchInfoSet, threadPara_t *threadPara)
{
	int32_t i, minValue, minRow, minNum, itemNum, minID;
	int32_t *zeroCovNumArray;
	int32_t maxValue, maxRow, maxNum;

	if((*matchItemNum)>1)
	{
		minValue = INT_MAX;
		minRow = -1;
		minNum = 0;
		for(i=0; i<(*matchItemNum); i++)
		{
			if(matchResultArray[i].mismatchNum<minValue)
			{
				minValue = matchResultArray[i].mismatchNum;
				minRow = i;
				minNum = 1;
			}else if(matchResultArray[i].mismatchNum==minValue)
			{
				minNum ++;
			}
		}

		if(minNum>1)
		{
			if(uniqueMapOpFlag==NO)
			{
				// select the item with maximum zero coverage
				zeroCovNumArray = (int32_t*) calloc (*matchItemNum, sizeof(int32_t));
				if(zeroCovNumArray==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}

				for(i=0; i<(*matchItemNum); i++)
				{
					if(matchResultArray[i].mismatchNum==minValue)
					{
						if(getZeroCovBaseNumSingleRead(zeroCovNumArray+i, matchResultArray+i, queryMatchInfoSet, threadPara)==FAILED)
						{
							printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
				}

				maxValue = INT_MIN;
				maxRow = -1;
				maxNum = 0;
				for(i=0; i<(*matchItemNum); i++)
				{
					if(matchResultArray[i].mismatchNum==minValue)
					{
						if(maxValue<zeroCovNumArray[i])
						{
							maxValue = zeroCovNumArray[i];
							maxRow = i;
							maxNum = 1;
						}else if(maxValue==zeroCovNumArray[i])
						{
							maxNum ++;
						}
					}
				}

				if(maxNum==1)
				{
					minRow = maxRow;
				}else if(maxNum>1)
				{ // select the item randomly
					minRow = -1;
					minID = rand() % maxNum;
					itemNum = 0;
					for(i=0; i<(*matchItemNum); i++)
					{
						if(matchResultArray[i].mismatchNum==minValue && zeroCovNumArray[i]==maxValue)
						{
							if(itemNum==minID)
							{
								minRow = i;
								break;
							}
							itemNum ++;
						}
					}
				}else
				{
					printf("line=%d, In %s(), maxNum=%d, error!\n", __LINE__, __func__, maxNum);
					return FAILED;
				}

				free(zeroCovNumArray);
			}else
			{
				minRow = -1;
			}
		}else if(minNum<=0)
		{
			printf("line=%d, In %s(), minNum=%d, error!\n", __LINE__, __func__, minNum);
			return FAILED;
		}

		if(minRow>0)
		{
			if(memcpy(matchResultArray, matchResultArray+minRow, sizeof(alignMatchItem_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			*matchItemNum = 1;
		}else if(minRow==0)
		{
			*matchItemNum = 1;
		}
	}

	return SUCCESSFUL;
}

/**
 * Get more align items from unmated paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMoreAlignItemsFromUnmatedPE(int32_t *matedNum, int64_t readID1, int64_t readID2, alignMatchItem_t *matchResultArray1, alignMatchItem_t *matchResultArray2, int32_t *matchItemNum1, int32_t *matchItemNum2, int32_t mismatchNumThreshold, queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet)
{
	read_t *pRead1, *pRead2;

	// get the two ends
	if(getReadByReadID(&pRead1, readID1, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the read %ld, error!\n", __LINE__, __func__, readID1);
		return FAILED;
	}
	if(getReadByReadID(&pRead2, readID2, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the read %ld, error!\n", __LINE__, __func__, readID2);
		return FAILED;
	}

	if(pRead1->validFlag==YES && pRead2->validFlag==YES)
	{
		// get aligned reads based on the second end
		if(getAlignItemsFromUnmatedEnd(matedNum, pRead1, matchResultArray1, matchItemNum1, pRead2, matchResultArray2, *matchItemNum2, mismatchNumThreshold, queryMatchInfoSet, readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot get align items from unmated end, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// get aligned reads based on second ends
		if(getAlignItemsFromUnmatedEnd(matedNum, pRead2, matchResultArray2, matchItemNum2, pRead1, matchResultArray1, *matchItemNum1, mismatchNumThreshold, queryMatchInfoSet, readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot get align items from unmated end, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Get more align items from unmated end.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getAlignItemsFromUnmatedEnd(int32_t *matedNum, read_t *pRead, alignMatchItem_t *matchResultArray, int32_t *matchItemNum, read_t *pReadBased, alignMatchItem_t *matchResultArrayBased, int32_t matchItemNumBased, int32_t mismatchNumThreshold, queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet)
{
	int32_t i;

	for(i=0; i<matchItemNumBased; i++)
	{
		if(matchResultArrayBased[i].pairRow==-1)
		{
			if(getAlignItemsFromSingleUnmatedPos(matedNum, pRead, matchResultArray, matchItemNum, pReadBased, matchResultArrayBased+i, i, mismatchNumThreshold, queryMatchInfoSet, readSet)==FAILED)
			{
				printf("line=%d, In %s(), cannot get align items from unmated position, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Get align items from unmated position.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getAlignItemsFromSingleUnmatedPos(int32_t *matedNum, read_t *pRead, alignMatchItem_t *matchResultArray, int32_t *matchItemNum, read_t *pReadBased, alignMatchItem_t *matchItemBased, int32_t itemRowBased, int32_t mismatchNumThreshold, queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet)
{
	int32_t startQueryPos, endQueryPos, queryIDBased, queryPosBased, orientBased, validFlag;
	query_t *queryItem;
	uint64_t *readSeqInt;
	int32_t softClipEnd, orient, clipNum, gapNumQuery, gapNumRead, misNum, startAlignQueryPos, endAlignQueryPos, startAlignReadPos, endAlignReadPos;
	double insertSize, standDev, fragSize;

	char *queryAlignSeq, *readAlignSeq;
	int32_t queryAlignSeqLen, readAlignSeqLen, maxSeqLen;
	char *alignResultArray[3];
	int32_t overlapLen, mismatchNum, queryLeftShiftLen, queryRightShiftLen, readLeftShiftLen, readRightShiftLen;

	insertSize = readSet->insertSize;
	standDev = readSet->standDev;

	queryIDBased = matchItemBased->queryID;
	queryPosBased = matchItemBased->queryPos;
	orientBased = matchItemBased->orientation;

	queryItem = queryMatchInfoSet->queryArray + queryIDBased - 1;
	readSeqInt = readSet->readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;

	validFlag = NO;
	if(orientBased==ORIENTATION_PLUS)
	{
		startQueryPos = queryPosBased;
		endQueryPos = startQueryPos + insertSize + standDev * 6;
		if(endQueryPos>queryItem->queryLen)
		{
			endQueryPos = queryItem->queryLen;
			validFlag = YES;
		}
		queryAlignSeqLen = endQueryPos - startQueryPos + 1;
	}else
	{
		endQueryPos = queryPosBased + pReadBased->seqlen - matchItemBased->startReadPos;
		startQueryPos = endQueryPos - insertSize - standDev * 6;
		if(startQueryPos<1)
		{
			startQueryPos = 1;
			validFlag = YES;
		}
		queryAlignSeqLen = endQueryPos - startQueryPos + 1;
	}
	readAlignSeqLen = pRead->seqlen;

	if(validFlag==NO)
		return SUCCESSFUL;

	if(queryAlignSeqLen<readAlignSeqLen)
		maxSeqLen = 2 * readAlignSeqLen;
	else
		maxSeqLen = 2 * queryAlignSeqLen;

	// initialize the alignment buffer
	if(initAlignBuf(alignResultArray, &queryAlignSeq, &readAlignSeq, maxSeqLen)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the alignment buffers, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the two sequences
	strncpy(queryAlignSeq, queryItem->querySeq+startQueryPos-1, queryAlignSeqLen);
	queryAlignSeq[queryAlignSeqLen] = '\0';

	if(orientBased==ORIENTATION_MINUS)
	{
		getReadBaseFromPosByInt(readAlignSeq, readSeqInt, pRead->seqlen, 0, pRead->seqlen);
		orient = ORIENTATION_PLUS;
	}else
	{
		getReverseReadBaseFromPosByInt(readAlignSeq, readSeqInt, pRead->seqlen, 0, pRead->seqlen);
		orient = ORIENTATION_MINUS;
	}

	// generate alignment
	if(computeSeqAlignment(alignResultArray, &overlapLen, &mismatchNum, &queryLeftShiftLen, &readLeftShiftLen, &queryRightShiftLen, &readRightShiftLen, queryAlignSeq, readAlignSeq, queryAlignSeqLen, readAlignSeqLen, NO)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the alignment, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(mismatchNum<=mismatchNumThreshold && overlapLen>=0.6*pRead->seqlen)
	{
		// compute alignment gap number
		if(computeAlignGapNum(&gapNumQuery, &gapNumRead, &misNum, alignResultArray, overlapLen)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the alignment gap count, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(gapNumQuery==0 && gapNumRead==0)
		{
			if(readLeftShiftLen>0 && startQueryPos-(readLeftShiftLen-queryLeftShiftLen)<1)
			{ // clip at 5' end
				softClipEnd = 1;
				clipNum = readLeftShiftLen - queryLeftShiftLen;
				mismatchNum += (readLeftShiftLen - clipNum) + readRightShiftLen;
				startAlignReadPos = clipNum + 1;
				endAlignReadPos = pRead->seqlen;
				startAlignQueryPos = 1;
				endAlignQueryPos = startAlignQueryPos + pRead->seqlen - clipNum - 1;
			}else if(readRightShiftLen>0 && endQueryPos+(readRightShiftLen-queryRightShiftLen)>queryItem->queryLen)
			{ // clip at 3' end
				softClipEnd = 2;
				clipNum = readRightShiftLen - queryRightShiftLen;
				mismatchNum += readLeftShiftLen + (readRightShiftLen - clipNum);
				startAlignReadPos = 1;
				endAlignReadPos = pRead->seqlen - clipNum;
				startAlignQueryPos = endQueryPos - queryRightShiftLen - (pRead->seqlen - readRightShiftLen) + 1;
				endAlignQueryPos = endQueryPos;
			}else
			{ // no clip
				softClipEnd = -1;
				clipNum = 0;
				mismatchNum += readLeftShiftLen + readRightShiftLen;
				startAlignReadPos = 1;
				endAlignReadPos = pRead->seqlen;
				startAlignQueryPos = startQueryPos + (queryLeftShiftLen - readLeftShiftLen);
				endAlignQueryPos = startAlignQueryPos + pRead->seqlen - 1;
			}

			if(mismatchNum<=mismatchNumThreshold)
			{
				matchResultArray[*matchItemNum].queryID = queryIDBased;
				matchResultArray[*matchItemNum].queryPos = startAlignQueryPos;
				matchResultArray[*matchItemNum].mismatchNum = mismatchNum;
				matchResultArray[*matchItemNum].startReadPos = startAlignReadPos;
				matchResultArray[*matchItemNum].alignSize = endAlignReadPos - startAlignReadPos + 1;
				matchResultArray[*matchItemNum].orientation = orient;
				matchResultArray[*matchItemNum].validFlag = YES;

				matchResultArray[*matchItemNum].pairRow = itemRowBased;
				matchItemBased->pairRow = *matchItemNum;
				matchItemBased->validFlag = YES;

				if(orient==ORIENTATION_PLUS)
					fragSize = queryPosBased + pReadBased->seqlen - startAlignQueryPos;
				else
					fragSize = startAlignQueryPos + pReadBased->seqlen - queryPosBased;

				matchResultArray[*matchItemNum].fragSize = fragSize;
				matchItemBased->fragSize = fragSize;

				(*matchItemNum) ++;
				(*matedNum) ++;
			}
		}
	}

	// free the buffer
	freeAlignBuf(alignResultArray, &queryAlignSeq, &readAlignSeq);

	return SUCCESSFUL;
}

/**
 * Update the query coverage flag.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateQueryCovFlag(queryMatchInfo_t *queryMatchInfoSet, alignMatchItem_t *alignMatchItem, threadPara_t *threadPara)
{
	int32_t i, arrRow, arrCol, queryPos, queryLen, alignSize;
	uint64_t *covFlagArray, flagInt;
	query_t *queryItem;
	ctrlLock_t *lockItem;

	queryItem = queryMatchInfoSet->queryArray + (alignMatchItem->queryID - 1);
	queryLen = queryItem->queryLen;
	covFlagArray = queryItem->covFlagArray;
	queryPos = alignMatchItem->queryPos;
	alignSize = alignMatchItem->alignSize;

	arrRow = (queryPos - 1) / 64;
	arrCol = (queryPos - 1) % 64;

	lockItem = threadPara->lockArray + alignMatchItem->queryID - 1;

	// increase the writer count
	pthread_mutex_lock(&lockItem->accessWriterCnt);
	{
		lockItem->writerCnt ++;
		if(lockItem->writerCnt==1)
		{
			pthread_mutex_lock(&lockItem->readerLock);
		}
	}
	pthread_mutex_unlock(&lockItem->accessWriterCnt);

	// write the data
	pthread_mutex_lock(&lockItem->writeLock);
	{
		for(i=0; i<alignSize; i++)
		{
			flagInt = 1LLU << (63 - arrCol);
			covFlagArray[arrRow] |= flagInt;
			arrCol ++;
			if(arrCol==64)
			{
				arrCol = 0;
				arrRow ++;
			}

			if(queryPos+i>queryLen)
			{
				printf("line=%d, In %s(), basePos=%d, queryLen=%d, error!\n", __LINE__, __func__, queryPos+i, queryLen);
				return FAILED;
			}
		}
	}
	pthread_mutex_unlock(&lockItem->writeLock);

	// decrease the writer count
	pthread_mutex_lock(&lockItem->accessWriterCnt);
	{
		lockItem->writerCnt --;
		if(lockItem->writerCnt==0)
		{
			pthread_mutex_unlock(&lockItem->readerLock);
		}
	}
	pthread_mutex_unlock(&lockItem->accessWriterCnt);

	return SUCCESSFUL;
}

/**
 * Generate the kmer integer sequence from a read specified by a startPos (>=0).
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateKmerSeqIntFromReadset(uint64_t *seqInt, uint64_t *readseq, int32_t startReadPos, int32_t kmerSize, int32_t entriesNumRead, int32_t baseNumLastEntry)
{
	int32_t i, j, entriesNumKmer, startEntryPos, remainedBaseNum, rightRemainedNum, entryRow, baseNumInEntry;
	uint64_t *readseqStart;

	entriesNumKmer = ((kmerSize - 1) >> 5) + 1;
	for(i=0; i<entriesNumKmer; i++) seqInt[i] = 0;

	entryRow = startReadPos >> 5;
	readseqStart = readseq + entryRow;
	startEntryPos = startReadPos % 32;

	i = 0;
	j = 0;
	remainedBaseNum = kmerSize;
	while(remainedBaseNum>0)
	{
		// process first entry
		if(entryRow!=entriesNumRead-1)
			baseNumInEntry = 32;
		else
			baseNumInEntry = baseNumLastEntry;

		rightRemainedNum = baseNumInEntry - remainedBaseNum - startEntryPos;
		if(rightRemainedNum<0)
			rightRemainedNum = 0;

		seqInt[i] = (readseqStart[j] << (2*startEntryPos)) >> (2*(startEntryPos+rightRemainedNum));
		remainedBaseNum -= baseNumInEntry - rightRemainedNum - startEntryPos;

		if(startEntryPos>remainedBaseNum)
			startEntryPos = remainedBaseNum;

		j++;
		entryRow ++;

		// process second entry
		if(startEntryPos>0)
		{
			if(entryRow!=entriesNumRead-1)
				baseNumInEntry = 32;
			else
				baseNumInEntry = baseNumLastEntry;

			seqInt[i] = (seqInt[i] << (2*startEntryPos)) | (readseqStart[j] >> (2*(baseNumInEntry-startEntryPos)));
			remainedBaseNum -= startEntryPos;
		}

		i++;
	}

	return SUCCESSFUL;
}

/**
 * Get the row range of the reads match information array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getRowRangeMatchArray(int32_t *startRow, int32_t *endRow, alignMatchItem_t *matchResultArray, int32_t arraySize, int32_t queryID)
{
	int32_t left, right, middle, existFlag;

	left = 0;
	right = arraySize - 1;
	existFlag = NO;
	while(left<=right)
	{
		middle = (left + right) / 2;
		if(matchResultArray[middle].queryID==queryID)
		{
			existFlag = YES;
			break;
		}
		if(matchResultArray[middle].queryID<queryID)
			left = middle + 1;
		else
			right = middle - 1;
	}

	if(existFlag==YES)
	{
		*startRow = middle;
		while((*startRow)>0 && matchResultArray[(*startRow)-1].queryID==queryID)
			(*startRow) --;
		*endRow = middle;
		while((*endRow)<arraySize-1 && matchResultArray[(*endRow)+1].queryID==queryID)
			(*endRow) ++;
	}else
	{
		*startRow = -1;
		*endRow = -1;
	}

	return SUCCESSFUL;
}

/**
 * Fill the reads match information to queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillReadMatchInfoQueries(queryMatchInfo_t *queryMatchInfoSet, readSetArr_t *readSetArray)
{
	int32_t i;

	// initialize the memory for query reads
	if(initQueryReadSets(queryMatchInfoSet, readSetArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory for query reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<readSetArray->readSetNum; i++)
	{
		if(fillQueryReadSet(queryMatchInfoSet, readSetArray->readSetArray+i)==FAILED)
		{
			printf("line=%d, In %s(), cannot fill query reads, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// sort the query reads
	if(sortQueryReadSets(queryMatchInfoSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot sort query reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Initialize the memory for query reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initQueryReadSets(queryMatchInfo_t *queryMatchInfoSet, readSetArr_t *readSetArray)
{
	int32_t i, j, k, queryReadSetNum;
	queryReadSet_t *queryReadSetArray;
	readSet_t *readSet;
	readBlock_t *readBlock;
	readMatchInfoBlock_t *readMatchInfoBlock;
	readMatchInfo_t *pReadMatchInfo;
	query_t *queryItem;

	queryReadSetNum = readSetArray->readSetNum;

	// compute the read count for each query
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		queryReadSetArray = (queryReadSet_t*) calloc (queryReadSetNum, sizeof(queryReadSet_t));
		if(queryReadSetArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		for(j=0; j<queryReadSetNum; j++)
		{
			queryReadSetArray[j].queryReadArray = NULL;
			queryReadSetArray[j].queryReadNum = 0;
			queryReadSetArray[j].setID = j + 1;
		}
		queryMatchInfoSet->queryArray[i].queryReadSetArray = queryReadSetArray;
		queryMatchInfoSet->queryArray[i].queryReadSetNum = queryReadSetNum;
	}

	for(i=0; i<queryReadSetNum; i++)
	{
		readSet = readSetArray->readSetArray + i;
		for(j=0; j<readSet->blocksNumRead; j++)
		{
			readBlock = readSet->readBlockArr + j;
			readMatchInfoBlock = readSet->readMatchInfoBlockArr + j;
			for(k=0; k<readBlock->itemNum; k++)
			{
				pReadMatchInfo = readMatchInfoBlock->readMatchInfoArr + k;
				if(pReadMatchInfo->queryID>0)
				{
					queryItem = queryMatchInfoSet->queryArray + pReadMatchInfo->queryID - 1;
					queryItem->queryReadSetArray[i].queryReadNum ++;
				}
			}
		}
	}

	// allocate memory
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		queryItem = queryMatchInfoSet->queryArray + i;
		for(j=0; j<queryItem->queryReadSetNum; j++)
		{
			if(queryItem->queryReadSetArray[j].queryReadNum>0)
			{
				queryItem->queryReadSetArray[j].queryReadArray = (queryRead_t *) calloc (queryItem->queryReadSetArray[j].queryReadNum, sizeof(queryRead_t));
				if(queryItem->queryReadSetArray[j].queryReadArray==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				queryItem->queryReadSetArray[j].queryReadNum = 0;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Fill the query reads set.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillQueryReadSet(queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet)
{
	int32_t i, j, k;
	readBlock_t *readBlock;
	readMatchInfoBlock_t *readMatchInfoBlock;
	readMatchInfo_t *pReadMatchInfo;
	queryReadSet_t *queryReadSet;
	int32_t queryReadSetID;

	queryReadSetID = readSet->setID;

	// fill the match information
	for(j=0; j<readSet->blocksNumRead; j++)
	{
		readBlock = readSet->readBlockArr + j;
		readMatchInfoBlock = readSet->readMatchInfoBlockArr + j;
		for(k=0; k<readBlock->itemNum; k++)
		{
			pReadMatchInfo = readMatchInfoBlock->readMatchInfoArr + k;
			if(pReadMatchInfo->queryID>0)
			{
				queryReadSet = queryMatchInfoSet->queryArray[pReadMatchInfo->queryID-1].queryReadSetArray + queryReadSetID - 1;
				queryReadSet->queryReadArray[queryReadSet->queryReadNum].readID = pReadMatchInfo->readID;
				queryReadSet->queryReadArray[queryReadSet->queryReadNum].queryPos = pReadMatchInfo->queryPos;
				queryReadSet->queryReadArray[queryReadSet->queryReadNum].seqlen = pReadMatchInfo->seqlen;
				queryReadSet->queryReadArray[queryReadSet->queryReadNum].orientation = pReadMatchInfo->readOrientation;
				queryReadSet->queryReadArray[queryReadSet->queryReadNum].startReadPos = pReadMatchInfo->startReadPos;
				queryReadSet->queryReadArray[queryReadSet->queryReadNum].alignSize = pReadMatchInfo->alignSize;
				queryReadSet->queryReadNum ++;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Sort the query reads sets.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short sortQueryReadSets(queryMatchInfo_t *queryMatchInfoSet)
{
	int32_t i, j, maxReadsNum;
	readBlock_t *readBlock;
	readMatchInfoBlock_t *readMatchInfoBlock;
	readMatchInfo_t *pReadMatchInfo;
	query_t *queryItem;
	queryReadSet_t *queryReadSetArray;
	queryRead_t *queryReadArrayBuf;
	int32_t queryReadSetNum;

	// get the maximal array size and allocate the buffer memory
	maxReadsNum = 0;
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		queryItem = queryMatchInfoSet->queryArray + i;
		queryReadSetArray = queryItem->queryReadSetArray;
		queryReadSetNum = queryItem->queryReadSetNum;
		for(j=0; j<queryReadSetNum; j++)
		{
			if(maxReadsNum<queryReadSetArray[j].queryReadNum)
				maxReadsNum = queryReadSetArray[j].queryReadNum;
		}
	}
	if(maxReadsNum<=0)
	{
		printf("line=%d, In %s(), invalid reads count %d, error!\n", __LINE__, __func__, maxReadsNum);
		return FAILED;
	}

	queryReadArrayBuf = (queryRead_t*) calloc(maxReadsNum, sizeof(queryRead_t));
	if(queryReadArrayBuf==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// radix sort the match information
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		queryItem = queryMatchInfoSet->queryArray + i;
		queryReadSetArray = queryItem->queryReadSetArray;
		queryReadSetNum = queryItem->queryReadSetNum;
		for(j=0; j<queryReadSetNum; j++)
		{
			if(queryReadSetArray[j].queryReadNum>0)
			{
				if(radixSortQueryReadArray(queryReadSetArray[j].queryReadArray, queryReadArrayBuf, queryReadSetArray[j].queryReadNum)==FAILED)
				{
					printf("line=%d, In %s(), cannot sort the query read array, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}

	// free queryReadArrayBuf
	free(queryReadArrayBuf);

	return SUCCESSFUL;
}

/**
 * Radix sort the reads match information on queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short radixSortQueryReadArray(queryRead_t *queryReadArray, queryRead_t *queryReadArrayBuf, int32_t itemNum)
{
	struct partNode
	{
		int curItemNum;
		int totalItemNum;
		int firstRow;
	};

	int64_t i, step, total;
	queryRead_t *data, *buf;
	struct partNode *part;
	int partArrSize, stepBits, maxStepLen;
	uint64_t bitMask;
	uint64_t hashcode, firstRow, curItemNum;

	stepBits = 16;
	maxStepLen = 32;
	partArrSize = 1 << stepBits;
	bitMask = (1 << stepBits) - 1;

	part = (struct partNode *) malloc(partArrSize * sizeof(struct partNode));
	if(part==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Begin to sort
	step = 0;
	while(step!=maxStepLen)
	{
		// set the data and buf
		if(step==stepBits)
		{
			buf = queryReadArray;
			data = queryReadArrayBuf;
		}else
		{
			data = queryReadArray;
			buf = queryReadArrayBuf;
		}

		// count the number
		if(memset(part, 0, partArrSize * sizeof(struct partNode))==NULL)
		{
			printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		for(i=0; i<itemNum; i++)
		{
			//part[ bitMask - ((data[i].queryPos >> step) & bitMask) ].totalItemNum ++;  // from big to small
			part[ (data[i].queryPos >> step) & bitMask ].totalItemNum ++;  // from small to big
		}

		// initialize the part array
		total = 0;
		for(i=0; i<partArrSize; i++)
		{
			part[i].firstRow = total;
			total += part[i].totalItemNum;
		}

		// copy the data to the right place
		for(i=0; i<itemNum; i++)
		{
			//hashcode = bitMask - ((data[i].queryPos >> step) & bitMask);  // from big to small
			hashcode = (data[i].queryPos >> step) & bitMask;  // from small to big
			firstRow = part[hashcode].firstRow;
			curItemNum = part[hashcode].curItemNum;
			if(memcpy(buf+firstRow+curItemNum, data+i, sizeof(queryRead_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			part[hashcode].curItemNum ++;
		}

		step += stepBits;
	}

	free(part);
	part = NULL;

	return SUCCESSFUL;
}

/**
 * Sort match information according to mismatch.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short sortMatchResults(alignMatchItem_t *matchResultArray, alignMatchItem_t *matchResultArrayBuf, int32_t matchItemNum)
{
	if(matchItemNum<20)
	{ // selection sort
		if(selectionSortMatchResults(matchResultArray, matchItemNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot do the selection sort, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{ // radix sort
		if(radixSortMatchResults(matchResultArray, matchResultArrayBuf, matchItemNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot do the radix sort, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Selection sort match information according to mismatch.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short selectionSortMatchResults(alignMatchItem_t *matchResultArray, int32_t matchItemNum)
{
	int32_t i, j, minValue, minRow;
	alignMatchItem_t tmp;

	for(i=0; i<matchItemNum; i++)
	{
		minValue = matchResultArray[i].mismatchNum;
		minRow = i;
		for(j=i+1; j<matchItemNum; j++)
		{
			if(minValue>matchResultArray[j].mismatchNum)
			{
				minValue = matchResultArray[j].mismatchNum;
				minRow = j;
			}
		}

		if(minRow!=i)
		{
			if(memcpy(&tmp, matchResultArray+i, sizeof(alignMatchItem_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(memcpy(matchResultArray+i, matchResultArray+minRow, sizeof(alignMatchItem_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(memcpy(matchResultArray+minRow, &tmp, sizeof(alignMatchItem_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Radix sort match information according to mismatch.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short radixSortMatchResults(alignMatchItem_t *matchResultArray, alignMatchItem_t *matchResultArrayBuf, int32_t itemNum)
{
	struct partNode
	{
		int curItemNum;
		int totalItemNum;
		int firstRow;
	};

	int64_t i, step, total;
	alignMatchItem_t *data, *buf;
	struct partNode *part;
	int partArrSize, stepBits, maxStepLen;
	uint64_t bitMask;
	uint64_t hashcode, firstRow, curItemNum;

	stepBits = 16;
	maxStepLen = 32;
	partArrSize = 1 << stepBits;
	bitMask = (1 << stepBits) - 1;

	part = (struct partNode *) malloc(partArrSize * sizeof(struct partNode));
	if(part==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Begin to sort
	step = 0;
	while(step!=maxStepLen)
	{
		// set the data and buf
		if(step==stepBits)
		{
			buf = matchResultArray;
			data = matchResultArrayBuf;
		}else
		{
			data = matchResultArray;
			buf = matchResultArrayBuf;
		}

		// count the number
		if(memset(part, 0, partArrSize * sizeof(struct partNode))==NULL)
		{
			printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		for(i=0; i<itemNum; i++)
		{
			//part[ bitMask - ((data[i].mismatchNum >> step) & bitMask) ].totalItemNum ++;  // from big to small
			part[ (data[i].mismatchNum >> step) & bitMask ].totalItemNum ++;  // from small to big
		}

		// initialize the part array
		total = 0;
		for(i=0; i<partArrSize; i++)
		{
			part[i].firstRow = total;
			total += part[i].totalItemNum;
		}

		// copy the data to the right place
		for(i=0; i<itemNum; i++)
		{
			//hashcode = bitMask - ((mismatchNum >> step) & bitMask);  // from big to small
			hashcode = (data[i].mismatchNum >> step) & bitMask;  // from small to big
			firstRow = part[hashcode].firstRow;
			curItemNum = part[hashcode].curItemNum;
			if(memcpy(buf+firstRow+curItemNum, data+i, sizeof(alignMatchItem_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			part[hashcode].curItemNum ++;
		}

		step += stepBits;
	}

	free(part);
	part = NULL;

	return SUCCESSFUL;
}

/**
 * Output the reads match information on queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputQueryReadArray(char *outfile, queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet)
{
	int32_t i, j, queryID, queryLen, queryPos, queryPos_paired, fragSize, itemNum, matedFlag, uniqueMapFlag;
	int64_t readID, readID_paired;
	double insertSize, standDev, fragSizeDif;
	char orientation, outFileName[256];
	query_t *queryItem;
	queryRead_t *queryReadArray;
	int32_t queryReadSetNum;
	FILE *fpOut;

	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo;
	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;

	queryReadSetNum = 1;

	for(i=0; i<queryReadSetNum; i++)
	{
		strcpy(outFileName, outfile);
		sprintf(outFileName+strlen(outFileName), "%d", i);

		fpOut = fopen(outFileName, "w");
		if(fpOut==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, outFileName);
			return FAILED;
		}

		readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
		maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

		insertSize = readSet->insertSize;
		standDev = readSet->standDev;

		for(j=0; j<queryMatchInfoSet->itemNumQueryArray; j++)
		{
			queryItem = queryMatchInfoSet->queryArray + j;
			queryID = queryItem->queryID;
			queryLen = queryItem->queryLen;

			if(queryItem->queryReadSetArray[i].queryReadNum>0)
			{
				queryReadArray = queryItem->queryReadSetArray[i].queryReadArray;
				itemNum = queryItem->queryReadSetArray[i].queryReadNum;
				fprintf(fpOut, ">%s\t%d\t0\t%d\n", queryItem->queryTitle, queryLen, itemNum);
				for(j=0; j<itemNum; j++)
				{
					readID = queryReadArray[j].readID;
					queryPos = queryReadArray[j].queryPos;

					if(queryReadArray[j].orientation==ORIENTATION_PLUS)
						orientation = '+';
					else
						orientation = '-';

					// get uniqueMap flag
					if(getUniqueMapFlag(&uniqueMapFlag, readID, readSet)==FAILED)
					{
						printf("line=%d, In %s(), cannot get unique map flag, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(readID%2==1)
						readID_paired = readID + 1;
					else
						readID_paired = readID - 1;

					readMatchInfoBlockID = (readID_paired - 1) / maxItemNumPerReadMatchInfoBlock;
					rowNumInReadMatchInfoBlock = (readID_paired - 1) % maxItemNumPerReadMatchInfoBlock;
					pReadMatchInfo = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;

					if(pReadMatchInfo->queryID==queryID && queryReadArray[j].orientation!=pReadMatchInfo->readOrientation)
					{
						queryPos_paired = pReadMatchInfo->queryPos;

						if(queryReadArray[j].orientation==ORIENTATION_PLUS)
							fragSize = queryPos_paired + pReadMatchInfo->seqlen - queryPos;
						else
							fragSize = queryPos + queryReadArray[j].seqlen - queryPos_paired;

						fragSizeDif = fragSize - insertSize;
						if(fragSizeDif>-5*standDev && fragSizeDif<5*standDev)
							matedFlag = YES;
						else
							matedFlag = NO;

						if(queryReadArray[j].orientation==ORIENTATION_MINUS)
							fragSize = - fragSize;
					}else
					{
						matedFlag = NO;
					}

					if(matedFlag==YES)
						fprintf(fpOut, "%d\t%ld\t%c\t%d\t%d\t%d\t%d\t%d\n", queryReadArray[j].queryPos, (int64_t)queryReadArray[j].readID, orientation, (int32_t)queryReadArray[j].startReadPos, (int32_t)queryReadArray[j].alignSize, uniqueMapFlag, (int32_t)queryReadArray[j].seqlen, fragSize);
					else
						fprintf(fpOut, "%d\t%ld\t%c\t%d\t%d\t%d\t%d\n", queryReadArray[j].queryPos, (int64_t)queryReadArray[j].readID, orientation, (int32_t)queryReadArray[j].startReadPos, (int32_t)queryReadArray[j].alignSize, uniqueMapFlag, (int32_t)queryReadArray[j].seqlen);
				}
			}
		}

		fclose(fpOut);
	}

	return SUCCESSFUL;
}

/**
 * Get the unique map flag of a read.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getUniqueMapFlag(int32_t *uniqueMapFlag, int64_t readID, readSet_t *readSet)
{
	int32_t readBlockID, rowNumInReadBlock, maxItemNumPerReadBlock;
	read_t *pRead;

	maxItemNumPerReadBlock = readSet->maxItemNumPerReadBlock;

	readBlockID = (readID - 1) / maxItemNumPerReadBlock;
	rowNumInReadBlock = (readID - 1) % maxItemNumPerReadBlock;
	pRead = readSet->readBlockArr[readBlockID].readArr + rowNumInReadBlock;

	if(pRead->validFlag==YES && pRead->successMapFlag==YES)
		*uniqueMapFlag = pRead->uniqueMapFlag;
	else
		*uniqueMapFlag = -1;

	return SUCCESSFUL;
}

/**
 * Get the paired read match information of a read.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getPairedMatchInfo(readMatchInfo_t **pReadMatchInfoPaired, int64_t readID, readSet_t *readSet)
{
	readMatchInfoBlock_t *readMatchInfoBlockArr;
	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;
	int64_t readID_paired;

	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

	if(readID%2==1)
		readID_paired = readID + 1;
	else
		readID_paired = readID - 1;

	readMatchInfoBlockID = (readID_paired - 1) / maxItemNumPerReadMatchInfoBlock;
	rowNumInReadMatchInfoBlock = (readID_paired - 1) % maxItemNumPerReadMatchInfoBlock;
	*pReadMatchInfoPaired = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;

	return SUCCESSFUL;
}
