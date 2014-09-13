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
short mapReads(queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet, queryIndex_t *queryIndex, double insertSize, double standDev)
{
	struct timeval tpstart, tpend;
	double timeused_align;
	gettimeofday(&tpstart, NULL);

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

	// map reads uniquely
	if(mapReadsOp(queryMatchInfoSet, readSet, queryIndex, insertSize, standDev, YES)==FAILED)
	{
		printf("line=%d, In %s(), cannot map reads uniquely, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// map reads non-uniquely
	if(mapReadsOp(queryMatchInfoSet, readSet, queryIndex, insertSize, standDev, NO)==FAILED)
	{
		printf("line=%d, In %s(), cannot map reads uniquely, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// free coverage flag array
	freeCovFlagArray(queryMatchInfoSet);

	// calculate aligning time
	gettimeofday(&tpend, NULL);
	timeused_align = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Aligning reads used time: %.2f seconds.\n", timeused_align);

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
 * Start the reads map.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mapReadsOp(queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet, queryIndex_t *queryIndex, double insertSize, double standDev, int32_t uniqueMapOpFlag)
{
	int32_t i, j, k, seqLen, maxArraySize, *matchItemNum, matchItemNums[2];
	readBlock_t *pReadBlockArray;
	read_t *pRead, *pReadArray;
	uint64_t rid, tmpRid, *readSeqInt, *readSeqIntRev;
	alignMatchItem_t *matchResultArray, *matchResultArrays[2], *matchResultArrayBuf1, *matchResultArrayBuf2;
	char readseqTmp[MAX_READ_LEN_IN_BUF+1];
	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo;
	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;

	if(uniqueMapOpFlag==YES)
		printf("Aligning reads uniquely ...\n");
	else
		printf("Aligning reads of multiple occurrences ...\n");

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

	matchResultArrays[0] = (alignMatchItem_t*) calloc (4*maxArraySize, sizeof(alignMatchItem_t));
	matchResultArrays[1] = (alignMatchItem_t*) calloc (4*maxArraySize, sizeof(alignMatchItem_t));
	matchResultArrayBuf1 = (alignMatchItem_t*) calloc (4*maxArraySize, sizeof(alignMatchItem_t));
	matchResultArrayBuf2 = (alignMatchItem_t*) calloc (4*maxArraySize, sizeof(alignMatchItem_t));
	if(matchResultArrays[0]==NULL || matchResultArrays[1]==NULL || matchResultArrayBuf1==NULL || matchResultArrayBuf2==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

	// map unique reads
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

//			if(rid==1320867)
//			{
//				printf("rid=%ld, seqLen=%d\n", rid, seqLen);
//			}

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
				if(mapSingleReadToQueries(rid, pRead, readSet, matchResultArray, matchResultArrayBuf1, matchResultArrayBuf2, matchItemNum, queryIndex, queryMatchInfoSet->queryArray, 5)==FAILED)
				{
					printf("line=%d, In %s(), cannot map the read %lu, error!\n", __LINE__, __func__, rid);
					return FAILED;
				}
			}

			if(rid%2==0)
			{
				// select the mated match information
				if(matchItemNums[0]>0 && matchItemNums[1]>0 && (matchItemNums[0]>=2 || matchItemNums[1]>=2))
				{
					if(selectMatedMatchPosPE(matchResultArrays[0], matchResultArrays[1], matchItemNums, matchItemNums+1, (pRead-1)->seqlen, pRead->seqlen, insertSize, standDev, uniqueMapOpFlag, queryMatchInfoSet)==FAILED)
					{
						printf("line=%d, In %s(), cannot adjust the match information of paired-ends, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}else if((matchItemNums[0]>=2 && matchItemNums[1]==0) || (matchItemNums[0]==0 && matchItemNums[1]>=2))
				{
					// select best alignment for single read
					if(selectBestAlignInfoSingleRead(matchResultArrays[0], matchItemNums, uniqueMapOpFlag, queryMatchInfoSet)==FAILED)
					{
						printf("line=%d, In %s(), cannot select the best alignment information for single read, error!\n", __LINE__, __func__);
						return FAILED;
					}

					// select best alignment for single read
					if(selectBestAlignInfoSingleRead(matchResultArrays[1], matchItemNums+1, uniqueMapOpFlag, queryMatchInfoSet)==FAILED)
					{
						printf("line=%d, In %s(), cannot select the best alignment information for single read, error!\n", __LINE__, __func__);
						return FAILED;
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
						readSet->totalValidItemNumReadMatchInfo ++;

						pReadMatchInfo->queryID = matchResultArray[0].queryID;
						pReadMatchInfo->queryPos = matchResultArray[0].queryPos;
						pReadMatchInfo->alignSize = matchResultArray[0].alignSize;
						pReadMatchInfo->startReadPos = matchResultArray[0].startReadPos;
						pReadMatchInfo->readID = tmpRid;
						pReadMatchInfo->seqlen = pRead->seqlen;
						pReadMatchInfo->readOrientation = matchResultArray[0].orientation;

						pRead->successMapFlag = YES;
						pRead->uniqueMapFlag = uniqueMapOpFlag;
						queryMatchInfoSet->queryArray[pReadMatchInfo->queryID-1].queryReadNum ++;

						// update coverage flag
						if(updateQueryCovFlag(queryMatchInfoSet, matchResultArray)==FAILED)
						{
							printf("line=%d, In %s(), cannot update query coverage flag, error!\n", __LINE__, __func__);
							return FAILED;
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

	return SUCCESSFUL;
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
short estimateInsertSize(double *insertSize, double *standDev, readSet_t *readSet, queryIndex_t *queryIndex)
{
	int32_t i, j, k, seqLen, maxArraySize, *matchItemNum, matchItemNums[2];
	readBlock_t *pReadBlockArray;
	read_t *pRead, *pReadArray;
	uint64_t rid;
	alignMatchItem_t *matchResultArray, *matchResultArrays[2], *matchResultArrayBuf1, *matchResultArrayBuf2;
	double insertSizeTmp, standDevTmp, fragSizeSum, fragSize, fragDif, fragDifSum;
	int32_t pairNum, pairNumTmp, pairNumEnoughFlag, uniqueMapOpFlag;

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

	matchResultArrays[0] = (alignMatchItem_t*) calloc (2*maxArraySize, sizeof(alignMatchItem_t));
	matchResultArrays[1] = (alignMatchItem_t*) calloc (2*maxArraySize, sizeof(alignMatchItem_t));
	matchResultArrayBuf1 = (alignMatchItem_t*) calloc (maxArraySize, sizeof(alignMatchItem_t));
	matchResultArrayBuf2 = (alignMatchItem_t*) calloc (maxArraySize, sizeof(alignMatchItem_t));
	if(matchResultArrays[0]==NULL || matchResultArrays[1]==NULL || matchResultArrayBuf1==NULL || matchResultArrayBuf2==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// estimate the insert size and standard deviation
	fragSizeSum = 0;
	pairNum = 0;
	uniqueMapOpFlag = YES;

	for(k=0; k<2; k++)
	{
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
								}else
								{
									fragDif = (fragSize - insertSizeTmp) * (fragSize - insertSizeTmp);
									fragDifSum += fragDif;
								}
								pairNumTmp ++;
								if(pairNumTmp>=2000)
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
			{
				pairNum = pairNumTmp;
				insertSizeTmp = fragSizeSum / pairNum;
			}else
			{
				printf("line=%d, In %s(), no valid paired reads for computing insert size of library, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else
		{
			if(pairNumTmp!=pairNum)
			{
				printf("line=%d, In %s(), cannot compute insert size of library, error!\n", __LINE__, __func__);
				return FAILED;
			}
			standDevTmp = sqrt((fragDifSum) / pairNum);
		}
	}

	*insertSize = insertSizeTmp;
	*standDev = standDevTmp;

	free(matchResultArrays[0]);
	free(matchResultArrays[1]);
	free(matchResultArrayBuf1);
	free(matchResultArrayBuf2);

	printf("insertSize=%.4f, SDev=%.4f\n", *insertSize, *standDev);

	return SUCCESSFUL;
}

/**
 * Perfectly map single read to queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mapSingleReadToQueries(int64_t rid, read_t *pRead, readSet_t *readSet, alignMatchItem_t *matchResultArray, alignMatchItem_t *matchResultArrayBuf1, alignMatchItem_t *matchResultArrayBuf2, int32_t *matchItemNum, queryIndex_t *queryIndex, query_t *queryArray, int32_t mismatchNumThreshold)
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
	int32_t matchDistance, matchFlag, startRow, endRow;

	entriesNumRead = ((seqLen - 1) >> 5) + 1;
	baseNumLastEntryRead = ((seqLen - 1) % 32) + 1;

	mapTimes = (seqLen - 1) / queryIndex->kmerSize + 1;
	matchItemNum1 = 0;
	newMatchItemNum1 = 0;
	if(seqLen>=queryIndex->kmerSize)
	{
		startBasePos = 0;

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
	int32_t i, seqLen, entriesNumRead, tmpItemNum;
	uint64_t *readSeqInt, *readSeqIntRev;

	readSeqInt = readSet->readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
	seqLen = pRead->seqlen;

	entriesNumRead = ((seqLen - 1) >> 5) + 1;

	*matchItemNum = 0;
	if(seqLen>=queryIndex->kmerSize)
	{
		if(getMatchedQueryPosWithMismatch(matchResultArray, &tmpItemNum, readSeqInt, seqLen, queryIndex, queryArray, mismatchNumThreshold)==FAILED)
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
		if(getMatchedQueryPosWithMismatch(matchResultArray+(*matchItemNum), &tmpItemNum, readSeqIntRev, seqLen, queryIndex, queryArray, mismatchNumThreshold)==FAILED)
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
short getMatchedQueryPosWithMismatch(alignMatchItem_t *matchResultArray, int32_t *matchItemNum, uint64_t *readSeqInt, int32_t seqLen, queryIndex_t *queryIndex, query_t *queryArray, int32_t mismatchNumThreshold)
{
	int32_t i, j, entriesNumRead, baseNumLastEntryRead, startBasePos, processedFlag;
	uint64_t hashcode, kmerSeqInt[queryIndex->entriesPerKmer];
	queryKmer_t *kmer;
	char *querySeq, readSeq[MAX_READ_LEN_IN_BUF+1];
	int32_t queryID, queryLen, queryPos, maxStartBasePos, remainBaseNum, misNum, startCmpRpos, startCmpQpos;
	int32_t startAlignReadPos, endAlignReadPos, startAlignQueryPos, queryEndIgnoreLen;

	entriesNumRead = ((seqLen - 1) >> 5) + 1;
	baseNumLastEntryRead = ((seqLen - 1) % 32) + 1;
	queryEndIgnoreLen = 15;

	*matchItemNum = 0;
	if(seqLen>=queryIndex->kmerSize)
	{
		// get the read base sequence
		getReadBaseFromPosByInt(readSeq, readSeqInt, seqLen, 0, seqLen);

		startBasePos = 0;
		maxStartBasePos = seqLen - queryIndex->kmerSize;
		while(startBasePos<=maxStartBasePos)
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
							// mismatch on right side
							misNum = 0;
							startCmpRpos = startBasePos + queryIndex->kmerSize;
							startCmpQpos = queryPos + queryIndex->kmerSize - 1;
							if(queryLen-startCmpQpos<seqLen-startCmpRpos)
								remainBaseNum = queryLen - startCmpQpos;
							else
								remainBaseNum = seqLen - startCmpRpos;
							for(j=0; j<remainBaseNum; j++)
							{
								if(readSeq[startCmpRpos+j]!=querySeq[startCmpQpos+j])
								{
									misNum ++;
									if(misNum>mismatchNumThreshold)
										break;
								}
							}

							// mismatch on left side
							startCmpRpos = startAlignReadPos - 1;
							if(startAlignReadPos==1)
							{
								startCmpQpos = queryPos - startBasePos - 1;
								remainBaseNum = startBasePos;
							}else
							{
								startCmpQpos = 0;
								remainBaseNum = startBasePos + 1 - startAlignReadPos;
							}

							for(j=0; j<remainBaseNum; j++)
							{
								if(readSeq[startCmpRpos+j]!=querySeq[startCmpQpos+j])
								{
									misNum ++;
									if(misNum>mismatchNumThreshold)
										break;
								}
							}

							if(misNum<=mismatchNumThreshold)
							{ // add the match information
								matchResultArray[*matchItemNum].queryID = queryID;
								matchResultArray[*matchItemNum].queryPos = startAlignQueryPos;
								matchResultArray[*matchItemNum].mismatchNum = misNum;
								matchResultArray[*matchItemNum].startReadPos = startAlignReadPos;
								matchResultArray[*matchItemNum].alignSize = endAlignReadPos - startAlignReadPos + 1;
								matchResultArray[*matchItemNum].validFlag = YES;
								(*matchItemNum) ++;
							}
						}
					}
				}
			}

			startBasePos += queryIndex->kmerSize;
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
 * Select the mated query position of paired reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short selectMatedMatchPosPE(alignMatchItem_t *matchResultArray1, alignMatchItem_t *matchResultArray2, int32_t *matchItemNum1, int32_t *matchItemNum2, int32_t seqLen1, int32_t seqLen2, double insertSize, double standDev, int32_t uniqueMapOpFlag, queryMatchInfo_t *queryMatchInfoSet)
{
	struct pairMismatch
	{
		int32_t arrRow1, arrRow2, mismatchNum1, mismatchNum2, totalMismatchNum, fragSize, fragDif;
	};

	int32_t i, j, queryID1, queryPos1, queryID2, queryPos2, orient1, orient2, matedNum, matedRow1, matedRow2, pairRow;
	double fragSize, fragDif, fragDifExist;
	int32_t *zeroCovNumArray, zeroNum, itemRow, maxValue, maxRow, maxNum;

	struct pairMismatch *pairMismatchArray;
	int32_t minValue, minRow, minNum, itemNum, arrRow1, arrRow2, minID;

	for(i=0; i<(*matchItemNum1); i++) matchResultArray1[i].pairRow = -1;
	for(i=0; i<(*matchItemNum2); i++) matchResultArray2[i].pairRow = -1;

	matedNum = 0;
	matedRow1 = matedRow2 = -1;
	for(i=0; i<(*matchItemNum1); i++)
	{
		queryID1 = matchResultArray1[i].queryID;
		queryPos1 = matchResultArray1[i].queryPos;
		orient1 = matchResultArray1[i].orientation;

		for(j=0; j<(*matchItemNum2); j++)
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

						matedNum ++;
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


	if(matedNum==1)
	{ // unique mated
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
						if(getZeroCovBaseNumSingleRead(&zeroNum, matchResultArray1+itemRow, queryMatchInfoSet)==FAILED)
						{
							printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
						zeroCovNumArray[i] += zeroNum;

						itemRow = pairMismatchArray[i].arrRow2;
						if(getZeroCovBaseNumSingleRead(&zeroNum, matchResultArray2+itemRow, queryMatchInfoSet)==FAILED)
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
		if(selectBestAlignInfoSingleRead(matchResultArray1, matchItemNum1, uniqueMapOpFlag, queryMatchInfoSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot select the best alignment information for single read, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// select best alignment for single read
		if(selectBestAlignInfoSingleRead(matchResultArray2, matchItemNum2, uniqueMapOpFlag, queryMatchInfoSet)==FAILED)
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
short getZeroCovBaseNumSingleRead(int32_t *zeroCovBaseNum, alignMatchItem_t *alignMatchItem, queryMatchInfo_t *queryMatchInfoSet)
{
	int32_t i, queryID, queryPos, queryLen, alignSize, arrRow, arrCol, baseCovFlag;
	uint64_t *covFlagArray;
	query_t *queryItem;

	queryID = alignMatchItem->queryID;
	queryPos = alignMatchItem->queryPos;
	queryItem = queryMatchInfoSet->queryArray + (queryID - 1);
	queryLen = queryItem->queryLen;
	covFlagArray = queryItem->covFlagArray;
	alignSize = alignMatchItem->alignSize;

	arrRow = (queryPos - 1) / 64;
	arrCol = (queryPos - 1) % 64;

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

	return SUCCESSFUL;
}

/**
 * Select best alignment for single read.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short selectBestAlignInfoSingleRead(alignMatchItem_t *matchResultArray, int32_t *matchItemNum, int32_t uniqueMapOpFlag, queryMatchInfo_t *queryMatchInfoSet)
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
						if(getZeroCovBaseNumSingleRead(zeroCovNumArray+i, matchResultArray+i, queryMatchInfoSet)==FAILED)
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
 * Update the query coverage flag.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateQueryCovFlag(queryMatchInfo_t *queryMatchInfoSet, alignMatchItem_t *alignMatchItem)
{
	int32_t i, arrRow, arrCol, queryPos, queryLen, alignSize;
	uint64_t *covFlagArray, flagInt;
	query_t *queryItem;

	queryItem = queryMatchInfoSet->queryArray + (alignMatchItem->queryID - 1);
	queryLen = queryItem->queryLen;
	covFlagArray = queryItem->covFlagArray;
	queryPos = alignMatchItem->queryPos;
	alignSize = alignMatchItem->alignSize;

	arrRow = (queryPos - 1) / 64;
	arrCol = (queryPos - 1) % 64;

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

	return SUCCESSFUL;
}

/**
 * Generate the kmer integer sequence from a read specified by a startPos (>=0).
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateKmerSeqIntFromReadset(uint64_t *seqInt, uint64_t *readseq, int32_t startReadPos, int32_t kmerSize, int32_t entriesNum, int32_t baseNumLastEntry)
{
	int32_t i, j, startEntryPos, remainedBaseNum, rightRemainedNum, entryRow, baseNumInEntry;
	uint64_t *readseqStart;

	for(i=0; i<entriesNum; i++) seqInt[i] = 0;

	entryRow = startReadPos >> 5;
	readseqStart = readseq + entryRow;
	startEntryPos = startReadPos % 32;

	i = 0;
	j = 0;
	remainedBaseNum = kmerSize;
	while(remainedBaseNum>0)
	{
		// process first entry
		if(entryRow!=entriesNum-1)
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
			if(entryRow!=entriesNum-1)
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
short fillReadMatchInfoQueries(queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet)
{
	int32_t i, j, maxReadsNum;
	readBlock_t *readBlock;
	readMatchInfoBlock_t *readMatchInfoBlock;
	readMatchInfo_t *pReadMatchInfo;
	query_t *queryItem;
	queryRead_t *queryReadArrayBuf;
	int32_t queryReadNumTmp;


	// allocate memory
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		queryItem = queryMatchInfoSet->queryArray + i;
		if(queryItem->queryReadNum>0)
		{
			queryItem->queryReadArray = (queryRead_t *) calloc (queryItem->queryReadNum, sizeof(queryRead_t));
			if(queryItem->queryReadArray==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
		queryItem->queryReadNum = 0;
	}

	// fill the match information
	queryReadNumTmp = 0;
	for(i=0; i<readSet->blocksNumRead; i++)
	{
		readBlock = readSet->readBlockArr + i;
		readMatchInfoBlock = readSet->readMatchInfoBlockArr + i;
		for(j=0; j<readBlock->itemNum; j++)
		{
			pReadMatchInfo = readMatchInfoBlock->readMatchInfoArr + j;
			if(pReadMatchInfo->queryID>0)
			{
				queryItem = queryMatchInfoSet->queryArray + pReadMatchInfo->queryID - 1;
				queryItem->queryReadArray[queryItem->queryReadNum].readID = pReadMatchInfo->readID;
				queryItem->queryReadArray[queryItem->queryReadNum].queryPos = pReadMatchInfo->queryPos;
				queryItem->queryReadArray[queryItem->queryReadNum].seqlen = pReadMatchInfo->seqlen;
				queryItem->queryReadArray[queryItem->queryReadNum].orientation = pReadMatchInfo->readOrientation;
				queryItem->queryReadArray[queryItem->queryReadNum].startReadPos = pReadMatchInfo->startReadPos;
				queryItem->queryReadArray[queryItem->queryReadNum].alignSize = pReadMatchInfo->alignSize;
				queryItem->queryReadNum ++;
			}
		}
	}

	// get the maximal array size and allocate the buffer memory
	maxReadsNum = 0;
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		queryItem = queryMatchInfoSet->queryArray + i;
		if(maxReadsNum<queryItem->queryReadNum)
			maxReadsNum = queryItem->queryReadNum;
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
		if(queryItem->queryReadNum>0)
		{
			if(radixSortQueryReadArray(queryItem->queryReadArray, queryReadArrayBuf, queryItem->queryReadNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot sort the query read array, error!\n", __LINE__, __func__);
				return FAILED;
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
	int32_t i, j, queryID, queryLen, queryPos, queryPos_paired, fragSize, fragSizeDif, itemNum, matedFlag;
	int64_t readID, readID_paired;
	char orientation;
	query_t *queryItem;
	queryRead_t *queryReadArray;
	FILE *fpOut;

	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo;
	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;


	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;


	fpOut = fopen(outfile, "w");
	if(fpOut==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, outfile);
		return FAILED;
	}

	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		queryItem = queryMatchInfoSet->queryArray + i;
		queryID = queryItem->queryID;
		queryLen = queryItem->queryLen;
		if(queryItem->queryReadNum>0)
		{
			queryReadArray = queryItem->queryReadArray;
			itemNum = queryItem->queryReadNum;
			fprintf(fpOut, ">%s\t%d\t0\t%d\n", queryItem->queryTitle, queryLen, itemNum);
			for(j=0; j<itemNum; j++)
			{
				readID = queryReadArray[j].readID;
				queryPos = queryReadArray[j].queryPos;

				if(queryReadArray[j].orientation==ORIENTATION_PLUS)
					orientation = '+';
				else
					orientation = '-';

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

					fragSizeDif = fragSize - 370;
					if(fragSizeDif>-300 && fragSizeDif<300)
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
					fprintf(fpOut, "%d\t%ld\t%c\t%d\t%d\n", queryReadArray[j].queryPos, queryReadArray[j].readID, orientation, queryReadArray[j].seqlen, fragSize);
				else
					fprintf(fpOut, "%d\t%ld\t%c\t%d\n", queryReadArray[j].queryPos, queryReadArray[j].readID, orientation, queryReadArray[j].seqlen);
			}
		}
	}

	fclose(fpOut);

	return SUCCESSFUL;
}

