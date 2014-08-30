/*
 * map.c
 *
 *  Created on: May 30, 2014
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Start the reads alignment.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mapReads(queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet, queryIndex_t *queryIndex)
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

	printf("Aligning the reads ...\n");

	// initialize read blocks
	if(initReadMatchInfoBlockInReadset(readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize read blocks, error!\n", __LINE__, __func__);
		return FAILED;
	}

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

	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

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

//			if(rid==1181367 || rid==1181368)
//			{
//				printf("rid=%ld\n", rid);
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
				// select the mated match information
				if(matchItemNums[0]>0 && matchItemNums[1]>0 && (matchItemNums[0]>=2 || matchItemNums[1]>=2))
				{
					if(selectMatedMatchPosPE(matchResultArrays[0], matchResultArrays[1], matchItemNums, matchItemNums+1)==FAILED)
					{
						printf("line=%d, In %s(), cannot adjust the match information of paired-ends, error!\n", __LINE__, __func__);
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
						pReadMatchInfo->matchlen = pRead->seqlen;
						pReadMatchInfo->readID = tmpRid;
						pReadMatchInfo->seqlen = pRead->seqlen;
						pReadMatchInfo->readOrientation = matchResultArray[0].orientation;

						queryMatchInfoSet->queryArray[pReadMatchInfo->queryID-1].queryReadNum ++;

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
 * Perfectly map single read to queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mapSingleReadToQueries(int64_t rid, read_t *pRead, readSet_t *readSet, alignMatchItem_t *matchResultArray, alignMatchItem_t *matchResultArrayBuf1, alignMatchItem_t *matchResultArrayBuf2, int32_t *matchItemNum, queryIndex_t *queryIndex, query_t *queryArray, int32_t mismatchNumThreshold)
{
	// perfectly map single read
	if(mapSingleReadToQueriesPerfect(rid, pRead, readSet, matchResultArray, matchResultArrayBuf1, matchResultArrayBuf2, matchItemNum, queryIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot map the read %lu, error!\n", __LINE__, __func__, rid);
		return FAILED;
	}

	if((*matchItemNum)<=0)
	{
		// perfectly map single read
		if(mapSingleReadToQueriesWithMismatch(rid, pRead, readSet, matchResultArray, matchResultArrayBuf1, matchItemNum, queryIndex, queryArray, mismatchNumThreshold)==FAILED)
		{
			printf("line=%d, In %s(), cannot map the read %lu, error!\n", __LINE__, __func__, rid);
			return FAILED;
		}
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
		if(tmpItemNum<=1)
		{
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
		if(tmpItemNum<=2)
		{
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
		}

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
	int32_t i, j, entriesNumRead, baseNumLastEntryRead, startBasePos;
	uint64_t hashcode, kmerSeqInt[queryIndex->entriesPerKmer];
	queryKmer_t *kmer;
	char *querySeq, readSeq[MAX_READ_LEN_IN_BUF+1];
	int32_t queryID, queryLen, queryPos, maxStartBasePos, remainBaseNum, misNum, startCmpRpos, startCmpQpos;

	entriesNumRead = ((seqLen - 1) >> 5) + 1;
	baseNumLastEntryRead = ((seqLen - 1) % 32) + 1;

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

					if((queryPos+seqLen-startBasePos-1<=queryLen) && (queryPos-1>=startBasePos))
					{
						misNum = 0;
						startCmpRpos = startBasePos + queryIndex->kmerSize;
						startCmpQpos = queryPos + queryIndex->kmerSize - 1;
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

						startCmpQpos = queryPos - startBasePos - 1;
						remainBaseNum = startBasePos;
						for(j=0; j<remainBaseNum; j++)
						{
							if(readSeq[j]!=querySeq[startCmpQpos+j])
							{
								misNum ++;
								if(misNum>mismatchNumThreshold)
									break;
							}
						}

						if(misNum<=mismatchNumThreshold)
						{ // add the match information
							matchResultArray[*matchItemNum].queryID = queryID;
							matchResultArray[*matchItemNum].queryPos = queryPos - startBasePos;
							matchResultArray[*matchItemNum].mismatchNum = misNum;
							matchResultArray[*matchItemNum].validFlag = YES;
							(*matchItemNum) ++;
						}
					}
				}
			}

			if((*matchItemNum)==0)
			{
				// set next start basePos
				startBasePos += queryIndex->kmerSize;
			}else
			{
				break;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Select the mated query position of paired reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short selectMatedMatchPosPE(alignMatchItem_t *matchResultArray1, alignMatchItem_t *matchResultArray2, int32_t *matchItemNum1, int32_t *matchItemNum2)
{
	int32_t i, j, queryID1, queryPos1, queryID2, queryPos2, orient1, orient2, matedNum, matedRow1, matedRow2, pairRow;

	for(i=0; i<(*matchItemNum1); i++) {matchResultArray1[i].validFlag = NO; matchResultArray1[i].pairRow = -1;}
	for(i=0; i<(*matchItemNum2); i++) {matchResultArray2[i].validFlag = NO; matchResultArray2[i].pairRow = -1;}

	matedNum = 0;
	matedRow1 = matedRow2 = -1;
	for(i=0; i<(*matchItemNum1); i++)
	{
		queryID1 = matchResultArray1[i].queryID;
		queryPos1 = matchResultArray1[i].queryPos;
		orient1 = matchResultArray1[i].orientation;

		if(matchResultArray1[i].validFlag==NO)
		{
			for(j=0; j<(*matchItemNum2); j++)
			{
				queryID2 = matchResultArray2[j].queryID;
				queryPos2 = matchResultArray2[j].queryPos;
				orient2 = matchResultArray2[j].orientation;

				if(matchResultArray2[j].validFlag==NO && queryID1==queryID2 && orient1!=orient2)
				{
					matchResultArray1[i].validFlag = YES;
					matchResultArray1[i].pairRow = j;
					matedRow1 = i;

					matchResultArray2[j].validFlag = YES;
					matchResultArray2[j].pairRow = i;
					matedRow2 = j;

					matedNum ++;

					break;
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
	{ // multiple mated, then select the first one
		if(matchResultArray1[0].validFlag==YES)
		{
			pairRow = matchResultArray1[i].pairRow;
			if(pairRow>1)
			{
				if(memcpy(matchResultArray2, matchResultArray2+pairRow, sizeof(alignMatchItem_t))==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				matchResultArray1[0].pairRow = 0;
			}
		}else
		{
			for(i=1; i<(*matchItemNum1); i++)
			{
				if(matchResultArray1[i].validFlag==YES)
				{
					if(memcpy(matchResultArray1, matchResultArray1+i, sizeof(alignMatchItem_t))==NULL)
					{
						printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
						return FAILED;
					}
					matchResultArray2[i].pairRow = 0;

					pairRow = matchResultArray1[0].pairRow;
					if(pairRow>=1)
					{
						if(memcpy(matchResultArray2, matchResultArray2+pairRow, sizeof(alignMatchItem_t))==NULL)
						{
							printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
						matchResultArray1[i].pairRow = 0;
					}

					break;
				}
			}
		}

		*matchItemNum1 = 1;
		*matchItemNum2 = 1;
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
			if(memcpy(matchResultArray+i, matchResultArray+minRow, sizeof(alignMatchItem_t))==NULL)
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
			//hashcode = bitMask - ((data[i].mismatchNum >> step) & bitMask);  // from big to small
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

