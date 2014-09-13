/*
 * misass.c
 *
 *  Created on: May 28, 2014
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Validate mis-assembled queries using paired-end reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short validateMisassQueries(char *queryMatchInfoFile, char *inputQueryFile, char *mergedSegFile, char **readsFileNames, int32_t readsFileNum, int32_t readsFileFormatType)
{
	// load the query match information from file
	if(loadQueryMatchInfoFromFile(&queryMatchInfoSet, queryMatchInfoFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the query match information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the queries
	if(fillQueries(queryMatchInfoSet, inputQueryFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the queries, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill subject sequences
	if(fillSubjectSeqs(queryMatchInfoSet, mergedSegFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill subject sequences, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// compute the number of potential mis-assembled queries
	if(computePotentMisassNum(queryMatchInfoSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the number of potential mis-assembled queries, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(queryMatchInfoSet->potentMisassNum>0)
	{

//		// #############################################
//		if(outputGapRegInQueries(queryMatchInfoSet)==FAILED)
//		{
//			printf("line=%d, In %s(), cannot get the gap regions, error!\n", __LINE__, __func__);
//			return FAILED;
//		}
//		// #############################################

		// generate the Query Index using the unmatched queries
		if(buildQueryIndex(&queryIndex, queryMatchInfoSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the query index, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// load the reads data to readSet
		if(constructReadset(&readSet, readsFileNames, readsFileNum, readsFileFormatType, reserveHashItemBlocksFlag)==FAILED)
		{
			printf("line=%d, In %s(), cannot construct read set, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(estimateInsertSize(&insertSize, &standDev, readSet, queryIndex)==FAILED)
		{
			printf("line=%d, In %s(), cannot estimate the insert size, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// align the reads to Query Index to generate reads alignment information
		if(mapReads(queryMatchInfoSet, readSet, queryIndex, insertSize, standDev)==FAILED)
		{
			printf("line=%d, In %s(), cannot align reads to queries, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// fill the reads match information to queries
		if(fillReadMatchInfoQueries(queryMatchInfoSet, readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot fill read match information to read set, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// ###################### Debug information ##################
//		if(outputQueryReadArray("../output_tmp/alignedReads.txt", queryMatchInfoSet, readSet)==FAILED)
//		{
//			printf("line=%d, In %s(), cannot output aligned reads, error!\n", __LINE__, __func__);
//			return FAILED;
//		}
		// ###################### Debug information ##################

		// compute the potential mis-assembled queries
		if(computeMisassQueries(queryMatchInfoSet, readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute potential mis-assembled queries, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// determine the structural variations in queries
		if(computeSVInQueries(queryMatchInfoSet, readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot determine structural variations in queries, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// correct the mis-assemblies
		if(correctMisassQueries(queryMatchInfoSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot correct the mis-assembled queries, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// extract mis-assembly regions
		if(extractMisReg(queryMatchInfoSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot extract mis-assembly regions, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// generate the data for Circos
//		if(generateCircosData(queryMatchInfoSet, readSet)==FAILED)
//		{
//			printf("line=%d, In %s(), cannot determine the query mis-assembly, error!\n", __LINE__, __func__);
//			return FAILED;
//		}

		// save mis new queries to file
		if(saveMisassQueries(errorsFile, svFile, misUncertainFile, gapFile, queryMatchInfoSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot save mis-assembled queries, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// save new queries to file
		if(saveNewQueries(newQueryFile, queryMatchInfoSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot save mis-assembled queries, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		printf("There are no potential mis-assembled queries.\n");
	}

	// free query match information, query index, readSet
	freeMisassMem();

	return SUCCESSFUL;
}

/**
 * Free the memory of mis-assembled queries validation.
 */
void freeMisassMem()
{
	// free the memory of query match information
	releaseQueryMatchInfo(&queryMatchInfoSet);

	// free queryIndex
	releaseQueryIndex(&queryIndex);

	// free readSet
	releaseReadset(&readSet);
}

/**
 * Compute the number of potential mis-assembled queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computePotentMisassNum(queryMatchInfo_t *queryMatchInfoSet)
{
	int32_t i, maxQueryID, secQueryID, maxQueryLen, secQueryLen, perfectNum;
	query_t *queryArray;

	queryArray = queryMatchInfoSet->queryArray;

	perfectNum = 0;
	queryMatchInfoSet->potentMisassNum = 0;
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{

		if(queryArray[i].globalMatchKind==PERFECT_MATCH_KIND)
			perfectNum ++;
		else
		{
			queryArray[i].misassFlag = POTENTIAL_MISASS;
			queryMatchInfoSet->potentMisassNum ++;
		}
	}

	if(queryMatchInfoSet->maxQueryID<=0)
	{
		maxQueryID = secQueryID = -1;
		maxQueryLen = secQueryLen = 0;
		for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
		{
			if(maxQueryLen<queryArray[i].queryLen)
			{
				secQueryLen = maxQueryLen;
				secQueryID = maxQueryID;
				maxQueryLen = queryArray[i].queryLen;
				maxQueryID = queryArray[i].queryID;
			}else if(secQueryLen<queryArray[i].queryLen)
			{
				secQueryLen = queryArray[i].queryLen;
				secQueryID = queryArray[i].queryID;
			}
		}
		queryMatchInfoSet->maxQueryID = maxQueryID;
		queryMatchInfoSet->secQueryID = secQueryID;
	}

	printf("Number of queries perfectly aligned to reference: %d\n", perfectNum);
	printf("Number of queries need to be further validated  : %d\n", queryMatchInfoSet->potentMisassNum);

	return SUCCESSFUL;
}

/**
 * Compute potential mis-assembled queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeMisassQueries(queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet)
{
	int32_t i, processedNum;
	query_t *queryArray;

	printf("Begin identifying assembly errors, please wait ...\n");

	queryArray = queryMatchInfoSet->queryArray;

	// compute insert size
//	if(computeInsertSize(&insertSize, &standDev, queryMatchInfoSet, readSet)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot compute insert size, error!\n", __LINE__, __func__);
//		return FAILED;
//	}

	// compute the S/P, S-/S, S+/S ratios in normal regions
	if(computeNormalRatios(&SP_ratio_Thres, &SMinus_ratio_Thres, &SPlus_ratio_Thres, &discorRatio_Thres, queryMatchInfoSet, readSet, insertSize, standDev)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the ratios of normal regions, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// validate the queries
	processedNum = 0;
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		if(queryArray[i].misassFlag==POTENTIAL_MISASS)
		{
			// ########################### Debug information ##############################
			//if(queryArray[i].queryID==4 || strcmp(queryArray[i].queryTitle, "ctg7180000002351")==0)
			//{
			//	printf("======= queryID=%d, queryTitle=%s, queryLen=%d, subjectNum=%d\n", queryMatchInfoSet->queryArray[i].queryID, queryMatchInfoSet->queryArray[i].queryTitle, queryMatchInfoSet->queryArray[i].queryLen, queryMatchInfoSet->queryArray[i].querySubjectNum);
			//}
			// ########################### Debug information ##############################

			if(computeSingleMisassQuery(queryArray+i, queryMatchInfoSet->subjectArray, readSet, SP_ratio_Thres, SMinus_ratio_Thres, SPlus_ratio_Thres, insertSize, standDev)==FAILED)
			{
				printf("line=%d, In %s(), cannot validate single potential mis-assembled query, error!\n", __LINE__, __func__);
				return FAILED;
			}

			processedNum ++;
			if(processedNum%100==0)
				printf("Queries processed: %d\n", processedNum);
		}
	}

	if(processedNum%100!=0)
		printf("Queries processed: %d\n", processedNum);

	return SUCCESSFUL;
}

/**
 * Compute insert size.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeInsertSize(double *insertSize, double *standDev, queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet)
{
	int64_t i, j, k, fragSizeSum, pairNum, pairNumTmp, fragSize, maxQueryID, secQueryID;
	int32_t queryPos1, queryPos2, orient1, orient2, seqLen1, seqLen2;
	double meanSize, sdev, fragDif, fragDifSum;

	readBlock_t *readBlock;
	readMatchInfoBlock_t *readMatchInfoBlock;
	readMatchInfo_t *pReadMatchInfo1, *pReadMatchInfo2;

	fragSizeSum = 0;
	pairNum = 0;

	if(queryMatchInfoSet->maxQueryID>0)
	{
		maxQueryID = queryMatchInfoSet->maxQueryID;
		secQueryID = queryMatchInfoSet->secQueryID;

		for(k=0; k<2; k++)
		{
			pairNumTmp = 0;
			fragDifSum = 0;
			for(i=0; i<readSet->blocksNumRead; i++)
			{
				readBlock = readSet->readBlockArr + i;
				readMatchInfoBlock = readSet->readMatchInfoBlockArr + i;
				for(j=0; j<readBlock->itemNum; j+=2)
				{
					pReadMatchInfo1 = readMatchInfoBlock->readMatchInfoArr + j;
					pReadMatchInfo2 = readMatchInfoBlock->readMatchInfoArr + j + 1;

					if(pReadMatchInfo1->queryID>0 && pReadMatchInfo1->queryID==pReadMatchInfo2->queryID && (pReadMatchInfo1->queryID==maxQueryID || pReadMatchInfo1->queryID==secQueryID))
					{
						if(pReadMatchInfo1->readID+1!=pReadMatchInfo2->readID)
						{
							printf("line=%d, In %s(), readID1=%ld, readID2=%ld, error!\n", __LINE__, __func__, (int64_t)pReadMatchInfo1->readID, (int64_t)pReadMatchInfo2->readID);
							return FAILED;
						}

						queryPos1 = pReadMatchInfo1->queryPos;
						seqLen1 = pReadMatchInfo1->seqlen;
						orient1 = pReadMatchInfo1->readOrientation;
						queryPos2 = pReadMatchInfo2->queryPos;
						seqLen2 = pReadMatchInfo2->seqlen;
						orient2 = pReadMatchInfo2->readOrientation;

						fragSize = 0;
						if(orient1==ORIENTATION_PLUS && orient2==ORIENTATION_MINUS)
							fragSize = queryPos2 + seqLen2 - queryPos1;
						else if(orient2==ORIENTATION_PLUS && orient1==ORIENTATION_MINUS)
							fragSize = queryPos1 + seqLen1 - queryPos2;

						if(fragSize>0)
						{
							if(k==0)
							{
								fragSizeSum += fragSize;
							}else
							{
								fragDif = (fragSize - meanSize) * (fragSize - meanSize);
								fragDifSum += fragDif;
							}
							pairNumTmp ++;
						}
					}
				}
			}

			if(k==0)
			{
				if(pairNumTmp>0)
				{
					pairNum = pairNumTmp;
					meanSize = (double)fragSizeSum / pairNum;
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

				sdev = sqrt((fragDifSum) / pairNumTmp);
			}
		}
	}

	if(pairNum>0)
	{
		*insertSize = meanSize;
		*standDev = sdev;

		printf("insertSize=%.4f, SDev=%.4f\n", *insertSize, *standDev);
	}else
	{
		printf("line=%d, In %s(), no valid paired reads for computing insert size of library, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Compute single potential mis-assembled query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeSingleMisassQuery(query_t *queryItem, subject_t *subjectArray, readSet_t *readSet, double SP_ratio_Thres, double SMinus_ratio_Thres, double SPlus_ratio_Thres, double insertSize, double standDev)
{
	// get the misInfo list
	if(getMisInfoList(queryItem, subjectArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the mis-assembly information list, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// determine mis-assembly information for each node in the list
	if(determineMisInfoSingleQuery(queryItem, subjectArray, readSet, SP_ratio_Thres, SMinus_ratio_Thres, SPlus_ratio_Thres, insertSize, standDev)==FAILED)
	{
		printf("line=%d, In %s(), cannot determine the mis-assembly information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ######################## Debug information ###################
	//outputMisInfoList(queryItem);
	// ######################## Debug information ###################

	return SUCCESSFUL;
}

/**
 * Compute base coverage of single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeBaseCovSingleQuery(baseCov_t *baseCovArray, query_t *queryItem, readSet_t *readSet)
{
	int32_t i, j, k, startQueryRow, endQueryRow, startReadRow;
	queryRead_t *queryRead;

	int32_t readBlockID, rowNumInReadBlock, maxItemNumPerReadBlock;
	readBlock_t *readBlockArray;
	read_t *pRead;
	readseqBlock_t *readseqBlockArray;
	uint64_t *pReadseq;
	char readseq[MAX_READ_LEN_IN_BUF+1];

	if(memset(baseCovArray, 0LU, queryItem->queryLen*sizeof(baseCov_t))==NULL)
	{
		printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	readBlockArray = readSet->readBlockArr;
	maxItemNumPerReadBlock = readSet->maxItemNumPerReadBlock;
	readseqBlockArray = readSet->readseqBlockArr;

	for(i=0; i<queryItem->queryReadNum; i++)
	{
		queryRead = queryItem->queryReadArray + i;

		readBlockID = (queryRead->readID - 1) / maxItemNumPerReadBlock;
		rowNumInReadBlock = (queryRead->readID - 1) % maxItemNumPerReadBlock;
		pRead = readBlockArray[readBlockID].readArr + rowNumInReadBlock;
		pReadseq = readseqBlockArray[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;

		if(queryRead->orientation==ORIENTATION_PLUS)
		{
			getReadBaseFromPosByInt(readseq, pReadseq, pRead->seqlen, 0, pRead->seqlen);
		}else
		{
			getReverseReadBaseFromPosByInt(readseq, pReadseq, pRead->seqlen, 0, pRead->seqlen);
		}

		startReadRow = queryRead->startReadPos - 1;
		startQueryRow = queryRead->queryPos - 1;
		endQueryRow = startQueryRow + queryRead->alignSize - 1;
		if(startQueryRow<0)
			startQueryRow = 0;
		if(endQueryRow>queryItem->queryLen-1)
		{
			endQueryRow = queryItem->queryLen - 1;

			printf("endQueryRow=%d, queryLen=%d\n", endQueryRow, queryItem->queryLen-1);
			return FAILED;
		}

		for(j=startQueryRow, k=startReadRow; j<=endQueryRow; j++, k++)
		{
			switch(readseq[k])
			{
				case 'A': baseCovArray[j].baseNumArray[0] ++; break;
				case 'C': baseCovArray[j].baseNumArray[1] ++; break;
				case 'G': baseCovArray[j].baseNumArray[2] ++; break;
				case 'T': baseCovArray[j].baseNumArray[3] ++; break;
				case 'N': baseCovArray[j].baseNumArray[4] ++; break;
				default: printf("line=%d, In %s(), unknown base %c, error!\n", __LINE__, __func__, readseq[k]); return FAILED;
			}
			baseCovArray[j].baseNumArray[5] ++;
		}
	}

	return SUCCESSFUL;
}

/**
 * Output base coverage of single query to file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputBaseCovSingleQueryToFile(char *covFileName, baseCov_t *baseCovArray, int32_t arraySize)
{
	FILE *fpCov;
	int32_t i, j, maxValue;

	fpCov = fopen(covFileName, "w");
	if(fpCov==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, covFileName);
		return FAILED;
	}

	fprintf(fpCov, "basePos\tcov\tcov_A\tcov_C\tcov_G\tcov_T\tcov_N\n");
	for(i=0; i<arraySize; i++)
	{
		if(baseCovArray[i].baseNumArray[5]>0)
		{
			maxValue = 0;
			for(j=0; j<5; j++)
			{
				if(maxValue<baseCovArray[i].baseNumArray[j])
					maxValue = baseCovArray[i].baseNumArray[j];
			}
			fprintf(fpCov, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\n", i+1, baseCovArray[i].baseNumArray[5], baseCovArray[i].baseNumArray[0], baseCovArray[i].baseNumArray[1], baseCovArray[i].baseNumArray[2], baseCovArray[i].baseNumArray[3], baseCovArray[i].baseNumArray[4], (double)maxValue/baseCovArray[i].baseNumArray[5]);
		}else
			fprintf(fpCov, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\n", i+1, baseCovArray[i].baseNumArray[5], baseCovArray[i].baseNumArray[0], baseCovArray[i].baseNumArray[1], baseCovArray[i].baseNumArray[2], baseCovArray[i].baseNumArray[3], baseCovArray[i].baseNumArray[4], 0.0f);
	}
	fclose(fpCov);

	return SUCCESSFUL;
}

/**
 * Output the reads covering a base.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputCovReadsToFile(char *covReadsFileName, int32_t basePos, char base, query_t *queryItem)
{
	FILE *fpCovReads;
	int32_t i, j, startQueryRow, endQueryRow, readRow;
	queryRead_t *queryRead;
	char orient;

	int32_t readBlockID, rowNumInReadBlock, maxItemNumPerReadBlock;
	readBlock_t *readBlockArray;
	read_t *pRead;
	readseqBlock_t *readseqBlockArray;
	uint64_t *pReadseq;
	char readseq[MAX_READ_LEN_IN_BUF+1];

	fpCovReads = fopen(covReadsFileName, "w");
	if(fpCovReads==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, covReadsFileName);
		return FAILED;
	}

	readBlockArray = readSet->readBlockArr;
	maxItemNumPerReadBlock = readSet->maxItemNumPerReadBlock;
	readseqBlockArray = readSet->readseqBlockArr;

	if(basePos<queryItem->queryLen)
	{
		for(i=0; i<queryItem->queryReadNum; i++)
		{
			queryRead = queryItem->queryReadArray + i;

			startQueryRow = queryRead->queryPos - 1;
			endQueryRow = startQueryRow + queryRead->seqlen - 1;
			if(startQueryRow<0)
				startQueryRow = 0;
			if(endQueryRow>queryItem->queryLen-1)
				endQueryRow = queryItem->queryLen - 1;

			if(startQueryRow+1<=basePos && endQueryRow+1>=basePos)
			{
				readBlockID = (queryRead->readID - 1) / maxItemNumPerReadBlock;
				rowNumInReadBlock = (queryRead->readID - 1) % maxItemNumPerReadBlock;
				pRead = readBlockArray[readBlockID].readArr + rowNumInReadBlock;
				pReadseq = readseqBlockArray[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;

				if(queryRead->orientation==ORIENTATION_PLUS)
					getReadBaseFromPosByInt(readseq, pReadseq, pRead->seqlen, 0, pRead->seqlen);
				else
					getReverseReadBaseFromPosByInt(readseq, pReadseq, pRead->seqlen, 0, pRead->seqlen);

				readRow = basePos - startQueryRow - 1;
				if(readseq[readRow]==base)
				{
					if(queryRead->orientation==ORIENTATION_PLUS)
						orient = '+';
					else
						orient = '-';
					fprintf(fpCovReads, "%ld\t%d\t%c\t%d\t%s\n", queryRead->readID, queryRead->queryPos, orient, readRow, readseq);
				}
			}
		}
	}else
	{
		printf("line=%d, In %s(), basePos=%d, queryLen=%d, error!\n", __LINE__, __func__, basePos, queryItem->queryLen);
		return FAILED;
	}

	fclose(fpCovReads);

	return SUCCESSFUL;
}

/**
 * Compute the disagreements of query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeDisagreements(int32_t *disagreeNum, int32_t *zeroCovNum, baseCov_t *baseCovArray, int32_t startRow, int32_t endRow, int32_t printFlag)
{
	int32_t i, j, maxID, maxValue, secID, secValue, total, disagreeFlag;
	double maxRatio;

	if(printFlag==YES)
	{
		printf("Disagreements for region [%d, %d]:\n", startRow+1, endRow+1);
		printf("basePos\tcov\tcov_A\tcov_C\tcov_G\tcov_T\tcov_N\tmaxRatio\n");
	}

	*disagreeNum = 0;
	*zeroCovNum = 0;
	for(i=startRow; i<=endRow; i++)
	{
		disagreeFlag = NO;
		if(baseCovArray[i].baseNumArray[5]>0)
		{
			maxID = secID = -1;
			maxValue = secValue = 0;
			total = baseCovArray[i].baseNumArray[5];
			for(j=0; j<5; j++)
			{
				if(baseCovArray[i].baseNumArray[j]>maxValue)
				{
					secValue = maxValue;
					secID = maxID;
					maxValue = baseCovArray[i].baseNumArray[j];
					maxID = j;
				}else if(baseCovArray[i].baseNumArray[j]>secValue)
				{
					secValue = baseCovArray[i].baseNumArray[j];
					secID = j;
				}
			}

			maxRatio = (double)maxValue / total;
			if(maxRatio<DISAGREE_RATIO_THRES)
			{
				disagreeFlag = YES;
				if(secValue==1 && maxValue>=4)
					disagreeFlag = NO;
			}
		}else
		{
			maxRatio = 0;
			disagreeFlag = YES;
		}

		if(disagreeFlag==YES)
		{
			if(printFlag==YES)
				printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\n", i+1, baseCovArray[i].baseNumArray[5], baseCovArray[i].baseNumArray[0], baseCovArray[i].baseNumArray[1], baseCovArray[i].baseNumArray[2], baseCovArray[i].baseNumArray[3], baseCovArray[i].baseNumArray[4], maxRatio);
			(*disagreeNum) ++;

			if(baseCovArray[i].baseNumArray[5]==0)
				(*zeroCovNum) ++;
		}
	}

	if(printFlag==YES)
	{
		printf("Disagreements: %d\n", *disagreeNum);
		printf("Zero coverage: %d\n", *zeroCovNum);
		printf("Region size  : %d\n", endRow-startRow+1);
	}

	return SUCCESSFUL;
}

/**
 * Compute the covRatio of a region.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeAbnormalCovRegNum(int32_t *highCovRegNum, int32_t *lowCovRegNum, baseCov_t *baseCovArray, int32_t startRow, int32_t endRow, int32_t arraySize, int32_t skipEndFlag)
{
	int64_t i, j, baseNumReg, baseNumTotal;
	double covReg, covTotal;
	int32_t subRegSize, newStartRow, newEndRow, tmp;

	if(startRow<500 || endRow>arraySize-500)
	{
		*highCovRegNum = *lowCovRegNum = 0;
		return SUCCESSFUL;
	}

	subRegSize = 50;

	if(skipEndFlag==1)
	{
		newStartRow = startRow + subRegSize;
		if(endRow>arraySize-subRegSize)
			newEndRow = endRow - subRegSize;
		else
			newEndRow = endRow;
	}else if(skipEndFlag==2)
	{
		if(startRow<subRegSize)
			newStartRow = startRow + subRegSize;
		else
			newStartRow = startRow;
		newEndRow = endRow - subRegSize;
	}else
	{
		newStartRow = startRow;
		newEndRow = endRow;
	}

	if(newStartRow>newEndRow)
	{
		tmp = newStartRow;
		newStartRow = newEndRow;
		newEndRow = tmp;
	}


	baseNumTotal = 0;
	for(i=0; i<arraySize; i++)
		baseNumTotal += baseCovArray[i].baseNumArray[5];
	covTotal = (double)baseNumTotal / arraySize;

	*highCovRegNum = *lowCovRegNum = 0;
	baseNumReg = 0;
	j = 0;
	for(i=newStartRow; i<=newEndRow; i++)
	{
		baseNumReg += baseCovArray[i].baseNumArray[5];

		j ++;
		if(j==subRegSize || i==endRow)
		{
			if(j>0)
			{
				covReg = (double)baseNumReg / j;
				//printf("[%d, %d], covReg=%.4f\n", i+2-j, i+1, covReg);
				if(covReg>1.5*covTotal)
					(*highCovRegNum) ++;
				else if(covReg<0.5*covTotal)
					(*lowCovRegNum) ++;
			}

			j = 0;
			baseNumReg = 0;
		}
	}

	//printf("highCovRegNum=%d, lowCovRegNum=%d, covTotal=%.4f\n", *highCovRegNum, *lowCovRegNum, covTotal);

	return SUCCESSFUL;
}

/**
 * Compute the align margin of query segments.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeMisassFlagAlignSeg(int32_t *misassFlag, globalValidSeg_t *leftAlignSeg, globalValidSeg_t *rightAlignSeg, int32_t subjectLen, int32_t circularFlag)
{
	int32_t distance;

	*misassFlag = NO;

	if(leftAlignSeg->subjectID!=rightAlignSeg->subjectID || leftAlignSeg->strand!=rightAlignSeg->strand)
	{
		*misassFlag = YES;
	}else
	{
		if(circularFlag==YES)
		{
			if(leftAlignSeg->strand==PLUS_STRAND)
			{
				if((leftAlignSeg->endSubPos>=subjectLen-varyEndLenThres && leftAlignSeg->endSubPos<=subjectLen)
					&& (rightAlignSeg->startSubPos>=1 && rightAlignSeg->startSubPos<=varyEndLenThres))
				{
					distance = (subjectLen - leftAlignSeg->endSubPos) + (rightAlignSeg->startSubPos - 1);
					if(distance<0)
						distance = -distance;

					if(distance > minDisjunctDistanceThres)
						*misassFlag = YES;
				}else
				{
					distance  = rightAlignSeg->startSubPos - leftAlignSeg->endSubPos;
					if(distance<0)
						distance = -distance;

					if(distance > minDisjunctDistanceThres)
						*misassFlag = YES;
				}
			}else
			{
				if((leftAlignSeg->endSubPos>=1 && leftAlignSeg->endSubPos<=varyEndLenThres)
					&& (rightAlignSeg->startSubPos>=subjectLen-varyEndLenThres && rightAlignSeg->startSubPos<=subjectLen))
				{
					distance = (subjectLen - rightAlignSeg->startSubPos) + (leftAlignSeg->endSubPos - 1);
					if(distance<0)
						distance = -distance;

					if(distance > minDisjunctDistanceThres)
						*misassFlag = YES;
				}else
				{
					distance  = rightAlignSeg->startSubPos - leftAlignSeg->endSubPos;
					if(distance<0)
						distance = -distance;

					if(distance > minDisjunctDistanceThres)
						*misassFlag = YES;
				}
			}
		}else
		{
			distance = leftAlignSeg->endSubPos - rightAlignSeg->startSubPos;
			if(distance<0)
				distance = -distance;
			if(distance>minDisjunctDistanceThres)
				*misassFlag = YES;
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the gap flag of misInfo node.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getGapFlagMisInfo(int32_t *gapFlag, globalValidSeg_t *leftAlignSeg, globalValidSeg_t *rightAlignSeg, query_t *queryItem)
{
	int32_t i, leftPos, rightPos;
	char *querySeq, base;

	querySeq = queryItem->querySeq;

	if(leftAlignSeg->endQueryPos<=rightAlignSeg->startQueryPos)
	{
		leftPos = leftAlignSeg->endQueryPos;
		rightPos = rightAlignSeg->startQueryPos;
	}else
	{
		leftPos = rightAlignSeg->startQueryPos;
		rightPos = leftAlignSeg->endQueryPos;
	}

	*gapFlag = NO;
	for(i=leftPos-1; i<rightPos; i++)
	{
		base = querySeq[i];
		if(base=='N' || base=='n' || base=='.')
		{
			*gapFlag = YES;
			break;
		}
	}

	return SUCCESSFUL;
}

/**
 * Compute the align margin of query segments.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeAlignSegMargin(int64_t *leftMargin, int64_t *rightMargin, globalValidSeg_t *leftAlignSeg, globalValidSeg_t *rightAlignSeg, int32_t difQuery, int32_t difSubject, query_t *queryItem, subject_t *subjectArray)
{
	int64_t i, leftQueryPos, rightQueryPos, leftSubjectPos, rightSubjectPos, leftMarginSubject, rightMarginSubject;
	char *queryAlignSeq, *subjectAlignSeq;
	int32_t queryAlignSeqLen, subjectAlignSeqLen, maxSeqLen;
	char *alignResultArray[3];
	int32_t overlapLen, mismatchNum, queryLeftShiftLen, queryRightShiftLen, subjectLeftShiftLen, subjectRightShiftLen;

	if(difQuery<-100)
	{
		// get the two base sequences for the first alignment
		if(getRightAlignSeqs(&queryAlignSeq, &subjectAlignSeq, &queryAlignSeqLen, &subjectAlignSeqLen, &leftQueryPos, &rightQueryPos, &leftSubjectPos, &rightSubjectPos, leftAlignSeg, rightAlignSeg, queryItem, subjectArray)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the right alignment sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// align the two sequences
		if(queryAlignSeqLen<subjectAlignSeqLen)
			maxSeqLen = 2 * subjectAlignSeqLen;
		else
			maxSeqLen = 2 * queryAlignSeqLen;
		for(i=0; i<3; i++)
		{
			alignResultArray[i] = (char*) calloc(maxSeqLen+1, sizeof(char));
			if(alignResultArray[i]==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		if(computeSeqAlignment(alignResultArray, &overlapLen, &mismatchNum, &queryLeftShiftLen, &subjectLeftShiftLen, &queryRightShiftLen, &subjectRightShiftLen, queryAlignSeq, subjectAlignSeq, queryAlignSeqLen, subjectAlignSeqLen, NO)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the alignment, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute the right margin
		if(computeRightMargin(rightMargin, &rightMarginSubject, alignResultArray, overlapLen, leftQueryPos, rightQueryPos, leftSubjectPos, rightSubjectPos, queryLeftShiftLen, subjectLeftShiftLen, leftAlignSeg)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the right margin, error!\n", __LINE__, __func__);
			return FAILED;
		}

		for(i=0; i<3; i++)
			free(alignResultArray[i]);
		free(queryAlignSeq);
		free(subjectAlignSeq);


		// get the two base sequences for the second alignment
		if(getLeftAlignSeqs(&queryAlignSeq, &subjectAlignSeq, &queryAlignSeqLen, &subjectAlignSeqLen, &leftQueryPos, &rightQueryPos, &leftSubjectPos, &rightSubjectPos, leftAlignSeg, rightAlignSeg, queryItem, subjectArray)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the right alignment sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// align the two sequences
		if(queryAlignSeqLen<subjectAlignSeqLen)
			maxSeqLen = 2 * subjectAlignSeqLen;
		else
			maxSeqLen = 2 * queryAlignSeqLen;
		for(i=0; i<3; i++)
		{
			alignResultArray[i] = (char*) calloc(maxSeqLen+1, sizeof(char));
			if(alignResultArray[i]==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		if(computeSeqAlignment(alignResultArray, &overlapLen, &mismatchNum, &queryLeftShiftLen, &subjectLeftShiftLen, &queryRightShiftLen, &subjectRightShiftLen, queryAlignSeq, subjectAlignSeq, queryAlignSeqLen, subjectAlignSeqLen, NO)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the alignment, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute the left margin
		if(computeLeftMargin(leftMargin, &leftMarginSubject, alignResultArray, overlapLen, leftQueryPos, rightQueryPos, leftSubjectPos, rightSubjectPos, queryRightShiftLen, subjectRightShiftLen, rightAlignSeg)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the right margin, error!\n", __LINE__, __func__);
			return FAILED;
		}

		for(i=0; i<3; i++)
			free(alignResultArray[i]);
		free(queryAlignSeq);
		free(subjectAlignSeq);
	}else
	{
		*rightMargin = leftAlignSeg->endQueryPos;
		*leftMargin = rightAlignSeg->startQueryPos;
	}

	return SUCCESSFUL;
}

/**
 * Get the right alignment base sequences including query sequence and the reference sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getRightAlignSeqs(char **queryAlignSeq, char **subjectAlignSeq, int32_t *queryAlignSeqLen, int32_t *subjectAlignSeqLen, int64_t *leftQueryPos, int64_t *rightQueryPos, int64_t *leftSubjectPos, int64_t *rightSubjectPos, globalValidSeg_t *leftAlignSeg, globalValidSeg_t *rightAlignSeg, query_t *queryItem, subject_t *subjectArray)
{
	char *subjectSeq;
	int32_t subjectID, subjectLen, baseNum, baseNum2;

	subjectID = leftAlignSeg->subjectID;
	subjectSeq = subjectArray[subjectID-1].subjectSeq;
	subjectLen = subjectArray[subjectID-1].subjectLen;

	// get the sequence for the first alignment
	if(leftAlignSeg->endQueryPos>rightAlignSeg->startQueryPos)
	{
		*leftQueryPos = rightAlignSeg->startQueryPos - 300;
		*rightQueryPos = leftAlignSeg->endQueryPos;
		if((*leftQueryPos)<1)
			*leftQueryPos = 1;
	}else
	{
		*leftQueryPos = leftAlignSeg->endQueryPos - 1000;
		*rightQueryPos = leftAlignSeg->endQueryPos;
		if((*leftQueryPos)<1)
			*leftQueryPos = 1;
	}

	if(queryItem->circularFlag==YES)
	{
		if(leftAlignSeg->strand==PLUS_STRAND)
		{
			*leftSubjectPos = leftAlignSeg->endSubPos - ((*rightQueryPos) - (*leftQueryPos));
			*rightSubjectPos = leftAlignSeg->endSubPos;
		}else
		{
			*leftSubjectPos = leftAlignSeg->endSubPos + ((*rightQueryPos) - (*leftQueryPos));
			*rightSubjectPos = leftAlignSeg->endSubPos;
		}
	}else
	{
		if(leftAlignSeg->strand==PLUS_STRAND)
		{
			*leftSubjectPos = leftAlignSeg->endSubPos - ((*rightQueryPos) - (*leftQueryPos));
			*rightSubjectPos = leftAlignSeg->endSubPos;
			if((*leftSubjectPos)<1)
				*leftSubjectPos = 1;
		}else
		{
			*leftSubjectPos = leftAlignSeg->endSubPos + ((*rightQueryPos) - (*leftQueryPos));
			*rightSubjectPos = leftAlignSeg->endSubPos;
			if((*leftSubjectPos)>subjectLen)
				*leftSubjectPos = subjectLen;
		}
	}

	// get the base sequences
	*queryAlignSeqLen = (*rightQueryPos) - (*leftQueryPos) + 1;
	if(leftAlignSeg->strand==PLUS_STRAND)
		*subjectAlignSeqLen = (*rightSubjectPos) - (*leftSubjectPos) + 1;
	else
		*subjectAlignSeqLen = (*leftSubjectPos) - (*rightSubjectPos) + 1;

	*queryAlignSeq = (char*)calloc((*queryAlignSeqLen)+1, sizeof(char));
	*subjectAlignSeq = (char*)calloc((*subjectAlignSeqLen)+1, sizeof(char));
	if((*queryAlignSeq)==NULL || (*subjectAlignSeq)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	strncpy(*queryAlignSeq, queryItem->querySeq+(*leftQueryPos)-1, *queryAlignSeqLen);
	(*queryAlignSeq)[*queryAlignSeqLen] = '\0';

	if(queryItem->circularFlag==YES)
	{
		if(leftAlignSeg->strand==PLUS_STRAND)
		{
			if((*leftSubjectPos)<1)
			{
				baseNum = 1 - (*leftSubjectPos);
				strncpy(*subjectAlignSeq, subjectSeq+subjectLen-baseNum, baseNum);
				(*subjectAlignSeq)[baseNum] = '\0';
				baseNum2 = (*subjectAlignSeqLen) - baseNum;
				strncpy((*subjectAlignSeq)+baseNum, subjectSeq, baseNum2);
				(*subjectAlignSeq)[*subjectAlignSeqLen] = '\0';
			}else
			{
				strncpy(*subjectAlignSeq, subjectSeq+(*leftSubjectPos)-1, *subjectAlignSeqLen);
				(*subjectAlignSeq)[*subjectAlignSeqLen] = '\0';
			}
		}else
		{
			if((*leftSubjectPos)>subjectLen)
			{
				baseNum = *subjectAlignSeqLen - ((*leftSubjectPos) - subjectLen);
				strncpy(*subjectAlignSeq, subjectSeq+(*rightSubjectPos)-1, baseNum);
				(*subjectAlignSeq)[baseNum] = '\0';
				baseNum2 = (*subjectAlignSeqLen) - baseNum;
				strncpy((*subjectAlignSeq)+baseNum, subjectSeq, baseNum2);
				(*subjectAlignSeq)[*subjectAlignSeqLen] = '\0';
			}else
			{
				strncpy(*subjectAlignSeq, subjectSeq+(*leftSubjectPos)-1, *subjectAlignSeqLen);
				(*subjectAlignSeq)[*subjectAlignSeqLen] = '\0';
			}

			if(reverseSeq(*subjectAlignSeq, *subjectAlignSeqLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the reverse sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}else
	{
		if(leftAlignSeg->strand==PLUS_STRAND)
		{
			strncpy(*subjectAlignSeq, subjectSeq+(*leftSubjectPos)-1, *subjectAlignSeqLen);
			(*subjectAlignSeq)[*subjectAlignSeqLen] = '\0';
		}else
		{
			strncpy(*subjectAlignSeq, subjectSeq+(*rightSubjectPos)-1, *subjectAlignSeqLen);
			(*subjectAlignSeq)[*subjectAlignSeqLen] = '\0';

			if(reverseSeq(*subjectAlignSeq, *subjectAlignSeqLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the reverse sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	//######################## Debug information ##########################
//	if(leftAlignSeg->strand==PLUS_STRAND)
//		orient = '+';
//	else
//		orient = '-';
//	printf("======== line=%d, In %s(), the information for the aligning base sequence are:\n", __LINE__, __func__);
//	printf("leftQueryPos=%ld, rightQueryPos=%ld, leftSubjectPos=%ld, rightSubjectPos=%ld\n", *leftQueryPos, *rightQueryPos, *leftSubjectPos, *rightSubjectPos);
//	printf("readSeq: seq=%s, len=%d\n", *queryAlignSeq, *queryAlignSeqLen);
//	printf(" refSeq: seq=%s, len=%d, orient=%c\n", *subjectAlignSeq, *subjectAlignSeqLen, orient);
	//######################## Debug information ##########################

	return SUCCESSFUL;
}

/**
 * Get the left alignment base sequences including query sequence and the reference sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getLeftAlignSeqs(char **queryAlignSeq, char **subjectAlignSeq, int32_t *queryAlignSeqLen, int32_t *subjectAlignSeqLen, int64_t *leftQueryPos, int64_t *rightQueryPos, int64_t *leftSubjectPos, int64_t *rightSubjectPos, globalValidSeg_t *leftAlignSeg, globalValidSeg_t *rightAlignSeg, query_t *queryItem, subject_t *subjectArray)
{
	char *subjectSeq;
	int32_t subjectID, subjectLen, baseNum, baseNum2;

	subjectID = rightAlignSeg->subjectID;
	subjectSeq = subjectArray[subjectID-1].subjectSeq;
	subjectLen = subjectArray[subjectID-1].subjectLen;

	if(leftAlignSeg->endQueryPos>rightAlignSeg->startQueryPos)
	{
		*leftQueryPos = rightAlignSeg->startQueryPos;
		*rightQueryPos = leftAlignSeg->endQueryPos + 300;
		if((*rightQueryPos)>queryItem->queryLen)
			*rightQueryPos = queryItem->queryLen;
	}else
	{
		*leftQueryPos = rightAlignSeg->startQueryPos;
		*rightQueryPos = rightAlignSeg->startQueryPos + 1000;
		if((*rightQueryPos)>queryItem->queryLen)
			*rightQueryPos = queryItem->queryLen;
	}

	if(queryItem->circularFlag==YES)
	{
		if(rightAlignSeg->strand==PLUS_STRAND)
		{
			*leftSubjectPos = rightAlignSeg->startSubPos;
			*rightSubjectPos = rightAlignSeg->startSubPos + ((*rightQueryPos) - (*leftQueryPos));
		}else
		{
			*leftSubjectPos = rightAlignSeg->startSubPos;
			*rightSubjectPos = rightAlignSeg->startSubPos - ((*rightQueryPos) - (*leftQueryPos));
		}
	}else
	{
		if(rightAlignSeg->strand==PLUS_STRAND)
		{
			*leftSubjectPos = rightAlignSeg->startSubPos;
			*rightSubjectPos = rightAlignSeg->startSubPos + ((*rightQueryPos) - (*leftQueryPos));
			if((*rightSubjectPos)>subjectLen)
				*rightSubjectPos = subjectLen;
		}else
		{
			*leftSubjectPos = rightAlignSeg->startSubPos;
			*rightSubjectPos = rightAlignSeg->startSubPos - ((*rightQueryPos) - (*leftQueryPos));
			if((*rightSubjectPos)<1)
				*rightSubjectPos = 1;
		}
	}

	// get the base sequences
	*queryAlignSeqLen = (*rightQueryPos) - (*leftQueryPos) + 1;
	if(rightAlignSeg->strand==PLUS_STRAND)
		*subjectAlignSeqLen = (*rightSubjectPos) - (*leftSubjectPos) + 1;
	else
		*subjectAlignSeqLen = (*leftSubjectPos) - (*rightSubjectPos) + 1;

	*queryAlignSeq = (char*)calloc((*queryAlignSeqLen)+1, sizeof(char));
	*subjectAlignSeq = (char*)calloc((*subjectAlignSeqLen)+1, sizeof(char));
	if((*queryAlignSeq)==NULL || (*subjectAlignSeq)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	strncpy(*queryAlignSeq, queryItem->querySeq+(*leftQueryPos)-1, *queryAlignSeqLen);
	(*queryAlignSeq)[*queryAlignSeqLen] = '\0';

	if(queryItem->circularFlag==YES)
	{
		if(rightAlignSeg->strand==PLUS_STRAND)
		{
			if((*rightSubjectPos)>subjectLen)
			{
				baseNum = subjectLen - (*leftSubjectPos) + 1;
				strncpy(*subjectAlignSeq, subjectSeq+(*leftSubjectPos)-1, baseNum);
				(*subjectAlignSeq)[baseNum] = '\0';
				baseNum2 = (*subjectAlignSeqLen) - baseNum;
				strncpy((*subjectAlignSeq)+baseNum, subjectSeq, baseNum2);
				(*subjectAlignSeq)[*subjectAlignSeqLen] = '\0';
			}else
			{
				strncpy(*subjectAlignSeq, subjectSeq+(*leftSubjectPos)-1, *subjectAlignSeqLen);
				(*subjectAlignSeq)[*subjectAlignSeqLen] = '\0';
			}
		}else
		{
			if((*rightSubjectPos)<1)
			{
				baseNum = 1 - (*rightSubjectPos);
				strncpy(*subjectAlignSeq, subjectSeq+subjectLen-baseNum, baseNum);
				(*subjectAlignSeq)[baseNum] = '\0';
				baseNum2 = (*subjectAlignSeqLen) - baseNum;
				strncpy(*subjectAlignSeq+baseNum, subjectSeq, baseNum2);
				(*subjectAlignSeq)[*subjectAlignSeqLen] = '\0';
			}else
			{
				strncpy(*subjectAlignSeq, subjectSeq+(*rightSubjectPos)-1, *subjectAlignSeqLen);
				(*subjectAlignSeq)[*subjectAlignSeqLen] = '\0';
			}

			if(reverseSeq(*subjectAlignSeq, *subjectAlignSeqLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the reverse sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}else
	{
		if(rightAlignSeg->strand==PLUS_STRAND)
		{
			strncpy(*subjectAlignSeq, subjectSeq+(*leftSubjectPos)-1, *subjectAlignSeqLen);
			(*subjectAlignSeq)[*subjectAlignSeqLen] = '\0';
		}else
		{
			strncpy(*subjectAlignSeq, subjectSeq+(*rightSubjectPos)-1, *subjectAlignSeqLen);
			(*subjectAlignSeq)[*subjectAlignSeqLen] = '\0';

			if(reverseSeq(*subjectAlignSeq, *subjectAlignSeqLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the reverse sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	//######################## Debug information ##########################
//	if(rightAlignSeg->strand==PLUS_STRAND)
//		orient = '+';
//	else
//		orient = '-';
//	printf("-------- line=%d, In %s(), the information for the aligning base sequence are:\n", __LINE__, __func__);
//	printf("leftQueryPos=%ld, rightQueryPos=%ld, leftSubjectPos=%ld, rightSubjectPos=%ld\n", *leftQueryPos, *rightQueryPos, *leftSubjectPos, *rightSubjectPos);
//	printf("  readSeq: seq=%s, len=%d\n", *queryAlignSeq, *queryAlignSeqLen);
//	printf("contigSeq: seq=%s, len=%d, orient=%c\n", *subjectAlignSeq, *subjectAlignSeqLen, orient);
	//######################## Debug information ##########################

	return SUCCESSFUL;
}

/**
 * Compute the reverse base sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short reverseSeq(char *seq, int32_t seqLen)
{
	int32_t i;
	char ch;

	for(i=0; i<seqLen/2; i++)
	{
		ch = seq[i];
		seq[i] = seq[seqLen-i-1];
		seq[seqLen-i-1] = ch;
	}

	for(i=0; i<seqLen; i++)
	{
		switch(seq[i])
		{
			case 'A': seq[i] = 'T'; break;
			case 'C': seq[i] = 'G'; break;
			case 'G': seq[i] = 'C'; break;
			case 'T': seq[i] = 'A'; break;
		}
	}

	return SUCCESSFUL;
}

/**
 * Compute sequence alignment.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short computeSeqAlignment(char **pAlignResultArray, int32_t *overlapLen, int32_t *mismatchNum,  int32_t *leftShiftLen1, int32_t *leftShiftLen2, int32_t *rightShiftLen1, int32_t *rightShiftLen2, char *seq1, char *seq2, int32_t seqLen1, int32_t seqLen2, int32_t printFlag)
{
	int32_t i, j, tmp;
	int32_t rowsNum, colsNum;
	int32_t maxValue, maxValueLastCol, maxValueLastRow, scoreIJ, maxRow, maxCol;
	int32_t itemNumInAlignArray, *scoreArray;
	int32_t matchScore;								// the match score
	int32_t mismatchScore;							// the mismatch score
	int32_t gapScore;								// the gap score

	// ####################### Debug information ########################
	if(printFlag==YES)
	{
		printf("seq1=%s, len=%d\nseq2=%s, len=%d\n", seq1, (int)strlen(seq1), seq2, (int)strlen(seq2));
	}
	// ####################### Debug information ########################

	matchScore = MATCH_SCORE;
	mismatchScore = MISMATCH_SCORE;
	gapScore = GAP_SCORE;


	rowsNum = seqLen1 + 1;
	colsNum = seqLen2 + 1;

	scoreArray = (int32_t*) calloc(rowsNum*colsNum, sizeof(int32_t));
	if(scoreArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initial values
	for(i=0; i<colsNum; i++)
		scoreArray[i] = 0;
	for(i=1; i<rowsNum; i++)
		scoreArray[colsNum*i] = 0;

	// compute the scores of each element
	for(i=1; i<rowsNum; i++)
	{
		for(j=1; j<colsNum; j++)
		{
			if(seq1[i-1]==seq2[j-1])
				scoreIJ = matchScore;
			else
				scoreIJ = mismatchScore;

			maxValue = INT_MIN;
			// compute the maximal score
			if(scoreArray[(i-1)*colsNum+j-1]+scoreIJ>maxValue)
			{ // from (i-1, j-1)
				maxValue = scoreArray[(i-1)*colsNum+j-1]+scoreIJ;
			}
			if(scoreArray[(i-1)*colsNum+j]+gapScore>maxValue)
			{ // from (i-1, j)
				maxValue = scoreArray[(i-1)*colsNum+j]+gapScore;
			}
			if(scoreArray[i*colsNum+j-1]+gapScore>maxValue)
			{ // from (i, j-1)
				maxValue = scoreArray[i*colsNum+j-1]+gapScore;
			}

			scoreArray[i*colsNum+j] = maxValue;
		}
	}

	// get the row and col of the maximal element in last row and column
	maxValueLastRow = INT_MIN;
	for(j=0; j<colsNum; j++)
		if(scoreArray[(rowsNum-1)*colsNum+j]>maxValueLastRow)
		{
			maxValueLastRow = scoreArray[(rowsNum-1)*colsNum+j];
			maxCol = j;
		}

	maxValueLastCol = INT_MIN;
	for(j=0; j<rowsNum; j++)
		if(scoreArray[j*colsNum+colsNum-1]>maxValueLastCol)
		{
			maxValueLastCol = scoreArray[j*colsNum+colsNum-1];
			maxRow = j;
		}

	if(maxValueLastRow>=maxValueLastCol)
	{
		maxValue = maxValueLastRow;
		maxRow = rowsNum - 1;
	}else
	{
		maxValue = maxValueLastCol;
		maxCol = colsNum - 1;
	}

	// ####################### Debug information ########################
	if(printFlag==YES)
	{
		printf("maxRow=%d, maxCol=%d, maxValue=%d\n", maxRow, maxCol, maxValue);
	}
	// ####################### Debug information ########################

	itemNumInAlignArray = 0;
	*mismatchNum = 0;
	i = maxRow;
	j = maxCol;
	while(i>0 && j>0)
	{
		if(seq1[i-1]==seq2[j-1])
			scoreIJ = matchScore;
		else
			scoreIJ = mismatchScore;

		if(scoreArray[(i-1)*colsNum+j-1]+scoreIJ==scoreArray[i*colsNum+j])
		{ // from (i-1, j-1)
			if(seq1[i-1]!=seq2[j-1])
			{
				(*mismatchNum) ++;
				//printf("0:(%d,%d)\n", i, j);

				pAlignResultArray[0][itemNumInAlignArray] = seq1[i-1];
				pAlignResultArray[1][itemNumInAlignArray] = ' ';
				pAlignResultArray[2][itemNumInAlignArray] = seq2[j-1];
				itemNumInAlignArray ++;
			}else
			{
				pAlignResultArray[0][itemNumInAlignArray] = seq1[i-1];
				pAlignResultArray[1][itemNumInAlignArray] = '|';
				pAlignResultArray[2][itemNumInAlignArray] = seq2[j-1];
				itemNumInAlignArray ++;
			}

			i --;
			j --;
		}else if(scoreArray[(i-1)*colsNum+j]+gapScore==scoreArray[i*colsNum+j])
		{ // from (i-1, j)
			(*mismatchNum) ++;
			//printf("1:(%d,%d)\n", i, j);

			pAlignResultArray[0][itemNumInAlignArray] = seq1[i-1];
			pAlignResultArray[1][itemNumInAlignArray] = ' ';
			pAlignResultArray[2][itemNumInAlignArray] = '-';
			itemNumInAlignArray ++;

			i --;
		}else
		{ // from (i, j-1)
			(*mismatchNum) ++;
			//printf("2:(%d,%d)\n", i, j);

			pAlignResultArray[0][itemNumInAlignArray] = '-';
			pAlignResultArray[1][itemNumInAlignArray] = ' ';
			pAlignResultArray[2][itemNumInAlignArray] = seq2[j-1];
			itemNumInAlignArray ++;

			j --;
		}
	}

	*leftShiftLen1 = i;
	*leftShiftLen2 = j;
	*rightShiftLen1 = rowsNum - 1 - maxRow;
	*rightShiftLen2 = colsNum - 1 - maxCol;
	*overlapLen = itemNumInAlignArray;

	// ####################### Debug information ########################
	if(printFlag==YES)
	{
		printf("i=%d, j=%d, mismatchNum=%d\n", i, j, *mismatchNum);
		printf("maxRow-i=%d, maxCol-j=%d, itemNumInAlignArray=%d\n", maxRow-i, maxCol-j, itemNumInAlignArray);
		printf("Alignment overlapLen=%d, leftShiftLen1=%d, leftShiftLen2=%d, rightShiftLen1=%d, rightShiftLen2=%d\n", *overlapLen, *leftShiftLen1, *leftShiftLen2, *rightShiftLen1, *rightShiftLen2);
	}
	// ####################### Debug information ########################

	// recover the alignment result
	for(i=0; i<3; i++)
	{
		for(j=0; j<itemNumInAlignArray/2; j++)
		{
			tmp = pAlignResultArray[i][j];
			pAlignResultArray[i][j] = pAlignResultArray[i][itemNumInAlignArray-1-j];
			pAlignResultArray[i][itemNumInAlignArray-1-j] = tmp;
		}
		pAlignResultArray[i][itemNumInAlignArray] = '\0';
	}

	// ####################### Debug information ########################
	if(printFlag==YES)
	{
		// print the alignment result
		for(i=0; i<3; i++)
			printf("%s\n", pAlignResultArray[i]);
	}
	// ####################### Debug information ########################

	free(scoreArray);

	return SUCCESSFUL;
}

/**
 * Compute right margin of query sequence alignment.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short computeRightMargin(int64_t *rightMargin, int64_t *rightMarginSubject, char **alignResultArray, int32_t overlapLen, int64_t leftQueryPos, int64_t rightQueryPos, int64_t leftSubjectPos, int64_t rightSubjectPos, int32_t queryLeftShiftLen, int32_t subjectLeftShiftLen, globalValidSeg_t *globalSeg)
{
	int64_t i, j;
	int32_t mismatchNum, totalMismatchNum, mismatchNum2, winSize, baseNum, baseNum2, marginBlockID, marginBlockID2, remainSize, startRow, endRow, targetRow, gapNum;
	double misDensity, misDensity2, misDensityThreshold;
	char *pMatchSeq, *pQueryMatchSeq, *pSubjectMatchSeq;
	int32_t newRightQueryPos, newRightSubjectPos, mismatchNumTrimmed, gapNumTrimmed, subjectPosShiftFactor;

	winSize = MISMATCH_WIN_SIZE;
	misDensityThreshold = MIS_DENSITY_THRES;

	newRightQueryPos = newRightSubjectPos = -1;
	pMatchSeq = alignResultArray[1];
	totalMismatchNum = mismatchNum = 0;
	marginBlockID = -1;
	baseNum = 0;
	for(i=0; i<overlapLen; i++)
	{
		if(pMatchSeq[i]!='|')
		{
			mismatchNum ++;
			totalMismatchNum ++;
		}

		baseNum ++;
		if(baseNum%winSize==0)
		{
			misDensity = (double)mismatchNum / winSize;
			if(misDensity>misDensityThreshold)
			{
				// check another margin block
				mismatchNum2 = 0;
				marginBlockID2 = -1;
				baseNum2 = baseNum;
				for(j=i+1; j<overlapLen; j++)
				{
					if(pMatchSeq[j]!='|')
						mismatchNum2 ++;

					baseNum2 ++;
					if(baseNum2%winSize==0)
					{
						misDensity2 = (double)mismatchNum2 / winSize;
						if(misDensity2>misDensityThreshold)
						{
							marginBlockID2 = (baseNum2 - 1) / winSize;
							break;
						}
						mismatchNum2 = 0;
					}
				}

				if(marginBlockID2>=0)
				{
					marginBlockID = (baseNum - 1) / winSize;
					break;
				}
			}
			mismatchNum = 0;
		}
	}

	if(marginBlockID==-1 && mismatchNum>0)
	{
		remainSize = overlapLen % winSize;
		misDensity = (double)mismatchNum / remainSize;
		if(misDensity>misDensityThreshold)
			marginBlockID = (overlapLen - 1) / winSize;
	}

	// compute the new right query position and right subject position
	if(globalSeg->strand==PLUS_STRAND)
		subjectPosShiftFactor = 1;
	else
		subjectPosShiftFactor = -1;

//	if(marginBlockID==-1)
//	{
//		startRow = overlapLen - winSize;
//		endRow = overlapLen - 1;
//		if(startRow<0)
//			startRow = 0;
//		targetRow = -1;
//		for(i=startRow; i<=endRow; i++)
//		{
//			if(pMatchSeq[i]!='|')
//			{
//				for(j=i-1; j>=0; j--)
//				{
//					if(pMatchSeq[j]=='|')
//					{
//						targetRow = j;
//						break;
//					}
//				}
//
//				if(targetRow>=0)
//					break;
//			}
//		}
//		if(targetRow<0)
//		{
//			if(i>endRow)
//				targetRow = endRow;
//			else
//			{
//				printf("line=%d, In %s(), targetRow=%d, error!\n", __LINE__, __func__, targetRow);
//				return FAILED;
//			}
//		}
//
//		// get the query position
//		gapNum = 0;
//		pQueryMatchSeq = alignResultArray[0];
//		for(i=0; i<=targetRow; i++)
//		{
//			if(pQueryMatchSeq[i]=='-')
//				gapNum ++;
//		}
//		newRightQueryPos = leftQueryPos + queryLeftShiftLen + targetRow - gapNum;
//
//		// get the subject position
//		gapNum = 0;
//		pSubjectMatchSeq = alignResultArray[2];
//		for(i=0; i<=targetRow; i++)
//		{
//			if(pSubjectMatchSeq[i]=='-')
//				gapNum ++;
//		}
//		newRightSubjectPos = leftSubjectPos + (subjectLeftShiftLen + targetRow - gapNum) * subjectPosShiftFactor;
//
//		// compute the mismatchNumTrimmed, gapNumTrimmed
//		mismatchNumTrimmed = gapNumTrimmed = 0;
//		for(i=targetRow+1; i<overlapLen; i++)
//		{
//			if(alignResultArray[1][i]!='|')
//			{
//				if(alignResultArray[0][i]=='-' || alignResultArray[2][i]=='-')
//					gapNumTrimmed ++;
//				mismatchNumTrimmed ++;
//			}
//		}
//
//		globalSeg->endQueryPos = newRightQueryPos;
//		globalSeg->endSubPos = newRightSubjectPos;
//		globalSeg->matchLen -= overlapLen - 1 - targetRow - mismatchNumTrimmed;
//		globalSeg->totalMatchLen -= overlapLen - 1 - targetRow;
//		if(globalSeg->matchLen>globalSeg->totalMatchLen)
//			globalSeg->matchLen = globalSeg->totalMatchLen;
//		globalSeg->matchPercent = (double)globalSeg->matchLen / globalSeg->totalMatchLen;
//		globalSeg->gapNum -= gapNumTrimmed;
//		if(globalSeg->gapNum>globalSeg->totalMatchLen-globalSeg->matchLen)
//			globalSeg->gapNum = globalSeg->totalMatchLen - globalSeg->matchLen;
//	}else
	if(marginBlockID>=0)
	{
		startRow = (marginBlockID - 1) * winSize;
		endRow = startRow + 2 * winSize - 1;
		if(startRow<0)
			startRow = 0;
		if(endRow>overlapLen-1)
			endRow = overlapLen - 1;
		targetRow = -1;
		for(i=startRow; i<=endRow; i++)
		{
			if(pMatchSeq[i]!='|')
			{
				for(j=i-1; j>=0; j--)
				{
					if(pMatchSeq[j]=='|')
					{
						targetRow = j;
						break;
					}
				}

				if(targetRow>=0)
					break;
			}
		}
		if(targetRow<0)
		{
			if(i>endRow)
				targetRow = endRow;
			else
			{
				printf("line=%d, In %s(), targetRow=%d, error!\n", __LINE__, __func__, targetRow);
				return FAILED;
			}
		}

		// get the query position
		gapNum = 0;
		pQueryMatchSeq = alignResultArray[0];
		for(i=0; i<=targetRow; i++)
		{
			if(pQueryMatchSeq[i]=='-')
				gapNum ++;
		}
		newRightQueryPos = leftQueryPos + queryLeftShiftLen + targetRow - gapNum;

		// get the subject position
		gapNum = 0;
		pSubjectMatchSeq = alignResultArray[2];
		for(i=0; i<=targetRow; i++)
		{
			if(pSubjectMatchSeq[i]=='-')
				gapNum ++;
		}
		newRightSubjectPos = leftSubjectPos + (subjectLeftShiftLen + targetRow - gapNum) * subjectPosShiftFactor;

		// compute the mismatchNumTrimmed, gapNumTrimmed
		mismatchNumTrimmed = gapNumTrimmed = 0;
		for(i=targetRow+1; i<overlapLen; i++)
		{
			if(alignResultArray[1][i]!='|')
			{
				if(alignResultArray[0][i]=='-' || alignResultArray[2][i]=='-')
					gapNumTrimmed ++;
				mismatchNumTrimmed ++;
			}
		}

		globalSeg->endQueryPos = newRightQueryPos;
		globalSeg->endSubPos = newRightSubjectPos;
		globalSeg->matchLen -= overlapLen - 1 - targetRow - mismatchNumTrimmed;
		globalSeg->totalMatchLen -= overlapLen - 1 - targetRow;
		if(globalSeg->matchLen>globalSeg->totalMatchLen)
			globalSeg->matchLen = globalSeg->totalMatchLen;
		globalSeg->matchPercent = (double)globalSeg->matchLen / globalSeg->totalMatchLen;
		globalSeg->gapNum -= gapNumTrimmed;
		if(globalSeg->gapNum>globalSeg->totalMatchLen-globalSeg->matchLen)
			globalSeg->gapNum = globalSeg->totalMatchLen - globalSeg->matchLen;
	}else
	{
		newRightQueryPos = globalSeg->endQueryPos;
		newRightSubjectPos = globalSeg->endSubPos;
	}

	if(newRightQueryPos==-1)
	{
		printf("line=%d, In %s(), cannot get the right margin of the query match information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	*rightMargin = newRightQueryPos;
	*rightMarginSubject = newRightSubjectPos;
	//printf("rightMargin=%ld, rightMarginSubject=%ld\n", *rightMargin, *rightMarginSubject);

	return SUCCESSFUL;
}

/**
 * Compute left margin of query sequence alignment.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short computeLeftMargin(int64_t *leftMargin, int64_t *leftMarginSubject, char **alignResultArray, int32_t overlapLen, int64_t leftQueryPos, int64_t rightQueryPos, int64_t leftSubjectPos, int64_t rightSubjectPos, int32_t queryRightShiftLen, int32_t subjectRightShiftLen, globalValidSeg_t *globalSeg)
{
	int64_t i, j;
	int32_t mismatchNum, totalMismatchNum, mismatchNum2, winSize, baseNum, baseNum2, marginBlockID, marginBlockID2, remainSize, startRow, endRow, targetRow, gapNum;
	double misDensity, misDensity2, misDensityThreshold;
	char *pMatchSeq, *pQueryMatchSeq, *pSubjectMatchSeq;
	int32_t newLeftQueryPos, newLeftSubjectPos, mismatchNumTrimmed, gapNumTrimmed, subjectPosShiftFactor;

	winSize = MISMATCH_WIN_SIZE;
	misDensityThreshold = MIS_DENSITY_THRES;

	newLeftQueryPos = newLeftSubjectPos = -1;
	pMatchSeq = alignResultArray[1];
	totalMismatchNum = mismatchNum = 0;
	marginBlockID = -1;
	baseNum = 0;
	for(i=overlapLen-1; i>=0; i--)
	{
		if(pMatchSeq[i]!='|')
		{
			mismatchNum ++;
			totalMismatchNum ++;
		}

		baseNum ++;
		if(baseNum%winSize==0)
		{
			misDensity = (double)mismatchNum / winSize;
			if(misDensity>misDensityThreshold)
			{
				// check another margin block
				mismatchNum2 = 0;
				marginBlockID2 = -1;
				baseNum2 = baseNum;
				for(j=i-1; j>=0; j--)
				{
					if(pMatchSeq[j]!='|')
						mismatchNum2 ++;

					baseNum2 ++;
					if(baseNum2%winSize==0)
					{
						misDensity2 = (double)mismatchNum2 / winSize;
						if(misDensity2>misDensityThreshold)
						{
							marginBlockID2 = (baseNum2 - 1) / winSize;
							break;
						}
						mismatchNum2 = 0;
					}
				}

				if(marginBlockID2>=0)
				{
					marginBlockID = (baseNum - 1) / winSize;
					break;
				}
				mismatchNum = 0;
			}
		}
	}

	if(marginBlockID==-1 && mismatchNum>0)
	{
		remainSize = overlapLen % winSize;
		misDensity = (double)mismatchNum / remainSize;
		if(misDensity>misDensityThreshold)
			marginBlockID = (overlapLen - 1) / winSize;
	}

	// compute the new left query position and left subject position
	if(globalSeg->strand==PLUS_STRAND)
		subjectPosShiftFactor = 1;
	else
		subjectPosShiftFactor = -1;

//	if(marginBlockID==-1)
//	{
//		startRow = winSize - 1;
//		endRow = 0;
//		if(startRow>overlapLen-1)
//			startRow = overlapLen - 1;
//		targetRow = -1;
//		for(i=startRow; i>=endRow; i--)
//		{
//			if(pMatchSeq[i]!='|')
//			{
//				for(j=i+1; j<overlapLen; j++)
//				{
//					if(pMatchSeq[j]=='|')
//					{
//						targetRow = j;
//						break;
//					}
//				}
//				if(targetRow>=0)
//					break;
//			}
//		}
//		if(targetRow<0)
//		{
//			if(i<endRow)
//				targetRow = endRow;
//			else
//			{
//				printf("line=%d, In %s(), targetRow=%d, error!\n", __LINE__, __func__, targetRow);
//				return FAILED;
//			}
//		}
//
//		// get the query position
//		gapNum = 0;
//		pQueryMatchSeq = alignResultArray[0];
//		for(i=overlapLen-1; i>=targetRow; i--)
//		{
//			if(pQueryMatchSeq[i]=='-')
//				gapNum ++;
//		}
//		newLeftQueryPos = rightQueryPos - (overlapLen-1+queryRightShiftLen-targetRow) + gapNum;
//
//		// get the subject position
//		gapNum = 0;
//		pSubjectMatchSeq = alignResultArray[2];
//		for(i=overlapLen-1; i>=targetRow; i--)
//		{
//			if(pSubjectMatchSeq[i]=='-')
//				gapNum ++;
//		}
//		newLeftSubjectPos = rightSubjectPos - (overlapLen-1 + subjectRightShiftLen - targetRow - gapNum) * subjectPosShiftFactor;
//
//		// compute the mismatchNumTrimmed, gapNumTrimmed
//		mismatchNumTrimmed = gapNumTrimmed = 0;
//		for(i=0; i<targetRow; i++)
//		{
//			if(alignResultArray[1][i]!='|')
//			{
//				if(alignResultArray[0][i]=='-' || alignResultArray[2][i]=='-')
//					gapNumTrimmed ++;
//				mismatchNumTrimmed ++;
//			}
//		}
//
//		globalSeg->startQueryPos = newLeftQueryPos;
//		globalSeg->startSubPos = newLeftSubjectPos;
//		globalSeg->matchLen -= targetRow - mismatchNumTrimmed;
//		globalSeg->totalMatchLen -= targetRow;
//		if(globalSeg->matchLen>globalSeg->totalMatchLen)
//			globalSeg->matchLen = globalSeg->totalMatchLen;
//		globalSeg->matchPercent = (double)globalSeg->matchLen / globalSeg->totalMatchLen;
//		globalSeg->gapNum -= gapNumTrimmed;
//		if(globalSeg->gapNum>globalSeg->totalMatchLen-globalSeg->matchLen)
//			globalSeg->gapNum = globalSeg->totalMatchLen - globalSeg->matchLen;
//	}else
	if(marginBlockID>=0)
	{
		startRow = overlapLen - 1 - (marginBlockID - 1) * winSize;
		endRow = startRow - 2 * winSize + 1;
		if(startRow>overlapLen-1)
			startRow = overlapLen - 1;
		if(endRow<0)
			endRow = 0;
		targetRow = -1;
		for(i=startRow; i>=endRow; i--)
		{
			if(pMatchSeq[i]!='|')
			{
				for(j=i+1; j<overlapLen; j++)
				{
					if(pMatchSeq[j]=='|')
					{
						targetRow = j;
						break;
					}
				}
				if(targetRow>=0)
					break;
			}
		}
		if(targetRow<0)
		{
			if(i<endRow)
				targetRow = endRow;
			else
			{
				printf("line=%d, In %s(), targetRow=%d, error!\n", __LINE__, __func__, targetRow);
				return FAILED;
			}
		}

		// get the query position
		gapNum = 0;
		pQueryMatchSeq = alignResultArray[0];
		for(i=overlapLen-1; i>=targetRow; i--)
		{
			if(pQueryMatchSeq[i]=='-')
				gapNum ++;
		}
		newLeftQueryPos = rightQueryPos - (overlapLen-1-targetRow) + gapNum;

		// get the subject position
		gapNum = 0;
		pSubjectMatchSeq = alignResultArray[2];
		for(i=overlapLen-1; i>=targetRow; i--)
		{
			if(pSubjectMatchSeq[i]=='-')
				gapNum ++;
		}
		newLeftSubjectPos = rightSubjectPos - (overlapLen-1 + subjectRightShiftLen - targetRow - gapNum) * subjectPosShiftFactor;

		// compute the mismatchNumTrimmed, gapNumTrimmed
		mismatchNumTrimmed = gapNumTrimmed = 0;
		for(i=0; i<targetRow; i++)
		{
			if(alignResultArray[1][i]!='|')
			{
				if(alignResultArray[0][i]=='-' || alignResultArray[2][i]=='-')
					gapNumTrimmed ++;
				mismatchNumTrimmed ++;
			}
		}

		globalSeg->startQueryPos = newLeftQueryPos;
		globalSeg->startSubPos = newLeftSubjectPos;
		globalSeg->matchLen -= targetRow - mismatchNumTrimmed;
		globalSeg->totalMatchLen -= targetRow;
		if(globalSeg->matchLen>globalSeg->totalMatchLen)
			globalSeg->matchLen = globalSeg->totalMatchLen;
		globalSeg->matchPercent = (double)globalSeg->matchLen / globalSeg->totalMatchLen;
		globalSeg->gapNum -= gapNumTrimmed;
		if(globalSeg->gapNum>globalSeg->totalMatchLen-globalSeg->matchLen)
			globalSeg->gapNum = globalSeg->totalMatchLen - globalSeg->matchLen;
	}else
	{
		newLeftQueryPos = globalSeg->startQueryPos;
		newLeftSubjectPos = globalSeg->startSubPos;
	}

	if(newLeftQueryPos==-1)
	{
		printf("line=%d, In %s(), cannot get the left margin of the query match information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	*leftMargin = newLeftQueryPos;
	*leftMarginSubject = newLeftSubjectPos;
	//printf("leftMargin=%ld, leftMarginSubject=%ld\n", *leftMargin, *leftMarginSubject);

	return SUCCESSFUL;
}

/**
 * Compute the ratios of normal regions.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short computeNormalRatios(double *SP_ratio, double *SMinus_ratio, double *SPlus_ratio, double *discorRatio, queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet, double insertSize, double standDev)
{
	ratioRegion_t *ratioRegionArray;
	int32_t i, maxRatioRegionNum, ratioRegionNum, subRegNum, subRegSize;

	subRegSize = 500;

	if(prepareRatioRegionArray(&ratioRegionArray, &maxRatioRegionNum, subRegSize, queryMatchInfoSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot prepare the ratio array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	ratioRegionNum = 0;
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		if(queryMatchInfoSet->queryArray[i].queryID==queryMatchInfoSet->maxQueryID || queryMatchInfoSet->queryArray[i].queryID==queryMatchInfoSet->secQueryID)
		{
			if(computeNormalRatiosSingleQuery(ratioRegionArray+ratioRegionNum, &subRegNum, subRegSize, queryMatchInfoSet->queryArray+i, readSet, insertSize, standDev)==FAILED)
			{
				printf("line=%d, In %s(), cannot compute the ratio array for single query, error!\n", __LINE__, __func__);
				return FAILED;
			}
			ratioRegionNum += subRegNum;
		}
	}

	if(ratioRegionNum>maxRatioRegionNum)
	{
		printf("line=%d, In %s(), ratioRegionNum=%d, maxRatioRegionNum=%d, error!\n", __LINE__, __func__, ratioRegionNum, maxRatioRegionNum);
		return FAILED;
	}

	// compute their average values
	if(computeAverRatios(SP_ratio, SMinus_ratio, SPlus_ratio, discorRatio, ratioRegionArray, ratioRegionNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the average ratios, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// free memory
	if(maxRatioRegionNum>0)
		free(ratioRegionArray);

	return SUCCESSFUL;
}

/**
 * Prepare the the memory of ratio array for normal regions of queries.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short prepareRatioRegionArray(ratioRegion_t **ratioRegionArray, int32_t *maxRatioRegionNum, int32_t subRegSize, queryMatchInfo_t *queryMatchInfoSet)
{
	int64_t sumLen;

	sumLen = 0;
	if(queryMatchInfoSet->secQueryID>0)
	{
		sumLen += queryMatchInfoSet->queryArray[queryMatchInfoSet->maxQueryID-1].queryLen;
		sumLen += queryMatchInfoSet->queryArray[queryMatchInfoSet->secQueryID-1].queryLen;
	}else if(queryMatchInfoSet->maxQueryID>0)
		sumLen += queryMatchInfoSet->queryArray[queryMatchInfoSet->maxQueryID-1].queryLen;
	else
	{
		printf("line=%d, In %s(), no valid queries, error!\n", __LINE__, __func__);
		return FAILED;
	}

	*maxRatioRegionNum = (sumLen - 2 *subRegSize) / subRegSize;
	if((*maxRatioRegionNum)>=1)
	{
		*ratioRegionArray = (ratioRegion_t*) calloc (*maxRatioRegionNum, sizeof(ratioRegion_t));
		if((*ratioRegionArray)==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		printf("line=%d, In %s(), no valid queries, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Compute the ratios of normal regions of single query.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short computeNormalRatiosSingleQuery(ratioRegion_t *ratioRegionArray, int32_t *ratioRegionNum, int32_t subRegSize, query_t *queryItem, readSet_t *readSet, double insertSize, double standDev)
{
	double SP_ratio, SMinus_ratio, SPlus_ratio, discorRatio;

	if(queryItem->querySeq)
	{
		// initialize the ratioRegion array
		if(initRatioRegionArrayNormal(ratioRegionArray, ratioRegionNum, subRegSize, queryItem)==FAILED)
		{
			printf("line=%d, In %s(), cannot initialize the ratioRegion array, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// fill the ratioRegion array
		if(fillRatioRegionArray(ratioRegionArray, *ratioRegionNum, queryItem, readSet, insertSize, standDev)==FAILED)
		{
			printf("line=%d, In %s(), cannot fill the ratioRegion array, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute their average values
		if(computeAverRatios(&SP_ratio, &SMinus_ratio, &SPlus_ratio, &discorRatio, ratioRegionArray, *ratioRegionNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the average ratios, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		printf("line=%d, In %s(), no query sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Compute the ratios of breakpoint regions.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short computeBreakpointRatios(query_t *queryItem, int32_t misjoinRegNum, readSet_t *readSet, double insertSize, double standDev)
{
	ratioRegion_t *ratioRegionArray;
	int32_t ratioRegionNum;

	if(misjoinRegNum>0)
	{
		ratioRegionNum = misjoinRegNum;

		// initialize the ratioRegion array
		if(initRatioRegionArrayBreakpoint(&ratioRegionArray, ratioRegionNum, queryItem)==FAILED)
		{
			printf("line=%d, In %s(), cannot initialize the ratioRegion array, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// fill the ratioRegion array
		if(fillRatioRegionArray(ratioRegionArray, ratioRegionNum, queryItem, readSet, insertSize, standDev)==FAILED)
		{
			printf("line=%d, In %s(), cannot fill the ratioRegion array, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute their values
		if(computeRatios(ratioRegionArray, ratioRegionNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the ratios, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// update ratios in query
		if(updateRatiosInQuery(queryItem, ratioRegionArray, ratioRegionNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the ratios, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// ###################### Debug information #########################
		//outputRatioRegionArray(ratioRegionArray, ratioRegionNum);
		// ###################### Debug information #########################

		// free memory
		if(ratioRegionNum>0)
			free(ratioRegionArray);
	}

	return SUCCESSFUL;
}

/**
 * Initialize the ratioRegion array of normal query regions.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short initRatioRegionArrayNormal(ratioRegion_t *ratioRegionArray, int32_t *ratioRegionNum, int32_t subRegSize, query_t *queryItem)
{
	int32_t i, j, globalSegNum, itemNum, initItemNum;
	globalValidSeg_t *globalSegArray;
	int64_t startQueryPos, endQueryPos, newStartQueryPos, newEndQueryPos, startRegPos, endRegPos, midPos;

	globalSegArray = queryItem->globalValidSegArray;
	globalSegNum = queryItem->globalValidSegNum;

	initItemNum = (queryItem->queryLen - 2 *subRegSize) / subRegSize;
	if(initItemNum<=0)
	{
		*ratioRegionNum = 0;
		return SUCCESSFUL;
	}

	itemNum = 0;
	for(i=0; i<globalSegNum; i++)
	{
		if(i<globalSegNum-1)
		{
			startQueryPos = globalSegArray[i].startQueryPos;
			if(globalSegArray[i].endQueryPos>globalSegArray[i+1].startQueryPos)
				endQueryPos = globalSegArray[i+1].startQueryPos;
			else
				endQueryPos = globalSegArray[i].endQueryPos;
		}else if(i>0)
		{
			if(globalSegArray[i].startQueryPos<globalSegArray[i-1].endQueryPos)
				startQueryPos = globalSegArray[i-1].endQueryPos;
			else
				startQueryPos = globalSegArray[i].startQueryPos;
			endQueryPos = globalSegArray[i].endQueryPos;
		}else
		{
			startQueryPos = globalSegArray[i].startQueryPos;
			endQueryPos = globalSegArray[i].endQueryPos;
		}

		// adjust the startQueryPos and endQueryPos
		newStartQueryPos = startQueryPos + 1000;
		newEndQueryPos = endQueryPos - 1000;
		if(newStartQueryPos<newEndQueryPos)
		{
			for(j=newStartQueryPos; j<=newEndQueryPos; j+=subRegSize)
			{
				startRegPos = j;
				endRegPos = startRegPos + subRegSize - 1;
				if(endRegPos>newEndQueryPos)
					break;

				midPos = (startRegPos + endRegPos) / 2;
				ratioRegionArray[itemNum].midQPos = midPos;
				ratioRegionArray[itemNum].startQPosLHalf = startRegPos;
				ratioRegionArray[itemNum].endQPosLHalf = midPos;
				ratioRegionArray[itemNum].startQPosRHalf = midPos + 1;
				ratioRegionArray[itemNum].endQPosRHalf = endRegPos;

				ratioRegionArray[itemNum].disagreeNum = 0;
				ratioRegionArray[itemNum].zeroCovNum = 0;
				ratioRegionArray[itemNum].discorNum = 0;

				ratioRegionArray[itemNum].SPRatio = -1;
				ratioRegionArray[itemNum].singleMinusRatio = -1;
				ratioRegionArray[itemNum].singlePlusRatio = -1;
				ratioRegionArray[itemNum].discorRatio = -1;
				itemNum ++;
			}
		}
	}
	*ratioRegionNum = itemNum;

	if(itemNum>initItemNum)
	{
		printf("line=%d, In %s(), itemNum=%d, initItemNum=%d, error!\n", __LINE__, __func__, itemNum, initItemNum);
		return FAILED;
	}

	// ###################### Debug information #########################
	//outputRatioRegionArray(ratioRegionArray, *ratioRegionNum);
	// ###################### Debug information #########################

	return SUCCESSFUL;
}

/**
 * Initialize the ratioRegion array of breakpoint query regions.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short initRatioRegionArrayBreakpoint(ratioRegion_t **ratioRegionArray, int32_t ratioRegionNum, query_t *queryItem)
{
	int32_t itemNum, midPos, misjoinNum, indelNum;
	misInfo_t *misInfo;
	queryMargin_t *queryMargin;

	if(ratioRegionNum>0)
	{
		// allocate memory
		*ratioRegionArray = (ratioRegion_t*) calloc (ratioRegionNum, sizeof(ratioRegion_t));
		if((*ratioRegionArray)==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		itemNum = 0;
		misInfo = queryItem->misInfoList;
		while(misInfo)
		{
			if(misInfo->misType==QUERY_MISJOIN_KIND)
			{
				queryMargin = misInfo->queryMargin;
				midPos = (queryMargin->leftMargin + queryMargin->rightMargin) / 2;
				(*ratioRegionArray)[itemNum].midQPos = midPos;

				if(queryMargin->leftMargin<=queryMargin->rightMargin)
				{
					(*ratioRegionArray)[itemNum].startQPosLHalf = queryMargin->leftMargin - 100;
					(*ratioRegionArray)[itemNum].endQPosLHalf = midPos;
					(*ratioRegionArray)[itemNum].startQPosRHalf = midPos + 1;
					(*ratioRegionArray)[itemNum].endQPosRHalf = queryMargin->rightMargin + 100;
					if((*ratioRegionArray)[itemNum].startQPosLHalf<1)
						(*ratioRegionArray)[itemNum].startQPosLHalf = 1;
					if((*ratioRegionArray)[itemNum].endQPosRHalf>queryItem->queryLen)
						(*ratioRegionArray)[itemNum].endQPosRHalf = queryItem->queryLen;

				}else
				{
					(*ratioRegionArray)[itemNum].startQPosLHalf = midPos + 1;
					(*ratioRegionArray)[itemNum].endQPosLHalf = queryMargin->leftMargin + 300;
					if((*ratioRegionArray)[itemNum].endQPosLHalf>queryItem->queryLen)
						(*ratioRegionArray)[itemNum].endQPosLHalf = queryItem->queryLen;

					(*ratioRegionArray)[itemNum].startQPosRHalf = queryMargin->rightMargin - 300;
					(*ratioRegionArray)[itemNum].endQPosRHalf = midPos;
					if((*ratioRegionArray)[itemNum].startQPosRHalf<1)
						(*ratioRegionArray)[itemNum].startQPosRHalf = 1;
				}

				(*ratioRegionArray)[itemNum].disagreeNum = 0;
				(*ratioRegionArray)[itemNum].zeroCovNum = 0;
				(*ratioRegionArray)[itemNum].discorNum = 0;

				(*ratioRegionArray)[itemNum].SPRatio = -1;
				(*ratioRegionArray)[itemNum].singleMinusRatio = -1;
				(*ratioRegionArray)[itemNum].singlePlusRatio = -1;
				(*ratioRegionArray)[itemNum].discorRatio = -1;
				itemNum ++;
			}

			misInfo = misInfo->next;
		}
	}else
	{
		*ratioRegionArray = NULL;
	}

	return SUCCESSFUL;
}

/**
 * Fill the ratioRegion array.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short fillRatioRegionArray(ratioRegion_t *ratioRegionArray, int32_t ratioRegionNum, query_t *queryItem, readSet_t *readSet, double insertSize, double standDev)
{
	int64_t i, j, queryID, queryID_paired, readID, readID_paired, queryPos, queryPos_paired, orient, orient_paired, seqLen, seqLen_paired;
	int64_t midPos, regRow, pairedFlag, sideFlag, discorFlag;
	double fragSize, difFragSize;
	queryRead_t *queryRead;
	int32_t *regIDArray, startRow, endRow;

	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;
	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo_paired;


	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;


	// fill the regIDArray
	regIDArray = (int32_t*) malloc (queryItem->queryLen * sizeof(int32_t));
	if(regIDArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(memset(regIDArray, -1, queryItem->queryLen*sizeof(int32_t))==NULL)
	{
		printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<ratioRegionNum; i++)
	{
		if(ratioRegionArray[i].startQPosLHalf<ratioRegionArray[i].startQPosRHalf)
			startRow = ratioRegionArray[i].startQPosLHalf - 1;
		else
			startRow = ratioRegionArray[i].startQPosRHalf - 1;

		if(ratioRegionArray[i].endQPosLHalf<ratioRegionArray[i].endQPosRHalf)
			endRow = ratioRegionArray[i].endQPosRHalf - 1;
		else
			endRow = ratioRegionArray[i].endQPosLHalf - 1;

		for(j=startRow; j<=endRow; j++) regIDArray[j] = i;
	}

	// fill the reads data
	queryID = queryItem->queryID;
	for(i=0; i<queryItem->queryReadNum; i++)
	{
		queryRead = queryItem->queryReadArray + i;

		readID = queryRead->readID;
		queryPos = queryRead->queryPos;
		orient = queryRead->orientation;
		seqLen = queryRead->seqlen;
		pairedFlag = -1;

		// get the regID
		midPos = queryPos + seqLen / 2 - 1;
		regRow = regIDArray[midPos-1];
		if(regRow>=0 && regRow<=ratioRegionNum-1)
		{
			if(ratioRegionArray[regRow].endQPosLHalf<=ratioRegionArray[regRow].midQPos)
			{
				if(midPos<=ratioRegionArray[regRow].midQPos)
					sideFlag = 1;
				else
					sideFlag = 2;
			}else
			{
				if(midPos<=ratioRegionArray[regRow].midQPos)
					sideFlag = 2;
				else
					sideFlag = 1;
			}

			if((readID&1)==1)
				readID_paired = readID + 1;
			else
				readID_paired = readID - 1;

			readMatchInfoBlockID = (readID_paired - 1) / maxItemNumPerReadMatchInfoBlock;
			rowNumInReadMatchInfoBlock = (readID_paired - 1) % maxItemNumPerReadMatchInfoBlock;
			pReadMatchInfo_paired = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;

			discorFlag = NO;
			if(pReadMatchInfo_paired)
			{
				queryID_paired = pReadMatchInfo_paired->queryID;
				queryPos_paired = pReadMatchInfo_paired->queryPos;
				orient_paired = pReadMatchInfo_paired->readOrientation;
				seqLen_paired = pReadMatchInfo_paired->seqlen;

				if(queryID_paired==queryID)
				{
					if(orient!=orient_paired)
					{
						if(orient==ORIENTATION_PLUS && orient_paired==ORIENTATION_MINUS)
						{
							if(queryPos<=queryPos_paired)
							{
								fragSize = queryPos_paired + seqLen_paired - queryPos;
								difFragSize = fragSize - insertSize;
								if(difFragSize<-3*standDev || difFragSize>3*standDev)
									discorFlag = YES;
							}else
								discorFlag = YES;
						}else if(orient==ORIENTATION_MINUS && orient_paired==ORIENTATION_PLUS)
						{
							if(queryPos>=queryPos_paired)
							{
								fragSize = queryPos + seqLen - queryPos_paired;
								difFragSize = fragSize - insertSize;
								if(difFragSize<-3*standDev || difFragSize>3*standDev)
									discorFlag = YES;
							}else
							{
								discorFlag = YES;
							}
						}
					}else
						discorFlag = YES;
					pairedFlag = YES;
				}else
					pairedFlag = NO;
			}else
				pairedFlag = NO;

			if(pairedFlag==NO)
			{
				if(sideFlag==1)
				{ // left half side
					if(orient==ORIENTATION_MINUS)
						ratioRegionArray[regRow].singleMinusNum ++;
					ratioRegionArray[regRow].singleNumLeftHalf ++;
				}else
				{ // right half side
					if(orient==ORIENTATION_PLUS)
						ratioRegionArray[regRow].singlePlusNum ++;
					ratioRegionArray[regRow].singleNumRightHalf ++;
				}
				ratioRegionArray[regRow].singleNum ++;
			}else
				ratioRegionArray[regRow].pairedNum ++;

			if(discorFlag==YES)
				ratioRegionArray[regRow].discorNum ++;
		}
	}

	free(regIDArray);

	return SUCCESSFUL;
}

/**
 * Compute the average ratios.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short computeAverRatios(double *SP_ratio, double *SMinus_ratio, double *SPlus_ratio, double *discorRatio, ratioRegion_t *ratioRegionArray, int32_t ratioRegionNum)
{
	int32_t i, itemNum_SP_ratio, itemNum_SMinus_ratio, itemNum_SPlus_ratio, itemNum_Discor_ratio;

	// compute their values
	if(computeRatios(ratioRegionArray, ratioRegionNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the ratios, error!\n", __LINE__, __func__);
		return FAILED;
	}

	*SP_ratio = *SMinus_ratio = *SPlus_ratio = *discorRatio = 0;
	itemNum_SP_ratio = itemNum_SMinus_ratio = itemNum_SPlus_ratio = itemNum_Discor_ratio = 0;

	for(i=0; i<ratioRegionNum; i++)
	{
		if(ratioRegionArray[i].pairedNum>0)
		{
			*SP_ratio += ratioRegionArray[i].SPRatio;
			itemNum_SP_ratio ++;
			*discorRatio += ratioRegionArray[i].discorRatio;
			itemNum_Discor_ratio ++;
		}

		if(ratioRegionArray[i].singleNumLeftHalf>0)
		{
			*SMinus_ratio += ratioRegionArray[i].singleMinusRatio;
			itemNum_SMinus_ratio ++;
		}

		if(ratioRegionArray[i].singleNumRightHalf>0)
		{
			*SPlus_ratio += ratioRegionArray[i].singlePlusRatio;
			itemNum_SPlus_ratio ++;
		}
	}

	if(itemNum_SP_ratio>0)
	{
		*SP_ratio /= itemNum_SP_ratio;
		*discorRatio /= itemNum_Discor_ratio;
	}else
		*SP_ratio = *discorRatio = -1;

	if(itemNum_SMinus_ratio>0)
	{
		*SMinus_ratio /= itemNum_SP_ratio;
		*SPlus_ratio /= itemNum_SPlus_ratio;
	}else
		*SMinus_ratio = *SPlus_ratio = -1;

	if((*SP_ratio)<0 || (*SMinus_ratio)<0 || (*SPlus_ratio)<0 || (*discorRatio)<0)
	{
		printf("line=%d, In %s(), SP_ratio=%.4f, SMinus_ratio=%.4f, SPlus_ratio=%.4f, discorRatio=%e, error!\n", __LINE__, __func__, *SP_ratio, *SMinus_ratio, *SPlus_ratio, *discorRatio);
		return FAILED;
	}

	//printf("SP_ratio=%.4f, SMinus_ratio=%.4f, SPlus_ratio=%.4f, discorRatio=%.4f, itemNum_SP_ratio=%d, itemNum_SMinus_ratio=%d, itemNum_SPlus_ratio=%d, itemNum_Discor_ratio=%d\n", *SP_ratio, *SMinus_ratio, *SPlus_ratio, *discorRatio, itemNum_SP_ratio, itemNum_SMinus_ratio, itemNum_SPlus_ratio, itemNum_Discor_ratio);

	return SUCCESSFUL;
}

/**
 * Compute the ratios.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short computeRatios(ratioRegion_t *ratioRegionArray, int32_t ratioRegionNum)
{
	int32_t i;

	for(i=0; i<ratioRegionNum; i++)
	{
		if(ratioRegionArray[i].pairedNum>0)
			ratioRegionArray[i].SPRatio = (double)ratioRegionArray[i].singleNum / ratioRegionArray[i].pairedNum;
		else
			ratioRegionArray[i].SPRatio = -1;

		if(ratioRegionArray[i].singleNumLeftHalf>0)
			ratioRegionArray[i].singleMinusRatio = (double)ratioRegionArray[i].singleMinusNum / ratioRegionArray[i].singleNumLeftHalf;
		else
			ratioRegionArray[i].singleMinusRatio = -1;

		if(ratioRegionArray[i].singleNumRightHalf>0)
			ratioRegionArray[i].singlePlusRatio = (double)ratioRegionArray[i].singlePlusNum / ratioRegionArray[i].singleNumRightHalf;
		else
			ratioRegionArray[i].singlePlusRatio = -1;

		if(ratioRegionArray[i].pairedNum>0)
		{
			ratioRegionArray[i].discorRatio = (double)ratioRegionArray[i].discorNum / ratioRegionArray[i].pairedNum;
		}else if(ratioRegionArray[i].discorNum>0)
		{
			printf("line=%d, In %s(), discorNum=%d, pairedNum=%d, error!\n", __LINE__, __func__, ratioRegionArray[i].discorNum, ratioRegionArray[i].pairedNum);
			return FAILED;
		}
	}

	// ###################### Debug information #########################
	//outputRatioRegionArray(ratioRegionArray, ratioRegionNum);
	// ###################### Debug information #########################

	return SUCCESSFUL;
}

/**
 * Update the ratios in query.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short updateRatiosInQuery(query_t *queryItem, ratioRegion_t *ratioRegionArray, int32_t ratioRegionNum)
{
	int32_t i;
	misInfo_t *misInfo;
	queryMargin_t *queryMargin;

	if(ratioRegionNum>0)
	{
		i = 0;
		misInfo = queryItem->misInfoList;
		while(misInfo)
		{
			if(misInfo->misType==QUERY_MISJOIN_KIND)
			{
				queryMargin = misInfo->queryMargin;
				queryMargin->SPRatio = ratioRegionArray[i].SPRatio;
				queryMargin->singleMinusRatio = ratioRegionArray[i].singleMinusRatio;
				queryMargin->singlePlusRatio = ratioRegionArray[i].singlePlusRatio;
				queryMargin->discorRatio = ratioRegionArray[i].discorRatio;
				queryMargin->discorNum = ratioRegionArray[i].discorNum;

				i ++;
			}
			misInfo = misInfo->next;
		}
	}

	return SUCCESSFUL;
}

/**
 * Output the ratios.
 */
void outputRatioRegionArray(ratioRegion_t *ratioRegionArray, int32_t ratioRegionNum)
{
	int32_t i;

	if(ratioRegionNum>0)
	{
		for(i=0; i<ratioRegionNum; i++)
			printf("ratioRegionArray[%d]: Left[%d, %d, %d], Right[%d, %d, %d]; SPRatio=%.4f, SMinusRatio=%.4f, SPlusRatio=%.4f, discorRatio=%.4f, pairedNum=%d, singleNum=%d, SMinusNum=%d, SNumLeft=%d, SPlusNum=%d, SNumRight=%d, discorNum=%d\n", i, ratioRegionArray[i].startQPosLHalf, ratioRegionArray[i].endQPosLHalf, ratioRegionArray[i].endQPosLHalf-ratioRegionArray[i].startQPosLHalf+1, ratioRegionArray[i].startQPosRHalf, ratioRegionArray[i].endQPosRHalf, ratioRegionArray[i].endQPosRHalf-ratioRegionArray[i].startQPosRHalf+1, ratioRegionArray[i].SPRatio, ratioRegionArray[i].singleMinusRatio, ratioRegionArray[i].singlePlusRatio, ratioRegionArray[i].discorRatio, ratioRegionArray[i].pairedNum, ratioRegionArray[i].singleNum, ratioRegionArray[i].singleMinusNum, ratioRegionArray[i].singleNumLeftHalf, ratioRegionArray[i].singlePlusNum, ratioRegionArray[i].singleNumRightHalf, ratioRegionArray[i].discorNum);
	}else
	{
		printf("There are no ratio regions.\n");
	}
}

/**
 * Determine the mis-assembly flag of query.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short determineMisassFlag(query_t *queryItem, int32_t misjoinRegNum, baseCov_t *baseCovArray, double SP_ratio_Thres, double SMinus_ratio_Thres, double SPlus_ratio_Thres)
{
	int32_t startQueryPos, endQueryPos;
	misInfo_t *misInfo;
	queryMargin_t *queryMargin;

	if(misjoinRegNum>0)
	{
		misInfo = queryItem->misInfoList;
		while(misInfo)
		{
			if(misInfo->misType==QUERY_MISJOIN_KIND)
			{
				queryMargin = misInfo->queryMargin;
				if(queryMargin->leftMargin<queryMargin->rightMargin)
				{
					startQueryPos = queryMargin->leftMargin - 100;
					endQueryPos = queryMargin->rightMargin + 100;
				}else
				{
					startQueryPos = queryMargin->rightMargin - 100;
					endQueryPos = queryMargin->leftMargin + 100;
				}

				if(startQueryPos<1)
					startQueryPos = 1;
				if(endQueryPos>queryItem->queryLen)
					endQueryPos = queryItem->queryLen;

				// compute the disagreements
				if(computeDisagreements(&queryMargin->disagreeNum, &queryMargin->zeroCovNum, baseCovArray, startQueryPos-1, endQueryPos-1, NO)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute the disagreements of single query, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// compute the highCovRegNum and lowCovRegNum of the region
				if(computeAbnormalCovRegNum(&queryMargin->highCovRegNum, &queryMargin->lowCovRegNum, baseCovArray, startQueryPos-1, endQueryPos-1, queryItem->queryLen, -1)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute the covRatio, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(queryMargin->zeroCovNum>3 || ((queryMargin->disagreeNum>=2 || queryMargin->zeroCovNum>0 || queryMargin->highCovRegNum>3 || queryMargin->lowCovRegNum>3) && (queryMargin->SPRatio>SP_ratio_Thres*3 || (queryMargin->singleMinusRatio>0.7 || queryMargin->singlePlusRatio>0.7) || (queryMargin->discorNum>=3 && queryMargin->discorRatio>0.1))))
					misInfo->misassFlag = queryMargin->misassFlag = TRUE_MISASS;
				else if(queryMargin->disagreeNum>=2 && (queryMargin->zeroCovNum>0 || queryMargin->highCovRegNum>3 || queryMargin->lowCovRegNum>3))
					misInfo->misassFlag = queryMargin->misassFlag = TRUE_MISASS;
				else
					misInfo->misassFlag = queryMargin->misassFlag = UNCERTAIN_MISASS;
			}

			misInfo = misInfo->next;
		}
	}

	return SUCCESSFUL;
}

/**
 * Save the mis-assembly queries.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short saveMisassQueries(char *errorsFile, char *svFile, char *misUncertainFile, char *gapFile, queryMatchInfo_t *queryMatchInfoSet)
{
	query_t *queryItem;
	misInfo_t *misInfo;
	misassSeq_t *misassSeq;
	char misKindStr[256], *subjectTitle, nullSubjectTitle[256];
	FILE *fpErrors, *fpSV, *fpUncertain, *fpGap;
	int32_t i, subjectID, errNum, svNum, gapNum, warningNum;
	subject_t *subjectArray;

	fpErrors = fopen(errorsFile, "w");
	if(fpErrors==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, errorsFile);
		return FAILED;
	}

	fpSV = fopen(svFile, "w");
	if(fpSV==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, svFile);
		return FAILED;
	}

	fpUncertain = fopen(misUncertainFile, "w");
	if(fpUncertain==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, misUncertainFile);
		return FAILED;
	}

	fpGap = fopen(gapFile, "w");
	if(fpGap==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, gapFile);
		return FAILED;
	}

	fprintf(fpErrors, "query\tmisKind\tqueryStart\tqueryEnd\tsubjectStart\tsubjectEnd\tSubject\n");
	fprintf(fpSV, "query\tmisKind\tqueryStart\tqueryEnd\tsubjectStart\tsubjectEnd\tSubject\n");
	fprintf(fpUncertain, "query\tmisKind\tqueryStart\tqueryEnd\tsubjectStart\tsubjectEnd\tSubject\n");
	fprintf(fpGap, "query\tmisKind\tqueryStart\tqueryEnd\tsubjectStart\tsubjectEnd\tSubject\n");

	subjectArray = queryMatchInfoSet->subjectArray;
	strcpy(nullSubjectTitle, "-");

	errNum = svNum = gapNum = warningNum = 0;
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		queryItem = queryMatchInfoSet->queryArray + i;
		if(queryItem->misInfoItemNum>0)
		{
			misInfo = queryItem->misInfoList;
			while(misInfo)
			{
				misassSeq = misInfo->misassSeqList;
				while(misassSeq)
				{
					if(misassSeq->subjectID>0)
						subjectTitle = subjectArray[misassSeq->subjectID-1].subjectTitle;
					else
						subjectTitle = nullSubjectTitle;

					if(misInfo->misassFlag==TRUE_MISASS)
					{
						if(misInfo->gapFlag==NO || misInfo->misType==QUERY_MISJOIN_KIND)
						{
							switch(misassSeq->misassKind)
							{
								case ERR_MISJOIN: strcpy(misKindStr, "err_misjoin"); break;
								case ERR_INSERT: strcpy(misKindStr, "err_insertion"); break;
								case ERR_DEL: strcpy(misKindStr, "err_deletion"); break;
								default: printf("line=%d, In %s(), invalid error kind=%d, error!\n", __LINE__, __func__, misassSeq->misassKind); return FAILED;
							}

							fprintf(fpErrors, "%s\t%s\t%d\t%d\t%d\t%d\t%s\n", queryItem->queryTitle, misKindStr, misassSeq->startQueryPos, misassSeq->endQueryPos, misassSeq->startSubjectPos, misassSeq->endSubjectPos, subjectTitle);
							errNum ++;
						}else if(misInfo->gapFlag==YES)
						{
							switch(misassSeq->misassKind)
							{
								case GAP_INSERT: strcpy(misKindStr, "gap_insertion"); break;
								case GAP_DEL: strcpy(misKindStr, "gap_deletion"); break;
								default: printf("line=%d, In %s(), invalid error kind=%d, error!\n", __LINE__, __func__, misassSeq->misassKind); return FAILED;
							}

							fprintf(fpGap, "%s\t%s\t%d\t%d\t%d\t%d\t%s\n", queryItem->queryTitle, misKindStr, misassSeq->startQueryPos, misassSeq->endQueryPos, misassSeq->startSubjectPos, misassSeq->endSubjectPos, subjectTitle);
							gapNum ++;
						}else
						{
							printf("line=%d, In %s(), invalid situation, error!\n", __LINE__, __func__);
							return FAILED;
						}

					}else if(misInfo->misassFlag==STRUCTURE_VARIATION)
					{
						switch(misassSeq->misassKind)
						{
							case SV_MISJOIN: strcpy(misKindStr, "sv_misjoin"); break;
							case SV_INSERT: strcpy(misKindStr, "sv_insertion"); break;
							case SV_DEL: strcpy(misKindStr, "sv_deletion"); break;
							default: printf("line=%d, In %s(), invalid SV kind=%d, error!\n", __LINE__, __func__, misassSeq->misassKind); return FAILED;
						}

						fprintf(fpSV, "%s\t%s\t%d\t%d\t%d\t%d\t%s\n", queryItem->queryTitle, misKindStr, misassSeq->startQueryPos, misassSeq->endQueryPos, misassSeq->startSubjectPos, misassSeq->endSubjectPos, subjectTitle);
						svNum ++;
					}else if(misInfo->misassFlag==UNCERTAIN_MISASS)
					{
						strcpy(misKindStr, "uncertain");

						fprintf(fpUncertain, "%s\t%s\t%d\t%d\t%d\t%d\t%s\n", queryItem->queryTitle, misKindStr, misassSeq->startQueryPos, misassSeq->endQueryPos, misassSeq->startSubjectPos, misassSeq->endSubjectPos, subjectTitle);
						warningNum ++;
					}else
					{
						printf("line=%d, In %s(), invalid misassFlag=%d, error!\n", __LINE__, __func__, misInfo->misassFlag);
						return FAILED;
					}

					misassSeq = misassSeq->next;
				}

				misInfo = misInfo->next;
			}
		}
	}

	printf("**** In summary, the statistics are as follows ****\n");
	printf("Number of mis-assemblies                : %d\n", errNum);
	printf("Number of correct assemblies due to SVs : %d\n", svNum);
	printf("Number of mis-assemblies in gap regions : %d\n", gapNum);
	printf("Number of warnings                      : %d\n", warningNum);

	fclose(fpErrors);
	fclose(fpSV);
	fclose(fpUncertain);
	fclose(fpGap);

	return SUCCESSFUL;
}

/**
 * Save the new queries after correction.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short saveNewQueries(char *newQueryFile, queryMatchInfo_t *queryMatchInfoSet)
{
	query_t *queryItem;
	querySeq_t *querySeqNode;
	char newTitle[256];
	char numCh[10];
	FILE *fpNew;
	int32_t i, newID;

	fpNew = fopen(newQueryFile, "w");
	if(fpNew==NULL)
	{
		printf("In %s(), cannot open file [ %s ], error!\n", __func__, newQueryFile);
		return FAILED;
	}

	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		queryItem = queryMatchInfoSet->queryArray + i;

		if(queryItem->misassFlag==TRUE_MISASS)
		{
			newID = 1;
			querySeqNode = queryItem->newQuerySeqList;
			while(querySeqNode)
			{
				strcpy(newTitle, queryItem->queryTitle);
				strcat(newTitle, "_");
				sprintf(numCh, "%d", newID);
				strcat(newTitle, numCh);
				fprintf(fpNew, ">%s_%d_%d\t%d\n%s\n", newTitle, querySeqNode->startQueryPos, querySeqNode->endQueryPos, querySeqNode->queryLen, querySeqNode->querySeq);

				querySeqNode = querySeqNode->next;
				newID ++;
			}
		}else
		{
			fprintf(fpNew, ">%s\t%d\n%s\n", queryItem->queryTitle, queryItem->queryLen, queryItem->querySeq);
		}
	}

	fclose(fpNew);

	return SUCCESSFUL;
}
