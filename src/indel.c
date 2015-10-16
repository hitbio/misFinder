/*
 * indel.c
 *
 *  Created on: Jul 6, 2014
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Determine the queryIndel information for single query.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short determineQueryIndel(query_t *queryItem, baseCov_t *baseCovArray, subject_t *subjectArray, readSet_t *readSet)
{
	misInfo_t *misInfo;
	queryIndel_t *queryIndel;
	int32_t indelRegNum;

	indelRegNum = 0;
	misInfo = queryItem->misInfoList;
	while(misInfo)
	{
		if(misInfo->misType==QUERY_INDEL_KIND)
		{
			indelRegNum ++;
			break;
		}
		misInfo = misInfo->next;
	}

	if(indelRegNum>0)
	{
		misInfo = queryItem->misInfoList;
		while(misInfo)
		{
			if(misInfo->misType==QUERY_INDEL_KIND)
			{
				queryIndel = misInfo->queryIndel;
				if(confirmQueryIndelKind(queryIndel, queryItem, misInfo->innerFlag, baseCovArray, subjectArray, readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot confirm query insert information, error!\n", __LINE__, __func__);
					return FAILED;
				}
				misInfo->misassFlag = queryIndel->misassFlag;
			}

			misInfo = misInfo->next;
		}
	}

	return SUCCESSFUL;
}

/**
 * Confirm the query insert information for single query.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short confirmQueryIndelKind(queryIndel_t *queryIndel, query_t *queryItem, int32_t innerFlag, baseCov_t *baseCovArray, subject_t *subjectArray, readSet_t *readSet)
{
	int32_t globalSegNum, leftMargin, rightMargin;
	globalValidSeg_t *globalSegArray, *leftSeg, *rightSeg;
	int32_t startQueryPosLeft, endQueryPosLeft, startQueryPosFragLeft, endQueryPosFragLeft;
	int32_t startQueryPosRight, endQueryPosRight, startQueryPosFragRight, endQueryPosFragRight;
	int32_t startRegPos, endRegPos, subRegSize, ratioRegionNum, leftPos, rightPos;
	ratioRegion_t *ratioRegionArray;
	double insertSize, standDev;

	insertSize = readSet->insertSize;
	standDev = readSet->standDev;
	globalSegArray = queryItem->globalValidSegArray;
	globalSegNum = queryItem->globalValidSegNum;

	if(queryIndel->leftSegRow<0) // head
		leftSeg = NULL;
	else
		leftSeg = globalSegArray + queryIndel->leftSegRow;

	if(queryIndel->rightSegRow<0) // tail
		rightSeg = NULL;
	else
		rightSeg = globalSegArray + queryIndel->rightSegRow;

	if(leftSeg==NULL && rightSeg)
	{ // head
		// get the right region
		rightMargin = queryIndel->rightMargin;
		startQueryPosRight = rightMargin;
		endQueryPosRight = rightMargin + insertSize;
		if(endQueryPosRight>queryItem->queryLen)
			endQueryPosRight = queryItem->queryLen;

		leftMargin = rightMargin - 1;
		endQueryPosLeft = leftMargin + END_BACK_CHECK_SIZE;
		startQueryPosLeft = leftMargin - insertSize;
		if(startQueryPosLeft<1)
			startQueryPosLeft = 1;
		if(endQueryPosLeft>queryItem->queryLen)
			endQueryPosLeft = queryItem->queryLen;

		endQueryPosFragLeft = leftMargin;
		startQueryPosFragLeft = endQueryPosFragLeft - insertSize;
		if(startQueryPosFragLeft<1)
			startQueryPosFragLeft = 1;

		// compute the disagreements
		if(computeDisagreements(&queryIndel->disagreeNum, &queryIndel->zeroCovNum, &queryIndel->disagreeRegSize, baseCovArray, startQueryPosLeft-1, endQueryPosRight-1, NO)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the disagreements, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute the highCovRegNum and lowCovRegNum of the region
		if(computeAbnormalCovRegNum(&queryIndel->highCovRegNum, &queryIndel->lowCovRegNum, baseCovArray, startQueryPosLeft-1, endQueryPosRight-1, queryItem, 1, NO)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the covRatio, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute the fragment size of paired-end reads in the left region
		if(computeFragSizeLeftRegQueryIndel(queryIndel, startQueryPosFragLeft, endQueryPosFragLeft, queryItem, readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the fragment size, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// determine the final queryIndel kind
		if(determineFinalQueryIndelKind(queryIndel, innerFlag, insertSize, standDev)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the fragment size, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else if(leftSeg && rightSeg==NULL)
	{ // tail
		// get the right region
		leftMargin = queryIndel->leftMargin;
		endQueryPosLeft = leftMargin;
		startQueryPosLeft = leftMargin - insertSize;
		if(startQueryPosLeft<1)
			startQueryPosLeft = 1;

		rightMargin = leftMargin + 1;
		startQueryPosRight = rightMargin - END_BACK_CHECK_SIZE;
		endQueryPosRight = rightMargin + insertSize;
		if(startQueryPosRight<1)
			startQueryPosRight = 1;
		if(endQueryPosRight>queryItem->queryLen)
			endQueryPosRight = queryItem->queryLen;

		startQueryPosFragRight = rightMargin;
		endQueryPosFragRight = startQueryPosFragRight + insertSize;
		if(endQueryPosFragRight>queryItem->queryLen)
			endQueryPosFragRight = queryItem->queryLen;

		// compute the disagreements
		if(computeDisagreements(&queryIndel->disagreeNum, &queryIndel->zeroCovNum, &queryIndel->disagreeRegSize, baseCovArray, startQueryPosLeft-1, endQueryPosRight-1, NO)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the disagreements, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute the highCovRegNum and lowCovRegNum of the region
		if(computeAbnormalCovRegNum(&queryIndel->highCovRegNum, &queryIndel->lowCovRegNum, baseCovArray, startQueryPosLeft-1, endQueryPosRight-1, queryItem, 2, NO)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the covRatio, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute the fragment size of paired-end reads in the right region
		if(computeFragSizeRightRegQueryIndel(queryIndel, startQueryPosFragRight, endQueryPosFragRight, queryItem, readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the fragment size, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// determine the final queryIndel kind
		if(determineFinalQueryIndelKind(queryIndel, innerFlag, insertSize, standDev)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the fragment size, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else if((leftSeg && rightSeg) || innerFlag==YES)
	{ // middle
		// get the left and right region
		leftMargin = queryIndel->leftMargin;
		endQueryPosLeft = leftMargin;
		startQueryPosLeft = leftMargin - insertSize;
		if(startQueryPosLeft<1)
			startQueryPosLeft = 1;

		rightMargin = queryIndel->rightMargin;
		startQueryPosRight = rightMargin;
		endQueryPosRight = startQueryPosRight + insertSize;
		if(endQueryPosRight>queryItem->queryLen)
			endQueryPosRight = queryItem->queryLen;

		// compute the disagreements
		if(computeDisagreements(&queryIndel->disagreeNum, &queryIndel->zeroCovNum, &queryIndel->disagreeRegSize, baseCovArray, startQueryPosLeft-1, endQueryPosRight-1, NO)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the disagreements, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute the highCovRegNum and lowCovRegNum of the region
		if(computeAbnormalCovRegNum(&queryIndel->highCovRegNum, &queryIndel->lowCovRegNum, baseCovArray, startQueryPosLeft-1, endQueryPosRight-1, queryItem, -1, NO)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the covRatio, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute the fragment size of paired-end reads between the two regions
		if(computeFragSizeBothRegQueryIndel(queryIndel, startQueryPosLeft, endQueryPosLeft, startQueryPosRight, endQueryPosRight, queryItem, readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the fragment size, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(computeStatisticsRatioRegionQueryIndel(queryIndel, baseCovArray, readSet, queryItem)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute statistics in the ratioRegion array for query indel, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// determine the final queryIndel kind
		if(determineFinalQueryIndelKind(queryIndel, innerFlag, insertSize, standDev)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the fragment size, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else // if(leftSeg==NULL && rightSeg==NULL)
	{
		// get the left and right region
		leftMargin = queryIndel->leftMargin;
		endQueryPosLeft = leftMargin;
		startQueryPosLeft = leftMargin - insertSize;
		if(startQueryPosLeft<1)
			startQueryPosLeft = 1;

		rightMargin = queryIndel->rightMargin;
		startQueryPosRight = rightMargin;
		endQueryPosRight = startQueryPosRight + insertSize;
		if(endQueryPosRight>queryItem->queryLen)
			endQueryPosRight = queryItem->queryLen;

		// compute the fragment size of paired-end reads between the two regions
		if(computeFragSizeBothRegQueryIndel(queryIndel, startQueryPosLeft, endQueryPosLeft, startQueryPosRight, endQueryPosRight, queryItem, readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the fragment size, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(computeStatisticsRatioRegionQueryIndel(queryIndel, baseCovArray, readSet, queryItem)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute statistics in the ratioRegion array for query indel, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// determine the final queryIndel kind
		if(determineFinalQueryIndelKindUnmatched(queryIndel)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the determine query indel mis-assembly kind, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Compute the fragment size of paired-end reads in the left region for single query.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeFragSizeLeftRegQueryIndel(queryIndel_t *queryIndel, int32_t startQueryPosLeft, int32_t endQueryPosLeft, query_t *queryItem, readSet_t *readSet)
{
	int64_t i, pairNumValid, queryID, queryID_paired, readID, readID_paired, queryPos, queryPos_paired, orient, orient_paired, seqLen, seqLen_paired, midPos;
	int32_t discorNum, setID, queryReadNum;
	queryRead_t *queryRead, *queryReadArray;
	double insertSize, standDev, fragSize, difFragSize, fragSizeSum;
	char orientCh1, orientCh2;

	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;
	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo_paired;

	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

	setID = readSet->setID;
	insertSize = readSet->insertSize;
	standDev = readSet->standDev;

	discorNum = 0;
	fragSizeSum = 0;
	pairNumValid = 0;
	queryID = queryItem->queryID;

	queryReadArray = queryItem->queryReadSetArray[setID-1].queryReadArray;
	queryReadNum = queryItem->queryReadSetArray[setID-1].queryReadNum;
	for(i=0; i<queryReadNum; i++)
	{
		queryRead = queryReadArray + i;

		readID = queryRead->readID;
		queryPos = queryRead->queryPos;
		orient = queryRead->orientation;
		seqLen = queryRead->seqlen;

		midPos = queryPos + (seqLen / 2) - 1;
		if(midPos>=startQueryPosLeft && midPos<=endQueryPosLeft)
		{
			if((readID&1)==1)
				readID_paired = readID + 1;
			else
				readID_paired = readID - 1;

			readMatchInfoBlockID = (readID_paired - 1) / maxItemNumPerReadMatchInfoBlock;
			rowNumInReadMatchInfoBlock = (readID_paired - 1) % maxItemNumPerReadMatchInfoBlock;
			pReadMatchInfo_paired = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;

			if(pReadMatchInfo_paired)
			{
				queryID_paired = pReadMatchInfo_paired->queryID;
				queryPos_paired = pReadMatchInfo_paired->queryPos;
				orient_paired = pReadMatchInfo_paired->readOrientation;
				seqLen_paired = pReadMatchInfo_paired->seqlen;

				if(queryID_paired==queryID)
				{
					if(orient==ORIENTATION_PLUS && orient_paired==ORIENTATION_MINUS)
					{
						if(queryPos<=queryPos_paired)
						{
							fragSize = queryPos_paired + seqLen_paired - queryPos;
							difFragSize = fragSize - insertSize;
							if(queryPos_paired+(seqLen_paired/2)>endQueryPosLeft)
							{
								fragSizeSum += fragSize;
								pairNumValid ++;
							}

							if(difFragSize<-3*standDev || difFragSize>3*standDev)
								discorNum ++;
						}else
						{
							discorNum ++;

							if(orient==ORIENTATION_PLUS)
								orientCh1 = '+';
							else
								orientCh1 = '-';
							if(orient_paired==ORIENTATION_PLUS)
								orientCh2 = '+';
							else
								orientCh2 = '-';
							//printf("line=%d, invalid read pair: (%ld, %ld, %ld, %c), (%ld, %ld, %ld, %c)\n", __LINE__, readID, queryPos, seqLen, orientCh1, readID_paired, queryPos_paired, seqLen_paired, orientCh2);
						}
					}else if(orient==ORIENTATION_MINUS && orient_paired==ORIENTATION_PLUS)
					{
						if(queryPos>=queryPos_paired)
						{
							fragSize = queryPos + seqLen - queryPos_paired;
							difFragSize = fragSize - insertSize;
							if(difFragSize<-3*standDev || difFragSize>3*standDev)
								discorNum ++;
						}else
						{
							discorNum ++;

							if(orient==ORIENTATION_PLUS)
								orientCh2 = '+';
							else
								orientCh2 = '-';
							if(orient_paired==ORIENTATION_PLUS)
								orientCh1 = '+';
							else
								orientCh1 = '-';
							//printf("line=%d, invalid read pair: (%ld, %ld, %ld, %c), (%ld, %ld, %ld, %c)\n", __LINE__, readID_paired, queryPos_paired, seqLen_paired, orientCh1, readID, queryPos, seqLen, orientCh2);
						}
					}
					else// if(orient==orient_paired)
						discorNum ++;
				}
			}
		}
	}

	queryIndel->pairNumLeft = pairNumValid;
	queryIndel->discorNumLeft = discorNum;
	if(pairNumValid>0)
	{
		queryIndel->averFragSizeLeft = fragSizeSum / pairNumValid;
		queryIndel->difFragSizeLeft = queryIndel->averFragSizeLeft - insertSize;

		//queryIndel->discorRatioLeft = (double)discorNum / pairNumTotal;
		queryIndel->discorRatioLeft = (double)discorNum / pairNumValid;
	}else
	{
		queryIndel->averFragSizeLeft = 0;
		queryIndel->difFragSizeLeft = 0;
		queryIndel->discorRatioLeft = 0;
	}

	return SUCCESSFUL;
}

/**
 * Compute the fragment size of paired-end reads in the right region for single query.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeFragSizeRightRegQueryIndel(queryIndel_t *queryIndel, int32_t startQueryPosRight, int32_t endQueryPosRight, query_t *queryItem, readSet_t *readSet)
{
	int64_t i, pairNumValid, queryID, queryID_paired, readID, readID_paired, queryPos, queryPos_paired, orient, orient_paired, seqLen, seqLen_paired, midPos;
	int32_t discorNum, setID, queryReadNum;
	queryRead_t *queryRead, *queryReadArray;
	double insertSize, standDev, fragSize, fragSizeSum;
	char orientCh1, orientCh2;

	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;
	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo_paired;

	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

	setID = readSet->setID;
	insertSize = readSet->insertSize;
	standDev = readSet->standDev;

	discorNum = 0;
	fragSizeSum = 0;
	pairNumValid = 0;
	queryID = queryItem->queryID;
	queryReadArray = queryItem->queryReadSetArray[setID-1].queryReadArray;
	queryReadNum = queryItem->queryReadSetArray[setID-1].queryReadNum;
	for(i=0; i<queryReadNum; i++)
	{
		queryRead = queryReadArray + i;

		readID = queryRead->readID;
		queryPos = queryRead->queryPos;
		orient = queryRead->orientation;
		seqLen = queryRead->seqlen;

		midPos = queryPos + (seqLen / 2) - 1;
		if(midPos>=startQueryPosRight && midPos<=endQueryPosRight)
		{
			if((readID&1)==1)
				readID_paired = readID + 1;
			else
				readID_paired = readID - 1;

			readMatchInfoBlockID = (readID_paired - 1) / maxItemNumPerReadMatchInfoBlock;
			rowNumInReadMatchInfoBlock = (readID_paired - 1) % maxItemNumPerReadMatchInfoBlock;
			pReadMatchInfo_paired = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;

			if(pReadMatchInfo_paired)
			{
				queryID_paired = pReadMatchInfo_paired->queryID;
				queryPos_paired = pReadMatchInfo_paired->queryPos;
				orient_paired = pReadMatchInfo_paired->readOrientation;
				seqLen_paired = pReadMatchInfo_paired->seqlen;

				if(queryID_paired==queryID)
				{
					if(orient==ORIENTATION_MINUS && orient_paired==ORIENTATION_PLUS)
					{
						if(queryPos>=queryPos_paired)
						{
							fragSize = queryPos + seqLen - queryPos_paired;
							if(queryPos_paired+(seqLen_paired/2)<startQueryPosRight)
							{
								fragSizeSum += fragSize;
								pairNumValid ++;
							}

							if(fragSize<insertSize-3*standDev || fragSize>insertSize+3*standDev)
								discorNum ++;
						}else
						{
							discorNum ++;

							if(orient==ORIENTATION_PLUS)
								orientCh2 = '+';
							else
								orientCh2 = '-';
							if(orient_paired==ORIENTATION_PLUS)
								orientCh1 = '+';
							else
								orientCh1 = '-';
							//printf("line=%d, invalid read pair: (%ld, %ld, %ld, %c), (%ld, %ld, %ld, %c)\n", __LINE__, readID_paired, queryPos_paired, seqLen_paired, orientCh1, readID, queryPos, seqLen, orientCh2);
						}
					}else if(orient==ORIENTATION_PLUS && orient_paired==ORIENTATION_MINUS)
					{
						if(queryPos<=queryPos_paired)
						{
							fragSize = queryPos_paired + seqLen_paired - queryPos;
							if(fragSize<insertSize-3*standDev || fragSize>insertSize+3*standDev)
								discorNum ++;
						}else
						{
							discorNum ++;

							if(orient==ORIENTATION_PLUS)
								orientCh1 = '+';
							else
								orientCh1 = '-';
							if(orient_paired==ORIENTATION_PLUS)
								orientCh2 = '+';
							else
								orientCh2 = '-';
							//printf("line=%d, invalid read pair: (%ld, %ld, %ld, %c), (%ld, %ld, %ld, %c)\n", __LINE__, readID, queryPos, seqLen, orientCh1, readID_paired, queryPos_paired, seqLen_paired, orientCh2);
						}
					}
					else// if(orient==orient_paired)
						discorNum ++;
				}
			}
		}
	}

	queryIndel->pairNumRight = pairNumValid;
	queryIndel->discorNumRight = discorNum;
	if(pairNumValid>0)
	{
		queryIndel->averFragSizeRight = fragSizeSum / pairNumValid;
		queryIndel->difFragSizeRight = queryIndel->averFragSizeRight - insertSize;

		//queryIndel->discorRatioRight = (double)discorNum / pairNumTotal;
		queryIndel->discorRatioRight = (double)discorNum / pairNumValid;
	}else
	{
		queryIndel->averFragSizeRight = 0;
		queryIndel->difFragSizeRight = 0;
		queryIndel->discorRatioRight = 0;
	}

	return SUCCESSFUL;
}

/**
 * Compute the fragment size of paired-end reads for the left and right regions for single query.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeFragSizeBothRegQueryIndel(queryIndel_t *queryIndel, int32_t startQueryPosLeft, int32_t endQueryPosLeft, int32_t startQueryPosRight, int32_t endQueryPosRight, query_t *queryItem, readSet_t *readSet)
{
	// compute the fragment size of paired-end reads in the left region
	if(computeFragSizeLeftRegQueryIndel(queryIndel, startQueryPosLeft, endQueryPosLeft, queryItem, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the fragment size, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// compute the fragment size of paired-end reads in the right region
	if(computeFragSizeRightRegQueryIndel(queryIndel, startQueryPosRight, endQueryPosRight, queryItem, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the fragment size, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Determine the final queryIndel kind.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short determineFinalQueryIndelKind(queryIndel_t *queryIndel, int32_t innerFlag, double insertSize, double standDev)
{
	double difFragSize, newDifFragSize, dif, difAll;

	if(queryIndel->leftSegRow<0 && queryIndel->rightSegRow>=0)
	{ // head
		if(queryIndel->disagreeNum>0 || queryIndel->zeroCovNum>0 || (queryIndel->lowCovRegNum>1 || queryIndel->highCovRegNum>1) || queryIndel->pairNumLeft==0)
			queryIndel->misassFlag = TRUE_MISASS;
		else
		{
			difFragSize = queryIndel->averFragSizeLeft - insertSize;
			if(difFragSize<0)
				difFragSize = -difFragSize;
			if(queryIndel->pairNumLeft>=3 && (difFragSize>2*standDev || difFragSize>0.2*insertSize))
				queryIndel->misassFlag = TRUE_MISASS;
			else
				queryIndel->misassFlag = UNCERTAIN_MISASS;
		}

		// ######################## Debug information #######################
		//printf("disagreeNum=%d, averFragSizeLeft=%.2f, validPairNumLeft=%d, discorRatioLeft=%.4f, discorNumLeft=%d, insertSize=%.2f, standDev=%.2f, misassFlag=%d\n", queryIndel->disagreeNum, queryIndel->averFragSizeLeft, queryIndel->pairNumLeft, queryIndel->discorRatioLeft, queryIndel->discorNumLeft, insertSize, standDev, queryIndel->misassFlag);
		// ######################## Debug information #######################

	}else if(queryIndel->leftSegRow>=0 && queryIndel->rightSegRow<0)
	{ // tail
		if(queryIndel->disagreeNum>0 || queryIndel->zeroCovNum>0 || (queryIndel->lowCovRegNum>1 || queryIndel->highCovRegNum>1) || queryIndel->pairNumRight==0)
			queryIndel->misassFlag = TRUE_MISASS;
		else
		{
			difFragSize = queryIndel->averFragSizeRight - insertSize;
			if(difFragSize<0)
				difFragSize = -difFragSize;
			if(queryIndel->pairNumRight>=3 && (difFragSize>2*standDev || difFragSize>0.2*insertSize))
				queryIndel->misassFlag = TRUE_MISASS;
			else
				queryIndel->misassFlag = UNCERTAIN_MISASS;
		}

		// ######################## Debug information #######################
		//printf("disagreeNum=%d, averFragSizeRight=%.2f, validPairNumRight=%d, discorRatioRight=%.4f, discorNumRight=%d, insertSize=%.2f, standDev=%.2f, misassFlag=%d\n", queryIndel->disagreeNum, queryIndel->averFragSizeRight, queryIndel->pairNumRight, queryIndel->discorRatioRight, queryIndel->discorNumRight, insertSize, standDev, queryIndel->misassFlag);
		// ######################## Debug information #######################
	}else if((queryIndel->leftSegRow>=0 && queryIndel->rightSegRow>=0) || innerFlag==YES)
	{ // mid

		if(queryIndel->zeroCovNum>0)
			queryIndel->misassFlag = TRUE_MISASS;
		else if(queryIndel->pairNumLeft>=3 && queryIndel->pairNumRight>=3)
		{
			difFragSize = (queryIndel->difFragSizeLeft + queryIndel->difFragSizeRight) / 2;

			if(queryIndel->queryIndelKind==QUERY_INSERT)
			{ // query insertion
				dif = queryIndel->difQuery - queryIndel->difSubject;
				difAll = dif - difFragSize;
				if(difAll<0)
					difAll = -difAll;

				if((difAll<2*standDev || difAll<0.2*insertSize) || (queryIndel->discorRatioLeft>0.1 || queryIndel->discorRatioRight>0.1))
				{
					if(queryIndel->disagreeNum>0 || queryIndel->zeroCovNum>0 || (queryIndel->lowCovRegNum>1 || queryIndel->highCovRegNum>1))
						queryIndel->misassFlag = TRUE_MISASS;
					else
					{
						if(queryIndel->discorRatioLeft>0.2 || queryIndel->discorRatioRight>0.2)
							queryIndel->misassFlag = TRUE_MISASS;
						else if((queryIndel->discorRatioLeft>0.1 || queryIndel->discorRatioRight>0.1) && (queryIndel->lowCovRegNum>0 || queryIndel->highCovRegNum>0 || ((difAll<standDev || difAll<0.1*insertSize) && (queryIndel->discorRegNumInRatioRegion>0 || queryIndel->multiReadsRegNumInRatioRegion>0))))
							queryIndel->misassFlag = TRUE_MISASS;
						else if((queryIndel->discorRatioLeft>0.05 && queryIndel->discorRatioRight>0.05) && queryIndel->multiReadsRegNumInRatioRegion>0)
							queryIndel->misassFlag = TRUE_MISASS;
						else
							queryIndel->misassFlag = UNCERTAIN_MISASS;
					}
				}else if(queryIndel->disagreeNum>0)
				{
					if((queryIndel->disagreeRegSize/queryIndel->disagreeNum<1000) || (queryIndel->lowCovRegNum>0 || queryIndel->highCovRegNum>0))
						queryIndel->misassFlag = TRUE_MISASS;
					else
						queryIndel->misassFlag = UNCERTAIN_MISASS;
				}else
				{
					if((queryIndel->lowCovRegNum>0 || queryIndel->highCovRegNum>0) && (queryIndel->zeroCovNumInRatioRegion>0 || queryIndel->discorRegNumInRatioRegion>0 || queryIndel->multiReadsRegNumInRatioRegion>0))
						queryIndel->misassFlag = TRUE_MISASS;
					else
						queryIndel->misassFlag = UNCERTAIN_MISASS;
				}

				//printf("dif=%.4f, difFragSize=%.4f, difAll=%.4f\n", dif, difFragSize, difAll);

			}else if(queryIndel->queryIndelKind==QUERY_DEL)
			{ // query deletion
				newDifFragSize = -difFragSize;
				dif = queryIndel->difSubject - queryIndel->difQuery;
				difAll = newDifFragSize - dif;
				if(difAll<0)
					difAll = -difAll;

				if((difAll<2*standDev || difAll<0.2*insertSize) || (queryIndel->discorRatioLeft>0.1 || queryIndel->discorRatioRight>0.1))
				{
					if(queryIndel->disagreeNum>0 || queryIndel->zeroCovNum>0 || queryIndel->lowCovRegNum>1 || queryIndel->highCovRegNum>1)
						queryIndel->misassFlag = TRUE_MISASS;
					else
					{
						if(queryIndel->discorRatioLeft>0.2 || queryIndel->discorRatioRight>0.2)
							queryIndel->misassFlag = TRUE_MISASS;
						else if((queryIndel->discorRatioLeft>0.1 || queryIndel->discorRatioRight>0.1) && (queryIndel->lowCovRegNum>0 || queryIndel->highCovRegNum>0 || ((difAll<standDev || difAll<0.1*insertSize) && (queryIndel->discorRegNumInRatioRegion>0 || queryIndel->multiReadsRegNumInRatioRegion>0))))
							queryIndel->misassFlag = TRUE_MISASS;
						else if((queryIndel->discorRatioLeft>0.05 && queryIndel->discorRatioRight>0.05) && queryIndel->multiReadsRegNumInRatioRegion>0)
							queryIndel->misassFlag = TRUE_MISASS;
						else
							queryIndel->misassFlag = UNCERTAIN_MISASS;
					}
				}else if(queryIndel->disagreeNum>0)
				{
					if((queryIndel->disagreeRegSize/queryIndel->disagreeNum<1000) || (queryIndel->lowCovRegNum>0 || queryIndel->highCovRegNum>0))
						queryIndel->misassFlag = TRUE_MISASS;
					else
						queryIndel->misassFlag = UNCERTAIN_MISASS;
				}else
				{
					if((queryIndel->lowCovRegNum>0 || queryIndel->highCovRegNum>0) && (queryIndel->zeroCovNumInRatioRegion>0 || queryIndel->discorRegNumInRatioRegion>0 || queryIndel->multiReadsRegNumInRatioRegion>0))
						queryIndel->misassFlag = TRUE_MISASS;
					else
						queryIndel->misassFlag = UNCERTAIN_MISASS;
				}

				//printf("dif=%.4f, difFragSize=%.4f, difAll=%.4f\n", dif, difFragSize, difAll);

			}else if(queryIndel->queryIndelKind==QUERY_GAP)
			{ // query gap
				dif = queryIndel->difQuery - queryIndel->difSubject;
				difAll = difFragSize - dif;
				if(difAll<0)
					difAll = -difAll;

				if((difAll<2*standDev || difAll<0.2*insertSize) || (queryIndel->discorRatioLeft>0.1 || queryIndel->discorRatioRight>0.1))
				{
					if(queryIndel->disagreeNum>0 || queryIndel->zeroCovNum>0 || queryIndel->lowCovRegNum>1 || queryIndel->highCovRegNum>1)
						queryIndel->misassFlag = TRUE_MISASS;
					else
					{
						if(queryIndel->discorRatioLeft>0.2 || queryIndel->discorRatioRight>0.2)
							queryIndel->misassFlag = TRUE_MISASS;
						else if((queryIndel->discorRatioLeft>0.1 || queryIndel->discorRatioRight>0.1) && (queryIndel->lowCovRegNum>0 || queryIndel->highCovRegNum>0 || ((difAll<standDev || difAll<0.1*insertSize) && (queryIndel->discorRegNumInRatioRegion>0 || queryIndel->multiReadsRegNumInRatioRegion>0))))
							queryIndel->misassFlag = TRUE_MISASS;
						else if((queryIndel->discorRatioLeft>0.05 && queryIndel->discorRatioRight>0.05) && queryIndel->multiReadsRegNumInRatioRegion>0)
							queryIndel->misassFlag = TRUE_MISASS;
						else
							queryIndel->misassFlag = UNCERTAIN_MISASS;
					}
				}else if(queryIndel->disagreeNum>0)
				{
					if((queryIndel->disagreeRegSize/queryIndel->disagreeNum<1000) || (queryIndel->lowCovRegNum>0 || queryIndel->highCovRegNum>0))
						queryIndel->misassFlag = TRUE_MISASS;
					else
						queryIndel->misassFlag = UNCERTAIN_MISASS;
				}else
				{
					if((queryIndel->lowCovRegNum>0 || queryIndel->highCovRegNum>0) && (queryIndel->zeroCovNumInRatioRegion>0 || queryIndel->discorRegNumInRatioRegion>0 || queryIndel->multiReadsRegNumInRatioRegion>0))
						queryIndel->misassFlag = TRUE_MISASS;
					else
						queryIndel->misassFlag = UNCERTAIN_MISASS;
				}

				//printf("dif=%.4f, difFragSize=%.4f, difAll=%.4f\n", dif, difFragSize, difAll);

			}else
			{ // uncertain
				queryIndel->misassFlag = UNCERTAIN_MISASS;
			}
		}else
		{
			queryIndel->misassFlag = TRUE_MISASS;
		}

		// ######################## Debug information #######################
		//printf("disagreeNum=%d, averFragSizeLeft=%.2f, averFragSizeRight=%.2f, validPairNumLeft=%d, validPairNumRight=%d, discorRatioLeft=%.4f, discorRatioRight=%.4f, discorNumLeft=%d, discorNumRight=%d, insertSize=%.2f, standDev=%.2f, misassFlag=%d\n", queryIndel->disagreeNum, queryIndel->averFragSizeLeft, queryIndel->averFragSizeRight, queryIndel->pairNumLeft, queryIndel->pairNumRight, queryIndel->discorRatioLeft, queryIndel->discorRatioRight, queryIndel->discorNumLeft, queryIndel->discorNumRight, insertSize, standDev, queryIndel->misassFlag);
		// ######################## Debug information #######################
	}else
	{
		printf("line=%d, In %s(), invalid kind, error.\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Compute the ratios of into sub-regions for queryIndel node.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeStatisticsRatioRegionQueryIndel(queryIndel_t *queryIndel, baseCov_t *baseCovArray, readSet_t *readSet, query_t *queryItem)
{
	int32_t startRegPos, endRegPos, subRegSize, ratioRegionNum, leftPos, rightPos;
	ratioRegion_t *ratioRegionArray;

	subRegSize = 500;
	leftPos = queryIndel->leftMargin;
	rightPos = queryIndel->rightMargin;

	if(rightPos-leftPos<readSet->insertSize)
	{
		startRegPos =  leftPos - readSet->insertSize;
		if(startRegPos<1)
			startRegPos = 1;

		endRegPos = rightPos + readSet->insertSize;
		if(endRegPos>queryItem->queryLen)
			endRegPos = queryItem->queryLen;
	}else
	{
		startRegPos = leftPos;
		endRegPos = rightPos;
	}

	// divide the region into sub-regions
	if(initRatioRegQueryIndel(&ratioRegionArray, &ratioRegionNum, startRegPos, endRegPos, subRegSize)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the ratio sub-regions for indel region, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// compute disagreements, discordantNum, insert size for each sub-region
	if(fillRatioRegionArray(ratioRegionArray, ratioRegionNum, queryItem, readSet, NO)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the ratioRegion array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// compute the disagreements of ratio regions
	if(computeDisagreeNumRatioRegs(ratioRegionArray, ratioRegionNum, baseCovArray, NO)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute disagreements for the ratioRegion array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// compute region ratios
	if(computeRatios(ratioRegionArray, ratioRegionNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute ratios for the ratioRegion array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// determine the final queryIndel kind
	if(statisticsSummaryInRatioRegQueryIndel(queryIndel, ratioRegionArray, ratioRegionNum, subRegSize, queryItem)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the determine query indel mis-assembly kind, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// free memory
	if(ratioRegionNum>0)
		free(ratioRegionArray);

	return SUCCESSFUL;
}

/**
 * Compute statistics summary in ratio regions of the  queryIndel.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short statisticsSummaryInRatioRegQueryIndel(queryIndel_t *queryIndel, ratioRegion_t *ratioRegionArray, int32_t ratioRegionNum, int32_t subRegSize, query_t *queryItem)
{
	int32_t totalDisagreeNum, totalZeroCovNum, discordantRegNum, headSkipRegNum, tailSkipRegNum;

//	if(queryIndel->leftSegRow==-1)
//		headSkipRegNum = 1;
//	else
//	{
//		if(ratioRegionArray[0].startQPosRHalf<subRegSize)
//			headSkipRegNum = 1;
//		else
//			headSkipRegNum = 0;
//	}
//
//	if(queryIndel->rightSegRow==-1)
//		tailSkipRegNum = 1;
//	else
//	{
//		if(ratioRegionArray[ratioRegionNum-1].endQPosLHalf>queryItem->queryLen-subRegSize)
//			tailSkipRegNum = 1;
//		else
//			tailSkipRegNum = 0;
//	}

	headSkipRegNum = 0;
	tailSkipRegNum = 0;

	// compute total disagreements and total zero coverage
	if(computeTotalDisagreeNum(&queryIndel->disagreeNumInRatioRegion, &queryIndel->zeroCovNumInRatioRegion, ratioRegionArray, ratioRegionNum, headSkipRegNum, tailSkipRegNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the total disagreements, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get discordant region count
	if(getDiscordantRegNum(&queryIndel->discorRegNumInRatioRegion, ratioRegionArray, ratioRegionNum, headSkipRegNum, tailSkipRegNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the discordant region count, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get multiple aligned reads region count
	if(getMultiReadsRegNum(&queryIndel->multiReadsRegNumInRatioRegion, ratioRegionArray, ratioRegionNum, headSkipRegNum, tailSkipRegNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the discordant region count, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// compute the paired reads count


	return SUCCESSFUL;
}

/**
 * Determine the final queryIndel kind of unmatched query.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short determineFinalQueryIndelKindUnmatched(queryIndel_t *queryIndel)
{
	if(queryIndel->pairNumLeft==0 || queryIndel->pairNumRight==0 || queryIndel->disagreeNumInRatioRegion>=3 || queryIndel->zeroCovNumInRatioRegion>0 || queryIndel->discorRegNumInRatioRegion>0)
		queryIndel->misassFlag = TRUE_MISASS;
	else
		queryIndel->misassFlag = UNCERTAIN_MISASS;

	return SUCCESSFUL;
}
