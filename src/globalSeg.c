/*
 * globalSeg.c
 *
 *  Created on: Jul 23, 2014
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Get global align segments for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short fillGlobalAlignSeg(queryMatchInfo_t *queryMatchInfoSet, segLinkSet_t *segLinkSet)
{
	int32_t i;

	if(segLinkSet==NULL || segLinkSet->linkArray==NULL)
	{
		printf("line=%d, In %s(), cannot fill global segments, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		// ########################### Debug information ##############################
//		if(queryMatchInfoSet->queryArray[i].queryID==74 || strcmp(queryMatchInfoSet->queryArray[i].queryTitle, "ctg7180000002429")==0)
//		{
//			printf("queryID=%d, queryTitle=%s, queryLen=%d, subjectNum=%d\n", queryMatchInfoSet->queryArray[i].queryID, queryMatchInfoSet->queryArray[i].queryTitle, queryMatchInfoSet->queryArray[i].queryLen, queryMatchInfoSet->queryArray[i].querySubjectNum);
//		}
		// ########################### Debug information ##############################

		// fill global segment array
		if(fillGlobalSegSingleQuery(queryMatchInfoSet->queryArray+i, queryMatchInfoSet->subjectArray, segLinkSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot fill global segments, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// trim alignment information of unmatched query ends
//		if(trimAlignInfoQueryEnds(queryMatchInfoSet->queryArray+i, queryMatchInfoSet->subjectArray)==FAILED)
//		{
//			printf("line=%d, In %s(), cannot trim the alignment information of unmatched query ends, error!\n", __LINE__, __func__);
//			return FAILED;
//		}

		// determine global matchKind
		if(determineGlobalMatchKindSingleQuery(queryMatchInfoSet->queryArray+i)==FAILED)
		{
			printf("line=%d, In %s(), cannot determine the global match kind for single query, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Fill global align segment for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short fillGlobalSegSingleQuery(query_t *queryItem, subject_t *subjectArray, segLinkSet_t *segLinkSet)
{
	int32_t querySubjectNum, bestSubjectID, bestSegNum, globalSegNum, maxGlobalSegNum;
	querySubject_t *querySubjectArray;
	validSegment_t *bestSegArray;
	globalValidSeg_t *wholeGlobalSegArray, *globalSegArrayBuf;


	if(queryItem->bestMatchRow==-1)
		return SUCCESSFUL;

	querySubjectArray = queryItem->querySubArray;
	querySubjectNum = queryItem->querySubjectNum;

	bestSubjectID = querySubjectArray[queryItem->bestMatchRow].subjectID;
	bestSegArray = querySubjectArray[queryItem->bestMatchRow].validSegArray;
	bestSegNum = querySubjectArray[queryItem->bestMatchRow].validSegmentNum;

	if(querySubjectArray[queryItem->bestMatchRow].matchKind==PERFECT_MATCH_KIND)
	{
		if(generateGlobalSegArrayFromValidSegArray(queryItem, bestSubjectID, bestSegArray, bestSegNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate global segment array, error!\n", __LINE__, __func__);
			return FAILED;
		}
		return SUCCESSFUL;
	}

	// fill the full segments
	if(initWholeGlobalSegArray(&wholeGlobalSegArray, &maxGlobalSegNum, queryItem)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize global segment array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate result buffer array
	globalSegArrayBuf = (globalValidSeg_t*) calloc (maxGlobalSegNum+bestSegNum, sizeof(globalValidSeg_t));
	if(globalSegArrayBuf==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the result buffer
	if(getGlobalSegSingleQuery(globalSegArrayBuf, &globalSegNum, wholeGlobalSegArray, maxGlobalSegNum, segLinkSet, queryItem, subjectArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot get global segments, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// remove the redundant items
	if(removeRedundantGlobalSegSingleQuery(globalSegArrayBuf, &globalSegNum, queryItem)==FAILED)
	{
		printf("line=%d, In %s(), cannot get global segments, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// adjust incorrect order of the segments
	if(adjustGlobalSegOrderSingleQuery(globalSegArrayBuf, globalSegNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust the order of align segments, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// generate global segment array from buffer array
	if(generateGlobalSegArraySingleQuery(queryItem, globalSegArrayBuf, globalSegNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate global segment array for single query, error!\n", __LINE__, __func__);
		return FAILED;
	}

	free(wholeGlobalSegArray);
	free(globalSegArrayBuf);

	return SUCCESSFUL;
}

/**
 * Generate global align segment array from valid segment array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short generateGlobalSegArrayFromValidSegArray(query_t *queryItem, int32_t bestSubjectID, validSegment_t *validSegArray, int32_t validSegNum)
{
	int32_t i;

	queryItem->globalValidSegNum = validSegNum;
	queryItem->globalValidSegArray = (globalValidSeg_t*) calloc (validSegNum, sizeof(globalValidSeg_t));
	if(queryItem->globalValidSegArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<validSegNum; i++)
	{
		queryItem->globalValidSegArray[i].subjectID = bestSubjectID;
		queryItem->globalValidSegArray[i].startSubPos = validSegArray[i].startSubPos;
		queryItem->globalValidSegArray[i].endSubPos = validSegArray[i].endSubPos;
		queryItem->globalValidSegArray[i].startQueryPos = validSegArray[i].startQueryPos;
		queryItem->globalValidSegArray[i].endQueryPos = validSegArray[i].endQueryPos;
		queryItem->globalValidSegArray[i].strand = validSegArray[i].strand;
		queryItem->globalValidSegArray[i].matchLen = validSegArray[i].matchLen;
		queryItem->globalValidSegArray[i].totalMatchLen = validSegArray[i].totalMatchLen;
		queryItem->globalValidSegArray[i].gapNum = validSegArray[i].gapNum;
		queryItem->globalValidSegArray[i].matchPercent = validSegArray[i].matchPercent;
	}

	return SUCCESSFUL;
}

/**
 * Initialize full global align segment array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short initWholeGlobalSegArray(globalValidSeg_t **wholeGlobalSegArray, int32_t *maxGlobalSegNum, query_t *queryItem)
{
	int32_t i, j, querySubjectNum, validSegNum, globalSegNum, subjectID;
	querySubject_t *querySubjectArray;
	validSegment_t *validSegArray;

	querySubjectArray = queryItem->querySubArray;
	querySubjectNum = queryItem->querySubjectNum;

	// fill the full segments
	*maxGlobalSegNum = 0;
	for(i=0; i<querySubjectNum; i++)
		if(i!=queryItem->bestMatchRow)
			(*maxGlobalSegNum) += querySubjectArray[i].validSegmentNum;

	*wholeGlobalSegArray = (globalValidSeg_t*) calloc(*maxGlobalSegNum, sizeof(globalValidSeg_t));
	if((*wholeGlobalSegArray)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	globalSegNum = 0;
	for(i=0; i<querySubjectNum; i++)
	{
		if(i!=queryItem->bestMatchRow)
		{
			validSegArray = querySubjectArray[i].validSegArray;
			validSegNum = querySubjectArray[i].validSegmentNum;
			subjectID = querySubjectArray[i].subjectID;
			for(j=0; j<validSegNum; j++)
			{
				(*wholeGlobalSegArray)[globalSegNum].subjectID = subjectID;
				(*wholeGlobalSegArray)[globalSegNum].startSubPos = validSegArray[j].startSubPos;
				(*wholeGlobalSegArray)[globalSegNum].endSubPos = validSegArray[j].endSubPos;
				(*wholeGlobalSegArray)[globalSegNum].startQueryPos = validSegArray[j].startQueryPos;
				(*wholeGlobalSegArray)[globalSegNum].endQueryPos = validSegArray[j].endQueryPos;
				(*wholeGlobalSegArray)[globalSegNum].strand = validSegArray[j].strand;
				(*wholeGlobalSegArray)[globalSegNum].matchLen = validSegArray[j].matchLen;
				(*wholeGlobalSegArray)[globalSegNum].totalMatchLen = validSegArray[j].totalMatchLen;
				(*wholeGlobalSegArray)[globalSegNum].gapNum = validSegArray[j].gapNum;
				(*wholeGlobalSegArray)[globalSegNum].matchPercent = validSegArray[j].matchPercent;
				globalSegNum ++;
			}
		}
	}

	// sort the whole global segments according to their lengths
	if(sortGlobalSegArray(*wholeGlobalSegArray, *maxGlobalSegNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot sort global segments, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Sort the global align segments according to their lengths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short sortGlobalSegArray(globalValidSeg_t *wholeGlobalSegArray, int32_t arraySize)
{
	if(arraySize<50)
	{ // selection sort
		if(selectionSortGlobalSegArray(wholeGlobalSegArray, arraySize)==FAILED)
		{
			printf("line=%d, In %s(), cannot sort global segments by selection sort, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{ // radix sort
		if(radixSortGlobalSegArray(wholeGlobalSegArray, arraySize)==FAILED)
		{
			printf("line=%d, In %s(), cannot sort global segments by selection sort, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Selection sort the global align segments according to their lengths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short selectionSortGlobalSegArray(globalValidSeg_t *wholeGlobalSegArray, int32_t arraySize)
{
	int32_t i, j, maxLen, maxRow;
	globalValidSeg_t tmp;

	for(i=0; i<arraySize-1; i++)
	{
		maxLen = wholeGlobalSegArray[i].matchLen;
		maxRow = -1;
		for(j=i+1; j<arraySize; j++)
		{
			if(maxLen<wholeGlobalSegArray[j].matchLen)
			{
				maxLen = wholeGlobalSegArray[j].matchLen;
				maxRow = j;
			}
		}

		if(maxRow>=0)
		{
			if(memcpy(&tmp, wholeGlobalSegArray+i, sizeof(globalValidSeg_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(memcpy(wholeGlobalSegArray+i, wholeGlobalSegArray+maxRow, sizeof(globalValidSeg_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(memcpy(wholeGlobalSegArray+maxRow, &tmp, sizeof(globalValidSeg_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Radix sort the global align segments according to their lengths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short radixSortGlobalSegArray(globalValidSeg_t *segArray, int32_t itemNum)
{
	struct partNode
	{
		int curItemNum;
		int totalItemNum;
		int firstRow;
	};

	int64_t i, step, total;
	globalValidSeg_t *segArrayBuf, *data, *buf;
	struct partNode *part;
	int32_t partArrSize, stepBits, maxStepLen;
	uint64_t bitMask, hashcode, firstRow, curItemNum;

	segArrayBuf = (globalValidSeg_t*) calloc (itemNum, sizeof(globalValidSeg_t));
	if(segArrayBuf==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

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
			buf = segArray;
			data = segArrayBuf;
		}else
		{
			data = segArray;
			buf = segArrayBuf;
		}

		// count the number
		if(memset(part, 0, partArrSize * sizeof(struct partNode))==NULL)
		{
			printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		for(i=0; i<itemNum; i++)
		{
			part[ bitMask - ((data[i].matchLen >> step) & bitMask) ].totalItemNum ++;  // from big to small
			//part[ (data[i].matchLen >> step) & bitMask ].totalItemNum ++;  // from small to big
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
			hashcode = bitMask - ((data[i].matchLen >> step) & bitMask);  // from big to small
			//hashcode = (data[i].matchLen >> step) & bitMask;  // from small to big
			firstRow = part[hashcode].firstRow;
			curItemNum = part[hashcode].curItemNum;
			if(memcpy(buf+firstRow+curItemNum, data+i, sizeof(globalValidSeg_t))==NULL)
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
	free(segArrayBuf);
	segArrayBuf = NULL;

	return SUCCESSFUL;
}

/**
 * Get global align segments for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short getGlobalSegSingleQuery(globalValidSeg_t *globalSegArrayBuf, int32_t *globalSegNum, globalValidSeg_t *wholeGlobalSegArray, int32_t maxGlobalSegNum, segLinkSet_t *segLinkSet, query_t *queryItem, subject_t *subjectArray)
{
	int32_t i, startQPos, endQPos, bestSubjectID, bestSegNum, queryLen;
	validSegment_t *bestSegArray;
	int8_t *usedFlagArray;

	usedFlagArray = (int8_t*) calloc(maxGlobalSegNum, sizeof(int8_t));
	if(usedFlagArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	for(i=0; i<maxGlobalSegNum; i++) usedFlagArray[i] = NO;

	queryLen = queryItem->queryLen;
	bestSubjectID = queryItem->querySubArray[queryItem->bestMatchRow].subjectID;
	bestSegArray = queryItem->querySubArray[queryItem->bestMatchRow].validSegArray;
	bestSegNum = queryItem->querySubArray[queryItem->bestMatchRow].validSegmentNum;

	*globalSegNum = 0;
	for(i=0; i<bestSegNum; i++)
	{
		if(i==0 && bestSegArray[i].startQueryPos>END_IGNORE_LEN)
		{ // check head
			startQPos = 1;
			endQPos = bestSegArray[i].endQueryPos - 1;
			if(fillLinkArrayGivenReg(segLinkSet, startQPos, endQPos, wholeGlobalSegArray, usedFlagArray, maxGlobalSegNum, subjectArray)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the link array given region, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// add new items
			if(copyNewItemsFromLinkArray(globalSegArrayBuf, globalSegNum, segLinkSet, wholeGlobalSegArray)==FAILED)
			{
				printf("line=%d, In %s(), cannot copy items of link array given region, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		// add the item
		if(addValidSegItemToGlobalArray(globalSegArrayBuf, globalSegNum, bestSegArray+i, bestSubjectID)==FAILED)
		{
			printf("line=%d, In %s(), cannot copy items of link array given region, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(i==bestSegNum-1)
		{ // check tail
			if(bestSegArray[i].endQueryPos<queryLen-END_IGNORE_LEN)
			{
				startQPos = bestSegArray[i].endQueryPos + 1;
				endQPos = queryLen;
				if(fillLinkArrayGivenReg(segLinkSet, startQPos, endQPos, wholeGlobalSegArray, usedFlagArray, maxGlobalSegNum, subjectArray)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the link array given region, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// add new items
				if(copyNewItemsFromLinkArray(globalSegArrayBuf, globalSegNum, segLinkSet, wholeGlobalSegArray)==FAILED)
				{
					printf("line=%d, In %s(), cannot copy items of link array given region, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}else
		{ // check middle
			if(bestSegArray[i+1].startQueryPos-bestSegArray[i].endQueryPos>varyEndLenThres)
			{
				startQPos = bestSegArray[i].endQueryPos + 1;
				endQPos = bestSegArray[i+1].startQueryPos - 1;
				if(fillLinkArrayGivenReg(segLinkSet, startQPos, endQPos, wholeGlobalSegArray, usedFlagArray, maxGlobalSegNum, subjectArray)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the link array given region, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// add new items
				if(copyNewItemsFromLinkArray(globalSegArrayBuf, globalSegNum, segLinkSet, wholeGlobalSegArray)==FAILED)
				{
					printf("line=%d, In %s(), cannot copy items of link array given region, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}

	free(usedFlagArray);

	return SUCCESSFUL;
}

/**
 * Fill link array for global align segments for given region.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short fillLinkArrayGivenReg(segLinkSet_t *segLinkSet, int32_t startQPos, int32_t endQPos, globalValidSeg_t *globalSegArray, int8_t *usedFlagArray, int32_t segArraySize, subject_t *subjectArray)
{
	int32_t i, startSegRow, adjacentRow, selectionRound, addOrder;

	segLinkSet->itemNum = 0;
	segLinkSet->headRow = segLinkSet->tailRow = -1;

	// get first unused item
	startSegRow = -1;
	for(i=0; i<segArraySize; i++)
	{
		if(usedFlagArray[i]==NO)
		{
			if(isValidGapSeg(globalSegArray+i, startQPos, endQPos)==YES)
			{
				startSegRow = i;
				break;
			}
		}
	}

	if(startSegRow>=0)
	{
		if(addNewItemSegLinkArray(segLinkSet, startSegRow, 0, 1)==FAILED)
		{
			printf("line=%d, In %s(), cannot add item to link array, error!\n", __LINE__, __func__);
			return FAILED;
		}

		usedFlagArray[startSegRow] = YES;
		adjacentRow = -1;

		if(globalSegArray[startSegRow].startQueryPos<=startQPos+varyEndLenThres && globalSegArray[startSegRow].endQueryPos>=startQPos-varyEndLenThres)
		{
			;
		}
		else
		{
			selectionRound = 1;
			while(selectionRound<=2)
			{
				// get adjacent row
				if(getAdjacentRowGlobalSeg(&adjacentRow, startSegRow, startQPos, endQPos, globalSegArray, usedFlagArray, segArraySize, subjectArray, selectionRound)==FAILED)
				{
					printf("line=%d, In %s(), cannot get adjacent row of global segments, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(selectionRound==1)
				{
					addOrder = segLinkSet->linkArray[segLinkSet->tailRow].addedOrder;
					if(adjacentRow>=0)
					{
						if(addNewItemSegLinkArray(segLinkSet, adjacentRow, addOrder, 1)==FAILED)
						{
							printf("line=%d, In %s(), cannot add item to link array, error!\n", __LINE__, __func__);
							return FAILED;
						}

						usedFlagArray[adjacentRow] = YES;
						startSegRow = adjacentRow;
					}else
					{
						startSegRow = segLinkSet->linkArray[segLinkSet->headRow].arrRow;
						selectionRound ++;
					}
				}else if(selectionRound==2)
				{
					addOrder =  segLinkSet->linkArray[segLinkSet->headRow].addedOrder;
					if(adjacentRow>=0)
					{
						if(addNewItemSegLinkArray(segLinkSet, adjacentRow, addOrder, 2)==FAILED)
						{
							printf("line=%d, In %s(), cannot add item to link array, error!\n", __LINE__, __func__);
							return FAILED;
						}

						usedFlagArray[adjacentRow] = YES;
						startSegRow = adjacentRow;
					}else
					{
						selectionRound ++;
					}
				}else
				{
					printf("line=%d, In %s(), invalid linkRound=%d, error!\n", __LINE__, __func__, selectionRound);
					return FAILED;
				}
			}
		}
	}


	return SUCCESSFUL;
}

/**
 * Check whether the global align segment is valid for given region.
 *  @return:
 *  	If valid, return YES; otherwise return NO.
 */
int32_t isValidGapSeg(globalValidSeg_t *globalSegNode, int32_t startQPos, int32_t endQPos)
{
	int32_t validFlag;
	if(globalSegNode->startQueryPos>=startQPos-varyEndLenThres && globalSegNode->endQueryPos<=endQPos+varyEndLenThres)
	{
		validFlag = YES;
	}else if(globalSegNode->startQueryPos<startQPos-varyEndLenThres && globalSegNode->endQueryPos>=startQPos+50)
	{
		validFlag = YES;
	}else if(globalSegNode->startQueryPos<=endQPos-50 && globalSegNode->endQueryPos>endQPos+varyEndLenThres)
	{
		validFlag = YES;
	}else
	{
		validFlag = NO;
	}

	return validFlag;
}

/**
 * Copy items in link array for global align segments for given region.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short copyNewItemsFromLinkArray(globalValidSeg_t *globalSegArrayBuf, int32_t *globalSegNum, segLinkSet_t *segLinkSet, globalValidSeg_t *wholeGlobalSegArray)
{
	int32_t itemRowLinkArray, dataRow;

	itemRowLinkArray = segLinkSet->headRow;
	while(itemRowLinkArray!=-1)
	{
		dataRow = segLinkSet->linkArray[itemRowLinkArray].arrRow;
		if(memcpy(globalSegArrayBuf+(*globalSegNum), wholeGlobalSegArray+dataRow, sizeof(globalValidSeg_t))==NULL)
		{
			printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		(*globalSegNum) ++;

		itemRowLinkArray = segLinkSet->linkArray[itemRowLinkArray].next;
	}

	return SUCCESSFUL;
}

/**
 * Copy items in valid segment array to global align segment array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short addValidSegItemToGlobalArray(globalValidSeg_t *globalValidSegArray, int32_t *globalSegNum, validSegment_t *validSegNode, int32_t subjectID)
{
	globalValidSegArray[*globalSegNum].subjectID = subjectID;
	globalValidSegArray[*globalSegNum].startSubPos = validSegNode->startSubPos;
	globalValidSegArray[*globalSegNum].endSubPos = validSegNode->endSubPos;
	globalValidSegArray[*globalSegNum].startQueryPos = validSegNode->startQueryPos;
	globalValidSegArray[*globalSegNum].endQueryPos = validSegNode->endQueryPos;
	globalValidSegArray[*globalSegNum].strand = validSegNode->strand;
	globalValidSegArray[*globalSegNum].matchLen = validSegNode->matchLen;
	globalValidSegArray[*globalSegNum].totalMatchLen = validSegNode->totalMatchLen;
	globalValidSegArray[*globalSegNum].gapNum = validSegNode->gapNum;
	globalValidSegArray[*globalSegNum].matchPercent = validSegNode->matchPercent;
	(*globalSegNum) ++;

	return SUCCESSFUL;
}

/**
 * Copy items in valid segment array to global align segment array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short getAdjacentRowGlobalSeg(int32_t *adjacentRow, int32_t startSegRow, int32_t startQPos, int32_t endQPos, globalValidSeg_t *globalSegArray, int8_t *pUsedArray, int32_t segArraySize, subject_t *subjectArray, int32_t selectionRound)
{
	int32_t i, gapLen, subjectID, subjectLen;
	int64_t distance, distanceQuery, distanceSubject, minDistQuery, minDistSubject;

	if((selectionRound==1 && globalSegArray[startSegRow].endQueryPos>=endQPos-50) || (selectionRound==2 && globalSegArray[startSegRow].startQueryPos<startQPos+50))
	{
		*adjacentRow = -1;
		return SUCCESSFUL;
	}

	gapLen = endQPos - startQPos + 1;
	*adjacentRow = -1;
	minDistQuery = minDistSubject = INT_MAX;
	for(i=0; i<segArraySize; i++)
	{
		subjectID = globalSegArray[i].subjectID;
		subjectLen = subjectArray[subjectID-1].subjectLen;

		if(pUsedArray[i]==NO && i!=startSegRow && (isValidGapSeg(globalSegArray+i, startQPos, endQPos)==YES) && (globalSegArray[i].totalMatchLen>minAlignedSegLenThres && gapLen*(1-matchPercentThres)<globalSegArray[i].totalMatchLen))
		{
			if(selectionRound==1)
			{ // the first selection round
				distanceQuery = globalSegArray[i].startQueryPos - globalSegArray[startSegRow].endQueryPos;
				if(distanceQuery<0)
					distanceQuery = -distanceQuery;
				if(distanceQuery<=varyEndLenThres)
				{ // distanceQuery <= varyEndLenThres
					if(globalSegArray[i].subjectID==globalSegArray[startSegRow].subjectID)
					{ // same subject
						if(globalSegArray[i].strand==globalSegArray[startSegRow].strand)
						{ // same strand
							if(globalSegArray[i].strand==PLUS_STRAND)
							{ // plus strand
								if(subjectArray[subjectID-1].circularFlag==YES)
								{ // circular subject
									if(globalSegArray[startSegRow].endSubPos>=subjectLen-varyEndLenThres && globalSegArray[startSegRow].endSubPos<=subjectLen)
									{ // subjectLen-varyEndLenThres <= endSubPos <= subjectLen
										if(globalSegArray[i].startSubPos>=1 && globalSegArray[i].startSubPos<=varyEndLenThres)
										{ // 1 <= [i].startSubPos <= varyEndLenThres
											distanceSubject = globalSegArray[i].startSubPos + subjectLen - globalSegArray[startSegRow].endSubPos;
										}else if(globalSegArray[i].startSubPos>=subjectLen-varyEndLenThres && globalSegArray[i].startSubPos<=subjectLen)
											distanceSubject = globalSegArray[i].startSubPos - globalSegArray[startSegRow].endSubPos;

										if(distanceSubject<0)
											distanceSubject = -distanceSubject;
										if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
										{
											*adjacentRow = i;
											minDistQuery = distanceQuery;
											minDistSubject = distanceSubject;
										}
									}else
									{
										distanceSubject = globalSegArray[i].startSubPos - globalSegArray[startSegRow].endSubPos;
										if(distanceSubject<0)
											distanceSubject = -distanceSubject;
										if(distanceSubject<=varyEndLenThres)
										{ // distanceSubject <= varyEndLenThres
											if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
											{
												*adjacentRow = i;
												minDistQuery = distanceQuery;
												minDistSubject = distanceSubject;
											}
										}
									}
								}else
								{ // linear subject
									distanceSubject = globalSegArray[i].startSubPos - globalSegArray[startSegRow].endSubPos;
									if(distanceSubject<0)
										distanceSubject = -distanceSubject;
									if(distanceSubject<=varyEndLenThres)
									{ // distanceSubject <= varyEndLenThres
										if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
										{
											*adjacentRow = i;
											minDistQuery = distanceQuery;
											minDistSubject = distanceSubject;
										}
									}
								}
							}else
							{ // minus strand
								if(subjectArray[subjectID-1].circularFlag==YES)
								{ // circular subject
									if(globalSegArray[startSegRow].startSubPos>=1 && globalSegArray[startSegRow].startSubPos<=varyEndLenThres)
									{ // 1 <= [startRow].startSubPos <= varyEndLenThres
										if(globalSegArray[i].startSubPos>=subjectLen-varyEndLenThres && globalSegArray[i].startSubPos<=subjectLen)
										{ // subjectLen-varyEndLenThres <= [i].startSubPos <= subjectLen
											distanceSubject = globalSegArray[startSegRow].endSubPos + subjectLen - globalSegArray[i].startSubPos;
										}else if(globalSegArray[i].startSubPos>=1 && globalSegArray[i].startSubPos<=varyEndLenThres)
											distanceSubject = globalSegArray[startSegRow].endSubPos - globalSegArray[i].startSubPos;

										if(distanceSubject<0)
											distanceSubject = -distanceSubject;
										if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
										{
											*adjacentRow = i;
											minDistQuery = distanceQuery;
											minDistSubject = distanceSubject;
										}
									}else
									{
										distanceSubject = globalSegArray[startSegRow].endSubPos - globalSegArray[i].startSubPos;
										if(distanceSubject<0)
											distanceSubject = -distanceSubject;
										if(distanceSubject<=varyEndLenThres)
										{ // distanceSubject <= varyEndLenThres
											if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
											{
												*adjacentRow = i;
												minDistQuery = distanceQuery;
												minDistSubject = distanceSubject;
											}
										}
									}
								}else
								{ // linear subject
									distanceSubject = globalSegArray[startSegRow].endSubPos - globalSegArray[i].startSubPos;
									if(distanceSubject<0)
										distanceSubject = -distanceSubject;
									if(distanceSubject<=varyEndLenThres)
									{ // distanceSubject <= varyEndLenThres
										if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
										{
											*adjacentRow = i;
											minDistQuery = distanceQuery;
											minDistSubject = distanceSubject;
										}
									}
								}
							}
						}
					}
				}
			}else
			{ // the second selection round
				distanceQuery = globalSegArray[i].endQueryPos - globalSegArray[startSegRow].startQueryPos;
				if(distanceQuery<0)
					distanceQuery = -distanceQuery;
				if(distanceQuery<=varyEndLenThres)
				{
					if(globalSegArray[i].subjectID==globalSegArray[startSegRow].subjectID)
					{ // same subject
						if(globalSegArray[i].strand==globalSegArray[startSegRow].strand)
						{ // same strand
							if(globalSegArray[i].strand==PLUS_STRAND)
							{ // plus strand
								if(subjectArray[subjectID-1].circularFlag==YES)
								{ // circular subject
									if(globalSegArray[startSegRow].startSubPos>=1 && globalSegArray[startSegRow].startSubPos<=varyEndLenThres)
									{ // 1 <= [startRow].startSubPos <= varyEndLenThres
										if(globalSegArray[i].endSubPos>=subjectLen-varyEndLenThres && globalSegArray[i].endSubPos<=subjectLen)
										{ // subjectLen-varyEndLenThres <= [i].endSubPos <= subjectLen
											distanceSubject = globalSegArray[startSegRow].startSubPos + subjectLen - globalSegArray[i].endSubPos;
										}else if(globalSegArray[i].endSubPos>=1 && globalSegArray[i].endSubPos<=varyEndLenThres)
											distanceSubject = globalSegArray[startSegRow].startSubPos - globalSegArray[i].endSubPos;

										if(distanceSubject<0)
											distanceSubject = -distanceSubject;
										if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
										{
											*adjacentRow = i;
											minDistQuery = distanceQuery;
											minDistSubject = distanceSubject;
										}
									}else
									{
										distanceSubject = globalSegArray[startSegRow].startSubPos - globalSegArray[i].endSubPos;
										if(distanceSubject<0)
											distanceSubject = -distanceSubject;
										if(distanceSubject<=varyEndLenThres)
										{ // distanceSubject <= varyEndLenThres
											if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
											{
												*adjacentRow = i;
												minDistQuery = distanceQuery;
												minDistSubject = distanceSubject;
											}
										}
									}
								}else
								{ // linear subject
									distanceSubject = globalSegArray[startSegRow].startSubPos - globalSegArray[i].endSubPos;
									if(distanceSubject<0)
										distanceSubject = -distanceSubject;
									if(distanceSubject<=+varyEndLenThres)
									{ // distanceSubject <= varyEndLenThres
										if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
										{
											*adjacentRow = i;
											minDistQuery = distanceQuery;
											minDistSubject = distanceSubject;
										}
									}
								}
							}else
							{ // minus strand
								if(subjectArray[subjectID-1].circularFlag==YES)
								{ // circular subject
									if(globalSegArray[startSegRow].startSubPos>=subjectLen-varyEndLenThres && globalSegArray[startSegRow].startSubPos<=subjectLen)
									{ // subjectLen-varyEndLenThres <= [startRow].startSubPos <= subjectLen
										if(globalSegArray[i].endSubPos>=1 && globalSegArray[i].endSubPos<=varyEndLenThres)
										{ // 1 <= [i].endSubPos <= varyEndLenThres
											distanceSubject = globalSegArray[i].endSubPos + subjectLen - globalSegArray[startSegRow].startSubPos;
										}else if(globalSegArray[i].endSubPos>=subjectLen-varyEndLenThres && globalSegArray[i].endSubPos<=subjectLen)
											distanceSubject = globalSegArray[i].endSubPos - globalSegArray[startSegRow].startSubPos;

										if(distanceSubject<0)
											distanceSubject = -distanceSubject;
										if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
										{
											*adjacentRow = i;
											minDistQuery = distanceQuery;
											minDistSubject = distanceSubject;
										}
									}else
									{
										distanceSubject = globalSegArray[i].endSubPos - globalSegArray[startSegRow].startSubPos;
										if(distanceSubject<0)
											distanceSubject = -distanceSubject;
										if(distanceSubject<=varyEndLenThres)
										{ // distanceSubject <= varyEndLenThres
											if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
											{
												*adjacentRow = i;
												minDistQuery = distanceQuery;
												minDistSubject = distanceSubject;
											}
										}
									}
								}else
								{ // linear subject
									distanceSubject = globalSegArray[i].endSubPos - globalSegArray[startSegRow].startSubPos;
									if(distanceSubject<0)
										distanceSubject = -distanceSubject;
									if(distanceSubject<=varyEndLenThres)
									{ // distanceSubject <= varyEndLenThres
										if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
										{
											*adjacentRow = i;
											minDistQuery = distanceQuery;
											minDistSubject = distanceSubject;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	// get the adjacent aligned segment according to same subject, strand, queryPos
	if((*adjacentRow)==-1)
	{
		for(i=0; i<segArraySize; i++)
		{
			if(pUsedArray[i]==NO && i!=startSegRow && (isValidGapSeg(globalSegArray+i, startQPos, endQPos)==YES) && globalSegArray[i].totalMatchLen>minAlignedSegLenThres)
			{
				if(selectionRound==1)
				{ // the first selection round
					if(globalSegArray[i].startQueryPos >= globalSegArray[startSegRow].endQueryPos-varyEndLenThres && globalSegArray[i].startQueryPos <= globalSegArray[startSegRow].endQueryPos+varyEndLenThres)
					{ // [startRow].endQueryPos-varyEndLenThres <= [i].startQueryPos <= [startRow].endQueryPos+varyEndLenThres
						if(globalSegArray[i].subjectID==globalSegArray[startSegRow].subjectID && globalSegArray[i].strand==globalSegArray[startSegRow].strand)
						{ // same strand
							*adjacentRow = i;
							break;
						}
					}
				}else
				{ // the second selection round
					if(globalSegArray[i].endQueryPos >= globalSegArray[startSegRow].startQueryPos-varyEndLenThres && globalSegArray[i].endQueryPos <= globalSegArray[startSegRow].startQueryPos+varyEndLenThres)
					{
						if(globalSegArray[i].subjectID==globalSegArray[startSegRow].subjectID && globalSegArray[i].strand==globalSegArray[startSegRow].strand)
						{ // same strand
							*adjacentRow = i;
							break;
						}
					}
				}
			}
		}
	}

	// get the adjacent aligned segment according to queryPos
	if((*adjacentRow)==-1)
	{
		for(i=0; i<segArraySize; i++)
		{
			if(pUsedArray[i]==NO && i!=startSegRow && (isValidGapSeg(globalSegArray+i, startQPos, endQPos)==YES) && globalSegArray[i].totalMatchLen>minAlignedSegLenThres)
			{
				if(selectionRound==1)
				{ // the first selection round
					if(globalSegArray[i].startQueryPos >= globalSegArray[startSegRow].endQueryPos-varyEndLenThres && globalSegArray[i].startQueryPos <= globalSegArray[startSegRow].endQueryPos+varyEndLenThres)
					{ // [startRow].endQueryPos-varyEndLenThres <= [i].startQueryPos <= [startRow].endQueryPos+varyEndLenThres
						*adjacentRow = i;
						break;
					}
				}else
				{ // the second selection round
					if(globalSegArray[i].endQueryPos >= globalSegArray[startSegRow].startQueryPos-varyEndLenThres && globalSegArray[i].endQueryPos <= globalSegArray[startSegRow].startQueryPos+varyEndLenThres)
					{
						*adjacentRow = i;
						break;
					}
				}
			}
		}
	}

	// get the most adjacent aligned segment according to queryPos
	if((*adjacentRow)==-1)
	{
		minDistQuery = INT64_MAX;
		for(i=0; i<segArraySize; i++)
		{
			if(pUsedArray[i]==NO && i!=startSegRow && (isValidGapSeg(globalSegArray+i, startQPos, endQPos)==YES) && globalSegArray[i].totalMatchLen>minAlignedSegLenThres)
			{
				if(selectionRound==1)
				{ // the first selection round
					if(globalSegArray[i].endQueryPos>globalSegArray[startSegRow].endQueryPos)
					{
						distance = globalSegArray[i].startQueryPos - globalSegArray[startSegRow].endQueryPos;
						if(distance<0)
							distance = -distance;

						if(distance<minDistQuery)
						{
							*adjacentRow = i;
							minDistQuery = distance;
						}
					}
				}else
				{ // the second selection round
					if(globalSegArray[i].startQueryPos<globalSegArray[startSegRow].startQueryPos)
					{
						distance = globalSegArray[startSegRow].startQueryPos - globalSegArray[i].endQueryPos;
						if(distance<0)
							distance = -distance;

						if(distance<minDistQuery)
						{
							*adjacentRow = i;
							minDistQuery = distance;
						}
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Remove redundant global align segment items.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short removeRedundantGlobalSegSingleQuery(globalValidSeg_t *globalSegArray, int32_t *globalSegNum, query_t *queryItem)
{
	int32_t i, j, itemNum, startItemRow, endItemRow, rowSeg1, rowSeg2, newItemNum;
	int64_t startSegPos, endSegPos, minPos, maxPos;
	int8_t *reduFlagArray;

	itemNum = *globalSegNum;
	if(itemNum>=2)
	{
		reduFlagArray = (int8_t *) calloc (itemNum, sizeof(int8_t));
		if(reduFlagArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		for(i=0; i<itemNum; i++) reduFlagArray[i] = NO;

		// get start and end alignment position, and item numbers
		startSegPos = endSegPos = -1;
		startItemRow = endItemRow = -1;
		for(i=1; i<itemNum; i++)
		{
			if(startSegPos==-1)
			{
				if(globalSegArray[i-1].endQueryPos>globalSegArray[i].startQueryPos)
				{
					startSegPos = globalSegArray[i].startQueryPos;
					endSegPos = globalSegArray[i-1].endQueryPos;
					startItemRow = i - 1;

					if(i==itemNum-1)
						endItemRow = i;
					else if(globalSegArray[i+1].startQueryPos>endSegPos)
					{
						if(globalSegArray[i+1].startSubPos!=globalSegArray[i].startSubPos && globalSegArray[i+1].endSubPos!=globalSegArray[i].endSubPos)
							endItemRow = i;
					}
				}
			}else
			{
				if(globalSegArray[i].startQueryPos<startSegPos)
					startSegPos = globalSegArray[i].startQueryPos;

				if(globalSegArray[i].endQueryPos>endSegPos)
				{
					if(globalSegArray[i].endQueryPos<endSegPos+200)
						endSegPos = globalSegArray[i].endQueryPos;
					else if(globalSegArray[i-1].endQueryPos>endSegPos)
						endSegPos = globalSegArray[i-1].endQueryPos;
				}

				if(i==itemNum-1)
					endItemRow = i;
				else if(globalSegArray[i+1].startQueryPos>endSegPos)
				{
					if(globalSegArray[i+1].startSubPos!=globalSegArray[i].startSubPos && globalSegArray[i+1].endSubPos!=globalSegArray[i].endSubPos)
						endItemRow = i;
				}
			}

			if(startItemRow>=0 && endItemRow>=0)
			{
				// get valid two large segments
				minPos = INT_MAX;
				maxPos = INT_MIN;
				for(j=startItemRow; j<=endItemRow; j++)
				{
					if(globalSegArray[j].startQueryPos<minPos)
						minPos = globalSegArray[j].startQueryPos;
					if(globalSegArray[j].endQueryPos>maxPos)
						maxPos = globalSegArray[j].endQueryPos;
				}

				rowSeg1 = rowSeg2 = -1;
				for(j=startItemRow; j<=endItemRow; j++)
				{
					if(globalSegArray[j].startQueryPos<startSegPos)
						rowSeg1 = j;
					else if(rowSeg1==-1 && globalSegArray[j].startQueryPos==startSegPos)
						rowSeg1 = j;
					else if(globalSegArray[j].startQueryPos==minPos && globalSegArray[j].endQueryPos>endSegPos)
						rowSeg1 = j;

					if(globalSegArray[j].endQueryPos>endSegPos)
						rowSeg2 = j;
					else if(rowSeg2==-1 && globalSegArray[j].endQueryPos==endSegPos)
						rowSeg2 = j;
					else if(globalSegArray[j].endQueryPos==maxPos && globalSegArray[j].startQueryPos<startSegPos)
						rowSeg2 = j;
				}

				if(rowSeg1>=0 && rowSeg2>=0)
				{
					if(rowSeg1!=rowSeg2 && globalSegArray[rowSeg1].startSubPos==globalSegArray[rowSeg2].startSubPos)
					{
						if(globalSegArray[rowSeg1].matchLen<globalSegArray[rowSeg2].matchLen)
							rowSeg1 = rowSeg2;
					}else if(rowSeg1!=rowSeg2 && globalSegArray[rowSeg1].endSubPos==globalSegArray[rowSeg2].endSubPos)
					{
						if(globalSegArray[rowSeg1].matchLen>globalSegArray[rowSeg2].matchLen)
							rowSeg2 = rowSeg1;
					}

					for(j=startItemRow; j<=endItemRow; j++)
					{
						if(j!=rowSeg1 && j!=rowSeg2)
							reduFlagArray[j] = YES;
					}
				}

				startSegPos = endSegPos = -1;
				startItemRow = endItemRow = -1;
			}
		}

		// remove short discordant segments at query ends
		if((globalSegArray[0].subjectID!=globalSegArray[1].subjectID || globalSegArray[0].strand!=globalSegArray[1].strand) && globalSegArray[0].matchLen<1.5*minAlignedSegLenThres)
			reduFlagArray[0] = YES;
		if((globalSegArray[itemNum-2].subjectID!=globalSegArray[itemNum-1].subjectID || globalSegArray[itemNum-2].strand!=globalSegArray[itemNum-1].strand) && globalSegArray[itemNum-1].matchLen<1.5*minAlignedSegLenThres)
			reduFlagArray[0] = YES;

		// remove redundant segments
		newItemNum = 0;
		for(i=0; i<itemNum; i++)
		{
			if(reduFlagArray[i]==NO)
			{
				if(i!=newItemNum)
				{
					if(memcpy(globalSegArray+newItemNum, globalSegArray+i, sizeof(globalValidSeg_t))==NULL)
					{
						printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
				newItemNum ++;
			}
		}

		*globalSegNum = newItemNum;

		free(reduFlagArray);
	}

	return SUCCESSFUL;
}

/**
 * Adjust the order of global align segments for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short adjustGlobalSegOrderSingleQuery(globalValidSeg_t *globalSegArray, int32_t globalSegNum)
{
	int32_t i;
	globalValidSeg_t tmpSeg;

	i = 1;
	while(i<globalSegNum)
	{
		if(globalSegArray[i-1].startQueryPos>globalSegArray[i].startQueryPos && globalSegArray[i-1].endQueryPos>globalSegArray[i].endQueryPos)
		{
			if(memcpy(&tmpSeg, globalSegArray+i-1, sizeof(globalValidSeg_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(memcpy(globalSegArray+i-1, globalSegArray+i, sizeof(globalValidSeg_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(memcpy(globalSegArray+i, &tmpSeg, sizeof(globalValidSeg_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(i>1)
				i --;
			else
				i ++;
		}else
		{
			i ++;
		}
	}

	return SUCCESSFUL;
}

/**
 * Generate global align segment array for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short generateGlobalSegArraySingleQuery(query_t *queryItem, globalValidSeg_t *globalSegArrayBuf, int32_t globalSegNum)
{
	queryItem->globalValidSegNum = 0;
	queryItem->globalValidSegArray = NULL;

	if(globalSegNum>0)
	{
		queryItem->globalValidSegArray = (globalValidSeg_t*) calloc (globalSegNum, sizeof(globalValidSeg_t));
		if(queryItem->globalValidSegArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(memcpy(queryItem->globalValidSegArray, globalSegArrayBuf, globalSegNum*sizeof(globalValidSeg_t))==NULL)
		{
			printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		queryItem->globalValidSegNum = globalSegNum;
	}

	return SUCCESSFUL;
}

/**
 * Trim alignment information of unmatched query ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short trimAlignInfoQueryEnds(query_t *queryItem, subject_t *subjectArray)
{
	int32_t i, globalSegNum, subjectID;
	globalValidSeg_t *globalSegArray;
	char *querySeq, *subjectSeq;

	globalSegArray = queryItem->globalValidSegArray;
	globalSegNum = queryItem->globalValidSegNum;
	querySeq = queryItem->querySeq;

	for(i=0; i<globalSegNum; i++)
	{
		subjectID = globalSegArray[i].subjectID;
		subjectSeq = subjectArray[subjectID-1].subjectSeq;

		// trim left end
		if(trimAlignLeftSegEnd(globalSegArray+i, querySeq, subjectSeq)==FAILED)
		{
			printf("line=%d, In %s(), cannot trim the alignment information at the left end of query segment, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// trim right end
		if(trimAlignRightSegEnd(globalSegArray+i, querySeq, subjectSeq)==FAILED)
		{
			printf("line=%d, In %s(), cannot trim the alignment information at the left end of query segment, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Trim left alignment information of unmatched query ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short trimAlignLeftSegEnd(globalValidSeg_t *globalSeg, char *querySeq, char *subjectSeq)
{
	int32_t winSize, alignSeqSize, increaseSize, remainLenQuery, remainLenSubject, validFlag;
	int32_t startQueryPos, endQueryPos, startSubjectPos, endSubjectPos;
	int64_t leftMargin, leftMarginSubject, leftQueryPos, rightQueryPos, leftSubjectPos, rightSubjectPos;

	char *queryAlignSeq, *subjectAlignSeq;
	int32_t queryAlignSeqLen, subjectAlignSeqLen, maxSeqLen, exactMatchFlag;
	char *alignResultArray[3];
	int32_t overlapLen, mismatchNum, queryLeftShiftLen, queryRightShiftLen, subjectLeftShiftLen, subjectRightShiftLen;
	int32_t validHeadAlignFlag, validTailAlignFlag;

	maxSeqLen = 10000;

	// initialize buffers
	if(initAlignBuf(alignResultArray, &queryAlignSeq, &subjectAlignSeq, maxSeqLen)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the alignment buffers, error!\n", __LINE__, __func__);
		return FAILED;
	}

	startQueryPos = globalSeg->startQueryPos;
	endQueryPos = globalSeg->endQueryPos;
	startSubjectPos = globalSeg->startSubPos;
	endSubjectPos = globalSeg->endSubPos;

	winSize = 500;
	increaseSize = 0;
	validFlag = NO;
	while(validFlag==NO)
	{
		alignSeqSize = winSize + increaseSize;

		if(alignSeqSize>endQueryPos-startQueryPos+1 && increaseSize>0)
			break;

		remainLenQuery = endQueryPos - startQueryPos + 1;
		if(remainLenQuery>=alignSeqSize)
			queryAlignSeqLen = alignSeqSize;
		else
			queryAlignSeqLen = remainLenQuery;

		leftQueryPos = startQueryPos;
		rightQueryPos = startQueryPos + queryAlignSeqLen - 1;


		if(globalSeg->strand==PLUS_STRAND)
		{ // plus strand
			remainLenSubject = endSubjectPos - startSubjectPos + 1;
			if(remainLenSubject>=alignSeqSize)
				subjectAlignSeqLen = alignSeqSize;
			else
				subjectAlignSeqLen = remainLenSubject;

			leftSubjectPos = startSubjectPos;
			rightSubjectPos = startSubjectPos + subjectAlignSeqLen - 1;
		}else
		{ // minus strand
			remainLenSubject = startSubjectPos - (endSubjectPos - 1);
			if(remainLenSubject>=alignSeqSize)
				subjectAlignSeqLen = alignSeqSize;
			else
				subjectAlignSeqLen = remainLenSubject;

			leftSubjectPos = startSubjectPos;
			rightSubjectPos = startSubjectPos - subjectAlignSeqLen + 1;
		}

		// control buffer size
		if(queryAlignSeqLen>maxSeqLen/2 || subjectAlignSeqLen>maxSeqLen/2)
		{
			if(increaseAlignBufSize(alignResultArray, &queryAlignSeq, &subjectAlignSeq, &maxSeqLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot reverse the base sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		// copy base sequences
		strncpy(queryAlignSeq, querySeq+startQueryPos-1, queryAlignSeqLen);
		queryAlignSeq[queryAlignSeqLen] = '\0';

		if(globalSeg->strand==PLUS_STRAND)
		{ // plus strand
			strncpy(subjectAlignSeq, subjectSeq+startSubjectPos-1, subjectAlignSeqLen);
			subjectAlignSeq[subjectAlignSeqLen] = '\0';
		}else
		{ // minus strand
			strncpy(subjectAlignSeq, subjectSeq+startSubjectPos-subjectAlignSeqLen, subjectAlignSeqLen);
			subjectAlignSeq[subjectAlignSeqLen] = '\0';
			if(reverseSeq(subjectAlignSeq, subjectAlignSeqLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot reverse the base sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		// exact match
		if(strcmp(queryAlignSeq, subjectAlignSeq)==0)
			exactMatchFlag = YES;
		else
			exactMatchFlag = NO;

		if(exactMatchFlag==NO)
		{
			printf("LLLLLLLLLL startQueryPos=%d, endQueryPos=%d, queryAlignSeqLen=%d, subjectAlignSeqLen=%d\n", startQueryPos, endQueryPos, queryAlignSeqLen, subjectAlignSeqLen);

			// generate alignment
			if(computeSeqAlignment(alignResultArray, &overlapLen, &mismatchNum, &queryLeftShiftLen, &subjectLeftShiftLen, &queryRightShiftLen, &subjectRightShiftLen, queryAlignSeq, subjectAlignSeq, queryAlignSeqLen, subjectAlignSeqLen, YES)==FAILED)
			{
				printf("line=%d, In %s(), cannot compute the alignment, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(isValidAlignment(&validHeadAlignFlag, &validTailAlignFlag, alignResultArray, overlapLen, queryLeftShiftLen, subjectLeftShiftLen, YES)==FAILED)
			{
				printf("line=%d, In %s(), cannot check the alignment gap, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(validTailAlignFlag==NO)
			{
				increaseSize += 300;
				continue;
			}

			// get the most left matched position
			if(computeLeftMargin(&leftMargin, &leftMarginSubject, alignResultArray, overlapLen, leftQueryPos, rightQueryPos, leftSubjectPos, rightSubjectPos, queryRightShiftLen, subjectRightShiftLen, globalSeg)==FAILED)
			{
				printf("line=%d, In %s(), cannot compute the right margin, error!\n", __LINE__, __func__);
				return FAILED;
			}

			validFlag = YES;
		}else
		{
			validFlag = YES;
		}
	}

	// free buffers
	freeAlignBuf(alignResultArray, &queryAlignSeq, &subjectAlignSeq);

	return SUCCESSFUL;
}

/**
 * Trim right alignment information of unmatched query ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short trimAlignRightSegEnd(globalValidSeg_t *globalSeg, char *querySeq, char *subjectSeq)
{
	int32_t winSize, alignSeqSize, increaseSize, remainLenQuery, remainLenSubject, validFlag;
	int32_t startQueryPos, endQueryPos, startSubjectPos, endSubjectPos;
	int64_t rightMargin, rightMarginSubject, leftQueryPos, rightQueryPos, leftSubjectPos, rightSubjectPos;

	char *queryAlignSeq, *subjectAlignSeq;
	int32_t queryAlignSeqLen, subjectAlignSeqLen, maxSeqLen, exactMatchFlag;
	char *alignResultArray[3];
	int32_t overlapLen, mismatchNum, queryLeftShiftLen, queryRightShiftLen, subjectLeftShiftLen, subjectRightShiftLen;
	int32_t validHeadAlignFlag, validTailAlignFlag;

	maxSeqLen = 10000;

	// initialize buffers
	if(initAlignBuf(alignResultArray, &queryAlignSeq, &subjectAlignSeq, maxSeqLen)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the alignment buffers, error!\n", __LINE__, __func__);
		return FAILED;
	}

	startQueryPos = globalSeg->startQueryPos;
	endQueryPos = globalSeg->endQueryPos;
	startSubjectPos = globalSeg->startSubPos;
	endSubjectPos = globalSeg->endSubPos;

	winSize = 500;
	increaseSize = 0;
	validFlag = NO;
	while(validFlag==NO)
	{
		alignSeqSize = winSize + increaseSize;

		if(alignSeqSize>endQueryPos-startQueryPos+1 && increaseSize>0)
			break;

		remainLenQuery = endQueryPos - startQueryPos + 1;
		if(remainLenQuery>=alignSeqSize)
			queryAlignSeqLen = alignSeqSize;
		else
			queryAlignSeqLen = remainLenQuery;

		leftQueryPos = endQueryPos - queryAlignSeqLen + 1;
		rightQueryPos = endQueryPos;

		if(globalSeg->strand==PLUS_STRAND)
		{ // plus strand
			remainLenSubject = endSubjectPos - startSubjectPos + 1;
			if(remainLenSubject>=alignSeqSize)
				subjectAlignSeqLen = alignSeqSize;
			else
				subjectAlignSeqLen = remainLenSubject;

			leftSubjectPos = endSubjectPos - subjectAlignSeqLen + 1;
			rightSubjectPos = endSubjectPos;
		}else
		{ // minus strand
			remainLenSubject = startSubjectPos - (endSubjectPos - 1);
			if(remainLenSubject>=alignSeqSize)
				subjectAlignSeqLen = alignSeqSize;
			else
				subjectAlignSeqLen = remainLenSubject;

			leftSubjectPos = endSubjectPos + subjectAlignSeqLen - 1;
			rightSubjectPos = endSubjectPos;
		}

		// control buffer size
		if(queryAlignSeqLen>maxSeqLen/2 || subjectAlignSeqLen>maxSeqLen/2)
		{
			if(increaseAlignBufSize(alignResultArray, &queryAlignSeq, &subjectAlignSeq, &maxSeqLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot reverse the base sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		// copy base sequences
		strncpy(queryAlignSeq, querySeq+endQueryPos-queryAlignSeqLen, queryAlignSeqLen);
		queryAlignSeq[queryAlignSeqLen] = '\0';

		if(globalSeg->strand==PLUS_STRAND)
		{ // plus strand
			strncpy(subjectAlignSeq, subjectSeq+endSubjectPos-subjectAlignSeqLen, subjectAlignSeqLen);
			subjectAlignSeq[subjectAlignSeqLen] = '\0';
		}else
		{ // minus strand
			strncpy(subjectAlignSeq, subjectSeq+endSubjectPos-1, subjectAlignSeqLen);
			subjectAlignSeq[subjectAlignSeqLen] = '\0';
			if(reverseSeq(subjectAlignSeq, subjectAlignSeqLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot reverse the base sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		// exact match
		if(strcmp(queryAlignSeq, subjectAlignSeq)==0)
			exactMatchFlag = YES;
		else
			exactMatchFlag = NO;

		if(exactMatchFlag==NO)
		{
			printf("RRRRRRRRRR startQueryPos=%d, endQueryPos=%d, queryAlignSeqLen=%d, subjectAlignSeqLen=%d\n", startQueryPos, endQueryPos, queryAlignSeqLen, subjectAlignSeqLen);

			// generate alignment
			if(computeSeqAlignment(alignResultArray, &overlapLen, &mismatchNum, &queryLeftShiftLen, &subjectLeftShiftLen, &queryRightShiftLen, &subjectRightShiftLen, queryAlignSeq, subjectAlignSeq, queryAlignSeqLen, subjectAlignSeqLen, YES)==FAILED)
			{
				printf("line=%d, In %s(), cannot compute the alignment, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(isValidAlignment(&validHeadAlignFlag, &validTailAlignFlag, alignResultArray, overlapLen, queryLeftShiftLen, subjectLeftShiftLen, YES)==FAILED)
			{
				printf("line=%d, In %s(), cannot check the alignment gap, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(validHeadAlignFlag==NO)
			{
				increaseSize += 300;
				continue;
			}

			// get the most right matched position
			if(computeRightMargin(&rightMargin, &rightMarginSubject, alignResultArray, overlapLen, leftQueryPos, rightQueryPos, leftSubjectPos, rightSubjectPos, queryLeftShiftLen, subjectLeftShiftLen, globalSeg)==FAILED)
			{
				printf("line=%d, In %s(), cannot compute the right margin, error!\n", __LINE__, __func__);
				return FAILED;
			}

			validFlag = YES;
		}else
		{
			validFlag = YES;
		}
	}

	// free buffers
	freeAlignBuf(alignResultArray, &queryAlignSeq, &subjectAlignSeq);

	return SUCCESSFUL;
}

/**
 * Determine global match kind for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short determineGlobalMatchKindSingleQuery(query_t *queryItem)
{
	int32_t i, globalSegNum, bestSubjectID, bestMatchKind, sameSubjectFlag, validItemNum, sumLen;
	globalValidSeg_t *globalSegArray;
	double autoMatchPercentThres, baseCovSingleQuery;

	globalSegArray = queryItem->globalValidSegArray;
	globalSegNum = queryItem->globalValidSegNum;

	// auto adjust the matchPercentThres
	if(queryItem->queryLen>shortQueryLenThres)
		autoMatchPercentThres = matchPercentThres;
	else
		autoMatchPercentThres = matchPercentThres * matchPercentFactor;

	if(globalSegNum==0)
	{
		queryItem->globalMatchKind = UNMATCHED_KIND;
	}else
	{
		sameSubjectFlag = YES;
		if(queryItem->bestMatchRow>=0)
		{
			bestMatchKind = queryItem->querySubArray[queryItem->bestMatchRow].matchKind;
			if(bestMatchKind==PERFECT_MATCH_KIND)
				queryItem->globalMatchKind = PERFECT_MATCH_KIND;
			else
			{
				bestSubjectID = queryItem->querySubArray[queryItem->bestMatchRow].subjectID;
				for(i=0; i<queryItem->globalValidSegNum; i++)
				{
					if(globalSegArray[i].subjectID!=bestSubjectID)
					{
						sameSubjectFlag = NO;
						break;
					}
				}

				if(sameSubjectFlag==YES)
					queryItem->globalMatchKind = bestMatchKind;
				else
				{
					if(bestMatchKind==MATCHED_KIND)
					{
						queryItem->globalMatchKind = DISJUNCT_MATCH_KIND;
					}else if(bestMatchKind==DISJUNCT_MATCH_KIND)
					{
						queryItem->globalMatchKind = DISJUNCT_MATCH_KIND;
					}else if(bestMatchKind==UNMATCHED_KIND)
					{
						// check single align segment, total length
						validItemNum = 0;
						sumLen = 0;
						for(i=0; i<globalSegNum; i++)
						{
							if((double)globalSegArray[i].matchLen/globalSegArray[i].totalMatchLen>autoMatchPercentThres)
								validItemNum ++;

							sumLen += globalSegArray[i].totalMatchLen;
						}

						if(validItemNum<=0 || sumLen<matchPercentThres*queryItem->queryLen)
							queryItem->globalMatchKind = UNMATCHED_KIND;
						else
						{
							if(getBaseCovPercentSingleQuery(&baseCovSingleQuery, queryItem)==FAILED)
							{
								printf("line=%d, In %s(), cannot get the base coverage of single query, error!\n", __LINE__, __func__);
								return FAILED;
							}

							if(baseCovSingleQuery>=matchPercentThres)
								queryItem->globalMatchKind = DISJUNCT_MATCH_KIND;
							else
								queryItem->globalMatchKind = UNMATCHED_KIND;
						}
					}else
					{
						printf("line=%d, In %s(), invalid matchKind=%d, error!\n", __LINE__, __func__, bestMatchKind);
						return FAILED;
					}
				}
			}
		}else
		{
			queryItem->globalMatchKind = UNMATCHED_KIND;
		}
	}

	return SUCCESSFUL;
}

/**
 * Get base coverage percentage for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short getBaseCovPercentSingleQuery(double *baseCovSingleQuery, query_t *queryItem)
{
	int32_t i, j, globalSegNum, queryLen, startQueryPos, endQueryPos, covNum;
	globalValidSeg_t *globalSegArray;
	int8_t *covArray;

	globalSegArray = queryItem->globalValidSegArray;
	globalSegNum = queryItem->globalValidSegNum;

	queryLen = queryItem->queryLen;
	covArray = (int8_t *) calloc (queryLen, sizeof(int8_t));
	if(covArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	for(i=0; i<queryLen; i++) covArray[i] = NO;

	for(i=0; i<globalSegNum; i++)
	{
		startQueryPos = globalSegArray[i].startQueryPos;
		endQueryPos = globalSegArray[i].endQueryPos;
		for(j=startQueryPos; j<=endQueryPos; j++)
			covArray[j-1] = YES;
	}

	covNum = 0;
	for(i=0; i<queryLen; i++)
		if(covArray[i]==YES)
			covNum ++;

	*baseCovSingleQuery = (double)covNum/queryLen;

	free(covArray);

	return SUCCESSFUL;
}
