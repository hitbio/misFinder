/*
 * misinfo.c
 *
 *  Created on: Aug 2, 2014
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Get the mis-assembly information list for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMisInfoList(query_t *queryItem, subject_t *subjectArray)
{
	// get the misInfo list from breakpoints
	if(getMisInfoListBreakpoint(queryItem, subjectArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the mis-assembly information list, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// adjust breakpoints
	if(adjustMisInfoList(queryItem)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust the mis-assembly information list, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the misInfo list from inner
	if(getMisInfoListInner(queryItem, subjectArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the mis-assembly information list, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get the mis-assembly information list from breakpoints for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMisInfoListBreakpoint(query_t *queryItem, subject_t *subjectArray)
{
	int32_t i, tmp, globalSegNum, subjectID, subjectLen, queryID, queryLen, startItemRow, endItemRow,  misType, misjoinFlag, circularFlag;
	int32_t difQuery, difSubject, indelKind, gapFlag, circularRegFlag;
	int64_t leftMargin, rightMargin, leftMarginSubject, rightMarginSubject, startSegPos, endSegPos;
	globalValidSeg_t *globalSegArray;
	querySubject_t *querySubItem;
	char *subjectSeq, *querySeq;

	queryItem->misInfoList = NULL;
	queryItem->misInfoItemNum = 0;

	globalSegArray = queryItem->globalValidSegArray;
	globalSegNum = queryItem->globalValidSegNum;

	if(queryItem->bestMatchRow>=0)
	{
		querySubItem = queryItem->querySubArray + queryItem->bestMatchRow;
		subjectID = querySubItem->subjectID;
		subjectLen = subjectArray[subjectID-1].subjectLen;
		subjectSeq = subjectArray[subjectID-1].subjectSeq;

		queryID = queryItem->queryID;
		querySeq = queryItem->querySeq;
		queryLen = queryItem->queryLen;
		circularFlag = queryItem->circularFlag;

		// check head segment
		if(globalSegArray[0].startQueryPos>=3)
		{
			rightMargin = globalSegArray[0].startQueryPos;
			leftMargin = rightMargin - 1;
			difQuery = difSubject = -1;

			if(addNewQueryMisInfoNode(queryItem, -1, 0, QUERY_INDEL_KIND, QUERY_INSERT, NO, NO, leftMargin, rightMargin, -1, -1, -1, difQuery, difSubject)==FAILED)
			{
				printf("line=%d, In %s(), cannot add new query mis-assembly information node, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		startSegPos = endSegPos = -1;
		startItemRow = endItemRow = -1;
		for(i=1; i<globalSegNum; i++)
		{
			misType = -1;

			// compute the overlapped align segments of query
			if(startSegPos==-1)
			{
				if(globalSegArray[i-1].subjectID!=globalSegArray[i].subjectID || globalSegArray[i-1].strand!=globalSegArray[i].strand || globalSegArray[i-1].endQueryPos>globalSegArray[i].startQueryPos || (globalSegArray[i-1].endSubPos-globalSegArray[i].startSubPos<-minDisjunctDistanceThres || globalSegArray[i-1].endSubPos-globalSegArray[i].startSubPos>minDisjunctDistanceThres))
				{
					startSegPos = globalSegArray[i].startQueryPos;
					endSegPos = globalSegArray[i-1].endQueryPos;
					startItemRow = i - 1;

					if(i==globalSegNum-1 || (globalSegArray[i-1].strand!=globalSegArray[i].strand || globalSegArray[i-1].subjectID!=globalSegArray[i].subjectID))
						endItemRow = i;
					else if(globalSegArray[i+1].startQueryPos>endSegPos)
						endItemRow = i;
				}
			}else
			{
				if(globalSegArray[i].startQueryPos<startSegPos)
					startSegPos = globalSegArray[i].startQueryPos;

				if(globalSegArray[i].endQueryPos>endSegPos && globalSegArray[i].endQueryPos<endSegPos+200)
					endSegPos = globalSegArray[i].endQueryPos;

				if(i==globalSegNum-1 || (globalSegArray[i-1].strand!=globalSegArray[i].strand || globalSegArray[i-1].subjectID!=globalSegArray[i].subjectID))
					endItemRow = i;
				else if(globalSegArray[i+1].startQueryPos>endSegPos)
					endItemRow = i;
			}

			// detect the misjoin
			if(startItemRow>=0 && endItemRow>=0)
			{
				// valid the segments
				if(computeMisassFlagAlignSeg(&misjoinFlag, globalSegArray+startItemRow, globalSegArray+endItemRow, subjectLen, circularFlag)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute the mis-assembly flag, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// compute the align segment margin
				if(misjoinFlag==YES)
				{
					difQuery = globalSegArray[i].startQueryPos - globalSegArray[i-1].endQueryPos;
					difSubject = -1;
					misType = QUERY_MISJOIN_KIND;
				}
			}

			// detect the indel
			if(misType==-1)
			{
				startItemRow = endItemRow = -1;
				indelKind = -1;
				circularRegFlag = NO;
				if(globalSegArray[i-1].subjectID==globalSegArray[i].subjectID && globalSegArray[i-1].strand==globalSegArray[i].strand)
				{ // same strand
					difQuery = globalSegArray[i].startQueryPos - globalSegArray[i-1].endQueryPos;
					if(circularFlag==YES)
					{
						if(globalSegArray[i].strand==PLUS_STRAND)
						{
							if(globalSegArray[i-1].endSubPos>=subjectLen-varyEndLenThres && globalSegArray[i].startSubPos<=varyEndLenThres)
							{
								difSubject = globalSegArray[i].startSubPos + subjectLen - globalSegArray[i-1].endSubPos;
								circularRegFlag = YES;
							}else
								difSubject = globalSegArray[i].startSubPos - globalSegArray[i-1].endSubPos;
						}else
						{
							if(globalSegArray[i-1].endSubPos<=varyEndLenThres && globalSegArray[i].startSubPos>=subjectLen-varyEndLenThres)
							{
								difSubject = globalSegArray[i-1].endSubPos + subjectLen - globalSegArray[i].startSubPos;
								circularRegFlag = YES;
							}else
								difSubject = globalSegArray[i-1].endSubPos - globalSegArray[i].startSubPos;
						}
					}else
					{
						if(globalSegArray[i].strand==PLUS_STRAND)
							difSubject = globalSegArray[i].startSubPos - globalSegArray[i-1].endSubPos;
						else
							difSubject = globalSegArray[i-1].endSubPos - globalSegArray[i].startSubPos;
					}

					startItemRow = i - 1;
					endItemRow = i;
					if(difQuery-difSubject>=indelSizeThres && difSubject<=indelSizeThres)
					{
						// further check misjoin
						misjoinFlag = NO;
						if(globalSegArray[i].strand==PLUS_STRAND && globalSegArray[i].endSubPos<=globalSegArray[i-1].startSubPos)
							misjoinFlag = YES;
						else if(globalSegArray[i].strand==MINUS_STRAND && globalSegArray[i].endSubPos>=globalSegArray[i-1].startSubPos)
							misjoinFlag = YES;

						if(misjoinFlag==YES)
							misType = QUERY_MISJOIN_KIND;
						else
							indelKind = QUERY_INSERT;
					}else if(difSubject-difQuery>=indelSizeThres && difQuery<=indelSizeThres)
					{
						indelKind = QUERY_DEL;
					}else if(difQuery>=indelSizeThres && difSubject>=indelSizeThres)
					{
						// further check misjoin
						misjoinFlag = NO;
						if(globalSegArray[i].strand==PLUS_STRAND && globalSegArray[i].endSubPos<=globalSegArray[i-1].startSubPos)
							misjoinFlag = YES;
						else if(globalSegArray[i].strand==MINUS_STRAND && globalSegArray[i].endSubPos>=globalSegArray[i-1].startSubPos)
							misjoinFlag = YES;

						if(misjoinFlag==YES)
							misType = QUERY_MISJOIN_KIND;
						else
							indelKind = QUERY_GAP;
					}else
					{
						if(circularRegFlag==YES)
						{
							startSegPos = endSegPos = -1;
							startItemRow = endItemRow = -1;
						}else
							printf("line=%d, In %s(), startItemRow=%d, endItemRow=%d, difQuery=%d, difSubject=%d\n", __LINE__, __func__, startItemRow, endItemRow, difQuery, difSubject);
					}
				}else
				{ // different strands
					printf("#+#+#+#+#+#+##+#+#+ line=%d, In %s(), different subjects or strands, need further checking ...\n", __LINE__, __func__);
					return FAILED;
				}

				if(indelKind>0)
					misType = QUERY_INDEL_KIND;
			}

			if(misType>=0 && (startItemRow>=0 && endItemRow>=0))
			{ // misjoin

				// get the gap flag
				if(getGapFlagMisInfo(&gapFlag, globalSegArray+startItemRow, globalSegArray+endItemRow, queryItem)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the gap flag of misInfo, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// compute the align segment margin
				if(misType==QUERY_MISJOIN_KIND)
				{ // misjoin
					if(computeAlignSegMargin(&leftMargin, &rightMargin, globalSegArray+startItemRow, globalSegArray+endItemRow, difQuery, difSubject, queryItem, subjectArray)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the initial overlap margins, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(addNewQueryMisInfoNode(queryItem, startItemRow, endItemRow, QUERY_MISJOIN_KIND, -1, gapFlag, NO, leftMargin, rightMargin, -1, -1, -1, -1, -1)==FAILED)
					{
						printf("line=%d, In %s(), cannot add new query mis-assembly information node, error!\n", __LINE__, __func__);
						return FAILED;
					}

					startSegPos = endSegPos = -1;
					startItemRow = endItemRow = -1;
				}else
				{ // indel
					leftMargin = globalSegArray[startItemRow].endQueryPos;
					rightMargin = globalSegArray[endItemRow].startQueryPos;
					leftMarginSubject = globalSegArray[startItemRow].endSubPos;
					rightMarginSubject = globalSegArray[endItemRow].startSubPos;
					subjectID = globalSegArray[startItemRow].subjectID;

					if(leftMargin>rightMargin)
					{
						tmp = leftMargin; leftMargin = rightMargin; rightMargin = tmp;
						tmp = leftMarginSubject; leftMarginSubject = rightMarginSubject; rightMarginSubject = tmp;
					}

					if(addNewQueryMisInfoNode(queryItem, startItemRow, endItemRow, QUERY_INDEL_KIND, indelKind, gapFlag, NO, leftMargin, rightMargin, leftMarginSubject, rightMarginSubject, subjectID, difQuery, difSubject)==FAILED)
					{
						printf("line=%d, In %s(), cannot add new query mis-assembly information node, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

				startSegPos = endSegPos = -1;
				startItemRow = endItemRow = -1;
			}
		}

		// check tail segment
		if(globalSegArray[globalSegNum-1].endQueryPos<=queryItem->queryLen-5)
		{
			leftMargin = globalSegArray[globalSegNum-1].endQueryPos;
			rightMargin = leftMargin + 1;

			if(addNewQueryMisInfoNode(queryItem, globalSegNum-1, -1, QUERY_INDEL_KIND, QUERY_INSERT, NO, NO, leftMargin, rightMargin, -1, -1, -1, -1, -1)==FAILED)
			{
				printf("line=%d, In %s(), cannot add new query mis-assembly information node, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}else
	{
		startItemRow = endItemRow = -1;
		leftMargin = 1;
		rightMargin = queryItem->queryLen;
		if(addNewQueryMisInfoNode(queryItem, startItemRow, endItemRow, QUERY_INDEL_KIND, QUERY_INSERT, NO, NO, leftMargin, rightMargin, -1, -1, -1, -1, -1)==FAILED)
		{
			printf("line=%d, In %s(), cannot add new query mis-assembly information node, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Adjust the mis-assembly information list.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short adjustMisInfoList(query_t *queryItem)
{
	misInfo_t *misInfo, *misInfoNext, *misInfoPre;
	queryIndel_t *queryIndel;
	int32_t queryLen, startQPos, endQPos, mismatchNum;

	queryLen = queryItem->queryLen;

	// merge adjacent misjoin nodes
	misInfo = queryItem->misInfoList;
	while(misInfo)
	{
		misInfoNext = misInfo->next;
		if(misInfo->misType==QUERY_MISJOIN_KIND && misInfoNext && misInfoNext->misType==QUERY_MISJOIN_KIND)
		{
			if(misInfoNext->queryMargin->rightMargin-misInfo->queryMargin->leftMargin<1000 || misInfoNext->queryMargin->leftMargin-misInfo->queryMargin->rightMargin<50)
			{
				if(misInfo->queryMargin->leftMargin>misInfoNext->queryMargin->leftMargin)
					misInfo->queryMargin->leftMargin = misInfoNext->queryMargin->leftMargin;

				if(misInfo->queryMargin->rightMargin<misInfoNext->queryMargin->rightMargin)
					misInfo->queryMargin->rightMargin = misInfoNext->queryMargin->rightMargin;

				misInfo->rightSegRow = misInfoNext->rightSegRow;

				free(misInfoNext->queryMargin);
				misInfo->next = misInfoNext->next;
				if(queryItem->tailMisInfo==misInfoNext)
					queryItem->tailMisInfo = misInfo;
				free(misInfoNext);
				queryItem->misInfoItemNum --;

				continue;
			}
		}

		misInfo = misInfo->next;
	}

/*
	// check the query ends
	misInfo = queryItem->misInfoList;
	misInfoPre = NULL;
	while(misInfo)
	{
		if(misInfo->misType==QUERY_MISJOIN_KIND && (misInfo->queryMargin->leftMargin<1000 || misInfo->queryMargin->rightMargin<1000 || misInfo->queryMargin->leftMargin>queryLen-1000 || misInfo->queryMargin->rightMargin>queryLen-1000))
		{
			misInfoNext = misInfo->next;
			free(misInfo->queryMargin);
			if(misInfoPre==NULL)
				queryItem->misInfoList = misInfoNext;
			else
				misInfoPre->next = misInfoNext;
			if(queryItem->tailMisInfo==misInfo)
				queryItem->tailMisInfo = misInfoPre;
			free(misInfo);
			misInfo = misInfoNext;
			queryItem->misInfoItemNum --;

			continue;
		}

		misInfoPre = misInfo;
		misInfo = misInfo->next;
	}
*/

	return SUCCESSFUL;
}

/**
 * Get the mis-assembly information list from inner for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMisInfoListInner(query_t *queryItem, subject_t *subjectArray)
{
	int32_t i, globalSegNum, subjectID, subjectLen, queryID, queryLen, validIndelFlag;
	int32_t winSize, stepSize, increaseSize, alignSeqSize, strand, startAlignRowQuery, startAlignRowSubject, remainLenQuery, remainLenSubject;
	int32_t startQueryPos, endQueryPos, startSubjectPos, endSubjectPos, rightMostQueryPos, rightMostSubjectPos;
	int32_t gapNumQuery, gapNumSubject, startAlignPosQuery, startAlignPosSubject;
	globalValidSeg_t *globalSegArray;
	char *subjectSeq, *querySeq;

	char *queryAlignSeq, *subjectAlignSeq;
	int32_t queryAlignSeqLen, subjectAlignSeqLen, maxSeqLen, exactMatchFlag, validHeadAlignFlag, validTailAlignFlag, continueFlag, recomputeStartAlignPosFlag;
	char *alignResultArray[3];
	int32_t overlapLen, mismatchNum, queryLeftShiftLen, queryRightShiftLen, subjectLeftShiftLen, subjectRightShiftLen;


	winSize = 300;
	stepSize = 200;
	globalSegArray = queryItem->globalValidSegArray;
	globalSegNum = queryItem->globalValidSegNum;

	if(globalSegNum>0)
	{
		maxSeqLen = 10000;

		// initialize buffers
		if(initAlignBuf(alignResultArray, &queryAlignSeq, &subjectAlignSeq, maxSeqLen)==FAILED)
		{
			printf("line=%d, In %s(), cannot initialize the alignment buffers, error!\n", __LINE__, __func__);
			return FAILED;
		}

		queryID = queryItem->queryID;
		querySeq = queryItem->querySeq;
		queryLen = queryItem->queryLen;

		for(i=0; i<globalSegNum; i++)
		{
			if(globalSegArray[i].gapNum>=indelSizeThres)
			{
				subjectID = globalSegArray[i].subjectID;
				subjectLen = subjectArray[subjectID-1].subjectLen;
				subjectSeq = subjectArray[subjectID-1].subjectSeq;
				strand = globalSegArray[i].strand;

				// compute start alignment positions for query and subject
				if(computeStartAlignPos(&startAlignPosQuery, &startAlignPosSubject, -1, -1, globalSegArray+i, querySeq, subjectSeq, 200)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute the first shift size, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(startAlignPosQuery==-1 || startAlignPosSubject==-1)
					continue;

				startQueryPos = startAlignPosQuery;
				startSubjectPos = startAlignPosSubject;
				endQueryPos = globalSegArray[i].endQueryPos - stepSize;
				if(globalSegArray[i].strand==PLUS_STRAND)
					endSubjectPos = globalSegArray[i].endSubPos - stepSize;
				else
					endSubjectPos = globalSegArray[i].endSubPos + stepSize;

				// loop to detect the indel regions
				increaseSize = 0;
				startAlignRowQuery = startQueryPos - 1;
				startAlignRowSubject = startSubjectPos - 1;
				rightMostQueryPos = startQueryPos;
				rightMostSubjectPos = startSubjectPos;
				while(startAlignRowQuery<endQueryPos)
				{
					//if(queryID==68 && startAlignRowQuery==120180)
					//{
					//	printf("queryID=%d, startAlignRowQuery=%d, startAlignRowSubject=%d\n", queryID, startAlignRowQuery, startAlignRowSubject);
					//}

					// get the align base sequences
					alignSeqSize = winSize + increaseSize;
					remainLenQuery = endQueryPos - startAlignRowQuery;
					if(remainLenQuery>=alignSeqSize)
						queryAlignSeqLen = alignSeqSize;
					else
						queryAlignSeqLen = remainLenQuery;

					if(strand==PLUS_STRAND)
					{ // plus strand
						remainLenSubject = endSubjectPos - startAlignRowSubject;
						if(remainLenSubject>=alignSeqSize)
							subjectAlignSeqLen = alignSeqSize;
						else
							subjectAlignSeqLen = remainLenSubject;
					}else
					{ // minus strand
						remainLenSubject = (startAlignRowSubject + 1) - (endSubjectPos - 1);
						if(remainLenSubject>=alignSeqSize)
							subjectAlignSeqLen = alignSeqSize;
						else
							subjectAlignSeqLen = remainLenSubject;
					}

					// control loop and buffer size
					if((queryAlignSeqLen<=0 || subjectAlignSeqLen<=0) || (alignSeqSize>remainLenQuery+winSize || alignSeqSize>remainLenSubject+winSize))
						break;
					else if(queryAlignSeqLen>maxSeqLen/2 || subjectAlignSeqLen>maxSeqLen/2)
					{
						if(increaseAlignBufSize(alignResultArray, &queryAlignSeq, &subjectAlignSeq, &maxSeqLen)==FAILED)
						{
							printf("line=%d, In %s(), cannot reverse the base sequence, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}

					// copy base sequences
					strncpy(queryAlignSeq, querySeq+startAlignRowQuery, queryAlignSeqLen);
					queryAlignSeq[queryAlignSeqLen] = '\0';

					if(strand==PLUS_STRAND)
					{ // plus strand
						strncpy(subjectAlignSeq, subjectSeq+startAlignRowSubject, subjectAlignSeqLen);
						subjectAlignSeq[subjectAlignSeqLen] = '\0';
					}else
					{ // minus strand
						strncpy(subjectAlignSeq, subjectSeq+(startAlignRowSubject+1)-subjectAlignSeqLen, subjectAlignSeqLen);
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
						// generate alignment
						if(computeSeqAlignment(alignResultArray, &overlapLen, &mismatchNum, &queryLeftShiftLen, &subjectLeftShiftLen, &queryRightShiftLen, &subjectRightShiftLen, queryAlignSeq, subjectAlignSeq, queryAlignSeqLen, subjectAlignSeqLen, NO)==FAILED)
						{
							printf("line=%d, In %s(), cannot compute the alignment, error!\n", __LINE__, __func__);
							return FAILED;
						}

						// adjust alignment gap
						if(adjustAlignmentGap(&gapNumQuery, &gapNumSubject, alignResultArray, overlapLen, NO)==FAILED)
						{
							printf("line=%d, In %s(), cannot adjust the alignment gap, error!\n", __LINE__, __func__);
							return FAILED;
						}

						if(isValidAlignment(&validHeadAlignFlag, &validTailAlignFlag, alignResultArray, overlapLen, queryLeftShiftLen, subjectLeftShiftLen, NO)==FAILED)
						{
							printf("line=%d, In %s(), cannot check the alignment gap, error!\n", __LINE__, __func__);
							return FAILED;
						}

						validIndelFlag = YES;
						if(validHeadAlignFlag==NO || validTailAlignFlag==NO)
						{
							validIndelFlag = NO;
							continueFlag = NO;
							recomputeStartAlignPosFlag = NO;
							if(validHeadAlignFlag==NO)
							{
								if(startAlignRowQuery>=startQueryPos-1)
								{
									startAlignRowQuery -= stepSize;
									if(strand==PLUS_STRAND)
										startAlignRowSubject -= stepSize;
									else
										startAlignRowSubject += stepSize;
									increaseSize += stepSize;
									continueFlag = YES;
								}

								if(rightMostQueryPos-startAlignRowQuery-1>=winSize*3)
								{
									recomputeStartAlignPosFlag = YES;
									continueFlag = YES;
								}
							}

							if(validTailAlignFlag==NO && increaseSize<=stepSize*10)
							{
								increaseSize += stepSize;
								continueFlag = YES;
							}

							if(continueFlag==YES)
							{
								if(recomputeStartAlignPosFlag==YES)
								{
									// compute start alignment positions for query and subject
									if(computeStartAlignPos(&startAlignPosQuery, &startAlignPosSubject, rightMostQueryPos, rightMostSubjectPos, globalSegArray+i, querySeq, subjectSeq, 200)==FAILED)
									{
										printf("line=%d, In %s(), cannot compute the first shift size, error!\n", __LINE__, __func__);
										return FAILED;
									}

									if(startAlignPosQuery==-1 || startAlignPosSubject==-1)
										break;

									startQueryPos = startAlignPosQuery;
									startSubjectPos = startAlignPosSubject;

									increaseSize = 0;
									startAlignRowQuery = startQueryPos - 1;
									startAlignRowSubject = startSubjectPos - 1;
									rightMostQueryPos = startQueryPos;
									rightMostSubjectPos = startSubjectPos;
								}

								continue;
							}
						}

						if(validIndelFlag==YES)
						{
							// add the inner indel to misInfo list
							if(gapNumQuery-gapNumSubject>=indelSizeThres || gapNumQuery-gapNumSubject<=-indelSizeThres)
							{
								if(addInnerIndelToMisInfoList(queryItem, startAlignRowQuery+1, startAlignRowSubject+1, subjectID, strand, alignResultArray, overlapLen, queryAlignSeqLen, subjectAlignSeqLen, queryLeftShiftLen, subjectLeftShiftLen, queryRightShiftLen, subjectRightShiftLen, gapNumQuery, gapNumSubject)==FAILED)
								{
									printf("line=%d, In %s(), cannot compute the indel margin, error!\n", __LINE__, __func__);
									return FAILED;
								}
							}

							if(increaseSize>0)
								increaseSize = 0;
						}
					}

					// compute new start row of query and subject
					if(exactMatchFlag==YES)
					{
						startAlignRowQuery += queryAlignSeqLen;

						if(strand==PLUS_STRAND)
							startAlignRowSubject += subjectAlignSeqLen;
						else
							startAlignRowSubject -= subjectAlignSeqLen;
					}else
					{
						if(queryAlignSeqLen-queryRightShiftLen==0 || subjectAlignSeqLen-subjectRightShiftLen==0)
							break;

						startAlignRowQuery += queryAlignSeqLen - queryRightShiftLen;

						if(strand==PLUS_STRAND)
							startAlignRowSubject += subjectAlignSeqLen - subjectRightShiftLen;
						else
							startAlignRowSubject -= subjectAlignSeqLen - subjectRightShiftLen;
					}

					if(rightMostQueryPos<startAlignRowQuery+1)
						rightMostQueryPos = startAlignRowQuery + 1;
					if((strand==PLUS_STRAND && rightMostSubjectPos<startAlignRowSubject+1) || (strand==MINUS_STRAND && rightMostSubjectPos>startAlignRowSubject+1))
						rightMostSubjectPos = startAlignRowSubject + 1;
				}
			}
		}

		// free buffers
		freeAlignBuf(alignResultArray, &queryAlignSeq, &subjectAlignSeq);
	}

	return SUCCESSFUL;
}

/**
 * Compute start alignment positions for query and subject.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeStartAlignPos(int32_t *startAlignPosQuery, int32_t *startAlignPosSubject, int32_t startQueryPos, int32_t startSubjectPos, globalValidSeg_t *globalSeg, char *querySeq, char *subjectSeq, int32_t initAlignSize)
{
	int32_t winSize, alignSeqSize, increaseSize, remainLenQuery, remainLenSubject, validFlag;
	int32_t endQueryPos, endSubjectPos;

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

	*startAlignPosQuery = -1;
	*startAlignPosSubject = -1;

	if(startQueryPos<=0)
		startQueryPos = globalSeg->startQueryPos;
	if(startSubjectPos<=0)
		startSubjectPos = globalSeg->startSubPos;
	endQueryPos = globalSeg->endQueryPos;
	endSubjectPos = globalSeg->endSubPos;

	winSize = initAlignSize;
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

		if(globalSeg->strand==PLUS_STRAND)
		{ // plus strand
			remainLenSubject = endSubjectPos - startSubjectPos + 1;
			if(remainLenSubject>=alignSeqSize)
				subjectAlignSeqLen = alignSeqSize;
			else
				subjectAlignSeqLen = remainLenSubject;
		}else
		{ // minus strand
			remainLenSubject = startSubjectPos - (endSubjectPos - 1);
			if(remainLenSubject>=alignSeqSize)
				subjectAlignSeqLen = alignSeqSize;
			else
				subjectAlignSeqLen = remainLenSubject;
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
			// generate alignment
			if(computeSeqAlignment(alignResultArray, &overlapLen, &mismatchNum, &queryLeftShiftLen, &subjectLeftShiftLen, &queryRightShiftLen, &subjectRightShiftLen, queryAlignSeq, subjectAlignSeq, queryAlignSeqLen, subjectAlignSeqLen, NO)==FAILED)
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
				increaseSize += 200;
				continue;
			}

			*startAlignPosQuery = startQueryPos + queryAlignSeqLen - queryRightShiftLen;
			if(globalSeg->strand==PLUS_STRAND)
				*startAlignPosSubject = startSubjectPos + subjectAlignSeqLen - subjectRightShiftLen;
			else
				*startAlignPosSubject = startSubjectPos - subjectAlignSeqLen + subjectRightShiftLen;

			validFlag = YES;
		}else
		{
			*startAlignPosQuery = startQueryPos + queryAlignSeqLen;
			if(globalSeg->strand==PLUS_STRAND)
				*startAlignPosSubject = startSubjectPos + subjectAlignSeqLen;
			else
				*startAlignPosSubject = startSubjectPos - subjectAlignSeqLen;

			validFlag = YES;
		}
	}

	// free buffers
	freeAlignBuf(alignResultArray, &queryAlignSeq, &subjectAlignSeq);

	return SUCCESSFUL;
}

/**
 * Initialize alignment buffers.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initAlignBuf(char **alignResultArray, char **alignSeq1, char **alignSeq2, int32_t maxSeqLen)
{
	int32_t i;

	for(i=0; i<3; i++)
	{
		alignResultArray[i] = (char*) calloc(maxSeqLen+1, sizeof(char));
		if(alignResultArray[i]==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}
	*alignSeq1 = (char *) calloc (maxSeqLen+1, sizeof(char));
	*alignSeq2 = (char *) calloc (maxSeqLen+1, sizeof(char));
	if((*alignSeq1)==NULL || (*alignSeq2)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Free alignment buffers.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
void freeAlignBuf(char **alignResultArray, char **alignSeq1, char **alignSeq2)
{
	int32_t i;

	for(i=0; i<3; i++)
	{
		free(alignResultArray[i]);
		alignResultArray[i] = NULL;
	}
	free(*alignSeq1);
	free(*alignSeq2);
	*alignSeq1 = NULL;
	*alignSeq2 = NULL;
}

/**
 * Increase alignment buffer size.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short increaseAlignBufSize(char **alignResultArray, char **alignSeq1, char **alignSeq2, int32_t *maxSeqLen)
{
	int32_t i;

	(*maxSeqLen) *= 2;

	for(i=0; i<3; i++)
	{
		alignResultArray[i] = (char*) realloc(alignResultArray[i], ((*maxSeqLen)+1) * sizeof(char));
		if(alignResultArray[i]==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	*alignSeq1 = (char *) realloc (*alignSeq1, ((*maxSeqLen)+1) * sizeof(char));
	*alignSeq2 = (char *) realloc (*alignSeq2, ((*maxSeqLen)+1) * sizeof(char));
	if((*alignSeq1)==NULL || (*alignSeq2)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Compute alignment gap count.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeAlignGapNum(int32_t *gapNum1, int32_t *gapNum2, int32_t *misNum, char **alignResultArray, int32_t overlapLen)
{
	int32_t i;

	*gapNum1 = *gapNum2 = 0;
	*misNum = 0;
	for(i=0; i<overlapLen; i++)
	{
		if(alignResultArray[1][i]!='|')
		{
			(*misNum) ++;

			if(alignResultArray[0][i]=='-' || alignResultArray[0][i]=='N' || alignResultArray[0][i]=='n' || alignResultArray[0][i]=='.')
				(*gapNum1) ++;
			else if(alignResultArray[2][i]=='-' || alignResultArray[2][i]=='N' || alignResultArray[2][i]=='n' || alignResultArray[2][i]=='.')
				(*gapNum2) ++;
		}
	}

	return SUCCESSFUL;
}

/**
 * Adjust alignment gap.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short adjustAlignmentGap(int32_t *gapNum1, int32_t *gapNum2, char **alignResultArray, int32_t overlapLen, int32_t printFlag)
{
	int32_t i, j, mismatchNum, tmpMismatchNum, tmpGapNum, matchRow, mostRightRow, mostLeftRow;
	char tmp;

	if(computeAlignGapNum(gapNum1, gapNum2, &mismatchNum, alignResultArray, overlapLen)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the alignment gap count, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if((*gapNum1)>0 || (*gapNum2)>0)
	{
		// check from left to right
		mostRightRow = -INT_MAX;
		for(i=0; i<overlapLen; i++)
		{
			if(alignResultArray[0][i]=='-')
			{ // gap in seq1
				matchRow = -1;
				tmpGapNum = tmpMismatchNum = 1;
				for(j=i+1; j<overlapLen; j++)
				{
					if(alignResultArray[0][j]!='-')
					{
						if(alignResultArray[1][j]=='|' && alignResultArray[0][j]==alignResultArray[2][i] && (tmpGapNum<(*gapNum1) || tmpMismatchNum<mismatchNum))
							matchRow = j;
						break;
					}else
					{
						tmpGapNum ++;
						tmpMismatchNum ++;
					}
				}

				// exchange the rows
				if(matchRow>=0)
				{
					tmp = alignResultArray[0][i];
					alignResultArray[0][i] = alignResultArray[0][matchRow];
					alignResultArray[0][matchRow] = tmp;

					tmp = alignResultArray[1][i];
					alignResultArray[1][i] = alignResultArray[1][matchRow];
					alignResultArray[1][matchRow] = tmp;

					if(mostRightRow<i)
						mostRightRow = i;
				}else
				{
					break;
				}
			}else if(alignResultArray[2][i]=='-')
			{ // gap in seq2
				matchRow = -1;
				tmpGapNum = tmpMismatchNum = 1;
				for(j=i+1; j<overlapLen; j++)
				{
					if(alignResultArray[2][j]!='-')
					{
						if(alignResultArray[1][j]=='|' && alignResultArray[2][j]==alignResultArray[0][i] && (tmpGapNum<(*gapNum2) || tmpMismatchNum<mismatchNum))
							matchRow = j;
						break;
					}else
					{
						tmpGapNum ++;
						tmpMismatchNum ++;
					}
				}

				// exchange the rows
				if(matchRow>=0)
				{
					tmp = alignResultArray[1][i];
					alignResultArray[1][i] = alignResultArray[1][matchRow];
					alignResultArray[1][matchRow] = tmp;

					tmp = alignResultArray[2][i];
					alignResultArray[2][i] = alignResultArray[2][matchRow];
					alignResultArray[2][matchRow] = tmp;

					if(mostRightRow<i)
						mostRightRow = i;
				}else
				{
					break;
				}
			}
		}

		// check from right to left
		mostLeftRow = INT_MAX;
		for(i=overlapLen-1; i>=0; i--)
		{
			if(alignResultArray[0][i]=='-')
			{ // gap in seq1
				matchRow = -1;
				tmpGapNum = tmpMismatchNum = 1;
				for(j=i-1; j>=0; j--)
				{
					if(alignResultArray[0][j]!='-')
					{
						if(j>mostRightRow && alignResultArray[1][j]=='|' && alignResultArray[0][j]==alignResultArray[2][i] && (tmpGapNum<(*gapNum1) || tmpMismatchNum<mismatchNum))
							matchRow = j;
						break;
					}else
					{
						tmpGapNum ++;
						tmpMismatchNum ++;
					}
				}

				// exchange the rows
				if(matchRow>=0)
				{
					tmp = alignResultArray[0][i];
					alignResultArray[0][i] = alignResultArray[0][matchRow];
					alignResultArray[0][matchRow] = tmp;

					tmp = alignResultArray[1][i];
					alignResultArray[1][i] = alignResultArray[1][matchRow];
					alignResultArray[1][matchRow] = tmp;

					if(mostLeftRow>i)
						mostLeftRow = i;
				}else
				{
					break;
				}
			}else if(alignResultArray[2][i]=='-')
			{ // gap in seq2
				matchRow = -1;
				tmpGapNum = tmpMismatchNum = 1;
				for(j=i-1; j>=0; j--)
				{
					if(alignResultArray[2][j]!='-')
					{
						if(j>mostRightRow && alignResultArray[1][j]=='|' && alignResultArray[2][j]==alignResultArray[0][i] && (tmpGapNum<(*gapNum2) || tmpMismatchNum<mismatchNum))
							matchRow = j;
						break;
					}else
					{
						tmpGapNum ++;
						tmpMismatchNum ++;
					}
				}

				// exchange the rows
				if(matchRow>=0)
				{
					tmp = alignResultArray[1][i];
					alignResultArray[1][i] = alignResultArray[1][matchRow];
					alignResultArray[1][matchRow] = tmp;

					tmp = alignResultArray[2][i];
					alignResultArray[2][i] = alignResultArray[2][matchRow];
					alignResultArray[2][matchRow] = tmp;

					if(mostLeftRow>i)
						mostLeftRow = i;
				}else
				{
					break;
				}
			}
		}

		// ####################### Debug information ########################
		if(printFlag==YES)
		{
			// print the alignment result
			for(i=0; i<3; i++)
				printf("%s\n", alignResultArray[i]);
		}
		// ####################### Debug information ########################
	}

	return SUCCESSFUL;
}

/**
 * Check whether the alignment is valid, i.e. no mismatch at ends (50 bp).
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short isValidAlignment(int32_t *validHeadAlignFlag, int32_t *validTailAlignFlag, char **alignResultArray, int32_t overlapLen, int32_t queryLeftShiftLen, int32_t subjectLeftShiftLen, int32_t allowLeftShiftFlag)
{
	int32_t i, endAlignedNum;

	// head
	endAlignedNum = 0;
	for(i=0; i<overlapLen; i++)
	{
		if(alignResultArray[1][i]=='|')
			endAlignedNum ++;
		else
			break;

		if(endAlignedNum>50)
			break;
	}

	if(allowLeftShiftFlag==NO)
	{
		if(endAlignedNum>=50 && queryLeftShiftLen==0 && subjectLeftShiftLen==0)
			*validHeadAlignFlag = YES;
		else
			*validHeadAlignFlag = NO;
	}else
	{
		if(endAlignedNum>=50)
			*validHeadAlignFlag = YES;
		else
			*validHeadAlignFlag = NO;
	}

	// tail
	endAlignedNum = 0;
	for(i=overlapLen-1; i>=0; i--)
	{
		if(alignResultArray[1][i]=='|')
			endAlignedNum ++;
		else
			break;

		if(endAlignedNum>50)
			break;
	}

	if(endAlignedNum>=50)
		*validTailAlignFlag = YES;
	else
		*validTailAlignFlag = NO;

	return SUCCESSFUL;
}

/**
 * Add new queryIndel node to mis-assembly information list for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addInnerIndelToMisInfoList(query_t *queryItem, int32_t startAlignQueryPos, int32_t startAlignSubjectPos, int32_t subjectID, int32_t strand, char **alignResultArray, int32_t overlapLen, int32_t queryAlignSeqLen, int32_t subjectAlignSeqLen, int32_t queryLeftShiftLen, int32_t subjectLeftShiftLen, int32_t queryRightShiftLen, int32_t subjectRightShiftLen, int32_t gapNumQuery, int32_t gapNumSubject)
{
	int32_t i, gapStartQueryPos, gapEndQueryPos, gapStartSubjectPos, gapEndSubjectPos, queryIndelKind, difQuery, difSubject, gapFlag, unknownBaseNum;
	char *querySeq;

	if(gapNumQuery>0 || gapNumSubject>0)
	{
		if(gapNumQuery<gapNumSubject)
			queryIndelKind = QUERY_INSERT;
		else if(gapNumQuery>gapNumSubject)
			queryIndelKind = QUERY_DEL;
		else
			queryIndelKind = QUERY_GAP;

		difQuery = gapNumSubject;
		difSubject = gapNumQuery;

		// compute the margins
		for(i=0; i<overlapLen; i++)
		{
			if((alignResultArray[0][i]=='-' || alignResultArray[0][i]=='N' || alignResultArray[0][i]=='n' || alignResultArray[0][i]=='.') || (alignResultArray[2][i]=='-' || alignResultArray[2][i]=='N' || alignResultArray[2][i]=='n' || alignResultArray[2][i]=='.'))
			{
				gapStartQueryPos = startAlignQueryPos + i + queryLeftShiftLen - 1;

				if(strand==PLUS_STRAND)
					gapStartSubjectPos =  startAlignSubjectPos + i + subjectLeftShiftLen - 1;
				else
					gapStartSubjectPos = startAlignSubjectPos - i - subjectLeftShiftLen + 1;

				break;
			}
		}

		for(i=overlapLen-1; i>=0; i--)
		{
			if((alignResultArray[0][i]=='-' || alignResultArray[0][i]=='N' || alignResultArray[0][i]=='n' || alignResultArray[0][i]=='.') || (alignResultArray[2][i]=='-' || alignResultArray[2][i]=='N' || alignResultArray[2][i]=='n' || alignResultArray[2][i]=='.'))
			{
				gapEndQueryPos = startAlignQueryPos + (queryAlignSeqLen - 1) - (overlapLen - 1 - i) - queryRightShiftLen + 1;

				if(strand==PLUS_STRAND)
					gapEndSubjectPos = startAlignSubjectPos + (subjectAlignSeqLen - 1) - (overlapLen - 1 - i) - subjectRightShiftLen + 1;
				else
					gapEndSubjectPos = startAlignSubjectPos - (subjectAlignSeqLen - 1) + (overlapLen - 1 - i) + subjectRightShiftLen - 1;

				break;
			}
		}

		unknownBaseNum = 0;
		querySeq = queryItem->querySeq;
		for(i=gapStartQueryPos-1; i<gapEndQueryPos; i++)
			if(querySeq[i]=='N' || querySeq[i]=='n' || querySeq[i]=='.')
				unknownBaseNum ++;
		if(unknownBaseNum>0)
			gapFlag = YES;
		else
			gapFlag = NO;

		if(addNewQueryMisInfoNode(queryItem, -1, -1, QUERY_INDEL_KIND, queryIndelKind, gapFlag, YES, gapStartQueryPos, gapEndQueryPos, gapStartSubjectPos, gapEndSubjectPos, subjectID, difQuery, difSubject)==FAILED)
		{
			printf("line=%d, In %s(), cannot add inner query mis-assembly information node, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Add new queryIndel node to mis-assembly information list for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewQueryMisInfoNode(query_t *queryItem, int32_t leftSegRow, int32_t rightSegRow, int32_t queryMisInfoType, int32_t queryIndelKind, int32_t gapFlag, int32_t innerFlag, int32_t leftMarginQueryPos, int32_t rightMarginQueryPos, int32_t leftMarginSubjectPos, int32_t rightMarginSubjectPos, int32_t subjectID, int32_t difQuery, int32_t difSubject)
{
	misInfo_t *misInfo, *misInfoTmp, *misInfoNext;

	misInfo = (misInfo_t*) calloc (1, sizeof(misInfo_t));
	if(misInfo==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	misInfo->misType = queryMisInfoType;
	misInfo->misassFlag = UNUSED_MISASS;
	misInfo->leftSegRow = leftSegRow;
	misInfo->rightSegRow = rightSegRow;
	misInfo->queryMargin = NULL;
	misInfo->queryIndel = NULL;
	misInfo->gapFlag = gapFlag;
	misInfo->innerFlag = innerFlag;
	misInfo->misassSeqList = misInfo->tailMisassSeq = NULL;
	misInfo->misassSeqNum = 0;
	misInfo->next = NULL;

	if(queryMisInfoType==QUERY_MISJOIN_KIND)
	{ // misjoin
		// allocate memory
		misInfo->queryMargin = (queryMargin_t*) calloc (1, sizeof(queryMargin_t));
		if(misInfo->queryMargin==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		misInfo->queryMargin->leftMargin = leftMarginQueryPos;
		misInfo->queryMargin->rightMargin = rightMarginQueryPos;
	}else
	{ // indel
		// allocate memory
		misInfo->queryIndel = (queryIndel_t*) calloc (1, sizeof(queryIndel_t));
		if(misInfo->queryIndel==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		misInfo->queryIndel->leftSegRow = leftSegRow;
		misInfo->queryIndel->rightSegRow = rightSegRow;
		misInfo->queryIndel->queryIndelKind = queryIndelKind;

		misInfo->queryIndel->leftMargin = leftMarginQueryPos;
		misInfo->queryIndel->rightMargin = rightMarginQueryPos;
		misInfo->queryIndel->leftMarginSubject = leftMarginSubjectPos;
		misInfo->queryIndel->rightMarginSubject = rightMarginSubjectPos;
		misInfo->queryIndel->subjectID = subjectID;
		misInfo->queryIndel->difQuery = difQuery;
		misInfo->queryIndel->difSubject = difSubject;
	}

	if(queryItem->misInfoItemNum==0)
	{
		queryItem->misInfoList = misInfo;
		queryItem->tailMisInfo = misInfo;
	}else
	{
		if(innerFlag==NO)
		{
			queryItem->tailMisInfo->next = misInfo;
			queryItem->tailMisInfo = misInfo;
		}else
		{
			misInfoTmp = queryItem->misInfoList;
			misInfoNext = NULL;
			while(misInfoTmp)
			{
				if(misInfoTmp->misType==QUERY_MISJOIN_KIND)
				{
					if(misInfoTmp->queryMargin->leftMargin>rightMarginQueryPos)
					{
						misInfoNext = misInfoTmp;
						break;
					}
				}else
				{
					if(misInfoTmp->queryIndel->leftMargin>rightMarginQueryPos)
					{
						misInfoNext = misInfoTmp;
						break;
					}
				}

				misInfoTmp = misInfoTmp->next;
			}

			// insert the node to the proper place
			misInfoTmp = queryItem->misInfoList;
			while(misInfoTmp)
			{
				if(misInfoTmp->next==misInfoNext)
				{
					misInfo->next = misInfoNext;
					misInfoTmp->next = misInfo;

					if(misInfoTmp==queryItem->tailMisInfo)
						queryItem->tailMisInfo = misInfo;

					break;
				}else if(misInfoNext==queryItem->misInfoList)
				{
					misInfo->next = misInfoNext;
					queryItem->misInfoList = misInfo;
					break;
				}

				misInfoTmp = misInfoTmp->next;
			}
		}
	}
	queryItem->misInfoItemNum ++;

	return SUCCESSFUL;
}

/**
 * Determine mis-assembly information for each node in the list.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short determineMisInfoSingleQuery(query_t *queryItem, subject_t *subjectArray, readSet_t *readSet)
{
	int32_t queryLen;
	baseCov_t *baseCovArray;
//	char covFileName[256], covReadsFileName[256], regCovFileName[256];

//	strcpy(covFileName, outputPathStr);
//	strcat(covFileName, queryItem->queryTitle);
//	strcat(covFileName, ".cov");

//	strcpy(regCovFileName, outputPathStr);
//	strcat(regCovFileName, queryItem->queryTitle);
//	strcat(regCovFileName, ".regcov");

//	strcpy(covReadsFileName, outputPathStr);
//	strcat(covReadsFileName, queryItem->queryTitle);
//	strcat(covReadsFileName, "_covReads");

	queryLen = queryItem->queryLen;
	baseCovArray = (baseCov_t *) calloc (queryLen, sizeof(baseCov_t));
	if(baseCovArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// compute base coverage
	if(computeBaseCovSingleQuery(baseCovArray, queryItem, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the base coverage of single query, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// compute the region coverage
	if(computeRegCovSingleQuery(queryItem, baseCovArray, queryLen)==FAILED)
	{
		printf("line=%d, In %s(), cannot output the base coverage of single query to file, error!\n", __LINE__, __func__);
		return FAILED;
	}

//	// output the region coverage to file
//	if(outputRegCovSingleQueryToFile(regCovFileName, baseCovArray, queryLen)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot output the base coverage of single query to file, error!\n", __LINE__, __func__);
//		return FAILED;
//	}

	// output the coverage to file
//	if(outputBaseCovSingleQueryToFile(covFileName, baseCovArray, queryLen)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot output the base coverage of single query to file, error!\n", __LINE__, __func__);
//		return FAILED;
//	}

//	if(queryItem->queryID==5 || strcmp(queryItem->queryTitle, "ctg7180000002352")==0)
//	{
//		if(outputCovReadsToFile(covReadsFileName, 10421, 'T', queryItem, readSet)==FAILED)
//		{
//			printf("line=%d, In %s(), cannot output the covered reads, error!\n", __LINE__, __func__);
//			return FAILED;
//		}
//	}

	// determine misjoin
	if(computeMisjoinSingleQuery(queryItem, baseCovArray, subjectArray, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the misjoin, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// determine indel
	if(determineQueryIndel(queryItem, baseCovArray, subjectArray, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot determine the queryIndel information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(determineQueryMis(queryItem)==FAILED)
	{
		printf("line=%d, In %s(), cannot determine the query mis-assembly, error!\n", __LINE__, __func__);
		return FAILED;
	}

	free(baseCovArray);

	return SUCCESSFUL;
}

/**
 * Determine misjoin for each node in the list.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeMisjoinSingleQuery(query_t *queryItem, baseCov_t *baseCovArray, subject_t *subjectArray, readSet_t *readSet)
{
	int32_t misjoinRegNum;
	misInfo_t *misInfo;

	misjoinRegNum = 0;
	misInfo = queryItem->misInfoList;
	while(misInfo)
	{
		if(misInfo->misType==QUERY_MISJOIN_KIND)
			misjoinRegNum ++;
		misInfo = misInfo->next;
	}

	if(misjoinRegNum>0)
	{
		// compute the discorRatio, multiRatio in breakpoint regions
		if(computeBreakpointRatios(queryItem, misjoinRegNum, readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the ratios of normal regions, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// determine the mis-assembly
		if(determineMisassFlag(queryItem, misjoinRegNum, baseCovArray, readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot determine the mis-assembly flag of query, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// adjust breakpoint margins


	}

	return SUCCESSFUL;
}

/**
 * Determine query mis-assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short determineQueryMis(query_t *queryItem)
{
	misInfo_t *misInfo;

	queryItem->misassFlag = UNCERTAIN_MISASS;

	misInfo = queryItem->misInfoList;
	while(misInfo)
	{
		if(misInfo->misassFlag==TRUE_MISASS)
		{
			if(misInfo->misType==QUERY_MISJOIN_KIND || (misInfo->misType==QUERY_INDEL_KIND && misInfo->gapFlag==NO))
			{
				queryItem->misassFlag = TRUE_MISASS;
				break;
			}
		}
		misInfo = misInfo->next;
	}

	return SUCCESSFUL;
}

/**
 * Output the mis-assembly information of single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputMisInfoList(query_t *queryItem)
{
	misInfo_t *tailMisInfo;
	queryMargin_t *queryMargin;
	queryIndel_t *queryIndel;
	int32_t i;

	i = 0;
	tailMisInfo = queryItem->misInfoList;
	while(tailMisInfo)
	{
		if(tailMisInfo->misType==QUERY_MISJOIN_KIND)
		{
			queryMargin = tailMisInfo->queryMargin;
			printf("reg[%d]: [%d, %d], leftRow=%d, rightRow=%d, misType=%d, misassFlag=%d, gapFlag=%d, innerFlag=%d, endFlag=%d, disagreeNum=%d, zeroCovNum=%d, discorNum=%d, highCovRegNum=%d, lowCovRegNum=%d, discorRatio=%.4f, multiRatio=%.4f\n", i, queryMargin->leftMargin, queryMargin->rightMargin, tailMisInfo->leftSegRow, tailMisInfo->rightSegRow, tailMisInfo->misType, tailMisInfo->misassFlag, tailMisInfo->gapFlag, tailMisInfo->innerFlag, queryMargin->queryEndFlag, queryMargin->disagreeNum, queryMargin->zeroCovNum, queryMargin->discorNum, queryMargin->highCovRegNum, queryMargin->lowCovRegNum, queryMargin->discorRatio, queryMargin->multiReadsRatio);
		}else
		{
			queryIndel = tailMisInfo->queryIndel;
			printf("reg[%d]: [%d, %d], leftRow=%d, rightRow=%d, misType=%d, misassFlag=%d, gapFlag=%d, innerFlag=%d, indelKind=%d, disagreeNum=%d, zeroCovNum=%d, highCovRegNum=%d, lowCovRegNum=%d, difQuery=%d, difSubject=%d, difFragSizeLeft=%.4f, difFragSizeRight=%.4f, averFragSizeLeft=%.4f, averFragSizeRight=%.4f, pairNumLeft=%d, pairNumRight=%d, discorRatioLeft=%.4f, discorRatioRight=%.4f, discorNumLeft=%d, discorNumRight=%d\n", i, queryIndel->leftMargin, queryIndel->rightMargin, tailMisInfo->leftSegRow, tailMisInfo->rightSegRow, tailMisInfo->misType, tailMisInfo->misassFlag, tailMisInfo->gapFlag, tailMisInfo->innerFlag, queryIndel->queryIndelKind, queryIndel->disagreeNum, queryIndel->zeroCovNum, queryIndel->highCovRegNum, queryIndel->lowCovRegNum, queryIndel->difQuery, queryIndel->difSubject, queryIndel->difFragSizeLeft, queryIndel->difFragSizeRight, queryIndel->averFragSizeLeft, queryIndel->averFragSizeRight, queryIndel->pairNumLeft, queryIndel->pairNumRight, queryIndel->discorRatioLeft, queryIndel->discorRatioRight, queryIndel->discorNumLeft, queryIndel->discorNumRight);
		}

		i ++;
		tailMisInfo = tailMisInfo->next;
	}

	return SUCCESSFUL;
}
