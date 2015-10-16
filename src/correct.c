/*
 * correct.c
 *
 *  Created on: Jul 20, 2014
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Correct mis-assembled queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short correctMisassQueries(queryMatchInfo_t *queryMatchInfoSet)
{
	int32_t i, processedNum;
	query_t *queryArray;

	printf("Begin correcting the mis-assemblies at their breakpoints ...\n");

	processedNum = 0;
	queryArray = queryMatchInfoSet->queryArray;
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		if(queryArray[i].misassFlag==TRUE_MISASS)
		{
			// ########################### Debug information ##############################
			//if(queryArray[i].queryID==71 || strcmp(queryArray[i].queryTitle, "ctg7180000002451")==0)
			//{
			//	printf("###### queryID=%d, queryTitle=%s, queryLen=%d, subjectNum=%d\n", queryArray[i].queryID, queryArray[i].queryTitle, queryArray[i].queryLen, queryArray[i].querySubjectNum);
			//}
			// ########################### Debug information ##############################

			// correct the mis-assemblies in queries
			if(correctSingleMisassQuery(queryArray+i, queryMatchInfoSet->subjectArray)==FAILED)
			{
				printf("line=%d, In %s(), cannot correct single potential mis-assembled query, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// ############# Debug information ###############
			//outputNewSeqInfo(queryArray+i);
			// ############# Debug information ###############
		}

		processedNum ++;
//		if(processedNum%100==0)
//			printf("Queries processed: %d\n", processedNum);
	}

//	if(processedNum%100!=0)
//		printf("Queries processed: %d\n", processedNum);

	return SUCCESSFUL;
}

/**
 * Correct single potential mis-assembled query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short correctSingleMisassQuery(query_t *queryItem, subject_t *subjectArray)
{
	int32_t newSeqLen, newSeqFlag;
	char *newSeq;
	misInfo_t *misInfo, *misInfoTmp, *startMisInfo, *endMisInfo, *processedMisInfo;

	newSeqFlag = NO;
	misInfo = queryItem->misInfoList;
	while(misInfo)
	{
		if(misInfo->misassFlag==TRUE_MISASS && misInfo->misType==QUERY_MISJOIN_KIND)
		{
			newSeqFlag = YES;
			break;
		}
		misInfo = misInfo->next;
	}

	if(newSeqFlag==NO)
	{
		queryItem->newQuerySeqList = NULL;
		queryItem->newQuerySeqNum = 0;
		return SUCCESSFUL;
	}

	processedMisInfo = NULL;
	misInfo = queryItem->misInfoList;
	while(misInfo)
	{
		// get start and end misInfo
		endMisInfo = NULL;
		startMisInfo = misInfo;
		misInfoTmp = startMisInfo;
		while(misInfoTmp)
		{
			if(misInfoTmp!=processedMisInfo && misInfoTmp->misassFlag==TRUE_MISASS && misInfoTmp->misType==QUERY_MISJOIN_KIND)
			{
				endMisInfo = misInfoTmp;
				break;
			}
			misInfoTmp = misInfoTmp->next;
		}

		// compute new sequence length
		if(computeNewSeqLen(&newSeqLen, startMisInfo, endMisInfo, processedMisInfo, queryItem)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute new sequence length, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// generate new sequence node
		if(generateNewSeqNode(queryItem, &newSeq, newSeqLen, startMisInfo, endMisInfo, processedMisInfo)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the new sequence nodes, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// fill the bases
		if(fillNewBases(queryItem, newSeq, newSeqLen, startMisInfo, endMisInfo, processedMisInfo, subjectArray)==FAILED)
		{
			printf("line=%d, In %s(), cannot fill the new sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}

		processedMisInfo = endMisInfo;
		misInfo = endMisInfo;
	}

	return SUCCESSFUL;
}

/**
 * Get start and end query positions.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getStartEndQueryPosMisCor(int32_t *startQueryPos, int32_t *endQneryPos, misInfo_t *startMisInfo, misInfo_t *endMisInfo, misInfo_t *processedMisInfo, query_t *queryItem)
{
	int32_t startQPos, endQPos, nullFlag, leftSegRow, rightSegRow, globalSegNum;
	globalValidSeg_t *globalSegArray;

	globalSegArray = queryItem->globalValidSegArray;
	globalSegNum = queryItem->globalValidSegNum;

	if(startMisInfo->misType==QUERY_MISJOIN_KIND)
	{
		if(startMisInfo==endMisInfo)
			startQPos = 1;
		else
		{
			if(startMisInfo==processedMisInfo)
				startQPos = startMisInfo->queryMargin->leftMargin;
			else
				//startQPos = globalSegArray[0].startQueryPos;
				startQPos = 1;
		}
	}else
	{
		leftSegRow = startMisInfo->leftSegRow;
		rightSegRow = startMisInfo->rightSegRow;
		if(leftSegRow==-1)
		{
			if(rightSegRow>=0 && startMisInfo->misassFlag==TRUE_MISASS)
				//startQPos = globalSegArray[rightSegRow].startQueryPos;
				startQPos = 1;
			else
				startQPos = 1;
		}else
		{
			if(startMisInfo==queryItem->misInfoList)
				//startQPos = globalSegArray[0].startQueryPos;
				startQPos = 1;
			else
				startQPos = globalSegArray[leftSegRow].startQueryPos;
		}
	}

	if(endMisInfo==NULL)
	{
		endMisInfo = queryItem->tailMisInfo;
		nullFlag = YES;
	}else
		nullFlag = NO;

	if(endMisInfo->misType==QUERY_MISJOIN_KIND)
	{
		if(nullFlag==NO)
			endQPos = endMisInfo->queryMargin->rightMargin;
		else
			endQPos = queryItem->queryLen;
	}else
	{
		leftSegRow = endMisInfo->leftSegRow;
		rightSegRow = endMisInfo->rightSegRow;
		if(rightSegRow==-1)
		{
			if(leftSegRow>=0 && endMisInfo->misassFlag==TRUE_MISASS)
				//endQPos = globalSegArray[leftSegRow].endQueryPos;
				endQPos = queryItem->queryLen;
			else
				endQPos = queryItem->queryLen;
		}else
		{
			if(endMisInfo->next)
				endQPos = globalSegArray[rightSegRow].endQueryPos;
			else
				//endQPos = globalSegArray[globalSegNum-1].endQueryPos;
				endQPos = queryItem->queryLen;
		}
	}

	*startQueryPos = startQPos;
	*endQneryPos = endQPos;

	return SUCCESSFUL;
}

/**
 * Compute new sequence length.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeNewSeqLen(int32_t *newSeqLen, misInfo_t *startMisInfo, misInfo_t *endMisInfo, misInfo_t *processedMisInfo, query_t *queryItem)
{
	int32_t startQPos, endQPos, nullFlag, leftSegRow, rightSegRow, globalSegNum, newLen, insLen, rmLen, newStartQPos, newStartSPos, distSubject, distQuery;
	misInfo_t *misInfo;
	queryIndel_t *queryIndel;
	globalValidSeg_t *globalSegArray;

	if(getStartEndQueryPosMisCor(&startQPos, &endQPos, startMisInfo, endMisInfo, processedMisInfo, queryItem)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the start and end position for new query sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}
	newLen = endQPos - startQPos + 1;

/*
	misInfo = startMisInfo;
	while(misInfo)
	{
		if(misInfo->misassFlag==TRUE_MISASS && misInfo->misType==QUERY_INDEL_KIND && misInfo->gapFlag==NO)
		{
			queryIndel = misInfo->queryIndel;
			leftSegRow = misInfo->leftSegRow;
			rightSegRow = misInfo->rightSegRow;

			if(leftSegRow>=0 && rightSegRow>=0)
			{
				// ????????????????????
				if(globalSegArray[leftSegRow].subjectID!=globalSegArray[rightSegRow].subjectID || globalSegArray[leftSegRow].strand!=globalSegArray[rightSegRow].strand)
				{
					printf("line=%d, In %s(), subjectID1=%d, subjectID2=%d, strand1=%d, strand2=%d, error!\n", __LINE__, __func__, globalSegArray[leftSegRow].subjectID, globalSegArray[rightSegRow].subjectID, globalSegArray[leftSegRow].strand, globalSegArray[rightSegRow].strand);
					return FAILED;
				}

				if(queryIndel->queryIndelKind==QUERY_INSERT)
				{
					if(globalSegArray[leftSegRow].strand==PLUS_STRAND)
					{ // plus strand
						if(globalSegArray[leftSegRow].endSubPos<globalSegArray[rightSegRow].startSubPos)
						{
							newStartQPos = globalSegArray[rightSegRow].startQueryPos + (globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos) + 1;
							rmLen = newStartQPos - globalSegArray[leftSegRow].endQueryPos -1;
						}else if(globalSegArray[leftSegRow].endSubPos<=globalSegArray[rightSegRow].endSubPos)
						{
							newStartQPos = globalSegArray[rightSegRow].startQueryPos + (globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos) + 1;
							rmLen = newStartQPos - globalSegArray[leftSegRow].endQueryPos -1;
						}else
						{
							//newStartQPos = globalSegArray[rightSegRow].startQueryPos + (globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos) + 1;
							//rmLen = globalSegArray[rightSegRow].endQueryPos - globalSegArray[leftSegRow].endQueryPos;

							printf("line=%d, In %s(), invalid query insertion, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}else
					{ // minus strand
						if(globalSegArray[leftSegRow].endSubPos>globalSegArray[rightSegRow].startSubPos)
						{
							newStartQPos = globalSegArray[rightSegRow].startQueryPos - (globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos) + 1;
							rmLen = newStartQPos - globalSegArray[leftSegRow].endQueryPos -1;
						}else if(globalSegArray[leftSegRow].endSubPos>=globalSegArray[rightSegRow].endSubPos)
						{
							newStartQPos = globalSegArray[rightSegRow].startQueryPos - (globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos) + 1;
							rmLen = newStartQPos - globalSegArray[leftSegRow].endQueryPos -1;
						}else
						{
							//newStartQPos = globalSegArray[rightSegRow].startQueryPos - (globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos) + 1;
							//rmLen = globalSegArray[rightSegRow].endQueryPos - globalSegArray[leftSegRow].endQueryPos;

							printf("line=%d, In %s(), invalid query insertion, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}

					newLen -= rmLen;

				}else if(queryIndel->queryIndelKind==QUERY_DEL)
				{
					if(globalSegArray[leftSegRow].strand==PLUS_STRAND)
					{ // plus strand
						if(globalSegArray[leftSegRow].endSubPos<globalSegArray[rightSegRow].startSubPos)
						{
							newStartSPos = globalSegArray[rightSegRow].startSubPos - (globalSegArray[rightSegRow].startQueryPos - globalSegArray[leftSegRow].endQueryPos) + 1;
							insLen = newStartSPos - globalSegArray[leftSegRow].endSubPos - 1;
						}else if(globalSegArray[leftSegRow].endSubPos<=globalSegArray[rightSegRow].endSubPos)
						{
							newStartSPos = globalSegArray[rightSegRow].startSubPos - (globalSegArray[rightSegRow].startQueryPos - globalSegArray[leftSegRow].endQueryPos) + 1;
							insLen = newStartSPos - globalSegArray[leftSegRow].endSubPos - 1;
							newStartQPos = globalSegArray[rightSegRow].startQueryPos + (globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos) + 1;
						}else
						{
							printf("line=%d, In %s(), invalid query deletion, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}else
					{ // minus strand
						if(globalSegArray[leftSegRow].endSubPos>globalSegArray[rightSegRow].startSubPos)
						{
							newStartSPos = globalSegArray[rightSegRow].startSubPos + (globalSegArray[rightSegRow].startQueryPos - globalSegArray[leftSegRow].endQueryPos) - 1;
							insLen = globalSegArray[leftSegRow].endSubPos - newStartSPos - 1;
						}else if(globalSegArray[leftSegRow].endSubPos>=globalSegArray[rightSegRow].endSubPos)
						{
							newStartSPos = globalSegArray[rightSegRow].startSubPos + (globalSegArray[rightSegRow].startQueryPos - globalSegArray[leftSegRow].endQueryPos) - 1;
							insLen = globalSegArray[leftSegRow].endSubPos - newStartSPos - 1;
							newStartQPos = globalSegArray[rightSegRow].startQueryPos + (globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos) + 1;
						}else
						{
							printf("line=%d, In %s(), invalid query deletion, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}

					newLen += insLen;
				}else if(queryIndel->queryIndelKind==QUERY_GAP)
				{
					if(globalSegArray[leftSegRow].strand==PLUS_STRAND)
					{ // plus strand
						if(globalSegArray[leftSegRow].endSubPos<globalSegArray[rightSegRow].startSubPos)
						{
							distSubject = globalSegArray[rightSegRow].startSubPos - globalSegArray[leftSegRow].endSubPos;
							distQuery = globalSegArray[rightSegRow].startQueryPos - globalSegArray[leftSegRow].endQueryPos;
							newLen += distSubject - distQuery;
						}else
						{
							if(globalSegArray[leftSegRow].endSubPos<=globalSegArray[rightSegRow].endSubPos)
							{
								newStartQPos = globalSegArray[rightSegRow].startQueryPos + (globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos) + 1;
								rmLen = newStartQPos - globalSegArray[leftSegRow].endQueryPos - 1;
							}else
							{
								//newStartQPos = globalSegArray[rightSegRow].startQueryPos + (globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos) + 1;
								//rmLen = globalSegArray[rightSegRow].endQueryPos - globalSegArray[leftSegRow].endQueryPos;

								printf("line=%d, In %s(), invalid query insertion, error!\n", __LINE__, __func__);
								return FAILED;
							}
							newLen -= rmLen;
						}
					}else
					{ // minus strand
						if(globalSegArray[leftSegRow].endSubPos>globalSegArray[rightSegRow].startSubPos)
						{
							distSubject = globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos;
							distQuery = globalSegArray[rightSegRow].startQueryPos - globalSegArray[leftSegRow].endQueryPos;
							newLen += distSubject - distQuery;
						}else
						{
							if(globalSegArray[leftSegRow].endSubPos>=globalSegArray[rightSegRow].endSubPos)
							{
								newStartQPos = globalSegArray[rightSegRow].startQueryPos - (globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos) + 1;
								rmLen = newStartQPos - globalSegArray[leftSegRow].endQueryPos - 1;
							}else
							{
								//newStartQPos = globalSegArray[rightSegRow].startQueryPos - (globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos) + 1;
								//rmLen = globalSegArray[rightSegRow].endQueryPos - globalSegArray[leftSegRow].endQueryPos;

								printf("line=%d, In %s(), invalid query insertion, error!\n", __LINE__, __func__);
								return FAILED;
							}
							newLen -= rmLen;
						}
					}
				}else
				{
					printf("line=%d, In %s(), invalid query indel kind, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}

		if(misInfo==endMisInfo)
			break;
		else
			misInfo = misInfo->next;
	}
*/

	*newSeqLen = newLen;
	if((*newSeqLen)<=0)
	{
		printf("queryTitle=%s, startQPos=%d, endQPos=%d, newSeqLen=%d\n", queryItem->queryTitle, startQPos, endQPos, *newSeqLen);
		return FAILED;
	}

	//printf("queryTitle=%s, startQPos=%d, endQPos=%d, oldQLen=%d, newSeqLen=%d\n", queryItem->queryTitle, startQPos, endQPos, endQPos-startQPos+1, *newSeqLen);

	return SUCCESSFUL;
}

/**
 * Generate new sequence node.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateNewSeqNode(query_t *queryItem, char **newSeq, int32_t newSeqLen, misInfo_t *startMisInfo, misInfo_t *endMisInfo, misInfo_t *processedMisInfo)
{
	querySeq_t *querySeqNode;
	globalValidSeg_t *globalSegArray;
	int32_t leftSegRow, rightSegRow, globalSegNum, nullFlag;

	*newSeq = (char *) calloc (newSeqLen+1, sizeof(char));
	if((*newSeq)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	querySeqNode = (querySeq_t *) calloc(1, sizeof(querySeq_t));
	if(querySeqNode==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	querySeqNode->queryLen = newSeqLen;
	querySeqNode->querySeq = *newSeq;
	querySeqNode->next = NULL;

	// set the start query position of new querySeq
	if(getStartEndQueryPosMisCor(&querySeqNode->startQueryPos, &querySeqNode->endQueryPos, startMisInfo, endMisInfo, processedMisInfo, queryItem)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the start and end position for new query sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// append the node
	if(queryItem->newQuerySeqList==NULL)
	{
		queryItem->newQuerySeqList = querySeqNode;
		queryItem->tailNewQuerySeqList = querySeqNode;
	}else
	{
		queryItem->tailNewQuerySeqList->next = querySeqNode;
		queryItem->tailNewQuerySeqList = querySeqNode;
	}
	queryItem->newQuerySeqNum ++;

	return SUCCESSFUL;
}

/**
 * Fill new bases for single query.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillNewBases(query_t *queryItem, char *newSeq, int32_t newSeqLen, misInfo_t *startMisInfo, misInfo_t *endMisInfo, misInfo_t *processedMisInfo, subject_t *subjectArray)
{
	int32_t startQPos, endQPos, nullFlag, leftSegRow, rightSegRow, newLen, preQueryPos, newStartSPos;
	int32_t subjectID, globalSegNum, queryLen, subLen;
	misInfo_t *misInfo;
	queryMargin_t *queryMargin;
	queryIndel_t *queryIndel;
	globalValidSeg_t *globalSegArray;
	char *querySeq, *subjectSeq;

	querySeq = queryItem->querySeq;
	queryLen = queryItem->queryLen;

	// set the start query position of new querySeq
	if(getStartEndQueryPosMisCor(&startQPos, &endQPos, startMisInfo, endMisInfo, processedMisInfo, queryItem)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the start and end position for new query sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	newLen = endQPos - startQPos + 1;
	strncpy(newSeq, querySeq+startQPos-1, newLen);
	newSeq[newLen] = '\0';

/*
	newLen = 0;
	preQueryPos = -1;
	misInfo = startMisInfo;
	while(misInfo)
	{
		if(misInfo->misassFlag==TRUE_MISASS)
		{
			if(misInfo->misType==QUERY_MISJOIN_KIND)
			{
				queryMargin = misInfo->queryMargin;

				if(preQueryPos==-1)
					preQueryPos = startQPos;

				if(nullFlag==NO)
				{
					subLen = queryMargin->rightMargin - preQueryPos + 1;
					if(subLen>0)
					{
						strncpy(newSeq+newLen, querySeq+preQueryPos-1, subLen);
						newSeq[newLen+subLen] = '\0';
						newLen += subLen;
					}else if(subLen<0)
					{
						printf("line=%d, In %s(), invalid query indel, error!\n", __LINE__, __func__);
						return FAILED;
					}

					preQueryPos = queryMargin->leftMargin;
				}
			}else if(misInfo->gapFlag==NO)
			{
				queryIndel = misInfo->queryIndel;
				leftSegRow = queryIndel->leftSegRow;
				rightSegRow = queryIndel->rightSegRow;

				if(leftSegRow==-1)
				{
					if(preQueryPos==-1)
						preQueryPos = startQPos;
					else
					{
						printf("line=%d, In %s(), preQueryPos=%d, error!\n", __LINE__, __func__, preQueryPos);
						return FAILED;
					}
				}else if(rightSegRow==-1)
				{
					if(preQueryPos==-1)
						preQueryPos = startQPos;

					subLen = globalSegArray[leftSegRow].endQueryPos - preQueryPos + 1;
					if(subLen>0)
					{
						strncpy(newSeq+newLen, querySeq+preQueryPos-1, subLen);
						newSeq[newLen+subLen] = '\0';
						newLen += subLen;
					}else if(subLen<0)
					{
						printf("line=%d, In %s(), invalid query indel, error!\n", __LINE__, __func__);
						return FAILED;
					}

					preQueryPos = globalSegArray[leftSegRow].endQueryPos + 1;
				}else
				{
					if(preQueryPos==-1)
						preQueryPos = startQPos;

					subLen = globalSegArray[leftSegRow].endQueryPos - preQueryPos + 1;
					if(subLen>0)
					{
						strncpy(newSeq+newLen, querySeq+preQueryPos-1, subLen);
						newSeq[newLen+subLen] = '\0';
						newLen += subLen;
					}else if(subLen<0)
					{
						printf("line=%d, In %s(), invalid query indel, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(queryIndel->queryIndelKind==QUERY_INSERT)
					{
						if(globalSegArray[leftSegRow].strand==PLUS_STRAND)
						{ // plus strand
							if(globalSegArray[leftSegRow].endSubPos<globalSegArray[rightSegRow].startSubPos)
							{
								subLen = globalSegArray[rightSegRow].startSubPos - globalSegArray[leftSegRow].endSubPos - 1;
								if(subLen>0)
								{
									subjectID = globalSegArray[leftSegRow].subjectID;
									subjectSeq = subjectArray[subjectID-1].subjectSeq;
									newStartSPos = globalSegArray[leftSegRow].endSubPos + 1;
									strncpy(newSeq+newLen, subjectSeq+newStartSPos-1, subLen);
									newSeq[newLen+subLen] = '\0';
									newLen += subLen;
								}

								preQueryPos = globalSegArray[rightSegRow].startQueryPos;
							}else if(globalSegArray[leftSegRow].endSubPos<=globalSegArray[rightSegRow].endSubPos)
							{
								preQueryPos = globalSegArray[rightSegRow].startQueryPos + (globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos) + 1;
							}else
							{
								printf("line=%d, In %s(), invalid query insertion, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}else
						{ // minus strand
							if(globalSegArray[leftSegRow].endSubPos>globalSegArray[rightSegRow].startSubPos)
							{
								subLen = globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos - 1;
								if(subLen>0)
								{
									subjectID = globalSegArray[rightSegRow].subjectID;
									subjectSeq = subjectArray[subjectID-1].subjectSeq;
									newStartSPos = globalSegArray[rightSegRow].startSubPos + 1;
									strncpy(newSeq+newLen, subjectSeq+newStartSPos-1, subLen);
									newSeq[newLen+subLen] = '\0';
									if(reverseSeq(newSeq+newLen, subLen)==FAILED)
									{
										printf("line=%d, In %s(), cannot reverse the base sequence, error!\n", __LINE__, __func__);
										return FAILED;
									}
									newLen += subLen;
								}

								preQueryPos = globalSegArray[rightSegRow].startQueryPos;
							}else if(globalSegArray[leftSegRow].endSubPos>=globalSegArray[rightSegRow].endSubPos)
							{
								preQueryPos = globalSegArray[rightSegRow].startQueryPos + (globalSegArray[rightSegRow].startSubPos - globalSegArray[leftSegRow].endSubPos) + 1;
							}else
							{
								printf("line=%d, In %s(), invalid query insertion, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}
					}else if(queryIndel->queryIndelKind==QUERY_DEL)
					{
						if(globalSegArray[leftSegRow].strand==PLUS_STRAND)
						{ // plus strand
							if(globalSegArray[leftSegRow].endSubPos<globalSegArray[rightSegRow].startSubPos)
							{
								subLen = globalSegArray[rightSegRow].startSubPos - globalSegArray[leftSegRow].endSubPos - 1;
								if(subLen>0)
								{
									subjectID = globalSegArray[leftSegRow].subjectID;
									subjectSeq = subjectArray[subjectID-1].subjectSeq;
									newStartSPos = globalSegArray[leftSegRow].endSubPos + 1;
									strncpy(newSeq+newLen, subjectSeq+newStartSPos-1, subLen);
									newSeq[newLen+subLen] = '\0';
									newLen += subLen;
								}
								preQueryPos = globalSegArray[rightSegRow].startQueryPos;
							}else if(globalSegArray[leftSegRow].endSubPos<=globalSegArray[rightSegRow].endSubPos)
							{
								preQueryPos = globalSegArray[rightSegRow].startQueryPos + (globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos) + 1;
							}else
							{
								printf("line=%d, In %s(), invalid query deletion, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}else
						{ // minus strand
							if(globalSegArray[leftSegRow].endSubPos>globalSegArray[rightSegRow].startSubPos)
							{
								subLen = globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos - 1;
								if(subLen>0)
								{
									subjectID = globalSegArray[rightSegRow].subjectID;
									subjectSeq = subjectArray[subjectID-1].subjectSeq;
									newStartSPos = globalSegArray[rightSegRow].startSubPos + 1;
									strncpy(newSeq+newLen, subjectSeq+newStartSPos-1, subLen);
									newSeq[newLen+subLen] = '\0';
									if(reverseSeq(newSeq+newLen, subLen)==FAILED)
									{
										printf("line=%d, In %s(), cannot reverse the base sequence, error!\n", __LINE__, __func__);
										return FAILED;
									}
									newLen += subLen;
								}
								preQueryPos = globalSegArray[rightSegRow].startQueryPos;
							}else if(globalSegArray[leftSegRow].endSubPos>=globalSegArray[rightSegRow].endSubPos)
							{
								preQueryPos = globalSegArray[rightSegRow].startQueryPos + (globalSegArray[rightSegRow].startSubPos - globalSegArray[leftSegRow].endSubPos) + 1;
							}else
							{
								printf("line=%d, In %s(), invalid query deletion, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}
					}else if(queryIndel->queryIndelKind==QUERY_GAP)
					{
						if(globalSegArray[leftSegRow].strand==PLUS_STRAND)
						{ // plus strand
							if(globalSegArray[leftSegRow].endSubPos<globalSegArray[rightSegRow].startSubPos)
							{
								subLen = globalSegArray[rightSegRow].startSubPos - globalSegArray[leftSegRow].endSubPos - 1;
								if(subLen>0)
								{
									subjectID = globalSegArray[leftSegRow].subjectID;
									subjectSeq = subjectArray[subjectID-1].subjectSeq;
									newStartSPos = globalSegArray[leftSegRow].endSubPos + 1;
									strncpy(newSeq+newLen, subjectSeq+newStartSPos-1, subLen);
									newSeq[newLen+subLen] = '\0';
									newLen += subLen;
								}
								preQueryPos = globalSegArray[rightSegRow].startQueryPos;
							}else if(globalSegArray[leftSegRow].endSubPos<=globalSegArray[rightSegRow].endSubPos)
							{
								preQueryPos = globalSegArray[rightSegRow].startQueryPos + (globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos) + 1;
							}else
							{
								printf("line=%d, In %s(), invalid query gap, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}else
						{ // minus strand
							if(globalSegArray[leftSegRow].endSubPos>globalSegArray[rightSegRow].startSubPos)
							{
								subLen = globalSegArray[leftSegRow].endSubPos - globalSegArray[rightSegRow].startSubPos - 1;
								if(subLen>0)
								{
									subjectID = globalSegArray[rightSegRow].subjectID;
									subjectSeq = subjectArray[subjectID-1].subjectSeq;
									newStartSPos = globalSegArray[rightSegRow].startSubPos + 1;
									strncpy(newSeq+newLen, subjectSeq+newStartSPos-1, subLen);
									newSeq[newLen+subLen] = '\0';
									if(reverseSeq(newSeq+newLen, subLen)==FAILED)
									{
										printf("line=%d, In %s(), cannot reverse the base sequence, error!\n", __LINE__, __func__);
										return FAILED;
									}
									newLen += subLen;
								}
								preQueryPos = globalSegArray[rightSegRow].startQueryPos;
							}else if(globalSegArray[leftSegRow].endSubPos>=globalSegArray[rightSegRow].endSubPos)
							{
								preQueryPos = globalSegArray[rightSegRow].startQueryPos + (globalSegArray[rightSegRow].startSubPos - globalSegArray[leftSegRow].endSubPos) + 1;
							}else
							{
								printf("line=%d, In %s(), invalid query gap, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}
					}else
					{
						printf("line=%d, In %s(), invalid query indel kind, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

			}
		}

		if(misInfo->next==NULL)
		{
			if(misInfo->misType==QUERY_MISJOIN_KIND)
			{
				queryMargin = misInfo->queryMargin;

				if(preQueryPos==-1)
					preQueryPos = startQPos;

				if(nullFlag==YES)
				{
					subLen = endQPos - preQueryPos + 1;
					if(subLen>0)
					{
						strncpy(newSeq+newLen, querySeq+preQueryPos-1, subLen);
						newSeq[newLen+subLen] = '\0';
						newLen += subLen;
					}else if(subLen<0)
					{
						printf("line=%d, In %s(), invalid query indel, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
			}else
			{
				queryIndel = misInfo->queryIndel;
				leftSegRow = queryIndel->leftSegRow;
				rightSegRow = queryIndel->rightSegRow;

				if(leftSegRow==-1)
				{
					subLen = globalSegArray[rightSegRow].endQueryPos - preQueryPos + 1;
					if(subLen>0)
					{
						strncpy(newSeq+newLen, querySeq+preQueryPos-1, subLen);
						newSeq[newLen+subLen] = '\0';
						newLen += subLen;
					}else if(subLen<0)
					{
						printf("line=%d, In %s(), invalid query indel, error!\n", __LINE__, __func__);
						return FAILED;
					}
					//endQueryPosNewSeq = globalSegArray[rightSegRow].endQueryPos;
				}else if(rightSegRow==-1)
				{ // do nothing
					if(preQueryPos==-1)
						preQueryPos = startQPos;

					if(queryIndel->misassFlag!=TRUE_MISASS)
					{
						subLen = queryLen - preQueryPos + 1;
						if(subLen>0)
						{
							strncpy(newSeq+newLen, querySeq+preQueryPos-1, subLen);
							newSeq[newLen+subLen] = '\0';
							newLen += subLen;
						}else if(subLen<0)
						{
							printf("line=%d, In %s(), invalid query indel, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
				}else
				{
					subLen = globalSegArray[globalSegNum-1].endQueryPos - preQueryPos + 1;
					if(subLen>0)
					{
						strncpy(newSeq+newLen, querySeq+preQueryPos-1, subLen);
						newSeq[newLen+subLen] = '\0';
						newLen += subLen;
					}else if(subLen<0)
					{
						printf("line=%d, In %s(), invalid query indel, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
			}
		}

		if(misInfo==endMisInfo)
			break;
		misInfo = misInfo->next;
	}
*/

	if(newLen!=newSeqLen)
	{
		printf("line=%d, In %s(), newLen=%d, newSeqLen=%d, error!\n", __LINE__, __func__, newLen, newSeqLen);
		return FAILED;
	}

	queryItem->tailNewQuerySeqList->queryLen = newLen;
	if(queryItem->tailNewQuerySeqList->startQueryPos!=startQPos)
	{
		printf("line=%d, In %s(), startQueryPos=%d, startQPos=%d, error!\n", __LINE__, __func__, queryItem->tailNewQuerySeqList->startQueryPos, startQPos);
		return FAILED;
	}
	if(queryItem->tailNewQuerySeqList->endQueryPos!=endQPos)
	{
		printf("line=%d, In %s(), endQueryPos=%d, endQPos=%d, error!\n", __LINE__, __func__, queryItem->tailNewQuerySeqList->endQueryPos, endQPos);
		return FAILED;
	}

	//printf("newLen=%d, startQueryPos=%d, endQueryPos=%d\n", newLen, startQPos, endQPos);

	return SUCCESSFUL;
}

/**
 * Output new sequence node information.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short outputNewSeqInfo(query_t *queryItem)
{
	querySeq_t *querySeqNode;
	int32_t i;

	if(queryItem->misassFlag==TRUE_MISASS)
	{
		i = 0;
		querySeqNode = queryItem->newQuerySeqList;
		while(querySeqNode)
		{
			printf("newSeq[%d]: startQPos=%d, endQPos=%d, len=%d\n", i, querySeqNode->startQueryPos, querySeqNode->endQueryPos, querySeqNode->queryLen);
			querySeqNode = querySeqNode->next;
			i ++;
		}
	}

	return SUCCESSFUL;
}
