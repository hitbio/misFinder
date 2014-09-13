/*
 * classfy.c
 *
 *  Created on: May 19, 2012
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Classify the queries into four match kinds: (1) perfect, (2) matched, (3) disjunct, (4) unmatched.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short classifyQueries(char *perfectQueryFile, char *matchedQueryFile, char *disjunctQueryFile, char *unmatchedQueryFile, char *queryMatchFile, char *queryMatchInfoFileNew, char *queryMatchInfoFile, char *inputQueryFile, char *mergedSegFile)
{
	int32_t subjectNum;
	querySubject_t *pQuerySubject;
	int64_t i, j;

	if(initMemClassification(perfectQueryFile, matchedQueryFile, disjunctQueryFile, unmatchedQueryFile, queryMatchFile, queryMatchInfoFile, inputQueryFile, mergedSegFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize memory for classification of queries, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		// ########################### Debug information ##############################
		//if(queryMatchInfoSet->queryArray[i].queryID==32 || strcmp(queryMatchInfoSet->queryArray[i].queryTitle, "ctg7180000002387")==0)
		//{
		//	printf("queryID=%d, queryTitle=%s, queryLen=%d, subjectNum=%d\n", queryMatchInfoSet->queryArray[i].queryID, queryMatchInfoSet->queryArray[i].queryTitle, queryMatchInfoSet->queryArray[i].queryLen, queryMatchInfoSet->queryArray[i].querySubjectNum);
		//}
		// ########################### Debug information ##############################

		if(queryMatchInfoSet->queryArray[i].queryLen<minQueryLenThres)
			continue;

		subjectNum = queryMatchInfoSet->queryArray[i].querySubjectNum;
		pQuerySubject = queryMatchInfoSet->queryArray[i].querySubArray;
		for(j=0; j<subjectNum; j++)
		{
			// determine the match kind of the querySubject item
			if(determineMatchKind(pQuerySubject+j, queryMatchInfoSet->queryArray[i].queryLen, queryMatchInfoSet->matchItemArray, queryMatchInfoSet->subjectArray)==FAILED)
			{
				printf("line=%d, In %s(), cannot determine the match kind of query %d, error!\n", __LINE__, __func__, queryMatchInfoSet->queryArray[i].queryID);
				return FAILED;
			}

			// save the valid segments
			if(saveValidSegments(pQuerySubject+j, queryMatchInfoSet->matchItemArray, segLinkSet)==FAILED)
			{
				printf("line=%d, In %s(), cannot save the valid segments, query=%s, subject=%s, error!\n", __LINE__, __func__, queryMatchInfoSet->queryArray[i].queryTitle, queryMatchInfoSet->subjectArray[pQuerySubject[j].subjectID-1].subjectTitle);
				return FAILED;
			}
		}

		// determine the best matched querySubject item
		if(determineBestQuerySubjectItem(queryMatchInfoSet->queryArray+i)==FAILED)
		{
			printf("line=%d, In %s(), cannot determine the best matched querySubject item for query %s, error!\n", __LINE__, __func__, queryMatchInfoSet->queryArray[i].queryTitle);
			return FAILED;
		}

		// remove redundant align segments
		if(removeRedundantAlignSegsSingleQuery(queryMatchInfoSet->queryArray+i)==FAILED)
		{
			printf("line=%d, In %s(), cannot remove redundant align segments, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// adjust incorrect order of the segments
		if(adjustValidSegOrderSingleQuery(queryMatchInfoSet->queryArray+i)==FAILED)
		{
			printf("line=%d, In %s(), cannot adjust the order of align segments, error!\n", __LINE__, __func__);
			return FAILED;
		}

	} // end for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)

	// ############################# Debug information ###########################
	//if(outputQueryMatchInfoText("../tmp", queryMatchInfoSet->queryArray, queryMatchInfoSet->matchItemArray, queryMatchInfoSet->itemNumQueryArray, queryMatchInfoSet->subjectArray)==FAILED)
	//{
	//	printf("line=%d, In %s(), cannot output the query match information to file, query=%s, subject=%s, error!\n", __LINE__, __func__, queryMatchInfoSet->queryArray[i].queryTitle, queryMatchInfoSet->subjectArray[pQuerySubject[j].subjectID-1].subjectTitle);
	//	return FAILED;
	//}
	// ############################# Debug information ###########################

	// fill global align segments
	if(fillGlobalAlignSeg(queryMatchInfoSet, segLinkSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot get global align segments, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// output the best match to file
	if(outputGlobalMatchResultToFile(queryMatchInfoSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot output the matched result of query %d to file, error!\n", __LINE__, __func__, queryMatchInfoSet->queryArray[i].queryID);
		return FAILED;
	}

	if(saveQueryMatchInfoToFile(queryMatchInfoFileNew, queryMatchInfoSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot save the valid segments, error!\n", __LINE__, __func__);
		return FAILED;
	}

	freeMemClassification();

	return SUCCESSFUL;
}

/**
 * Initialize the memory for classification of queries.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemClassification(char *perfectQueryFile, char *matchedQueryFile, char *disjunctQueryFile, char *unmatchedQueryFile, char *queryMatchFile, char *queryMatchInfoFile, char *inputQueryFile, char *mergedSegFile)
{
	// load the query match information from the binary file
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

	segLinkSet = (segLinkSet_t*) calloc (1, sizeof(segLinkSet_t));
	if(segLinkSet==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	segLinkSet->maxArraySize = MAX_ARR_SIZE_SEG_LINK;
	segLinkSet->linkArray = (segmentLink_t*) malloc (segLinkSet->maxArraySize * sizeof(segmentLink_t));
	if(segLinkSet->linkArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	fpPerfectQuery = fopen(perfectQueryFile, "w");
	if(fpPerfectQuery==NULL)
	{
		printf("In %s(), cannot open file [ %s ], error!\n", __func__, perfectQueryFile);
		return FAILED;
	}

	fpMatchedQuery = fopen(matchedQueryFile, "w");
	if(fpMatchedQuery==NULL)
	{
		printf("In %s(), cannot open file [ %s ], error!\n", __func__, matchedQueryFile);
		return FAILED;
	}

	fpDisjunctQuery = fopen(disjunctQueryFile, "w");
	if(fpDisjunctQuery==NULL)
	{
		printf("In %s(), cannot open file [ %s ], error!\n", __func__, disjunctQueryFile);
		return FAILED;
	}

	fpUnmatchedQuery = fopen(unmatchedQueryFile, "w");
	if(fpUnmatchedQuery==NULL)
	{
		printf("In %s(), cannot open file [ %s ], error!\n", __func__, unmatchedQueryFile);
		return FAILED;
	}

	fpQueryMatch = fopen(queryMatchFile, "w");
	if(fpQueryMatch==NULL)
	{
		printf("In %s(), cannot open file [ %s ], error!\n", __func__, queryMatchFile);
		return FAILED;
	}

	return SUCCESSFUL;
}

void freeMemClassification()
{
	releaseQueryMatchInfo(&queryMatchInfoSet);

	free(segLinkSet->linkArray);
	segLinkSet->linkArray = NULL;
	free(segLinkSet);
	segLinkSet = NULL;

	fclose(fpPerfectQuery);
	fpPerfectQuery = NULL;
	fclose(fpMatchedQuery);
	fpMatchedQuery = NULL;
	fclose(fpDisjunctQuery);
	fpDisjunctQuery = NULL;
	fclose(fpUnmatchedQuery);
	fpUnmatchedQuery = NULL;
	fclose(fpQueryMatch);
	fpQueryMatch = NULL;
}

/**
 * Determine the match kind of a querySubject item.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short determineMatchKind(querySubject_t *pQuerySubject, int64_t queryLen, matchItem_t *matchItemArray, subject_t *subjectArray)
{
	int64_t i, matchItemNum, headItemRow, tailItemRow;
	int64_t validTotalMatchLen, subjectLen;
	char usedArray[pQuerySubject->matchItemNum];
	int64_t selectionRound, adjacentRowID, startRowID, addOrder, totalAlignedSegLen;
	int32_t perfectFlag, totalLenEqualToQueryFlag, circularFlag, subjectCircularFlag, disjunctFlag;
	matchItem_t *pMatchItemArray;
	float autoMatchPercentThres;

	for(i=0; i<pQuerySubject->matchItemNum; i++) usedArray[i] = NO;

	pMatchItemArray = matchItemArray + pQuerySubject->firstRow;
	matchItemNum = pQuerySubject->matchItemNum;
	subjectCircularFlag = subjectArray[pQuerySubject->subjectID-1].circularFlag;
	subjectLen = subjectArray[pQuerySubject->subjectID-1].subjectLen;

	if(matchItemNum==0)
	{
		pQuerySubject->circularFlag = NO;
		pQuerySubject->matchKind = UNMATCHED_KIND;
		return SUCCESSFUL;
	}

	// get the valid total match length
	validTotalMatchLen = 0;
	if(pMatchItemArray[0].totalMatchLen>=queryLen)
	{
		validTotalMatchLen = pMatchItemArray[0].totalMatchLen;
	}else
	{
		validTotalMatchLen = queryLen;
	}

	// auto adjust the matchPercentThres
	if(queryLen>shortQueryLenThres)
		autoMatchPercentThres = matchPercentThres;
	else
		autoMatchPercentThres = matchPercentThres * matchPercentFactor;


	segLinkSet->itemNum = 0;
	segLinkSet->headRow = segLinkSet->tailRow = -1;

	// determine the match kind
	if(pMatchItemArray[0].matchLen==validTotalMatchLen)
	{ // perfect match

		// get more items ===============================
		for(i=0; i<matchItemNum; i++)
		{
			if(pMatchItemArray[i].matchLen>=validTotalMatchLen*autoMatchPercentThres)
			{
				if(addNewItemSegLinkArray(segLinkSet, i, 0, 1)==FAILED)
				{
					printf("line=%d, In %s(), cannot add item to link array, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}

		pQuerySubject->circularFlag = NO;
		pQuerySubject->matchKind = PERFECT_MATCH_KIND;

	//}else if(pMatchItemArray[0].matchLen>=validTotalMatchLen*autoMatchPercentThres && validTotalMatchLen-pMatchItemArray[0].totalMatchLen<END_IGNORE_LEN)
	}else if(pMatchItemArray[0].matchLen>=validTotalMatchLen*autoMatchPercentThres && (pMatchItemArray[0].endQueryPos>queryLen-END_IGNORE_LEN && pMatchItemArray[0].startQueryPos<END_IGNORE_LEN))
	{ // match kind with the end length can be ignored

		// get more items =================================
		for(i=0; i<matchItemNum; i++)
		{
			if(pMatchItemArray[i].matchLen>=validTotalMatchLen*autoMatchPercentThres)
			{
				if(addNewItemSegLinkArray(segLinkSet, i, 0, 1)==FAILED)
				{
					printf("line=%d, In %s(), cannot add item to link array, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}

		pQuerySubject->circularFlag = NO;
		pQuerySubject->matchKind = MATCHED_KIND;

	}else
	{ // further to determine the match kind
		// get the first segment
		if(pMatchItemArray[0].matchLen<pMatchItemArray[0].totalMatchLen*autoMatchPercentThres)
		{ // unmatched segment
			// get all the unmatched items =================================
			for(i=0; i<matchItemNum; i++)
			{
				if(pMatchItemArray[i].totalMatchLen>=minAlignedSegLenThres)
				{
					if(addNewItemSegLinkArray(segLinkSet, i, 0, 1)==FAILED)
					{
						printf("line=%d, In %s(), cannot add item to link array, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
			}

			pQuerySubject->circularFlag = NO;
			pQuerySubject->matchKind = UNMATCHED_KIND;

		}else
		{
			// set first segment and the start selection round
			//if(pMatchItemArray[0].endQueryPos<=queryLen-END_IGNORE_LEN)
			//{
				selectionRound = 1;
			//}else// if(pMatchItemArray[0].startQueryPos>=END_IGNORE_LEN)
			//{
			//	selectionRound = 2;
			//}


			if(addNewItemSegLinkArray(segLinkSet, 0, 0, 1)==FAILED)
			{
				printf("line=%d, In %s(), cannot add item to link array, error!\n", __LINE__, __func__);
				return FAILED;
			}

			startRowID = 0;
			usedArray[startRowID] = YES;

			// get the other segments
			while(selectionRound<=2)
			{
				// get the adjacent segments
				if(getAdjacentSegment(&adjacentRowID, startRowID, pMatchItemArray, matchItemNum, usedArray, queryLen, subjectLen, subjectCircularFlag, selectionRound)==FAILED)
				{
					printf("line=%d, In %s(), cannot get adjacent segment, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(selectionRound==1)
				{ // the first selection round
					addOrder = segLinkSet->linkArray[segLinkSet->tailRow].addedOrder;
					if(adjacentRowID>=0)
					{ // valid adjacent segment
						if(addNewItemSegLinkArray(segLinkSet, adjacentRowID, addOrder, 1)==FAILED)
						{
							printf("line=%d, In %s(), cannot add item to link array, error!\n", __LINE__, __func__);
							return FAILED;
						}

						usedArray[adjacentRowID] = YES;
						startRowID = adjacentRowID;

					}else
					{
						startRowID = segLinkSet->linkArray[segLinkSet->headRow].arrRow;
						selectionRound ++;
					}
				}else
				{ // the second selection round
					addOrder =  segLinkSet->linkArray[segLinkSet->headRow].addedOrder;
					if(adjacentRowID>=0)
					{
						if(addNewItemSegLinkArray(segLinkSet, adjacentRowID, addOrder, 2)==FAILED)
						{
							printf("line=%d, In %s(), cannot add item to link array, error!\n", __LINE__, __func__);
							return FAILED;
						}

						usedArray[adjacentRowID] = YES;
						startRowID = adjacentRowID;
					}else
					{
						selectionRound ++;
					}
				}
			}

			// check the head and tail of the linked segments
			headItemRow = segLinkSet->linkArray[segLinkSet->headRow].arrRow;
			tailItemRow = segLinkSet->linkArray[segLinkSet->tailRow].arrRow;
			if(headItemRow==tailItemRow)
			{ // only one segment
				if(pMatchItemArray[headItemRow].matchLen >= validTotalMatchLen*autoMatchPercentThres)  //==============================================
				{
					pQuerySubject->circularFlag = NO;
					pQuerySubject->matchKind = MATCHED_KIND;
				}else                              //==============================================
				{
					pQuerySubject->circularFlag = NO;
					pQuerySubject->matchKind = UNMATCHED_KIND;
				}
			}else
			{ // more than one segment

				// remove redundant segments
				if(removeRedundantSegments(segLinkSet, pMatchItemArray, queryLen, subjectLen)==FAILED)
				{
					printf("line=%d, In %s(), cannot remove the perfect redundant segments, error!\n", __LINE__, __func__);
					return FAILED;
				}

				//if(pQuerySubject->matchKind!=DISJUNCT_MATCH_KIND)
				{
//					if(determineSameOrder(&orderSameFlag, segmentLinkArr, headRowSegmentLinkArr)==FAILED)
//					{
//						printf("line=%d, In %s(), cannot determine the same order, error!\n", __LINE__, __func__);
//						return FAILED;
//					}

					if(determinePerfectSegmentFlag(&perfectFlag, segLinkSet, pMatchItemArray)==FAILED)
					{
						printf("line=%d, In %s(), cannot determine the perfect flag of all segments, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(determineTotalLenEqualToQueryFlag(&totalLenEqualToQueryFlag, segLinkSet, queryLen, pMatchItemArray)==FAILED)
					{
						printf("line=%d, In %s(), cannot determine the perfect flag of all segments, error!\n", __LINE__, __func__);
						return FAILED;
					}

					circularFlag = NO;
					if(subjectArray[pQuerySubject->subjectID-1].circularFlag==YES)
					{
						if(determineCircularQueryFlag(&circularFlag, segLinkSet, pMatchItemArray, subjectLen)==FAILED)
						{
							printf("line=%d, In %s(), cannot determine the perfect flag of all segments, error!\n", __LINE__, __func__);
							return FAILED;
						}
						pQuerySubject->circularFlag = circularFlag;
					}else
					{
						pQuerySubject->circularFlag = circularFlag = NO;
					}

					// determine the match kind
					if(perfectFlag==YES && totalLenEqualToQueryFlag==YES && circularFlag==YES)
					{
						pQuerySubject->matchKind = PERFECT_MATCH_KIND;
					}else
					{
						// compute the aligned segment length
						if(computeAlignedSegLen(&totalAlignedSegLen, segLinkSet, pMatchItemArray)==FAILED)
						{
							printf("line=%d, In %s(), cannot compute the aligned segments length, error!\n", __LINE__, __func__);
							return FAILED;
						}

						// get the disjunct flag
						if(checkDisjunctFlag(&disjunctFlag, circularFlag, subjectLen, segLinkSet, pMatchItemArray)==FAILED)
						{
							printf("line=%d, In %s(), cannot compute the aligned segments length, error!\n", __LINE__, __func__);
							return FAILED;
						}

						// determine the match kind
						if(totalAlignedSegLen<validTotalMatchLen*autoMatchPercentThres)
						{
							pQuerySubject->matchKind = UNMATCHED_KIND;
						}else if(disjunctFlag==YES)
						{
							pQuerySubject->matchKind = DISJUNCT_MATCH_KIND;
						}else
						{
							pQuerySubject->matchKind = MATCHED_KIND;
						}
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the adjacent segment given the start match item.
 * 	If the segment is valid, get the valid segment by adjacentRowID; otherwise, set the adjacentRowID to -1.
 *
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short getAdjacentSegment(int64_t *adjacentRowID, int64_t startRowID, matchItem_t *matchItemArray, int64_t matchItemNum, char *pUsedArray, int64_t queryLen, int64_t subjectLen, int subjectCircularFlag, int selectionRound)
{
	int32_t i;
	int64_t distance, distanceQuery, distanceSubject, minDistQuery, minDistSubject;

	if((selectionRound==1 && matchItemArray[startRowID].endQueryPos>=queryLen-50) || (selectionRound==2 && matchItemArray[startRowID].startQueryPos<50))
	{
		*adjacentRowID = -1;
		return SUCCESSFUL;
	}
	//if((selectionRound==1 && matchItemArray[startRowID].endQueryPos>queryLen-END_IGNORE_LEN) || (selectionRound==2 && matchItemArray[startRowID].startQueryPos<END_IGNORE_LEN))
	//else if(queryLen>2*END_IGNORE_LEN && ((selectionRound==1 && matchItemArray[startRowID].endQueryPos>queryLen-END_IGNORE_LEN) || (selectionRound==2 && matchItemArray[startRowID].startQueryPos<END_IGNORE_LEN)))
	//{
	//	*adjacentRowID = -1;
		//return SUCCESSFUL;
	//}

	// get the adjacent aligned segment according to subjectPos of same strand, and queryPos
	*adjacentRowID = -1;
	minDistQuery = minDistSubject = INT_MAX;
	for(i=0; i<matchItemNum; i++)
	{
		//if(pUsedArray[i]==NO && i!=startRowID && matchItemArray[i].totalMatchLen>minAlignedSegLenThres)
		//if(pUsedArray[i]==NO && i!=startRowID && ((queryLen>END_IGNORE_LEN && matchItemArray[i].totalMatchLen>minAlignedSegLenThres) || (queryLen<=END_IGNORE_LEN && matchItemArray[i].totalMatchLen>0.5*minAlignedSegLenThres)))
		if(pUsedArray[i]==NO && i!=startRowID && ((queryLen>END_IGNORE_LEN && matchItemArray[i].totalMatchLen>minAlignedSegLenThres) || (queryLen*(1-matchPercentThres)<matchItemArray[i].totalMatchLen && matchItemArray[i].totalMatchLen>0.5*minAlignedSegLenThres)))
		{
			if(selectionRound==1)
			{ // the first selection round
				distanceQuery = matchItemArray[i].startQueryPos - matchItemArray[startRowID].endQueryPos;
				if(distanceQuery<0)
					distanceQuery = -distanceQuery;
				if(distanceQuery<=varyEndLenThres)
				{ // distanceQuery <= varyEndLenThres
					if(matchItemArray[i].strand==matchItemArray[startRowID].strand)
					{ // same strand
						if(matchItemArray[i].strand==PLUS_STRAND)
						{ // plus strand
							if(subjectCircularFlag==YES)
							{ // circular subject
								if(matchItemArray[startRowID].endSubPos>=subjectLen-varyEndLenThres && matchItemArray[startRowID].endSubPos<=subjectLen)
								{ // subjectLen-varyEndLenThres <= endSubPos <= subjectLen
									if(matchItemArray[i].startSubPos>=1 && matchItemArray[i].startSubPos<=varyEndLenThres)
									{ // 1 <= [i].startSubPos <= varyEndLenThres
										distanceSubject = matchItemArray[i].startSubPos + subjectLen - matchItemArray[startRowID].endSubPos;
									}else // if(matchItemArray[i].startSubPos>=subjectLen-varyEndLenThres && matchItemArray[i].startSubPos<=subjectLen)
										distanceSubject = matchItemArray[i].startSubPos - matchItemArray[startRowID].endSubPos;

									if(distanceSubject<0)
										distanceSubject = -distanceSubject;
									if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
									{
										*adjacentRowID = i;
										minDistQuery = distanceQuery;
										minDistSubject = distanceSubject;
									}
								}else
								{
									distanceSubject = matchItemArray[i].startSubPos - matchItemArray[startRowID].endSubPos;
									if(distanceSubject<0)
										distanceSubject = -distanceSubject;
									if(distanceSubject<=varyEndLenThres)
									{ // distanceSubject <= varyEndLenThres
										if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
										{
											*adjacentRowID = i;
											minDistQuery = distanceQuery;
											minDistSubject = distanceSubject;
										}
									}
								}
							}else
							{ // linear subject
								distanceSubject = matchItemArray[i].startSubPos - matchItemArray[startRowID].endSubPos;
								if(distanceSubject<0)
									distanceSubject = -distanceSubject;
								if(distanceSubject<=varyEndLenThres)
								{ // distanceSubject <= varyEndLenThres
									if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
									{
										*adjacentRowID = i;
										minDistQuery = distanceQuery;
										minDistSubject = distanceSubject;
									}
								}
							}
						}else
						{ // minus strand
							if(subjectCircularFlag==YES)
							{ // circular subject
								if(matchItemArray[startRowID].startSubPos>=1 && matchItemArray[startRowID].startSubPos<=varyEndLenThres)
								{ // 1 <= [startRow].startSubPos <= varyEndLenThres
									if(matchItemArray[i].startSubPos>=subjectLen-varyEndLenThres && matchItemArray[i].startSubPos<=subjectLen)
									{ // subjectLen-varyEndLenThres <= [i].startSubPos <= subjectLen
										distanceSubject = matchItemArray[startRowID].endSubPos + subjectLen - matchItemArray[i].startSubPos;
									}else// if(matchItemArray[i].startSubPos>=1 && matchItemArray[i].startSubPos<=varyEndLenThres)
										distanceSubject = matchItemArray[startRowID].endSubPos - matchItemArray[i].startSubPos;

									if(distanceSubject<0)
										distanceSubject = -distanceSubject;
									if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
									{
										*adjacentRowID = i;
										minDistQuery = distanceQuery;
										minDistSubject = distanceSubject;
									}
								}else
								{
									distanceSubject = matchItemArray[startRowID].endSubPos - matchItemArray[i].startSubPos;
									if(distanceSubject<0)
										distanceSubject = -distanceSubject;
									if(distanceSubject<=varyEndLenThres)
									{ // distanceSubject <= varyEndLenThres
										if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
										{
											*adjacentRowID = i;
											minDistQuery = distanceQuery;
											minDistSubject = distanceSubject;
										}
									}
								}
							}else
							{ // linear subject
								distanceSubject = matchItemArray[startRowID].endSubPos - matchItemArray[i].startSubPos;
								if(distanceSubject<0)
									distanceSubject = -distanceSubject;
								if(distanceSubject<=varyEndLenThres)
								{ // distanceSubject <= varyEndLenThres
									if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
									{
										*adjacentRowID = i;
										minDistQuery = distanceQuery;
										minDistSubject = distanceSubject;
									}
								}
							}
						}
					}
				}
			}else
			{ // the second selection round
				distanceQuery = matchItemArray[i].endQueryPos - matchItemArray[startRowID].startQueryPos;
				if(distanceQuery<0)
					distanceQuery = -distanceQuery;
				if(distanceQuery<=varyEndLenThres)
				{
					if(matchItemArray[i].strand==matchItemArray[startRowID].strand)
					{ // same strand
						if(matchItemArray[i].strand==PLUS_STRAND)
						{ // plus strand
							if(subjectCircularFlag==YES)
							{ // circular subject
								if(matchItemArray[startRowID].startSubPos>=1 && matchItemArray[startRowID].startSubPos<=varyEndLenThres)
								{ // 1 <= [startRow].startSubPos <= varyEndLenThres
									if(matchItemArray[i].endSubPos>=subjectLen-varyEndLenThres && matchItemArray[i].endSubPos<=subjectLen)
									{ // subjectLen-varyEndLenThres <= [i].endSubPos <= subjectLen
										distanceSubject = matchItemArray[startRowID].startSubPos + subjectLen - matchItemArray[i].endSubPos;
									}else// if(matchItemArray[i].endSubPos>=1 && matchItemArray[i].endSubPos<=varyEndLenThres)
										distanceSubject = matchItemArray[startRowID].startSubPos - matchItemArray[i].endSubPos;

									if(distanceSubject<0)
										distanceSubject = -distanceSubject;
									if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
									{
										*adjacentRowID = i;
										minDistQuery = distanceQuery;
										minDistSubject = distanceSubject;
									}
								}else
								{
									distanceSubject = matchItemArray[startRowID].startSubPos - matchItemArray[i].endSubPos;
									if(distanceSubject<0)
										distanceSubject = -distanceSubject;
									if(distanceSubject<=varyEndLenThres)
									{ // distanceSubject <= varyEndLenThres
										if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
										{
											*adjacentRowID = i;
											minDistQuery = distanceQuery;
											minDistSubject = distanceSubject;
										}
									}
								}
							}else
							{ // linear subject
								distanceSubject = matchItemArray[startRowID].startSubPos - matchItemArray[i].endSubPos;
								if(distanceSubject<0)
									distanceSubject = -distanceSubject;
								if(distanceSubject<=+varyEndLenThres)
								{ // distanceSubject <= varyEndLenThres
									if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
									{
										*adjacentRowID = i;
										minDistQuery = distanceQuery;
										minDistSubject = distanceSubject;
									}
								}
							}
						}else
						{ // minus strand
							if(subjectCircularFlag==YES)
							{ // circular subject
								if(matchItemArray[startRowID].startSubPos>=subjectLen-varyEndLenThres && matchItemArray[startRowID].startSubPos<=subjectLen)
								{ // subjectLen-varyEndLenThres <= [startRow].startSubPos <= subjectLen
									if(matchItemArray[i].endSubPos>=1 && matchItemArray[i].endSubPos<=varyEndLenThres)
									{ // 1 <= [i].endSubPos <= varyEndLenThres
										distanceSubject = matchItemArray[i].endSubPos + subjectLen - matchItemArray[startRowID].startSubPos;
									}else// if(matchItemArray[i].endSubPos>=subjectLen-varyEndLenThres && matchItemArray[i].endSubPos<=subjectLen)
										distanceSubject = matchItemArray[i].endSubPos - matchItemArray[startRowID].startSubPos;

									if(distanceSubject<0)
										distanceSubject = -distanceSubject;
									if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
									{
										*adjacentRowID = i;
										minDistQuery = distanceQuery;
										minDistSubject = distanceSubject;
									}
								}else
								{
									distanceSubject = matchItemArray[i].endSubPos - matchItemArray[startRowID].startSubPos;
									if(distanceSubject<0)
										distanceSubject = -distanceSubject;
									if(distanceSubject<=varyEndLenThres)
									{ // distanceSubject <= varyEndLenThres
										if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
										{
											*adjacentRowID = i;
											minDistQuery = distanceQuery;
											minDistSubject = distanceSubject;
										}
									}
								}
							}else
							{ // linear subject
								distanceSubject = matchItemArray[i].endSubPos - matchItemArray[startRowID].startSubPos;
								if(distanceSubject<0)
									distanceSubject = -distanceSubject;
								if(distanceSubject<=varyEndLenThres)
								{ // distanceSubject <= varyEndLenThres
									if(distanceQuery<minDistQuery && distanceSubject<minDistSubject)
									{
										*adjacentRowID = i;
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

	// get the adjacent aligned segment according to same strand, queryPos
	if((*adjacentRowID)==-1)
	{
		for(i=0; i<matchItemNum; i++)
		{
			if(pUsedArray[i]==NO && i!=startRowID && matchItemArray[i].totalMatchLen>minAlignedSegLenThres)
			{
				if(selectionRound==1)
				{ // the first selection round
					if(matchItemArray[i].startQueryPos >= matchItemArray[startRowID].endQueryPos-varyEndLenThres && matchItemArray[i].startQueryPos <= matchItemArray[startRowID].endQueryPos+varyEndLenThres)
					{ // [startRow].endQueryPos-varyEndLenThres <= [i].startQueryPos <= [startRow].endQueryPos+varyEndLenThres
						if(matchItemArray[i].strand==matchItemArray[startRowID].strand)
						{ // same strand
							*adjacentRowID = i;
							break;
						}
					}
				}else
				{ // the second selection round
					if(matchItemArray[i].endQueryPos >= matchItemArray[startRowID].startQueryPos-varyEndLenThres && matchItemArray[i].endQueryPos <= matchItemArray[startRowID].startQueryPos+varyEndLenThres)
					{
						if(matchItemArray[i].strand==matchItemArray[startRowID].strand)
						{ // same strand
							*adjacentRowID = i;
							break;
						}
					}
				}
			}
		}
	}

	// get the adjacent aligned segment according to queryPos
	if((*adjacentRowID)==-1)
	{
		for(i=0; i<matchItemNum; i++)
		{
			if(pUsedArray[i]==NO && i!=startRowID && matchItemArray[i].totalMatchLen>minAlignedSegLenThres)
			{
				if(selectionRound==1)
				{ // the first selection round
					if(matchItemArray[i].startQueryPos >= matchItemArray[startRowID].endQueryPos-varyEndLenThres && matchItemArray[i].startQueryPos <= matchItemArray[startRowID].endQueryPos+varyEndLenThres)
					{ // [startRow].endQueryPos-varyEndLenThres <= [i].startQueryPos <= [startRow].endQueryPos+varyEndLenThres
						*adjacentRowID = i;
						break;
					}
				}else
				{ // the second selection round
					if(matchItemArray[i].endQueryPos >= matchItemArray[startRowID].startQueryPos-varyEndLenThres && matchItemArray[i].endQueryPos <= matchItemArray[startRowID].startQueryPos+varyEndLenThres)
					{
						*adjacentRowID = i;
						break;
					}
				}
			}
		}
	}

	// get the most adjacent aligned segment according to queryPos
	if((*adjacentRowID)==-1)
	{
		minDistQuery = INT64_MAX;
		for(i=0; i<matchItemNum; i++)
		{
			if(pUsedArray[i]==NO && i!=startRowID && matchItemArray[i].totalMatchLen>minAlignedSegLenThres)
			{
				if(selectionRound==1)
				{ // the first selection round
					if(matchItemArray[i].endQueryPos>matchItemArray[startRowID].endQueryPos)
					{
						distance = matchItemArray[i].startQueryPos - matchItemArray[startRowID].endQueryPos;
						if(distance<0)
							distance = -distance;

						if(distance<minDistQuery)
						{
							*adjacentRowID = i;
							minDistQuery = distance;
						}
					}
				}else
				{ // the second selection round
					if(matchItemArray[i].startQueryPos<matchItemArray[startRowID].startQueryPos)
					{
						distance = matchItemArray[startRowID].startQueryPos - matchItemArray[i].endQueryPos;
						if(distance<0)
							distance = -distance;

						if(distance<minDistQuery)
						{
							*adjacentRowID = i;
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
 * Add item to link array.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewItemSegLinkArray(segLinkSet_t *segLinkSet, int32_t dataRow, int32_t addOrder, int32_t headTailFlag)
{
	if(headTailFlag==1)
	{ // tail
		segLinkSet->linkArray[segLinkSet->itemNum].arrRow = dataRow;
		segLinkSet->linkArray[segLinkSet->itemNum].addedOrder = addOrder + 1;
		segLinkSet->linkArray[segLinkSet->itemNum].validFlag = YES;
		segLinkSet->linkArray[segLinkSet->itemNum].previous = segLinkSet->tailRow;
		segLinkSet->linkArray[segLinkSet->itemNum].next = -1;

		if(segLinkSet->itemNum==0)
			segLinkSet->headRow = segLinkSet->tailRow = segLinkSet->itemNum;
		else
		{
			segLinkSet->linkArray[segLinkSet->tailRow].next = segLinkSet->itemNum;
			segLinkSet->tailRow = segLinkSet->itemNum;
		}
		segLinkSet->itemNum ++;
	}else if(headTailFlag==2)
	{ // head
		segLinkSet->linkArray[segLinkSet->itemNum].arrRow = dataRow;
		segLinkSet->linkArray[segLinkSet->itemNum].addedOrder = addOrder - 1;
		segLinkSet->linkArray[segLinkSet->itemNum].validFlag = YES;
		segLinkSet->linkArray[segLinkSet->itemNum].previous = -1;
		segLinkSet->linkArray[segLinkSet->itemNum].next = segLinkSet->headRow;

		if(segLinkSet->itemNum==0)
			segLinkSet->headRow = segLinkSet->tailRow = segLinkSet->itemNum;
		else
		{
			segLinkSet->linkArray[segLinkSet->headRow].previous = segLinkSet->itemNum;
			segLinkSet->headRow = segLinkSet->itemNum;
		}
		segLinkSet->itemNum ++;
	}else
	{
		printf("line=%d, In %s(), cannot add item to linkArray, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(segLinkSet->itemNum>=segLinkSet->maxArraySize)
	{
		segLinkSet->linkArray = (segmentLink_t *) realloc(segLinkSet->linkArray, 2 * segLinkSet->maxArraySize);
		if(segLinkSet->linkArray==NULL)
		{
			printf("line=%d, In %s(), cannot reallocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		segLinkSet->maxArraySize *= 2;
	}

	return SUCCESSFUL;
}

/**
 * Determine the same order of segments.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
//short determineSameOrder(int *orderSameFlag, segmentLink_t *segmentLinkArray, int headRowSegmentLinkArray)
//{
//	int itemRowSegmentLinkArray, nextRowSegmentLinkArray;
//
//	*orderSameFlag = NO;
//	itemRowSegmentLinkArray = headRowSegmentLinkArray;
//	nextRowSegmentLinkArray = segmentLinkArray[itemRowSegmentLinkArray].next;
//	while(nextRowSegmentLinkArray!=-1)
//	{
//		if(segmentLinkArray[itemRowSegmentLinkArray].addedOrder==segmentLinkArray[nextRowSegmentLinkArray].addedOrder)
//		{
//			*orderSameFlag = YES;
//			break;
//		}
//
//		itemRowSegmentLinkArray = nextRowSegmentLinkArray;
//		nextRowSegmentLinkArray = segmentLinkArray[nextRowSegmentLinkArray].next;
//	}
//
//	return SUCCESSFUL;
//}

/**
 * Determine the perfect segment flag.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeRedundantSegments(segLinkSet_t *segLinkSet, matchItem_t *matchItemArray, int64_t queryLen, int64_t subjectLen)
{
	int32_t itemRow, tmpRow, nextRow, itemRowSegmentLinkArray, nextRowSegmentLinkArray;
	int32_t strand, startQueryPos, endQueryPos, startSubjectPos, endSubjectPos;
	int32_t queryPosMatchFlag, subjectPosMatchFlag, removeFlag;

	// reset the used flag
	itemRowSegmentLinkArray = segLinkSet->headRow;
	while(itemRowSegmentLinkArray!=-1)
	{
		segLinkSet->linkArray[itemRowSegmentLinkArray].usedFlag = NO;
		itemRowSegmentLinkArray = segLinkSet->linkArray[itemRowSegmentLinkArray].next;
	}

	removeFlag = NO;

	itemRowSegmentLinkArray = segLinkSet->headRow;
	itemRow = segLinkSet->linkArray[itemRowSegmentLinkArray].arrRow;
	segLinkSet->linkArray[itemRowSegmentLinkArray].usedFlag = YES;
	strand = matchItemArray[itemRow].strand;
	startQueryPos = matchItemArray[itemRow].startQueryPos;
	endQueryPos = matchItemArray[itemRow].endQueryPos;
	startSubjectPos = matchItemArray[itemRow].startSubPos;
	endSubjectPos = matchItemArray[itemRow].endSubPos;

	nextRowSegmentLinkArray = segLinkSet->linkArray[itemRowSegmentLinkArray].next;
	while(nextRowSegmentLinkArray!=-1)
	{
		nextRow = segLinkSet->linkArray[nextRowSegmentLinkArray].arrRow;

		queryPosMatchFlag = NO;
		subjectPosMatchFlag = NO;

		if(matchItemArray[nextRow].strand==strand)
		{
			if(endQueryPos+1==matchItemArray[nextRow].startQueryPos)
			{
				queryPosMatchFlag = YES;
			}

			if((strand==PLUS_STRAND && (endSubjectPos+1)%subjectLen==matchItemArray[nextRow].startSubPos) ||
				(strand==MINUS_STRAND && endSubjectPos==(matchItemArray[nextRow].startSubPos+1)%subjectLen))
			{
				subjectPosMatchFlag = YES;
			}


			if(queryPosMatchFlag==YES && subjectPosMatchFlag==YES)
			{
				startQueryPos = matchItemArray[nextRow].startQueryPos;
				endQueryPos = matchItemArray[nextRow].endQueryPos;
				startSubjectPos = matchItemArray[nextRow].startSubPos;
				endSubjectPos = matchItemArray[nextRow].endSubPos;
				segLinkSet->linkArray[nextRowSegmentLinkArray].usedFlag = YES;
			}
		}

		if(endQueryPos==queryLen)
		{
			removeFlag = YES;
			break;
		}

		nextRowSegmentLinkArray = segLinkSet->linkArray[nextRowSegmentLinkArray].next;
	}

	// remove the segments
	if(removeFlag==YES)
	{
		itemRowSegmentLinkArray = segLinkSet->headRow;
		while(itemRowSegmentLinkArray!=-1)
		{
			if(segLinkSet->linkArray[itemRowSegmentLinkArray].usedFlag==NO)
			{ // remove this segment
				if(segLinkSet->linkArray[itemRowSegmentLinkArray].previous!=-1)
				{
					tmpRow = segLinkSet->linkArray[itemRowSegmentLinkArray].previous;
					segLinkSet->linkArray[tmpRow].next = segLinkSet->linkArray[itemRowSegmentLinkArray].next;
				}

				if(segLinkSet->linkArray[itemRowSegmentLinkArray].next!=-1)
				{
					tmpRow = segLinkSet->linkArray[itemRowSegmentLinkArray].next;
					segLinkSet->linkArray[tmpRow].previous = segLinkSet->linkArray[itemRowSegmentLinkArray].previous;
				}

				segLinkSet->itemNum --;

				// print this node
				//itemRow = segmentLinkArray[itemRowSegmentLinkArray].arrRow;
				//printf("startSubPos=%d, endSubPos=%d, startQueryPos=%d, endQueryPos=%d, orientation=%d\n", matchItemArray[itemRow].startSubPos, matchItemArray[itemRow].endSubPos, matchItemArray[itemRow].startQueryPos, matchItemArray[itemRow].endQueryPos, matchItemArray[itemRow].strand);

			}else
			{
				segLinkSet->tailRow = itemRowSegmentLinkArray;
			}

			itemRowSegmentLinkArray = segLinkSet->linkArray[itemRowSegmentLinkArray].next;
		}
	}

	return SUCCESSFUL;
}

/**
 * Determine the perfect segment flag.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short determinePerfectSegmentFlag(int32_t *perfectFlag, segLinkSet_t *segLinkSet, matchItem_t *matchItemArray)
{
	int32_t itemRow, itemRowsegmentLinkArray;

	*perfectFlag = YES;
	itemRowsegmentLinkArray = segLinkSet->headRow;
	while(itemRowsegmentLinkArray!=-1)
	{
		itemRow = segLinkSet->linkArray[itemRowsegmentLinkArray].arrRow;
		if(matchItemArray[itemRow].totalMatchLen!=matchItemArray[itemRow].matchLen)
		{
			*perfectFlag = NO;
			break;
		}
		itemRowsegmentLinkArray = segLinkSet->linkArray[itemRowsegmentLinkArray].next;
	}

	return SUCCESSFUL;
}

/**
 * Determine the total length Equal flag.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short determineTotalLenEqualToQueryFlag(int32_t *totalLenEqualToQueryFlag, segLinkSet_t *segLinkSet, int32_t queryLen, matchItem_t *matchItemArray)
{
	int32_t itemRow, itemRowSegmentLinkArray, totalLen;

	totalLen = 0;
	itemRowSegmentLinkArray = segLinkSet->headRow;
	while(itemRowSegmentLinkArray!=-1)
	{
		itemRow = segLinkSet->linkArray[itemRowSegmentLinkArray].arrRow;
		totalLen += matchItemArray[itemRow].matchLen;
		itemRowSegmentLinkArray = segLinkSet->linkArray[itemRowSegmentLinkArray].next;
	}

	if(totalLen==queryLen)
	{
		*totalLenEqualToQueryFlag = YES;
	}else
	{
		*totalLenEqualToQueryFlag = NO;
	}

	return SUCCESSFUL;
}

/**
 * Determine the circular flag of a query.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short determineCircularQueryFlag(int32_t *circularFlag, segLinkSet_t *segLinkSet, matchItem_t *matchItemArray, int64_t subjectLen)
{
	int32_t itemRow, nextRow, itemRowSegmentLinkArray, nextRowSegmentLinkArray;

	*circularFlag = NO;
	itemRowSegmentLinkArray = segLinkSet->headRow;
	nextRowSegmentLinkArray = segLinkSet->linkArray[itemRowSegmentLinkArray].next;
	while(nextRowSegmentLinkArray!=-1)
	{
		itemRow = segLinkSet->linkArray[itemRowSegmentLinkArray].arrRow;
		nextRow = segLinkSet->linkArray[nextRowSegmentLinkArray].arrRow;

		if(matchItemArray[itemRow].strand==matchItemArray[nextRow].strand)
		{
			if(matchItemArray[itemRow].strand==PLUS_STRAND)
			{ // plus strand
				if(matchItemArray[itemRow].endSubPos>=subjectLen-varyEndLenThres && matchItemArray[nextRow].startSubPos<=varyEndLenThres)
				{
					*circularFlag = YES;
					break;
				}
			}else
			{ // minus strand
				if(matchItemArray[itemRow].endSubPos<=varyEndLenThres && matchItemArray[nextRow].startSubPos>=subjectLen-varyEndLenThres)
				{
					*circularFlag = YES;
					break;
				}
			}
		}

		itemRowSegmentLinkArray = nextRowSegmentLinkArray;
		nextRowSegmentLinkArray = segLinkSet->linkArray[nextRowSegmentLinkArray].next;
	}

	return SUCCESSFUL;
}

/**
 * Compute the aligned segment length of queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short computeAlignedSegLen(int64_t *totalAlignedSegLen, segLinkSet_t *segLinkSet, matchItem_t *pMatchItemArray)
{
	int32_t itemRow, itemRowSegmentLinkArray;

	*totalAlignedSegLen = 0;
	itemRowSegmentLinkArray = segLinkSet->headRow;
	while(itemRowSegmentLinkArray!=-1)
	{
		itemRow = segLinkSet->linkArray[itemRowSegmentLinkArray].arrRow;
		*totalAlignedSegLen += pMatchItemArray[itemRow].totalMatchLen;
		itemRowSegmentLinkArray = segLinkSet->linkArray[itemRowSegmentLinkArray].next;
	}

	return SUCCESSFUL;
}

/**
 * Check whether the query is a disjunct one.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short checkDisjunctFlag(int32_t *disjunctFlag, int32_t circularFlag, int64_t subjectLen, segLinkSet_t *segLinkSet, matchItem_t *pMatchItemArray)
{
	int32_t itemRow, nextRow, itemRowSegmentLinkArray, nextRowSegmentLinkArray;
	int64_t distance;

	*disjunctFlag = NO;
	itemRowSegmentLinkArray = segLinkSet->headRow;
	nextRowSegmentLinkArray = segLinkSet->linkArray[itemRowSegmentLinkArray].next;
	while(nextRowSegmentLinkArray!=-1)
	{
		itemRow = segLinkSet->linkArray[itemRowSegmentLinkArray].arrRow;
		nextRow = segLinkSet->linkArray[nextRowSegmentLinkArray].arrRow;

		if(pMatchItemArray[nextRow].totalMatchLen>=shortQueryLenThres)
		{
			if(pMatchItemArray[itemRow].strand!=pMatchItemArray[nextRow].strand)
			{
				*disjunctFlag = YES;
			}else
			{
				if(circularFlag==YES)
				{
					if(pMatchItemArray[itemRow].strand==PLUS_STRAND)
					{
						if((pMatchItemArray[itemRow].endSubPos>=subjectLen-varyEndLenThres && pMatchItemArray[itemRow].endSubPos<=subjectLen)
							&& (pMatchItemArray[nextRow].startSubPos>=1 && pMatchItemArray[nextRow].startSubPos<=varyEndLenThres))
						{
							distance = (subjectLen - pMatchItemArray[itemRow].endSubPos) + (pMatchItemArray[nextRow].startSubPos - 1);
							if(distance<0)
								distance = -distance;

							if(distance > minDisjunctDistanceThres)
							{
								*disjunctFlag = YES;
							}
						}else
						{
							distance  = pMatchItemArray[nextRow].startSubPos - pMatchItemArray[itemRow].endSubPos;
							if(distance<0)
								distance = -distance;

							if(distance > minDisjunctDistanceThres)
							{
								*disjunctFlag = YES;
							}
						}
					}else
					{
						if((pMatchItemArray[itemRow].endSubPos>=1 && pMatchItemArray[itemRow].endSubPos<=varyEndLenThres)
							&& (pMatchItemArray[nextRow].startSubPos>=subjectLen-varyEndLenThres && pMatchItemArray[nextRow].startSubPos<=subjectLen))
						{
							distance = (subjectLen - pMatchItemArray[nextRow].startSubPos) + (pMatchItemArray[itemRow].endSubPos - 1);
							if(distance<0)
								distance = -distance;

							if(distance > minDisjunctDistanceThres)
							{
								*disjunctFlag = YES;
							}
						}
					}
				}else
				{
					distance  = pMatchItemArray[nextRow].startSubPos - pMatchItemArray[itemRow].endSubPos;
					if(distance<0)
						distance = -distance;

					if(distance > minDisjunctDistanceThres)
					{
						*disjunctFlag = YES;
					}
				}
			}
		}

		itemRowSegmentLinkArray = nextRowSegmentLinkArray;
		nextRowSegmentLinkArray = segLinkSet->linkArray[nextRowSegmentLinkArray].next;
	}

	return SUCCESSFUL;
}

/**
 * Output the best match information of queries to file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short outputGlobalMatchResultToFile(queryMatchInfo_t *queryMatchInfoSet)
{
	int32_t i, j, globalSegNum, subjectLen;
	query_t *queryItem;
	subject_t *subjectArray;
	globalValidSeg_t *globalSegArray;
	char *subjectTitle;
	FILE *fpCategoryQuery;

	subjectArray = queryMatchInfoSet->subjectArray;

	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		queryItem = queryMatchInfoSet->queryArray + i;
		globalSegArray = queryItem->globalValidSegArray;
		globalSegNum = queryItem->globalValidSegNum;

		switch(queryItem->globalMatchKind)
		{
			case PERFECT_MATCH_KIND: fpCategoryQuery = fpPerfectQuery; break;
			case MATCHED_KIND: fpCategoryQuery = fpMatchedQuery; break;
			case DISJUNCT_MATCH_KIND: fpCategoryQuery = fpDisjunctQuery; break;
			case UNMATCHED_KIND: fpCategoryQuery = fpUnmatchedQuery; break;
			default: printf("line=%d, In %s(), matchKind=%d, error!\n", __LINE__, __func__, queryItem->globalMatchKind); return FAILED;
		}

		if(globalSegNum>0)
		{
			fprintf(fpCategoryQuery, ">%s\t%d\t%d\t%d\t%d\n", queryItem->queryTitle, queryItem->queryLen, globalSegNum, queryItem->globalMatchKind, queryItem->circularFlag);
			fprintf(fpQueryMatch, ">%s\t%d\t%d\t%d\t%d\n", queryItem->queryTitle, queryItem->queryLen, globalSegNum, queryItem->globalMatchKind, queryItem->circularFlag);

			for(j=0; j<globalSegNum; j++)
			{
				subjectTitle = subjectArray[ globalSegArray[j].subjectID-1 ].subjectTitle;
				subjectLen = subjectArray[ globalSegArray[j].subjectID-1 ].subjectLen;

				fprintf(fpCategoryQuery, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%s\t%d\n", globalSegArray[j].startSubPos, globalSegArray[j].endSubPos, globalSegArray[j].startQueryPos, globalSegArray[j].endQueryPos, globalSegArray[j].strand, globalSegArray[j].matchLen, globalSegArray[j].totalMatchLen, globalSegArray[j].gapNum, globalSegArray[j].matchPercent, subjectTitle, subjectLen);
				fprintf(fpQueryMatch, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%s\t%d\n", globalSegArray[j].startSubPos, globalSegArray[j].endSubPos, globalSegArray[j].startQueryPos, globalSegArray[j].endQueryPos, globalSegArray[j].strand, globalSegArray[j].matchLen, globalSegArray[j].totalMatchLen, globalSegArray[j].gapNum, globalSegArray[j].matchPercent, subjectTitle, subjectLen);
			}
		}else
		{
			fprintf(fpCategoryQuery, ">%s\t%d\t%d\t%d\t%d\n", queryItem->queryTitle, queryItem->queryLen, 0, queryItem->globalMatchKind, queryItem->circularFlag);
			fprintf(fpQueryMatch, ">%s\t%d\t%d\t%d\t%d\n", queryItem->queryTitle, queryItem->queryLen, 0, queryItem->globalMatchKind, queryItem->circularFlag);
		}
	}

	return SUCCESSFUL;
}

/**
 * save the valid segments of querySubject item.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short saveValidSegments(querySubject_t *pQuerySubject, matchItem_t *matchItemArray, segLinkSet_t *segLinkSet)
{
	int32_t newMatchItemNum, headItemRow, headRow;
	matchItem_t *pMatchItemArray;

	pMatchItemArray = matchItemArray + pQuerySubject->firstRow;

	// allocate the memory for valid segments
	if(segLinkSet->itemNum>0)
	{
		pQuerySubject->validSegArray = (validSegment_t *) malloc(segLinkSet->itemNum * sizeof(validSegment_t));
		if(pQuerySubject->validSegArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		pQuerySubject->validSegmentNum = segLinkSet->itemNum;

		// fill the linked segments
		newMatchItemNum = 0;
		headRow = segLinkSet->headRow;
		while(headRow!=-1)
		{
			headItemRow = segLinkSet->linkArray[headRow].arrRow;

			if(memcpy(pQuerySubject->validSegArray+newMatchItemNum, pMatchItemArray+headItemRow, sizeof(validSegment_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			newMatchItemNum ++;

			headRow = segLinkSet->linkArray[headRow].next;
		}

		// ############################# Debug information ##########################
		if(newMatchItemNum!=segLinkSet->itemNum)
		{
			printf("line=%d, In %s(), newMatchItemNum=%d != itemNumSegmentLinkArray=%d, error!\n", __LINE__, __func__, newMatchItemNum, segLinkSet->itemNum);
			return FAILED;
		}
		// ############################# Debug information ##########################
	}else
	{
		pQuerySubject->validSegArray = NULL;
		pQuerySubject->validSegmentNum = 0;
	}

	return SUCCESSFUL;
}

/**
 * Determine the best matched querySubject item.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short determineBestQuerySubjectItem(query_t *pQuery)
{
	int32_t i, j;
	int32_t matchKindArray[4] = {PERFECT_MATCH_KIND, MATCHED_KIND, DISJUNCT_MATCH_KIND, UNMATCHED_KIND};
	int32_t subjectNum, validSegNum, queryMatchKind, tmpMatchKind, maxSubjectIndex;
	int64_t total, maxTotalLen;
	querySubject_t *pQuerySubjectArray;
	validSegment_t *validSegmentArray;

	pQuerySubjectArray = pQuery->querySubArray;
	subjectNum = pQuery->querySubjectNum;

	if(subjectNum==0)
	{
		pQuery->bestMatchRow = -1;
		pQuery->circularFlag = NO;

		return SUCCESSFUL;
	}

	// determine the match kind of the query
	queryMatchKind = -1;
	for(i=0; i<4; i++)
	{
		tmpMatchKind = matchKindArray[i];
		for(j=0; j<subjectNum; j++)
		{
			if(pQuerySubjectArray[j].matchKind==tmpMatchKind)
			{
				queryMatchKind = tmpMatchKind;
				break;
			}
		}

		if(queryMatchKind!=-1)
			break;
	}

	// determine the best querySubject item
	if(queryMatchKind==PERFECT_MATCH_KIND)
	{ // choose the first perfect item
		for(i=0; i<subjectNum; i++)
		{
			if(pQuerySubjectArray[i].matchKind==queryMatchKind)
			{
				pQuery->bestMatchRow = i;
				pQuery->circularFlag = pQuerySubjectArray[i].circularFlag;
				break;
			}
		}

		// ############################ Debug information ############################
		if(i==subjectNum)
		{
			printf("line=%d, In %s(), i=%d != subjectNum=%d, error!\n", __LINE__, __func__, i, subjectNum);
			return FAILED;
		}
		// ############################ Debug information ############################
	}else
	{ // choose the segment having the largest summed matchLen;

		maxTotalLen = 0;
		maxSubjectIndex = -1;
		for(i=0; i<subjectNum; i++)
		{
			if(pQuerySubjectArray[i].matchKind==queryMatchKind)
			{
				total = 0;
				validSegNum = pQuerySubjectArray[i].validSegmentNum;
				validSegmentArray = pQuerySubjectArray[i].validSegArray;
				for(j=0; j<validSegNum; j++)
				{
					total += validSegmentArray[j].matchLen;
				}

				if(total>maxTotalLen)
				{
					maxSubjectIndex = i;
					maxTotalLen = total;
				}
			}
		}

		if(maxSubjectIndex>=0)
		{
			pQuery->bestMatchRow = maxSubjectIndex;
			pQuery->circularFlag = pQuerySubjectArray[maxSubjectIndex].circularFlag;
		}else
		{
			pQuery->bestMatchRow = -1;
			pQuery->circularFlag = NO;
		}
	}

	return SUCCESSFUL;
}

/**
 * Remove redundant align segments for single querySubject item.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short removeRedundantAlignSegsSingleQuery(query_t *pQuery)
{
	int32_t i, j, validSegNum, startItemRow, endItemRow, rowValidSeg1, rowValidSeg2, *reduFlagArray, newItemNum;
	validSegment_t *validSegArray;
	int64_t startSegPos, endSegPos;

	if(pQuery->querySubArray==NULL || pQuery->bestMatchRow<0)
	{
		return SUCCESSFUL;
	}

	validSegArray = (pQuery->querySubArray + pQuery->bestMatchRow)->validSegArray;
	validSegNum = (pQuery->querySubArray + pQuery->bestMatchRow)->validSegmentNum;

	if(validSegNum>=2)
	{
		reduFlagArray = (int32_t *) calloc (validSegNum, sizeof(int32_t));
		if(reduFlagArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		for(i=0; i<validSegNum; i++) reduFlagArray[i] = NO;

		// get start and end alignment position, and item numbers
		startSegPos = endSegPos = -1;
		startItemRow = endItemRow = -1;
		for(i=1; i<validSegNum; i++)
		{
			if(startSegPos==-1)
			{
				if(validSegArray[i-1].endQueryPos>validSegArray[i].startQueryPos)
				{
					startSegPos = validSegArray[i].startQueryPos;
					endSegPos = validSegArray[i-1].endQueryPos;
					startItemRow = i - 1;

					if(i==validSegNum-1)
						endItemRow = i;
					else if(validSegArray[i+1].startQueryPos>endSegPos)
						endItemRow = i;
				}
			}else
			{
				if(validSegArray[i].startQueryPos<startSegPos)
					startSegPos = validSegArray[i].startQueryPos;

				if(validSegArray[i].endQueryPos>endSegPos && validSegArray[i].endQueryPos<endSegPos+200)
					endSegPos = validSegArray[i].endQueryPos;

				if(i==validSegNum-1)
					endItemRow = i;
				else if(validSegArray[i+1].startQueryPos>endSegPos)
					endItemRow = i;
			}

			if(startItemRow>=0 && endItemRow>=0)
			{
				// get valid two large segments
				rowValidSeg1 = rowValidSeg2 = -1;
				for(j=startItemRow; j<=endItemRow; j++)
				{
					if(validSegArray[j].startQueryPos<startSegPos)
						rowValidSeg1 = j;
					else if(rowValidSeg1==-1 && validSegArray[j].startQueryPos==startSegPos)
						rowValidSeg1 = j;

					if(validSegArray[j].endQueryPos>endSegPos)
						rowValidSeg2 = j;
					else if(rowValidSeg2==-1 && validSegArray[j].endQueryPos==endSegPos)
						rowValidSeg2 = j;
				}

				if(rowValidSeg1>=0 && rowValidSeg2>=0)
				{
					for(j=startItemRow; j<=endItemRow; j++)
					{
						if(j!=rowValidSeg1 && j!=rowValidSeg2)
							reduFlagArray[j] = YES;
					}
				}

				startSegPos = endSegPos = -1;
				startItemRow = endItemRow = -1;
			}
		}

		// remove redundant segments
		newItemNum = 0;
		for(i=0; i<validSegNum; i++)
		{
			if(reduFlagArray[i]==NO)
			{
				if(i!=newItemNum)
				{
					if(memcpy(validSegArray+newItemNum, validSegArray+i, sizeof(validSegment_t))==NULL)
					{
						printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
				newItemNum ++;
			}
		}
		(pQuery->querySubArray + pQuery->bestMatchRow)->validSegmentNum = newItemNum;

		free(reduFlagArray);
	}

	return SUCCESSFUL;
}

/**
 * Adjust the order of valid align segments for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short adjustValidSegOrderSingleQuery(query_t *queryItem)
{
	int32_t i, validSegNum;
	validSegment_t *validSegArray, tmpSeg;

	if(queryItem->bestMatchRow>=0)
	{
		validSegArray = queryItem->querySubArray[queryItem->bestMatchRow].validSegArray;
		validSegNum = queryItem->querySubArray[queryItem->bestMatchRow].validSegmentNum;

		i = 1;
		while(i<validSegNum)
		{
			if(validSegArray[i-1].startQueryPos>validSegArray[i].startQueryPos && validSegArray[i-1].endQueryPos>validSegArray[i].endQueryPos)
			{
				if(memcpy(&tmpSeg, validSegArray+i-1, sizeof(validSegment_t))==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				if(memcpy(validSegArray+i-1, validSegArray+i, sizeof(validSegment_t))==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				if(memcpy(validSegArray+i, &tmpSeg, sizeof(validSegment_t))==NULL)
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
	}

	return SUCCESSFUL;
}
