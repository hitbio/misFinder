/*
 * statistics.c
 *
 *  Created on: May 19, 2012
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Get the length statistics for queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short queryLenStatistics(metrics_t *queryMetrics, char *sortedQueryFile)
{
	// initialize the memory
	if(initMemLenStatistics()==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory for query statistics, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the data for queryLenStatisticArr
	if(fillQueryDataLenStatistic(&lenStatisticArr, &itemNumLenStatisticArr, queryMatchInfoSet->queryArray, queryMatchInfoSet->itemNumQueryArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the query data for statistics, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// sort the queries according to their lengths
	if(radixSortOfLengths(lenStatisticArr, lenStatisitcBufArr, itemNumLenStatisticArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot sort the queries according to their lenghts, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(saveSortedQueries(sortedQueryFile, lenStatisticArr, itemNumLenStatisticArr, queryMatchInfoSet->queryArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot output the sorted queries to file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the maxSize, N50, mean size, median size
	if(computeLenStatistics(queryMetrics, lenStatisticArr, itemNumLenStatisticArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the statistics of queries, error!\n", __LINE__, __func__);
		return FAILED;
	}

	freeMemLenStatistics();

	return SUCCESSFUL;
}

/**
 * Initialize the query statistics memory.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemLenStatistics()
{
	maxItemNumLenStatisticArr = queryMatchInfoSet->itemNumQueryArray;
	lenStatisticArr = (queryLenStatistic_t *) calloc (maxItemNumLenStatisticArr, sizeof(queryLenStatistic_t));
	if(lenStatisticArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	lenStatisitcBufArr = (queryLenStatistic_t *) calloc (maxItemNumLenStatisticArr, sizeof(queryLenStatistic_t));
	if(lenStatisitcBufArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Free the statistics memory.
 */
void freeMemLenStatistics()
{
	itemNumLenStatisticArr =0;
	maxItemNumLenStatisticArr = 0;

	free(lenStatisticArr);
	lenStatisticArr = NULL;
	free(lenStatisitcBufArr);
	lenStatisitcBufArr = NULL;
}

/**
 * Fill the data for query statistics.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillQueryDataLenStatistic(queryLenStatistic_t **lenStatisticArray, int64_t *itemNumLenStatisticArray, query_t *queryArray, int64_t itemNumQueryArray)
{
	int64_t i;

	*itemNumLenStatisticArray = 0;
	for(i=0; i<itemNumQueryArray; i++)
	{
		lenStatisticArr[*itemNumLenStatisticArray].queryID = queryArray[i].queryID;
		lenStatisticArr[*itemNumLenStatisticArray].queryLen = queryArray[i].queryLen;
		(*itemNumLenStatisticArray) ++;
	}

	return SUCCESSFUL;
}

/**
 * Radix sort of the query lengths in descending order.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short radixSortOfLengths(queryLenStatistic_t *lenStatisticArray, queryLenStatistic_t *lenStatisitcBufArray, int64_t itemNumLenStatisticArray)
{
	struct partNode
	{
		int32_t curItemNum;
		int32_t totalItemNum;
		int32_t firstRow;
	};

	int64_t i, step, total;
	queryLenStatistic_t *data, *buf;
	struct partNode *part;
	int32_t partArrSize, stepBits, maxStepLen;
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
			buf = lenStatisticArray;
			data = lenStatisitcBufArray;
		}else
		{
			data = lenStatisticArray;
			buf = lenStatisitcBufArray;
		}

		// count the number
		if(memset(part, 0, partArrSize * sizeof(struct partNode))==NULL)
		{
			printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
			free(part);
			return FAILED;
		}
		for(i=0; i<itemNumLenStatisticArray; i++)
			//part[ (data[i].queryLen >> step) & bitMask ].totalItemNum ++;
			part[ bitMask - ((data[i].queryLen >> step) & bitMask) ].totalItemNum ++;

		// initialize the part array
		total = 0;
		for(i=0; i<partArrSize; i++)
		{
			part[i].firstRow = total;
			total += part[i].totalItemNum;
		}

		// copy the data to the right place
		for(i=0; i<itemNumLenStatisticArray; i++)
		{
			hashcode = bitMask - ((data[i].queryLen >> step) & bitMask);
			firstRow = part[hashcode].firstRow;
			curItemNum = part[hashcode].curItemNum;
			if(memcpy(buf+firstRow+curItemNum, data+i, sizeof(queryLenStatistic_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				free(part);
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
 * Save the sorted queries to file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short saveSortedQueries(char *sortedQueriesFile, queryLenStatistic_t *lenStatisticArray, int64_t itemNumLenStatisticArray, query_t *queryArray)
{
	FILE *fpQueries;
	int64_t i;

	fpQueries = fopen(sortedQueriesFile, "w");
	if(fpQueries==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, sortedQueriesFile);
		return FAILED;
	}

	// output the queries
	fprintf(fpQueries, "Rank\tQuery\tLength\n");
	for(i=0; i<itemNumLenStatisticArray; i++)
	{
		fprintf(fpQueries, "%ld\t%s\t%d\n", i+1, queryArray[lenStatisticArray[i].queryID-1].queryTitle, lenStatisticArray[i].queryLen);
	}

	fclose(fpQueries);
	fpQueries = NULL;

	return SUCCESSFUL;
}

/**
 * Compute the statistics of queries, including maximal size, N50 size, mean size, median size.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short computeLenStatistics(metrics_t *queryMetrics, queryLenStatistic_t *lenStatisticArray, int64_t itemNumLenStatisticArray)
{
	int i;
	int64_t N50_ArrIndex, totalRefBaseNum;
	double totalQueryLen, tmpTotalLen, N50_TotalLen;
	int validQueryNum;

	totalRefBaseNum = 0;
	for(i=0; i<queryMetrics->subjectNum; i++)
	{
		totalRefBaseNum += queryMetrics->baseBitArrays[i].totalRefBaseNum;
	}

	totalQueryLen = 0;
	validQueryNum = 0;
	for(i=0; i<itemNumLenStatisticArray; i++)
	{
		if(lenStatisticArray[i].queryLen>=minQueryLenThres)
		{
			totalQueryLen += lenStatisticArray[i].queryLen;
			validQueryNum ++;
		}
	}

	N50_ArrIndex = -1;
	N50_TotalLen = totalQueryLen * 0.5;
	tmpTotalLen = 0;
	for(i=0; i<itemNumLenStatisticArray; i++)
	{
		if(lenStatisticArray[i].queryLen>=minQueryLenThres)
		{
			tmpTotalLen += lenStatisticArray[i].queryLen;
			if(tmpTotalLen>=N50_TotalLen)
			{
				N50_ArrIndex = i;
				break;
			}
		}
	}

	if(N50_ArrIndex==-1)
	{
		printf("line=%d, In %s(), N50_ArrIndex=%lu, error!\n", __LINE__, __func__, N50_ArrIndex);
		return FAILED;
	}

	// save the statistics
	queryMetrics->lengthMetrics.totalNum = validQueryNum;
	queryMetrics->lengthMetrics.totalLen = (uint64_t)totalQueryLen;
	queryMetrics->lengthMetrics.totalRefBaseNum = totalRefBaseNum;
	if(totalRefBaseNum>0)
		queryMetrics->lengthMetrics.lengthRatio = (double)totalQueryLen / totalRefBaseNum;
	else
	{
		printf("line=%d, In %s(), totalRefBaseNum=%ld, error!\n", __LINE__, __func__, totalRefBaseNum);
		return FAILED;
	}
	queryMetrics->lengthMetrics.maxSize = lenStatisticArray[0].queryLen;
	queryMetrics->lengthMetrics.N50 = lenStatisticArray[N50_ArrIndex].queryLen;
	if(validQueryNum>0)
		queryMetrics->lengthMetrics.meanSize = ((double)totalQueryLen)/validQueryNum;
	else
	{
		printf("line=%d, In %s(), validQueryNum=%d, error!\n", __LINE__, __func__, validQueryNum);
		return FAILED;
	}
	queryMetrics->lengthMetrics.medianSize = lenStatisticArray[validQueryNum/2].queryLen;

	return SUCCESSFUL;
}

/**
 * Compute the genome covered ratio.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short computeReferenceCovRatio(metrics_t *queryMetrics, subject_t *subjectArray, int64_t itemNumSubjectArray, query_t *queryArray, int64_t itemNumQueryArray, char *blastnResultFile)
{
	if(fillReferenceBaseFlag(queryMetrics, refDeletionFile, subjectArray, itemNumSubjectArray, queryArray, itemNumQueryArray, blastnResultFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the reference base flags, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// compute the reference covered ratio
	if(computeRefCovRatioByBaseFlag(queryMetrics)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the reference covered ratio, error!\n", __LINE__, __func__);
		return FAILED;
	}


	return SUCCESSFUL;
}

/**
 * Fill the reference base flags.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short fillReferenceBaseFlag(metrics_t *queryMetrics, char *refDeletionFile, subject_t *subjectArray, int64_t itemNumSubjectArray, query_t *queryArray, int64_t itemNumQueryArray, char *blastnResultFile)
{
	int j, len, fileStatus, matchFlag;
	char line[LINE_CHAR_MAX+1], *pch;
	FILE *fpBlastnRe, *fpRefDel;
	char patStr[20][256], strandStr[2][20];
	int stage;
	int step;
	int stepMatchFlag;
	char queryTitle[LINE_CHAR_MAX+1], subjectTitle[LINE_CHAR_MAX+1];
	int32_t queryLen, subjectLen, matchLen, totalMatchLen, gapNum, strand;
	int32_t startQueryPos, endQueryPos, startSubjectPos, endSubjectPos;
	float matchPercent, gapPercent;
	char matchLen_ch[20], totalMatchLen_ch[20], matchPercent_ch[20], gapNum_ch[20], gapPercent_ch[20];
	char queryPos_ch[20], subjectPos_ch[20];
	char queryBases[LINE_CHAR_MAX+1], middleMatches[LINE_CHAR_MAX+1], subjectBases[LINE_CHAR_MAX+1];
	char *alignBaseBlock[3];
	int alignBlockLen, middleMatchesSkipLen;
	int64_t querySubjectNum;
	int validQuerySubjectFlag, validSegmentFlag;


	fpBlastnRe = fopen(blastnResultFile, "r");
	if(fpBlastnRe==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, blastnResultFile);
		return FAILED;
	}

	fpRefDel = fopen(refDeletionFile, "w");
	if(fpRefDel==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, refDeletionFile);
		return FAILED;
	}

	strcpy(patStr[0], "Query=");
	strcpy(patStr[1], "Length=");
	strcpy(patStr[2], "Subject=");
	strcpy(patStr[3], " Score =");
	strcpy(patStr[4], " Identities =");
	strcpy(patStr[5], " Strand=");

	strcpy(patStr[6], "Query");
	strcpy(patStr[7], "Sbjct");
	strcpy(patStr[8], "Lambda");
	strcpy(patStr[9], "***** No hits found *****");

	strcpy(strandStr[0], "Plus");
	strcpy(strandStr[1], "Minus");


	alignBaseBlock[0] = queryBases;
	alignBaseBlock[1] = middleMatches;
	alignBaseBlock[2] = subjectBases;

	querySubjectNum = 0;
	stage = MATCH_START_STAGE;
	step = MATCH_INFO_START_STEP;
	stepMatchFlag = MATCH_INFO_START_STEP;
	while(1)
	{
		if(readLine(&fileStatus, line, &len, LINE_CHAR_MAX, fpBlastnRe)==FAILED)
		{
			printf("line=%d, In %s(), cannot read a line, error!\n", __LINE__, __func__);
			return FAILED;
		}

		//printf("#### %s\n", line);

		if(len==0)
		{
			if(fileStatus==EOF_STATUS)
				break;
			else
				continue;
		}

		// parse the blastn result
		if(stage==MATCH_START_STAGE || stage==MATCH_FINISHED_STAGE)
		{  //0
			if(len>QUERY_ID_SKIP_NUM)
			{
				matchFlag = YES;
				for(j=0; j<QUERY_ID_SKIP_NUM; j++) // check for "QUERY_TITLE_STAGE"
				{
					if(patStr[0][j]!=line[j])
					{
						matchFlag = NO;
						break;
					}
				}

				if(matchFlag==YES)
				{
					pch = line + QUERY_ID_SKIP_NUM;
					while(*pch==' ')
					{
						pch ++;
					}
					//printf("queryTitle=%d\n", pch);
					strcpy(queryTitle, pch);

					querySubjectNum ++;

					// check the valid querySubject item
//					if(checkValidQuerySubject(&validQuerySubjectFlag, querySubjectNum, queryArr, itemNumQueryArr, subjectArr, itemNumSubjectArr)==FAILED)
//					{
//						printf("line=%d, In %s(), cannot check the valid querySubject item, error!\n", __LINE__, __func__);
//						return FAILED;
//					}

					stage = QUERY_TITLE_STAGE;

					continue;
				}
			}
		}else if(stage==QUERY_TITLE_STAGE)
		{ // 1
			if(len>QUERY_LEN_SKIP_NUM)
			{
				matchFlag = YES;
				for(j=0; j<QUERY_LEN_SKIP_NUM; j++) // check for "QUERY_LEN_STAGE"
				{
					if(patStr[1][j]!=line[j])
					{
						matchFlag = NO;
						break;
					}
				}

				if(matchFlag==YES)
				{
					pch = line + QUERY_LEN_SKIP_NUM;
					while(*pch==' ')
					{
						pch ++;
					}
					queryLen = atoi(pch);
					//printf("queryLen=%d\n", queryLen);

					stage = QUERY_LEN_STAGE;

					continue;
				}
				else
				{
					strcat(queryTitle, line);
					continue;
				}
			}
		}else if(stage==QUERY_LEN_STAGE)
		{ // 2
			if(len>SUBJECT_HEAD_SKIP_NUM)
			{
				matchFlag = YES;
				for(j=0; j<SUBJECT_HEAD_SKIP_NUM; j++) // check for "SUBJECT_HEAD_STAGE"
				{
					if(patStr[2][j]!=line[j])
					{
						matchFlag = NO;
						break;
					}
				}

				if(matchFlag==YES)
				{
					pch = line + SUBJECT_HEAD_SKIP_NUM;
					while(*pch==' ')
					{
						pch ++;
					}
					strcpy(subjectTitle, pch);

					stage = SUBJECT_HEAD_STAGE;
					continue;
				}
			}
		}else if(stage==SUBJECT_HEAD_STAGE)
		{ // 3
			if(len>SUBJECT_LEN_SKIP_NUM)
			{
				matchFlag = YES;
				for(j=0; j<SUBJECT_LEN_SKIP_NUM; j++) // check for "SUBJECT_LEN_STAGE"
				{
					if(patStr[1][j]!=line[j])
					{
						matchFlag = NO;
						break;
					}
				}

				if(matchFlag==YES)
				{
					pch = line + SUBJECT_LEN_SKIP_NUM;
					while(*pch==' ')
					{
						pch ++;
					}
					subjectLen = atoi(pch);

					//fprintf(fpParseResult, ">%s\t%d\t%s\t%d\n", queryTitle, queryLen, subjectTitle, subjectLen);// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

					stage = SUBJECT_LEN_STAGE;
					continue;
				}
				else
				{
					strcat(subjectTitle, line);
					continue;
				}
			}
		}else if(stage==SUBJECT_LEN_STAGE)
		{ // 4
			if(len>SCORE_SKIP_NUM)
			{
				matchFlag = YES;
				for(j=0; j<SCORE_SKIP_NUM; j++) // check for "SCORE_STEP"
				{
					if(patStr[3][j]!=line[j])
					{
						matchFlag = NO;
						break;
					}
				}

				if(matchFlag==YES)
				{
					stage = MATCH_INFO_STAGE;
					step = SCORE_STEP;

					continue;
				}
			}

			if(strcmp(line, patStr[9])==0)
			{
				stage = UNMATCH_INFO_STAGE;
			}

		}else if(stage==MATCH_INFO_STAGE)
		{ // 5
			if(step==SCORE_STEP)
			{
				if(len>IDENTITY_SKIP_NUM)
				{
					matchFlag = YES;
					for(j=0; j<IDENTITY_SKIP_NUM; j++) // check for "IDENTITY_GAP_STEP"
					{
						if(patStr[4][j]!=line[j])
						{
							matchFlag = NO;
							break;
						}
					}

					if(matchFlag==YES)
					{
						pch = line + IDENTITY_SKIP_NUM;
						while(*pch==' ') // skip the blanks
						{
							pch ++;
						}

						j = 0;
						while(*pch!='/') // get the matchLen
						{
							matchLen_ch[j++] = *pch;
							pch ++;
						}
						matchLen_ch[j] = '\0';
						matchLen = atoi(matchLen_ch);
						//printf("matchLen=%d\n", matchLen);

						pch ++;
						j = 0;
						while(*pch!=' ') // get the totalMatchLen
						{
							totalMatchLen_ch[j++] = *pch;
							pch ++;
						}
						totalMatchLen_ch[j] = '\0';
						totalMatchLen = atoi(totalMatchLen_ch);
						//printf("totalMatchLen=%d\n", totalMatchLen);

						while(*pch!='(')
						{
							pch ++;
						}

						pch ++;
						j = 0;
						while(*pch!='%') // get the matchPercent
						{
							matchPercent_ch[j++] = *pch;
							pch ++;
						}
						matchPercent_ch[j] = '\0';
						matchPercent = atof(matchPercent_ch) / 100.0;
						//printf("matchPercent=%f\n", matchPercent);

						while(*pch!='=')
						{
							pch ++;
						}
						pch ++;
						while(*pch==' ')
						{
							pch ++;
						}

						//pch ++;
						j = 0;
						while(*pch!='/') // get the gapNum
						{
							gapNum_ch[j++] = *pch;
							pch ++;
						}
						gapNum_ch[j] = '\0';
						gapNum = atoi(gapNum_ch);
						//printf("gapNum=%d\n", gapNum);

						while(*pch!='(')
						{
							pch ++;
						}

						pch ++;
						j = 0;
						while(*pch!='%') // get the matchPercent
						{
							gapPercent_ch[j++] = *pch;
							pch ++;
						}
						gapPercent_ch[j] = '\0';
						gapPercent = atof(gapPercent_ch) / 100.0;
						//printf("gapPercent=%f\n", gapPercent);

						if(totalMatchLen>=minQueryLenThres*2)
							validQuerySubjectFlag = YES;
						else
							validQuerySubjectFlag = NO;

						step = IDENTITY_GAP_STEP;

						continue;
					}
				}
			}else if(step==IDENTITY_GAP_STEP)
			{
				if(len>STRAND_SKIP_NUM)
				{
					matchFlag = YES;
					for(j=0; j<STRAND_SKIP_NUM; j++) // check for "IDENTITY_GAP_STEP"
					{
						if(patStr[5][j]!=line[j])
						{
							matchFlag = NO;
							break;
						}
					}

					if(matchFlag==YES)
					{
						pch = line + len - 1;
						while(*pch!='/') // skip the characters
						{
							pch --;
						}
						pch ++;

						if(strcmp(pch, strandStr[0])==0)
							strand = PLUS_STRAND;
						else
							strand = MINUS_STRAND;

						//printf("strand=%d\n", strand);

						step = STRAND_STEP;
						stepMatchFlag = STEP_START_FLAG;
						continue;
					}
				}
			}else if(step==STRAND_STEP)
			{
				if(len>0)
				{
					if(stepMatchFlag==STEP_START_FLAG)
					{
						pch = line;
						while(*pch!=' ') // skip the "Query"
							pch ++;
						while(*pch==' ') // skip the blanks
							pch ++;

						j = 0;
						while(*pch!=' ')
						{
							queryPos_ch[j++] = *pch;
							pch ++;
						}
						queryPos_ch[j] = '\0';
						startQueryPos = atoi(queryPos_ch);
						//printf("startQueryPos=%d\n", startQueryPos);

						// get the query bases
						while(*pch==' ') // skip the blanks
							pch ++;

						middleMatchesSkipLen = pch - line;

						alignBlockLen = 0;
						while(*pch!=' ')
						{
							alignBaseBlock[0][alignBlockLen++] = *pch;
							pch ++;
						}
						alignBaseBlock[0][alignBlockLen] = '\0';
						//printf("QueryBases: %s\n", alignBaseBlock[0]);

						// get the end query pos
						pch  = line + len - 1;
						while(*pch!=' ')
						{
							pch --;
						}
						pch ++;

						endQueryPos = atoi(pch);
						//printf("endQueryPos=%d\n", endQueryPos);

						stepMatchFlag = STEP_QUERY_FLAG;

						continue;

					}else if(stepMatchFlag==STEP_QUERY_FLAG)
					{
						// get the match flags for each base
						if(middleMatchesSkipLen>=len)
						{
							printf("line=%d, In %s(), middleMatchesSkipLen=%d >= len=%d, error!\n", __LINE__, __func__, middleMatchesSkipLen, len);
							return FAILED;
						}

						pch = line + middleMatchesSkipLen;
						j = 0;
						while(*pch!='\n' && *pch!='\0')
						{
							alignBaseBlock[1][j++] = *pch;
							pch ++;
						}
						alignBaseBlock[1][j] = '\0';
						//printf("            %s\n", alignBaseBlock[1]);

						stepMatchFlag = STEP_MIDDLE_FLAG;
						continue;
					}else if(stepMatchFlag==STEP_MIDDLE_FLAG)
					{
						pch = line;
						while(*pch!=' ')
						{
							pch ++;
						}
						while(*pch==' ')
						{
							pch ++;
						}

						j = 0;
						while(*pch!=' ')
						{
							subjectPos_ch[j++] = *pch;
							pch ++;
						}
						subjectPos_ch[j] = '\0';
						startSubjectPos = atoi(subjectPos_ch);
						//printf("startSubjectPos=%d\n", startSubjectPos);

						// get the subject bases
						while(*pch==' ') // skip the blanks
						{
							pch ++;
						}

						j = 0;
						while(*pch!=' ')
						{
							alignBaseBlock[2][j++] = *pch;
							pch ++;
						}
						alignBaseBlock[2][j] = '\0';
						//printf("SbjctBases: %s\n\n", alignBaseBlock[2]);

						// get the end subject pos
						pch = line + len - 1;
						while(*pch!=' ')
						{
							pch --;
						}
						pch ++;

						endSubjectPos = atoi(pch);
						//printf("endSubjectPos=%d\n", endSubjectPos);

						//printf("QueryBases: %s  %d  %d\n", alignBaseBlock[0], startQueryPos, endQueryPos);
						//printf("            %s\n", alignBaseBlock[1]);
						//printf("SbjctBases: %s  %d  %d\n\n", alignBaseBlock[2], startSubjectPos, endSubjectPos);


						// check the valid segments
//						if(validQuerySubjectFlag==YES)
//						{
//							if(checkValidSegment(&validSegmentFlag, querySubjectNum, strand, matchLen, totalMatchLen, gapNum, startQueryPos, startSubjectPos, queryArr, itemNumQueryArr, subjectArr, itemNumSubjectArr)==FAILED)
//							{
//								printf("line=%d, In %s(), cannot check the valid segments, error!\n", __LINE__, __func__);
//								return FAILED;
//							}
//						}
//
//						// update the base flags of the reference
//						if(validQuerySubjectFlag==YES && validSegmentFlag==YES)
						if(validQuerySubjectFlag==YES)
						{
							if(updateGenomeBaseFlag(queryMetrics, fpRefDel, alignBaseBlock, alignBlockLen, querySubjectNum, strand, startSubjectPos, endSubjectPos, startQueryPos, endQueryPos, subjectArray, itemNumSubjectArray, queryArray, itemNumQueryArray)==FAILED)
							{
								printf("line=%d, In %s(), cannot update the genome base flag, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}

						step = STEP_INFO_STEP;
						stepMatchFlag = STEP_SUBJECT_FLAG;
						continue;
					}
				}
			}else if(step==STEP_INFO_STEP)
			{
				if(stepMatchFlag==STEP_SUBJECT_FLAG)
				{
					matchFlag = YES;
					for(j=0; j<STEP_QUERY_SKIP_NUM; j++) // check for "STEP_QUERY_FLAG"
					{
						if(patStr[6][j]!=line[j])
						{
							matchFlag = NO;
							break;
						}
					}

					if(matchFlag==YES)
					{
						pch = line;
						while(*pch!=' ') // skip "Query"
						{
							pch ++;
						}
						while(*pch==' ') // skip the blanks
						{
							pch ++;
						}

						j = 0;
						while(*pch!=' ')
						{
							queryPos_ch[j++] = *pch;
							pch ++;
						}
						queryPos_ch[j] = '\0';
						startQueryPos = atoi(queryPos_ch);
						//printf("startQueryPos=%d\n", startQueryPos);

						// get the query bases
						while(*pch==' ') // skip the blanks
							pch ++;

						middleMatchesSkipLen = pch - line;

						alignBlockLen = 0;
						while(*pch!=' ')
						{
							alignBaseBlock[0][alignBlockLen++] = *pch;
							pch ++;
						}
						alignBaseBlock[0][alignBlockLen] = '\0';
						//printf("QueryBases: %s\n", alignBaseBlock[0]);

						pch = line + len - 1;
						while(*pch!=' ')
						{
							pch --;
						}
						pch ++;
						endQueryPos = atoi(pch);
						//printf("endQueryPos=%d\n", endQueryPos);

						stepMatchFlag = STEP_QUERY_FLAG;
						continue;
					}

				}else if(stepMatchFlag==STEP_QUERY_FLAG)
				{
					// get the match flags for each base
					if(middleMatchesSkipLen>=len)
					{
						printf("line=%d, In %s(), middleMatchesSkipLen=%d >= len=%d, error!\n", __LINE__, __func__, middleMatchesSkipLen, len);
						return FAILED;
					}

					pch = line + middleMatchesSkipLen;
					j = 0;
					while(*pch!='\n' && *pch!='\0')
					{
						alignBaseBlock[1][j++] = *pch;
						pch ++;
					}
					alignBaseBlock[1][j] = '\0';
					//printf("            %s\n", alignBaseBlock[1]);

					stepMatchFlag = STEP_MIDDLE_FLAG;
					continue;
				}else if(stepMatchFlag==STEP_MIDDLE_FLAG)
				{
					matchFlag = YES;
					for(j=0; j<STEP_SUBJECT_SKIP_NUM; j++) // check for "STEP_SUBJECT_FLAG"
					{
						if(patStr[7][j]!=line[j])
						{
							matchFlag = NO;
							break;
						}
					}

					if(matchFlag==YES)
					{
						pch = line;
						while(*pch!=' ')
						{
							pch ++;
						}
						while(*pch==' ')
						{
							pch ++;
						}

						j = 0;
						while(*pch!=' ')
						{
							subjectPos_ch[j++] = *pch;
							pch ++;
						}
						subjectPos_ch[j] = '\0';
						startSubjectPos = atoi(subjectPos_ch);
						//printf("startSubjectPos=%d\n", startSubjectPos);

						// get the subject bases
						while(*pch==' ') // skip the blanks
						{
							pch ++;
						}

						j = 0;
						while(*pch!=' ')
						{
							alignBaseBlock[2][j++] = *pch;
							pch ++;
						}
						alignBaseBlock[2][j] = '\0';
						//printf("SbjctBases: %s\n\n", alignBaseBlock[2]);

						pch = line + len - 1;
						while(*pch!=' ')
						{
							pch --;
						}
						pch ++;
						endSubjectPos = atoi(pch);
						//printf("endSubjectPos=%d\n", endSubjectPos);

						//printf("QueryBases: %s  %d  %d\n", alignBaseBlock[0], startQueryPos, endQueryPos);
						//printf("            %s\n", alignBaseBlock[1]);
						//printf("SbjctBases: %s  %d  %d\n\n", alignBaseBlock[2], startSubjectPos, endSubjectPos);

						// update the base flags of the reference
//						if(validQuerySubjectFlag==YES && validSegmentFlag==YES)
						if(validQuerySubjectFlag==YES)
						{
							if(updateGenomeBaseFlag(queryMetrics, fpRefDel, alignBaseBlock, alignBlockLen, querySubjectNum, strand, startSubjectPos, endSubjectPos, startQueryPos, endQueryPos, subjectArray, itemNumSubjectArray, queryArray, itemNumQueryArray)==FAILED)
							{
								printf("line=%d, In %s(), cannot update the genome base flag, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}

						stepMatchFlag = STEP_SUBJECT_FLAG;
						continue;
					}
				}

				//...
				//.... Lambda
				matchFlag = YES;
				for(j=0; j<LAMBDA_SKIP_NUM; j++) // check for "Lambda"
				{
					if(patStr[8][j]!=line[j])
					{
						matchFlag = NO;
						break;
					}
				}

				if(matchFlag==YES)
				{
					//stage = MATCH_FINISHED_STAGE;
					stage = LAMBDA_STAGE;

					// print the result !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					//fprintf(fpParseResult, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\n", startSubjectPos, endSubjectPos, startQueryPos, endQueryPos, strand, matchLen, totalMatchLen, gapNum, matchPercent);


					// save the match result to files


					continue;
				}

				// probe the finish of a step
				matchFlag = YES;
				for(j=0; j<SCORE_SKIP_NUM; j++) // check for "SCORE_STEP"
				{
					if(patStr[3][j]!=line[j])
					{
						matchFlag = NO;
						break;
					}
				}

				if(matchFlag==YES)
				{
					// print the result !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					//fprintf(fpParseResult, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\n", startSubjectPos, endSubjectPos, startQueryPos, endQueryPos, strand, matchLen, totalMatchLen, gapNum, matchPercent);

					// reset the parameters
					stepMatchFlag = STEP_FINISH_FLAG;
					step = SCORE_STEP;

					continue;
				}
			}
		}else if(stage==UNMATCH_INFO_STAGE)
		{
			if(len>LAMBDA_SKIP_NUM)
			{
				matchFlag = YES;
				for(j=0; j<LAMBDA_SKIP_NUM; j++) // check for "Lambda"
				{
					if(patStr[8][j]!=line[j])
					{
						matchFlag = NO;
						break;
					}
				}

				if(matchFlag==YES)
				{
					stage = LAMBDA_STAGE;
				}
			}
		}else if(stage==LAMBDA_STAGE)
		{ // 6
			if(len>LAMBDA_SKIP_NUM)
			{
				matchFlag = YES;
				for(j=0; j<LAMBDA_SKIP_NUM; j++) // check for "Lambda"
				{
					if(patStr[8][j]!=line[j])
					{
						matchFlag = NO;
						break;
					}
				}

				if(matchFlag==YES)
				{
					stage = MATCH_FINISHED_STAGE;
				}
			}
		}
	}

	fclose(fpRefDel);
	fpRefDel = NULL;
	fclose(fpBlastnRe);
	fpBlastnRe = NULL;

	return SUCCESSFUL;
}

/**
 * Check the valid querySubject item.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short checkValidQuerySubject(int *validQuerySubjectFlag, int64_t querySubjectNum, query_t *queryArray, int64_t itemNumQueryArray, subject_t *subjectArray, int64_t itemNumSubjectArray)
{
	int64_t queryID, subjectID;

	queryID = (querySubjectNum - 1) / itemNumSubjectArray + 1;
	subjectID = (querySubjectNum - 1) % itemNumSubjectArray + 1;

	if(queryArray[queryID-1].queryLen>=minQueryLenThres && queryArray[queryID-1].querySubjectNum>0)
	{
		if(queryArray[queryID-1].querySubArray[queryArray[queryID-1].bestMatchRow].subjectID == subjectID)
			*validQuerySubjectFlag = YES;
		else
			*validQuerySubjectFlag = NO;
	}else
	{
		*validQuerySubjectFlag = NO;
	}

	return SUCCESSFUL;
}

/**
 * Check the valid segment item.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short checkValidSegment(int *validSegmentFlag, int64_t querySubjectNum, int strand, int32_t matchLen, int32_t totalMatchLen, int32_t gapNum, int32_t startQueryPos, int32_t startSubjectPos, query_t *queryArray, int64_t itemNumQueryArray, subject_t *subjectArray, int64_t itemNumSubjectArray)
{
	int64_t i;
	int64_t queryID, subjectID;
	validSegment_t *pValidSegArray;
	int32_t validSegmentNum;

	queryID = (querySubjectNum - 1) / itemNumSubjectArray + 1;
	subjectID = (querySubjectNum - 1) % itemNumSubjectArray + 1;

	validSegmentNum = queryArray[queryID-1].querySubArray[queryArray[queryID-1].bestMatchRow].validSegmentNum;
	pValidSegArray = queryArray[queryID-1].querySubArray[queryArray[queryID-1].bestMatchRow].validSegArray;

	*validSegmentFlag = NO;
	for(i=0; i<validSegmentNum; i++)
	{
		if(pValidSegArray[i].strand==strand && pValidSegArray[i].matchLen==matchLen && pValidSegArray[i].totalMatchLen==totalMatchLen && pValidSegArray[i].gapNum==gapNum && pValidSegArray[i].startQueryPos==startQueryPos && pValidSegArray[i].startSubPos==startSubjectPos)
		{
			*validSegmentFlag = YES;
			break;
		}
	}

	return SUCCESSFUL;
}

/**
 * Update the genome base flag.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short updateGenomeBaseFlag(metrics_t *queryMetrics, FILE *fpRefDel, char **alignBaseBlock, int alignBlockLen, int querySubjectNum, int strand, int startSubjectPos, int endSubjectPos, int startQueryPos, int endQueryPos, subject_t *subjectArray, int64_t itemNumSubjectArray, query_t *queryArray, int64_t itemNumQueryArray)
{
	int i, rowID, colID;
	int queryID, subjectID, queryPos, subjectPos, baseStepSize, baseFlag;
	int originalBaseFlag;
	char queryBase, middleFlag, subjectBase;
	uint8_t *pBaseBitArray;


	queryID = (querySubjectNum - 1) / queryMetrics->subjectNum + 1;
	subjectID = (querySubjectNum - 1) % queryMetrics->subjectNum + 1;
	pBaseBitArray = queryMetrics->baseBitArrays[subjectID-1].subjectBitArray;


	// only the matched and perfect matched kinds are valid for coverage ratio
	if(ONLY_MATCHED_FOR_COV==YES)	// added 2013-02-22
	{
		if(queryArray[queryID-1].globalMatchKind==DISJUNCT_MATCH_KIND || queryArray[queryID-1].globalMatchKind==UNMATCHED_KIND)
		{
			return SUCCESSFUL;
		}
	}

	if(strand==PLUS_STRAND)
		baseStepSize = 1;
	else
		baseStepSize = -1;

	// initialize the queryPos and subjectPos
	queryPos = startQueryPos;
	subjectPos = startSubjectPos;

	for(i=0; i<alignBlockLen; i++)
	{
		queryBase = alignBaseBlock[0][i];
		middleFlag = alignBaseBlock[1][i];
		subjectBase = alignBaseBlock[2][i];

		if(queryBase=='-')
		{ // deletion in query, then the queryPos--
			queryPos --;
		}
		if(subjectBase=='-')
		{ // deletion in subject, then the subjectPos--
			subjectPos -= baseStepSize;
		}

		originalBaseFlag = pBaseBitArray[subjectPos-1];

		// change the lower case to upper case
		if(queryBase>=97) // lower case
			queryBase -= 32;
		if(subjectBase>=97) // lower case
			subjectBase -= 32;

		if(middleFlag=='|')
		{ // matched base, then assign the flag '1'
			baseFlag = 1;
		}else
		{
			switch(queryBase)
			{
				case 'A': rowID = 0; break;
				case 'C': rowID = 1; break;
				case 'G': rowID = 2; break;
				case 'T': rowID = 3; break;
				case 'N': rowID = 4; break;
				case '.': rowID = 5; break;
				case '-': rowID = 6; break;
				default: rowID = 4; break;
				//default: printf("line=%d, In %s(), unknown base %c (%d), error!\n", __LINE__, __func__, queryBase, queryBase); return FAILED;
			}

			switch(subjectBase)
			{
				case 'A': colID = 0; break;
				case 'C': colID = 1; break;
				case 'G': colID = 2; break;
				case 'T': colID = 3; break;
				case 'N': colID = 4; break;
				case '.': colID = 5; break;
				case '-': colID = 6; break;
				default: rowID = 4; break;
				//default: printf("line=%d, In %s(), unknown base %c (%d), error!\n", __LINE__, __func__, subjectBase, subjectBase); return FAILED;
			}

			baseFlag = baseFlagArr[rowID][colID];
		}

		// determine the cases of changing the reference base flag
		if(baseFlag<MIN_BASE_FLAG_DELETION)
		{
			if(originalBaseFlag==0 && baseFlag>=1)
			{ // apply the change
				pBaseBitArray[subjectPos-1] = baseFlag;
			}else if(originalBaseFlag>1 && baseFlag==1)
			{ // apply the change
				pBaseBitArray[subjectPos-1] = baseFlag;
			}
		}else
		{
			fprintf(fpRefDel, "%d\t%d\t%d\t%d\t%d\n", queryID, queryPos, subjectID, subjectPos, baseFlag);
		}

		queryPos ++;
		subjectPos += baseStepSize;
	}

	// ################################# Debug information ##############################
	queryPos --;
	subjectPos -= baseStepSize;
	if(queryPos!=endQueryPos)
	{
		printf("line=%d, In %s(), queryPos=%d != endQueryPos=%d, error!\n", __LINE__, __func__, queryPos, endQueryPos);
		return FAILED;
	}
	if(subjectPos!=endSubjectPos)
	{
		printf("line=%d, In %s(), subjectPos=%d != endSubjectPos=%d, error!\n", __LINE__, __func__, subjectPos, endSubjectPos);
		return FAILED;
	}
	// ################################# Debug information ##############################

	return SUCCESSFUL;
}

/**
 * Compute the reference covered ratio given filled reference base flag array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short computeRefCovRatioByBaseFlag(metrics_t *queryMetrics)
{
	int64_t i, j;
	int64_t coveredBaseNum, totalRefBaseNum;
	uint8_t *pSubjectBitArray;

	queryMetrics->lengthMetrics.coveredBaseNum = 0;
	for(i=0; i<queryMetrics->subjectNum; i++)
	{
		coveredBaseNum = 0;
		totalRefBaseNum = queryMetrics->baseBitArrays[i].totalRefBaseNum;
		pSubjectBitArray = queryMetrics->baseBitArrays[i].subjectBitArray;
		for(j=0; j<totalRefBaseNum; j++)
		{
			if(pSubjectBitArray[j]>0)
				coveredBaseNum ++;
		}
		queryMetrics->baseBitArrays[i].coveredBaseNum = coveredBaseNum;
		queryMetrics->baseBitArrays[i].coveredRatio = (double) coveredBaseNum / totalRefBaseNum;

		queryMetrics->lengthMetrics.coveredBaseNum += coveredBaseNum;
	}

	queryMetrics->lengthMetrics.coveredRatio = (double) queryMetrics->lengthMetrics.coveredBaseNum / queryMetrics->lengthMetrics.totalRefBaseNum;


	return SUCCESSFUL;
}
