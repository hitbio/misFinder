/*
 * incoming.c
 *
 *  Created on: Apr 27, 2012
 *      Author: david
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Convert matchInfo to a binary file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short convertMatchInfo(char *queryMatchInfoFile, char *parseResultFile)
{
	tmpQuery_t *tmpQueryArray;
	int64_t itemNumTmpQueryArray;
	queryMatchInfo_t *queryMatchInfoSet;

	if(initMemConvertMatchInfo(&queryMatchInfoSet, &tmpQueryArray, &itemNumTmpQueryArray, parseResultFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory for converting the match information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the data for tmpQuery array
	if(fillDataTmpQueryArray(tmpQueryArray, itemNumTmpQueryArray, queryMatchInfoSet, parseResultFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the data for tmpQuery array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// converting
	if(convertQueryMatchInfo(queryMatchInfoSet, tmpQueryArray, itemNumTmpQueryArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot convert the tmpQuery array to query match information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// output the result to binary file
	if(saveQueryMatchInfoToFile(queryMatchInfoFile, queryMatchInfoSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot save query match information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	freeMemConvertMatchInfo(&queryMatchInfoSet, &tmpQueryArray, &itemNumTmpQueryArray);

	return SUCCESSFUL;
}

/**
 * Initialize the memory for converting match information.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short initMemConvertMatchInfo(queryMatchInfo_t **queryMatchInfoSet, tmpQuery_t **tmpQueryArray, int64_t *itemNumTmpQueryArray, char *parseResultFile)
{
	*queryMatchInfoSet = (queryMatchInfo_t *) calloc (1, sizeof(queryMatchInfo_t));
	if((*queryMatchInfoSet)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the numbers for converting
	if(getItemNums(&((*queryMatchInfoSet)->itemNumSubjectArray), &((*queryMatchInfoSet)->itemNumQueryArray), itemNumTmpQueryArray, &((*queryMatchInfoSet)->itemNumMatchItemArray), parseResultFile)==FAILED)
	{
		printf("In %s(), cannot get the arrayNums from the parseResultFile, error!\n", __func__);
		return 1;
	}

	// allocate the memory for converting
	(*queryMatchInfoSet)->queryArray = (query_t *) calloc ((*queryMatchInfoSet)->itemNumQueryArray, sizeof(query_t));
	if((*queryMatchInfoSet)->queryArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	(*queryMatchInfoSet)->subjectArray = (subject_t *) calloc ((*queryMatchInfoSet)->itemNumSubjectArray, sizeof(subject_t));
	if((*queryMatchInfoSet)->subjectArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	(*queryMatchInfoSet)->matchItemArray = (matchItem_t *) calloc ((*queryMatchInfoSet)->itemNumMatchItemArray, sizeof(matchItem_t));
	if((*queryMatchInfoSet)->matchItemArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	(*queryMatchInfoSet)->potentMisassNum = 0;
	(*queryMatchInfoSet)->trueMisassNum = 0;
	(*queryMatchInfoSet)->SVNum = 0;


	*tmpQueryArray = (tmpQuery_t *) calloc (*itemNumTmpQueryArray, sizeof(tmpQuery_t));
	if((*tmpQueryArray)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Free the memory of converting match information.
 */
void freeMemConvertMatchInfo(queryMatchInfo_t **queryMatchInfoSet, tmpQuery_t **tmpQueryArray, int64_t *itemNumTmpQueryArray)
{
	int64_t i;

	for(i=0; i<*itemNumTmpQueryArray; i++) free((*tmpQueryArray)[i].queryTitle);
	free(*tmpQueryArray);
	*tmpQueryArray = NULL;

	*itemNumTmpQueryArray = 0;

	releaseQueryMatchInfo(queryMatchInfoSet);
}

/**
 * Get the itemNums	itemNumSubjectArr, itemNumQueryArr, itemNumQueryArr.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short getItemNums(int64_t *itemNumSubjectArray, int64_t *itemNumQueryArray, int64_t *itemNumTmpQueryArray, int64_t *itemNumMatchItemArray, char *parseResultFile)
{
	if(getSubjectNum(itemNumSubjectArray, parseResultFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the subject numbers, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(getQuerySubjectNumAndMatchItemNum(itemNumTmpQueryArray, itemNumMatchItemArray, parseResultFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the total querySubject numbers, error!\n", __LINE__, __func__);
		return FAILED;
	}

	*itemNumQueryArray = (*itemNumTmpQueryArray) / (*itemNumSubjectArray);

	return SUCCESSFUL;
}

/**
 * Get the subject number.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short getSubjectNum(int64_t *itemNumSubjectArray, char *parseResultFile)
{
	int32_t len, fileStatus, validFirstHeadFlag;
	char line[LINE_CHAR_MAX+1], *pch, tmpQueryTitle[LINE_CHAR_MAX+1];
	FILE *fpParseResult;


	if(getFirstQueryTitle(tmpQueryTitle, parseResultFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, parseResultFile);
		return FAILED;
	}


	fpParseResult = fopen(parseResultFile, "r");
	if(fpParseResult==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, parseResultFile);
		return FAILED;
	}

	*itemNumSubjectArray = 0;
	validFirstHeadFlag = NO;
	while(1)
	{
		if(readLine(&fileStatus, line, &len, LINE_CHAR_MAX, fpParseResult)==FAILED)
		{
			printf("In %s(), cannot read a line, error!\n", __func__);
			return FAILED;
		}

		if(len==0)
		{
			if(fileStatus==EOF_STATUS)
				break;
			else
				continue;
		}


		pch = line;
		if(*pch == '>')
		{
			pch = strtok(line+1, "\t");
			if(pch)
			{
				if(strcmp(tmpQueryTitle, pch)==0)
				{
					(*itemNumSubjectArray) ++;
				}else
				{
					break;
				}
			}else
			{
				printf("In %s(), cannot get queryTitle, error!\n", __func__);
				return FAILED;
			}

		}
	}

	fclose(fpParseResult);
	fpParseResult = NULL;

	return SUCCESSFUL;
}

/**
 * Get first query title.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short getFirstQueryTitle(char *queryTitle, char *parseResultFile)
{
	int len, fileStatus;
	char line[LINE_CHAR_MAX+1], *pch;
	FILE *fpParseResult;

	fpParseResult = fopen(parseResultFile, "r");
	if(fpParseResult==NULL)
	{
		printf("In %s(), cannot open file [ %s ], error!\n", __func__, parseResultFile);
		return FAILED;
	}

	while(1)
	{
		if(readLine(&fileStatus, line, &len, LINE_CHAR_MAX, fpParseResult)==FAILED)
		{
			printf("In %s(), cannot read a line, error!\n", __func__);
			return FAILED;
		}

		if(len==0)
		{
			if(fileStatus==EOF_STATUS)
				break;
			else
				continue;
		}

		pch = line;
		if((*pch) == '>')
		{
			pch = strtok(line+1, "\t");
			if(pch)
			{
				strcpy(queryTitle, pch);
				break;
			}else
			{
				printf("In %s(), cannot get queryTitle, error!\n", __func__);
				return FAILED;
			}
		}else
		{
			printf("In %s(), the  File format is wrong\n, error!\n", __func__);
			return FAILED;
		}
	}

	fclose(fpParseResult);
	fpParseResult = NULL;

	return SUCCESSFUL;
}

/**
 * Get the querySubject item number and the match item number.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short getQuerySubjectNumAndMatchItemNum(int64_t *tmpQuerySubjectNum, int64_t *itemNumMatchItemArray, char *parseResultFile)
{
	int len, fileStatus;
	char line[LINE_CHAR_MAX+1];
	FILE *fpParseResult;

	fpParseResult = fopen(parseResultFile, "r");
	if(fpParseResult==NULL)
	{
		printf("In %s(), cannot open file [ %s ], error!\n", __func__, parseResultFile);
		return FAILED;
	}

	*tmpQuerySubjectNum = 0;
	*itemNumMatchItemArray = 0;
	while(1)
	{
		if(readLine(&fileStatus, line, &len, LINE_CHAR_MAX, fpParseResult)==FAILED)
		{
			printf("In %s(), cannot read a line, error!\n", __func__);
			return FAILED;
		}

		if(len==0)
		{
			if(fileStatus==EOF_STATUS)
				break;
			else
				continue;
		}

		if(line[0]!='>')
			(*itemNumMatchItemArray) ++;
		else
			(*tmpQuerySubjectNum) ++;
	}

	fclose(fpParseResult);
	fpParseResult = NULL;

	return SUCCESSFUL;
}

/**
 * Fill data of tmpQuery array from the parseResultFile.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short fillDataTmpQueryArray(tmpQuery_t *tmpQueryArray, int64_t itemNumTmpQueryArray, queryMatchInfo_t *queryMatchInfoSet, char *parseResultFile)
{
	int64_t returnCode, tmpItemNum, firstRow;
	FILE *fpParsedResult;

	// fill subject array
	if(fillSubjectArray(queryMatchInfoSet->subjectArray, queryMatchInfoSet->itemNumSubjectArray, parseResultFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill subject array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	fpParsedResult = fopen(parseResultFile, "r");
	if(fpParsedResult==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, parseResultFile);
		return FAILED;
	}


	// fill the tmpQuery data
	tmpItemNum = 0;
	firstRow = 0;
	while(1)
	{
		// fill the single tmpQuery
		returnCode = fillSingleTmpQueryMatchInfo(tmpQueryArray+tmpItemNum, queryMatchInfoSet->matchItemArray+firstRow, queryMatchInfoSet->itemNumMatchItemArray, tmpItemNum+1, queryMatchInfoSet->subjectArray, queryMatchInfoSet->itemNumSubjectArray, fpParsedResult);
		if(returnCode==FAILED)
		{ // end of file
			break;
		}else if(returnCode==ERROR)
		{
			printf("line=%d, In %s(), cannot fill single tmpQuery data, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(tmpQueryArray[tmpItemNum].matchItemNum>0)
		{
			tmpQueryArray[tmpItemNum].firstRow = firstRow;
			firstRow += tmpQueryArray[tmpItemNum].matchItemNum;
		}

		tmpItemNum ++;
	}

	fclose(fpParsedResult);
	fpParsedResult = NULL;

	// #################################### Debug information ########################################
	if(tmpItemNum!=itemNumTmpQueryArray)
	{
		printf("line=%d, In %s(), tmpItemNum=%ld != itemNumTmpQueryArray=%ld, error!\n", __LINE__, __func__, tmpItemNum, itemNumTmpQueryArray);
		return FAILED;
	}
	// #################################### Debug information ########################################

	return SUCCESSFUL;
}

/**
 * Fill the subject array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short fillSubjectArray(subject_t *subjectArray, int64_t itemNumSubjectArray, char *parseResultFile)
{
	int subjectID, tmpItemNumSubjectArray, len, fileStatus;
	char line[LINE_CHAR_MAX+1], *pch;
	FILE *fpParseResult;

	fpParseResult = fopen(parseResultFile, "r");
	if(fpParseResult==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, parseResultFile);
		return FAILED;
	}

	tmpItemNumSubjectArray = 0;
	subjectID = 0;
	while(1)
	{
		if(readLine(&fileStatus, line, &len, LINE_CHAR_MAX, fpParseResult)==FAILED)
		{
			printf("In %s(), cannot read a line, error!\n", __func__);
			return FAILED;
		}

		if(len==0)
		{
			if(fileStatus==EOF_STATUS)
				break;
			else
				continue;
		}


		pch = line;
		if(*pch == '>')
		{
			pch = strtok(line+1, "\t");
			if(pch==NULL) // queryTitle
			{
				printf("line=%d, In %s(), cannot get queryTitle, error!\n", __LINE__, __func__);
				return FAILED;
			}

			pch = strtok(NULL, "\t");
			if(pch==NULL) // queryLen
			{
				printf("line=%d, In %s(), cannot get queryLen, error!\n", __LINE__, __func__);
				return FAILED;
			}

			pch = strtok(NULL, "\t");
			if(pch) // subjectTitle
			{
				subjectArray[tmpItemNumSubjectArray].subjectID = ++subjectID;
				subjectArray[tmpItemNumSubjectArray].subjectTitle = (char*) calloc (strlen(pch)+1, sizeof(char));
				if(subjectArray[tmpItemNumSubjectArray].subjectTitle==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				strcpy(subjectArray[tmpItemNumSubjectArray].subjectTitle, pch);
				subjectArray[tmpItemNumSubjectArray].subjectTitleLen = strlen(pch);
			}else
			{
				printf("line=%d, In %s(), cannot get subjectTitle, error!\n", __LINE__, __func__);
				return FAILED;
			}

			pch = strtok(NULL, "\t");
			if(pch) // subjectLen
			{
				subjectArray[tmpItemNumSubjectArray].subjectLen = atoi(pch);
			}else
			{
				printf("line=%d, In %s(), cannot get subjectLen, error!\n", __LINE__, __func__);
				return FAILED;
			}
			subjectArray[tmpItemNumSubjectArray].circularFlag = YES;

			tmpItemNumSubjectArray ++;
			if(tmpItemNumSubjectArray>=itemNumSubjectArray)
				break;
		}
	}

	fclose(fpParseResult);
	fpParseResult = NULL;

	return SUCCESSFUL;
}

/**
 * Fill data of the single tmpQuery item from parsedResult file in fasta file.
 *  @return:
 *  	If succeed, return SUCCESSFUL; if end of file, return FAILED; else, return ERROR.
 */
short fillSingleTmpQueryMatchInfo(tmpQuery_t *pTmpQuery, matchItem_t *pMatchItemArray, int64_t itemNumMatchItemArray, int64_t tmpQueryID, subject_t *subjectArray, int64_t itemNumSubjectArray, FILE *fpParsedResult)
{
	char ch, *pch;
	char line[LINE_CHAR_MAX+1];
	int32_t len;
	int64_t matchItemNum, subjectID;
	int32_t startSubPos, endSubPos, startQueryPos, endQueryPos, strand, matchLen, totalMatchLen, gapNum;
	double matchPercent;

	ch = fgetc(fpParsedResult);
	if(ch=='>')//start of a new contig
	{
		len = 0;
		while((ch=fgetc(fpParsedResult))!='\n') line[len++] = ch;
		line[len] = '\0';

		pch = line;
		pch = strtok(pch, "\t");
		if(pch)
		{ // query title
			pTmpQuery->queryTitle = (char *) calloc (strlen(pch)+1, sizeof(char));
			if(pTmpQuery->queryTitle==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return ERROR;
			}
			strcpy(pTmpQuery->queryTitle, pch);
		}else
		{
			printf("line=%d, In %s(), cannot get queryTitle, error!\n", __LINE__, __func__);
			return ERROR;
		}

		pch = strtok(NULL, "\t");
		if(pch)
		{ // query length
			pTmpQuery->queryLen = atoi(pch);
		}else
		{
			printf("line=%d, In %s(), cannot get queryTitle, error!\n", __LINE__, __func__);
			return ERROR;
		}

		pch = strtok(NULL, "\t");
		if(pch)
		{ // subject title
			subjectID = ((tmpQueryID-1) % itemNumSubjectArray) + 1;
			if(strcmp(pch, subjectArray[subjectID-1].subjectTitle)==0)
			{
				pTmpQuery->subjectID = subjectID;
			}else
			{
				printf("line=%d, In %s(), incorrect subjectTitle, error!\n", __LINE__, __func__);
				return ERROR;
			}
		}else
		{
			printf("line=%d, In %s(), cannot get queryTitle, error!\n", __LINE__, __func__);
			return ERROR;
		}

		pch = strtok(NULL, "\t");
		if(pch)
		{ // subject length
			subjectArray[subjectID-1].subjectLen = atoi(pch);
		}else
		{
			printf("line=%d, In %s(), cannot get queryTitle, error!\n", __LINE__, __func__);
			return ERROR;
		}
	}
	else // end of file, return FAILED
	{
		return FAILED;
	}

	matchItemNum = 0;
	while(1)
	{
		// check the first character of a line
		ch = fgetc(fpParsedResult);
		fseek(fpParsedResult, -1, SEEK_CUR);
		if(ch=='>' || ch==EOF)
		{
			break;
		}

		// get the match item information
		if(fscanf(fpParsedResult, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\n", &startSubPos, &endSubPos, &startQueryPos, &endQueryPos, &strand, &matchLen, &totalMatchLen, &gapNum, &matchPercent)!=COL_NUM)
		{
			printf("In %s(), cannot read the match item data, error!\n", __func__);
			return ERROR;
		}

		pMatchItemArray[matchItemNum].startSubPos = startSubPos;
		pMatchItemArray[matchItemNum].endSubPos = endSubPos;
		pMatchItemArray[matchItemNum].startQueryPos = startQueryPos;
		pMatchItemArray[matchItemNum].endQueryPos = endQueryPos;
		pMatchItemArray[matchItemNum].matchLen = matchLen;
		pMatchItemArray[matchItemNum].strand = strand;
		pMatchItemArray[matchItemNum].totalMatchLen = totalMatchLen;
		pMatchItemArray[matchItemNum].gapNum = gapNum;
		pMatchItemArray[matchItemNum].matchPercent = matchPercent;
		matchItemNum ++;

		if(matchItemNum>itemNumMatchItemArray)
		{
			printf("line=%d, In %s(), matchItemNum=%ld > itemNumMatchItemArray=%ld, error!\n", __LINE__, __func__, matchItemNum, itemNumMatchItemArray);
			return ERROR;
		}
	}

	pTmpQuery->matchItemNum = matchItemNum;

	return SUCCESSFUL;
}

/**
 * Convert the tmpQuery array to query match information.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short convertQueryMatchInfo(queryMatchInfo_t *queryMatchInfoSet, tmpQuery_t *tmpQueryArray, int64_t itemNumTmpQueryArray)
{
	int64_t i, j, queryID, tmpQuerySubjectNum, tmpSubjectRow;

	queryID = 0;
	for(i=0; i<itemNumTmpQueryArray; i++)
	{
		if(i%queryMatchInfoSet->itemNumSubjectArray==0)
		{ // add the item into the query array
			queryID ++;
			queryMatchInfoSet->queryArray[queryID-1].queryID = queryID;
			queryMatchInfoSet->queryArray[queryID-1].queryTitle = (char *) calloc (strlen(tmpQueryArray[i].queryTitle)+1, sizeof(char));
			if(queryMatchInfoSet->queryArray[queryID-1].queryTitle==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			strcpy(queryMatchInfoSet->queryArray[queryID-1].queryTitle, tmpQueryArray[i].queryTitle);
			queryMatchInfoSet->queryArray[queryID-1].queryTitleLen = strlen(queryMatchInfoSet->queryArray[queryID-1].queryTitle);
			queryMatchInfoSet->queryArray[queryID-1].queryLen = tmpQueryArray[i].queryLen;
			queryMatchInfoSet->queryArray[queryID-1].querySeq = NULL;
			queryMatchInfoSet->queryArray[queryID-1].circularFlag = -1;
			queryMatchInfoSet->queryArray[queryID-1].bestMatchRow = -1;
			queryMatchInfoSet->queryArray[queryID-1].globalMatchKind = -1;
			queryMatchInfoSet->queryArray[queryID-1].misassFlag = UNUSED_MISASS;

			// get the subjectNum
			tmpQuerySubjectNum = 0;
			for(j=0; j<queryMatchInfoSet->itemNumSubjectArray; j++)
			{
				if(tmpQueryArray[i+j].matchItemNum>0)
					tmpQuerySubjectNum ++;
			}
			queryMatchInfoSet->queryArray[queryID-1].querySubjectNum = tmpQuerySubjectNum;
			// allocate the querySubject item array
			if(tmpQuerySubjectNum>0)
			{
				queryMatchInfoSet->queryArray[queryID-1].querySubArray = (querySubject_t *) calloc (tmpQuerySubjectNum, sizeof(querySubject_t));
				if(queryMatchInfoSet->queryArray[queryID-1].querySubArray==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}else
			{
				queryMatchInfoSet->queryArray[queryID-1].querySubArray = NULL;
			}

			tmpSubjectRow = 0;
		}

		// fill the querySubject item
		if(tmpQueryArray[i].matchItemNum>0)
		{
			queryMatchInfoSet->queryArray[queryID-1].querySubArray[tmpSubjectRow].subjectID = tmpQueryArray[i].subjectID;
			queryMatchInfoSet->queryArray[queryID-1].querySubArray[tmpSubjectRow].matchItemNum = tmpQueryArray[i].matchItemNum;
			queryMatchInfoSet->queryArray[queryID-1].querySubArray[tmpSubjectRow].firstRow = tmpQueryArray[i].firstRow;
			queryMatchInfoSet->queryArray[queryID-1].querySubArray[tmpSubjectRow].matchKind = -1;
			queryMatchInfoSet->queryArray[queryID-1].querySubArray[tmpSubjectRow].circularFlag = -1;
			queryMatchInfoSet->queryArray[queryID-1].querySubArray[tmpSubjectRow].validSegmentNum = 0;
			queryMatchInfoSet->queryArray[queryID-1].querySubArray[tmpSubjectRow].validSegArray = NULL;
			tmpSubjectRow ++;
		}
	}

	return SUCCESSFUL;
}

/**
 * Release the queryMatchInfo structure.
 */
void releaseQueryMatchInfo(queryMatchInfo_t **queryMatchInfoSet)
{
	int64_t i, j, querySubjectNum;
	querySubject_t *pQuerySubArray;
	misInfo_t *misInfo, *misInfo_tmp;
	querySeq_t *newQuerySeqNode, *newQuerySeqNode_tmp;
	misassSeq_t *misassSeq, *misassSeq_tmp;

	if((*queryMatchInfoSet)==NULL)
		return;

	// release the subjectArray
	for(i=0; i<(*queryMatchInfoSet)->itemNumSubjectArray; i++)
	{
		free((*queryMatchInfoSet)->subjectArray[i].subjectTitle);
		(*queryMatchInfoSet)->subjectArray[i].subjectTitle = NULL;
		if((*queryMatchInfoSet)->subjectArray[i].subjectSeq)
		{
			free((*queryMatchInfoSet)->subjectArray[i].subjectSeq);
			(*queryMatchInfoSet)->subjectArray[i].subjectSeq = NULL;
		}
	}
	free((*queryMatchInfoSet)->subjectArray);
	(*queryMatchInfoSet)->subjectArray = NULL;
	(*queryMatchInfoSet)->itemNumSubjectArray = 0;

	// release the queryArray
	for(i=0; i<(*queryMatchInfoSet)->itemNumQueryArray; i++)
	{
		pQuerySubArray = (*queryMatchInfoSet)->queryArray[i].querySubArray;
		querySubjectNum = (*queryMatchInfoSet)->queryArray[i].querySubjectNum;
		for(j=0; j<querySubjectNum; j++)
		{
			if(pQuerySubArray[j].validSegmentNum>0)
			{
				free(pQuerySubArray[j].validSegArray);
				pQuerySubArray[j].validSegArray = NULL;
				pQuerySubArray[j].validSegmentNum = 0;
			}
		}

		free((*queryMatchInfoSet)->queryArray[i].queryTitle);
		(*queryMatchInfoSet)->queryArray[i].queryTitle = NULL;

		if((*queryMatchInfoSet)->queryArray[i].querySeq)
		{
			free((*queryMatchInfoSet)->queryArray[i].querySeq);
			(*queryMatchInfoSet)->queryArray[i].querySeq = NULL;
		}

		// queryReadSet array
		if((*queryMatchInfoSet)->queryArray[i].queryReadSetNum>0)
		{
			for(j=0; j<(*queryMatchInfoSet)->queryArray[i].queryReadSetNum; j++)
			{
				if((*queryMatchInfoSet)->queryArray[i].queryReadSetArray[j].queryReadNum>0)
				{
					free((*queryMatchInfoSet)->queryArray[i].queryReadSetArray[j].queryReadArray);
					(*queryMatchInfoSet)->queryArray[i].queryReadSetArray[j].queryReadArray = NULL;
					(*queryMatchInfoSet)->queryArray[i].queryReadSetArray[j].queryReadNum = 0;
				}
			}
			free((*queryMatchInfoSet)->queryArray[i].queryReadSetArray);
			(*queryMatchInfoSet)->queryArray[i].queryReadSetArray = NULL;
			(*queryMatchInfoSet)->queryArray[i].queryReadSetNum = 0;
		}

		if((*queryMatchInfoSet)->queryArray[i].querySubjectNum>0)
		{
			free((*queryMatchInfoSet)->queryArray[i].querySubArray);
			(*queryMatchInfoSet)->queryArray[i].querySubArray = NULL;
			(*queryMatchInfoSet)->queryArray[i].querySubjectNum = 0;
		}

		if((*queryMatchInfoSet)->queryArray[i].newQuerySeqNum>0)
		{
			newQuerySeqNode = (*queryMatchInfoSet)->queryArray[i].newQuerySeqList;
			while(newQuerySeqNode)
			{
				newQuerySeqNode_tmp = newQuerySeqNode->next;
				free(newQuerySeqNode->querySeq);
				free(newQuerySeqNode);
				newQuerySeqNode = newQuerySeqNode_tmp;
			}
			(*queryMatchInfoSet)->queryArray[i].newQuerySeqNum = 0;
			(*queryMatchInfoSet)->queryArray[i].newQuerySeqList = NULL;
			(*queryMatchInfoSet)->queryArray[i].tailNewQuerySeqList = NULL;
		}

		if((*queryMatchInfoSet)->queryArray[i].globalValidSegNum>0)
		{
			free((*queryMatchInfoSet)->queryArray[i].globalValidSegArray);
			(*queryMatchInfoSet)->queryArray[i].globalValidSegArray = NULL;
			(*queryMatchInfoSet)->queryArray[i].globalValidSegNum = 0;
		}

		if((*queryMatchInfoSet)->queryArray[i].misInfoItemNum>0)
		{
			misInfo = (*queryMatchInfoSet)->queryArray[i].misInfoList;
			while(misInfo)
			{
				misInfo_tmp = misInfo->next;
				if(misInfo->misType==QUERY_MISJOIN_KIND)
					free(misInfo->queryMargin);
				else
					free(misInfo->queryIndel);
				misassSeq = misInfo->misassSeqList;
				while(misassSeq)
				{
					misassSeq_tmp = misassSeq->next;
					free(misassSeq);
					misassSeq = misassSeq_tmp;
				}
				free(misInfo);
				misInfo = misInfo_tmp;
			}
			(*queryMatchInfoSet)->queryArray[i].misInfoList = NULL;
			(*queryMatchInfoSet)->queryArray[i].tailMisInfo = NULL;
			(*queryMatchInfoSet)->queryArray[i].misInfoItemNum = 0;
		}
	}

	free((*queryMatchInfoSet)->queryArray);
	(*queryMatchInfoSet)->queryArray = NULL;
	(*queryMatchInfoSet)->itemNumQueryArray = 0;

	free((*queryMatchInfoSet)->matchItemArray);
	(*queryMatchInfoSet)->matchItemArray = NULL;
	(*queryMatchInfoSet)->itemNumMatchItemArray = 0;

	free(*queryMatchInfoSet);
	*queryMatchInfoSet = NULL;
}


/**
 * Save the query match information to a binary file.
 *  File format:
 *  	(1) queryMatchInfo node;
 *  	(2) subjectArray and their titleArray;
 *  	(3) queryArray and their titleArray, querySubjectArray and the validSegArray;
 *  	(4) queryReadSetArray;
 *  	(5) globalSegArray;
 *  	(6) matchItemArray;
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short saveQueryMatchInfoToFile(char *queryMatchInfoFile, queryMatchInfo_t *queryMatchInfoSet)
{
	FILE *fpQueryMatchInfo;
	int64_t i, j;

	fpQueryMatchInfo = fopen(queryMatchInfoFile, "wb");
	if(fpQueryMatchInfo==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, queryMatchInfoFile);
		return FAILED;
	}

	// write the numbers of head items
	if(fwrite(queryMatchInfoSet, sizeof(queryMatchInfo_t), 1, fpQueryMatchInfo) != 1)
	{
		printf("line=%d, In %s(), cannot write information to  binary file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// write the subject array and the subject titles
	if(fwrite(queryMatchInfoSet->subjectArray, sizeof(subject_t), queryMatchInfoSet->itemNumSubjectArray, fpQueryMatchInfo) != queryMatchInfoSet->itemNumSubjectArray)
	{
		printf("line=%d, In %s(), cannot write information to  binary file, error!\n", __LINE__, __func__);
		return FAILED;
	}
	for(i=0; i < queryMatchInfoSet->itemNumSubjectArray; i++)
	{
		if(fwrite(queryMatchInfoSet->subjectArray[i].subjectTitle, sizeof(char), queryMatchInfoSet->subjectArray[i].subjectTitleLen+1, fpQueryMatchInfo)!=queryMatchInfoSet->subjectArray[i].subjectTitleLen+1)
		{
			printf("line=%d, In %s(), cannot write information to  binary file, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// write the query array, query titles, querySubArray
	if(fwrite(queryMatchInfoSet->queryArray, sizeof(query_t), queryMatchInfoSet->itemNumQueryArray, fpQueryMatchInfo) != queryMatchInfoSet->itemNumQueryArray)
	{
		printf("line=%d, In %s(), cannot write information to  binary file, error!\n", __LINE__, __func__);
		return FAILED;
	}
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		if(fwrite(queryMatchInfoSet->queryArray[i].queryTitle, sizeof(char), queryMatchInfoSet->queryArray[i].queryTitleLen+1, fpQueryMatchInfo)!=queryMatchInfoSet->queryArray[i].queryTitleLen+1)
		{
			printf("line=%d, In %s(), cannot write information to  binary file, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(queryMatchInfoSet->queryArray[i].querySubjectNum>0)
		{
			// save the querySubject items
			if(fwrite(queryMatchInfoSet->queryArray[i].querySubArray, sizeof(querySubject_t), queryMatchInfoSet->queryArray[i].querySubjectNum, fpQueryMatchInfo)!=queryMatchInfoSet->queryArray[i].querySubjectNum)
			{
				printf("line=%d, In %s(), cannot write information to  binary file, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// write the valid segment array
			for(j=0; j<queryMatchInfoSet->queryArray[i].querySubjectNum; j++)
			{
				if(queryMatchInfoSet->queryArray[i].querySubArray[j].validSegmentNum>0)
				{
					if(fwrite(queryMatchInfoSet->queryArray[i].querySubArray[j].validSegArray, sizeof(validSegment_t), queryMatchInfoSet->queryArray[i].querySubArray[j].validSegmentNum, fpQueryMatchInfo)!=queryMatchInfoSet->queryArray[i].querySubArray[j].validSegmentNum)
					{
						printf("line=%d, In %s(), cannot write information to  binary file, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
			}
		}

		// queryReadSet array
		if(queryMatchInfoSet->queryArray[i].queryReadSetNum>0)
		{
			if(fwrite(queryMatchInfoSet->queryArray[i].queryReadSetArray, sizeof(queryReadSet_t), queryMatchInfoSet->queryArray[i].queryReadSetNum, fpQueryMatchInfo)!=queryMatchInfoSet->queryArray[i].queryReadSetNum)
			{
				printf("line=%d, In %s(), cannot write information to  binary file, error!\n", __LINE__, __func__);
				return FAILED;
			}

			for(j=0; j<queryMatchInfoSet->queryArray[i].queryReadSetNum; j++)
			{
				if(queryMatchInfoSet->queryArray[i].queryReadSetArray[j].queryReadNum>0)
				{
					if(fwrite(queryMatchInfoSet->queryArray[i].queryReadSetArray[j].queryReadArray, sizeof(queryRead_t), queryMatchInfoSet->queryArray[i].queryReadSetArray[j].queryReadNum, fpQueryMatchInfo)!=queryMatchInfoSet->queryArray[i].queryReadSetArray[j].queryReadNum)
					{
						printf("line=%d, In %s(), cannot write information to  binary file, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
			}
		}

		// global segments
		if(queryMatchInfoSet->queryArray[i].globalValidSegNum>0)
		{
			if(fwrite(queryMatchInfoSet->queryArray[i].globalValidSegArray, sizeof(globalValidSeg_t), queryMatchInfoSet->queryArray[i].globalValidSegNum, fpQueryMatchInfo)!=queryMatchInfoSet->queryArray[i].globalValidSegNum)
			{
				printf("line=%d, In %s(), cannot write information to  binary file, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	// write the match item array
	if(fwrite(queryMatchInfoSet->matchItemArray, sizeof(matchItem_t), queryMatchInfoSet->itemNumMatchItemArray, fpQueryMatchInfo) != queryMatchInfoSet->itemNumMatchItemArray)
	{
		printf("line=%d, In %s(), cannot write information to  binary file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	fclose(fpQueryMatchInfo);
	fpQueryMatchInfo = NULL;

	return SUCCESSFUL;
}

/**
 * Load the binary file of query match information to memory.
 *  File format:
 *  	(1) queryMatchInfo node;
 *  	(2) subjectArray and their titleArray;
 *  	(3) queryArray and their titleArray, querySubjectArray and the validSegArray;
 *  	(4) queryReadSetArray;
 *  	(5) globalSegArray;
 *  	(6) matchItemArray.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short loadQueryMatchInfoFromFile(queryMatchInfo_t **queryMatchInfoSet, char *queryMatchInfoFile)
{
	FILE *fpQueryMatchInfo;
	int64_t i, j;

	*queryMatchInfoSet = (queryMatchInfo_t *) calloc (1, sizeof(queryMatchInfo_t));
	if(*queryMatchInfoSet==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	fpQueryMatchInfo = fopen(queryMatchInfoFile, "rb");
	if(fpQueryMatchInfo==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, queryMatchInfoFile);
		return FAILED;
	}

	// read the item numbers of head
	if(fread(*queryMatchInfoSet, sizeof(queryMatchInfo_t), 1, fpQueryMatchInfo) != 1)
	{
		printf("line=%d, In %s(), cannot read information to  binary file, error!\n", __LINE__, __func__);
		return FAILED;
	}


	// allocate the memory
	(*queryMatchInfoSet)->subjectArray = (subject_t *) calloc ((*queryMatchInfoSet)->itemNumSubjectArray, sizeof(subject_t));
	if((*queryMatchInfoSet)->subjectArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	(*queryMatchInfoSet)->queryArray = (query_t *) calloc ((*queryMatchInfoSet)->itemNumQueryArray, sizeof(query_t));
	if((*queryMatchInfoSet)->queryArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	(*queryMatchInfoSet)->matchItemArray = (matchItem_t *) calloc ((*queryMatchInfoSet)->itemNumMatchItemArray, sizeof(matchItem_t));
	if((*queryMatchInfoSet)->matchItemArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// read the subject array and the query titles
	if(fread((*queryMatchInfoSet)->subjectArray, sizeof(subject_t), (*queryMatchInfoSet)->itemNumSubjectArray, fpQueryMatchInfo) != (*queryMatchInfoSet)->itemNumSubjectArray)
	{
		printf("line=%d, In %s(), cannot read information to  binary file, error!\n", __LINE__, __func__);
		return FAILED;
	}
	for(i=0; i < (*queryMatchInfoSet)->itemNumSubjectArray; i++)
	{
		(*queryMatchInfoSet)->subjectArray[i].subjectTitle = (char *) calloc ((*queryMatchInfoSet)->subjectArray[i].subjectTitleLen+1, sizeof(char));
		if((*queryMatchInfoSet)->subjectArray[i].subjectTitle==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(fread((*queryMatchInfoSet)->subjectArray[i].subjectTitle, sizeof(char), (*queryMatchInfoSet)->subjectArray[i].subjectTitleLen+1, fpQueryMatchInfo)!=(*queryMatchInfoSet)->subjectArray[i].subjectTitleLen+1)
		{
			printf("line=%d, In %s(), cannot read information to  binary file, error!\n", __LINE__, __func__);
			return FAILED;
		}

		(*queryMatchInfoSet)->subjectArray[i].subjectSeq = NULL;
	}

	// read the queryArray and their titleArray, querySubjectArray and the validSegArray
	if(fread((*queryMatchInfoSet)->queryArray, sizeof(query_t), (*queryMatchInfoSet)->itemNumQueryArray, fpQueryMatchInfo) != (*queryMatchInfoSet)->itemNumQueryArray)
	{
		printf("line=%d, In %s(), cannot read information to  binary file, error!\n", __LINE__, __func__);
		return FAILED;
	}
	for(i=0; i<(*queryMatchInfoSet)->itemNumQueryArray; i++)
	{
		// allocate memory for queryTitle
		(*queryMatchInfoSet)->queryArray[i].queryTitle = (char *) calloc ((*queryMatchInfoSet)->queryArray[i].queryTitleLen+1, sizeof(char));
		if((*queryMatchInfoSet)->queryArray[i].queryTitle==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// fill queryTitle
		if(fread((*queryMatchInfoSet)->queryArray[i].queryTitle, sizeof(char), (*queryMatchInfoSet)->queryArray[i].queryTitleLen+1, fpQueryMatchInfo)!=(*queryMatchInfoSet)->queryArray[i].queryTitleLen+1)
		{
			printf("line=%d, In %s(), cannot read information to  binary file, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if((*queryMatchInfoSet)->queryArray[i].querySubjectNum>0)
		{
			(*queryMatchInfoSet)->queryArray[i].querySubArray = (querySubject_t *) calloc ((*queryMatchInfoSet)->queryArray[i].querySubjectNum, sizeof(querySubject_t));
			if((*queryMatchInfoSet)->queryArray[i].querySubArray==NULL)
			{
				printf("line=%d, In %s(), hello , cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(fread((*queryMatchInfoSet)->queryArray[i].querySubArray, sizeof(querySubject_t), (*queryMatchInfoSet)->queryArray[i].querySubjectNum, fpQueryMatchInfo)!=(*queryMatchInfoSet)->queryArray[i].querySubjectNum)
			{
				printf("line=%d, In %s(), cannot read information to  binary file, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// read the valid segments
			for(j=0; j<(*queryMatchInfoSet)->queryArray[i].querySubjectNum; j++)
			{
				if((*queryMatchInfoSet)->queryArray[i].querySubArray[j].validSegmentNum>0)
				{
					(*queryMatchInfoSet)->queryArray[i].querySubArray[j].validSegArray = (validSegment_t *) malloc ((*queryMatchInfoSet)->queryArray[i].querySubArray[j].validSegmentNum * sizeof(validSegment_t));
					if((*queryMatchInfoSet)->queryArray[i].querySubArray[j].validSegArray==NULL)
					{
						printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
						return FAILED;
					}
					if(fread((*queryMatchInfoSet)->queryArray[i].querySubArray[j].validSegArray, sizeof(validSegment_t), (*queryMatchInfoSet)->queryArray[i].querySubArray[j].validSegmentNum, fpQueryMatchInfo)!=(*queryMatchInfoSet)->queryArray[i].querySubArray[j].validSegmentNum)
					{
						printf("line=%d, In %s(), cannot read information from  binary file, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
			}
		}

		// queryReadSetArray
		if((*queryMatchInfoSet)->queryArray[i].queryReadSetNum>0)
		{
			(*queryMatchInfoSet)->queryArray[i].queryReadSetArray = (queryReadSet_t *) calloc ((*queryMatchInfoSet)->queryArray[i].queryReadSetNum, sizeof(queryReadSet_t));
			if((*queryMatchInfoSet)->queryArray[i].queryReadSetArray==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(fread((*queryMatchInfoSet)->queryArray[i].queryReadSetArray, sizeof(queryReadSet_t), (*queryMatchInfoSet)->queryArray[i].queryReadSetNum, fpQueryMatchInfo)!=(*queryMatchInfoSet)->queryArray[i].queryReadSetNum)
			{
				printf("line=%d, In %s(), cannot read information to  binary file, error!\n", __LINE__, __func__);
				return FAILED;
			}

			for(j=0; j<(*queryMatchInfoSet)->queryArray[i].queryReadSetNum; j++)
			{
				if((*queryMatchInfoSet)->queryArray[i].queryReadSetArray[j].queryReadNum>0)
				{
					(*queryMatchInfoSet)->queryArray[i].queryReadSetArray[j].queryReadArray = (queryRead_t *) calloc ((*queryMatchInfoSet)->queryArray[i].queryReadSetArray[j].queryReadNum, sizeof(queryRead_t));
					if((*queryMatchInfoSet)->queryArray[i].queryReadSetArray[j].queryReadArray==NULL)
					{
						printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
						return FAILED;
					}
					if(fread((*queryMatchInfoSet)->queryArray[i].queryReadSetArray[j].queryReadArray, sizeof(queryRead_t), (*queryMatchInfoSet)->queryArray[i].queryReadSetArray[j].queryReadNum, fpQueryMatchInfo)!=(*queryMatchInfoSet)->queryArray[i].queryReadSetArray[j].queryReadNum)
					{
						printf("line=%d, In %s(), cannot read information to  binary file, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
			}
		}

		// global segments
		if((*queryMatchInfoSet)->queryArray[i].globalValidSegNum>0)
		{
			(*queryMatchInfoSet)->queryArray[i].globalValidSegArray = (globalValidSeg_t *) calloc ((*queryMatchInfoSet)->queryArray[i].globalValidSegNum, sizeof(globalValidSeg_t));
			if((*queryMatchInfoSet)->queryArray[i].globalValidSegArray==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(fread((*queryMatchInfoSet)->queryArray[i].globalValidSegArray, sizeof(globalValidSeg_t), (*queryMatchInfoSet)->queryArray[i].globalValidSegNum, fpQueryMatchInfo)!=(*queryMatchInfoSet)->queryArray[i].globalValidSegNum)
			{
				printf("line=%d, In %s(), cannot read information to  binary file, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		(*queryMatchInfoSet)->queryArray[i].querySeq = NULL;
		(*queryMatchInfoSet)->queryArray[i].subRegSize = REG_COV_SIZE_THRES;
	}

	// read the match items
	if(fread((*queryMatchInfoSet)->matchItemArray, sizeof(matchItem_t), (*queryMatchInfoSet)->itemNumMatchItemArray, fpQueryMatchInfo) != (*queryMatchInfoSet)->itemNumMatchItemArray)
	{
		printf("line=%d, In %s(), cannot read information to  binary file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	fclose(fpQueryMatchInfo);
	fpQueryMatchInfo = NULL;

	return SUCCESSFUL;
}
