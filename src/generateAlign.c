/*
 * generateAlign.c
 *
 *  Created on: May 29, 2012
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Generate the blastn alignment result information.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short generateBlastnResult(const char *outputPathStr, const char *inputBlastnFile, const char *inputQueryFile, const char *mergedSegFile, int64_t threadNum)
{
	if(threadNum==1)
	{
		if(generateBlastnResultNoThread(inputBlastnFile, inputQueryFile, mergedSegFile)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the alignment result, error.\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		// initialize the thread parameters
		if(initThreadParas(&threadArr, &threadParaArr, threadNum, outputPathStr, mergedSegFile)==FAILED)
		{
			printf("line=%d, In %s(), cannot initialize the thread parameters, error.\n", __LINE__, __func__);
			return FAILED;
		}

		// divide the query files
		if(divideQueryFiles(threadParaArr, threadNum, inputQueryFile)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the sub-query files, error.\n", __LINE__, __func__);
			return FAILED;
		}

		// create threads to do the alignment
		if(createThreadsBlastn(threadArr, threadParaArr, threadNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot create threads, error.\n", __LINE__, __func__);
			return FAILED;
		}

		if(waitThreadsBlastn(threadArr, threadParaArr, threadNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot wait threads, error.\n", __LINE__, __func__);
			return FAILED;
		}

		// merge the alignment file into a whole one
		if(mergeBlastnResults(inputBlastnFile, threadParaArr, threadNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot merge blastn results, error.\n", __LINE__, __func__);
			return FAILED;
		}

		// free memory for thread parameters
		if(freeThreadParas(&threadArr, &threadParaArr, threadNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot free the thread parameters, error.\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Initialize the thread parameters for blastn alignment.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short initThreadParas(pthread_t **threadArray, threadPara_t **threadParaArray, int64_t threadNum, const char *outputPathStr, const char *mergedSegFile)
{
	int32_t i;

	*threadArray = (pthread_t*) calloc (threadNum, sizeof(pthread_t));
	if((*threadArray)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error.\n", __LINE__, __func__);
		return FAILED;
	}


	*threadParaArray = (threadPara_t*) calloc (threadNum, sizeof(threadPara_t));
	if((*threadParaArray)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error.\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<threadNum; i++)
	{
		strcpy((*threadParaArray)[i].queryFileName, outputPathStr);
		strcat((*threadParaArray)[i].queryFileName, "subquery_");
		sprintf((*threadParaArray)[i].queryFileName+strlen((*threadParaArray)[i].queryFileName), "%d", i);
		strcat((*threadParaArray)[i].queryFileName, ".fa");

		strcpy((*threadParaArray)[i].subjectFileName, mergedSegFile);

		strcpy((*threadParaArray)[i].blastnFileName, outputPathStr);
		strcat((*threadParaArray)[i].blastnFileName, "subBlastnResult_");
		sprintf((*threadParaArray)[i].blastnFileName+strlen((*threadParaArray)[i].blastnFileName), "%d", i);

		(*threadParaArray)[i].validFlag = NO;
	}

	return SUCCESSFUL;
}

/**
 * Free the thread parameters for blastn alignment.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short freeThreadParas(pthread_t **threadArray, threadPara_t **threadParaArray, int64_t threadNum)
{
	int32_t i;

	// remove sub-query files
	for(i=0; i<threadNum; i++)
	{
		if((*threadParaArray)[i].validFlag==YES)
		{
			if(remove((*threadParaArray)[i].queryFileName)!=0)
			{
				printf("line=%d, In %s(), cannot remove file [ %s ], error!\n", __LINE__, __func__, (*threadParaArray)[i].queryFileName);
				return FAILED;
			}

			if(remove((*threadParaArray)[i].blastnFileName)!=0)
			{
				printf("line=%d, In %s(), cannot remove file [ %s ], error!\n", __LINE__, __func__, (*threadParaArray)[i].blastnFileName);
				return FAILED;
			}
		}
	}

	free(*threadArray);
	*threadArray = NULL;

	free(*threadParaArray);
	*threadParaArray = NULL;

	return SUCCESSFUL;
}

/**
 * Divide query file into sub-query files for alignment.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short divideQueryFiles(threadPara_t *threadParaArray, int64_t threadNum, const char *inputQueryFile)
{
	int32_t sizeOrderFlag; // 1: descending order; 2: ascending order; 3: equal order

	// get the order of the query size
	if(getQueryOrderFlag(&sizeOrderFlag, inputQueryFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the order of queries, error.\n", __LINE__, __func__);
		return FAILED;
	}

	// compute the sub-query size for each file
	if(computeSizeSubQueries(threadParaArray, threadNum, sizeOrderFlag, inputQueryFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the queries, error.\n", __LINE__, __func__);
		return FAILED;
	}

	if(generateSubQueries(threadParaArray, threadNum, inputQueryFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the queries, error.\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get the size order of queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short getQueryOrderFlag(int32_t *sizeOrderFlag, const char *inputQueryFile)
{
	int64_t queryNum, sumQueryLen;
	int32_t *sizeArray;

	// get the query item number
	if(getSumLengthQuery(&sumQueryLen, &queryNum, inputQueryFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the queries, error.\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate the size array
	sizeArray = (int32_t *) calloc (queryNum, sizeof(int32_t));
	if(sizeArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error.\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the size array
	if(fillQuerySizeArray(sizeArray, queryNum, inputQueryFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the query size array, error.\n", __LINE__, __func__);
		return FAILED;
	}

	// determine the size order
	if(determineQuerySizeOrder(sizeOrderFlag, sizeArray, queryNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot determine the query size order, error.\n", __LINE__, __func__);
		return FAILED;
	}

	free(sizeArray);

	return SUCCESSFUL;
}

/**
 * Fill the query size array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short fillQuerySizeArray(int32_t *sizeArray, int64_t queryNum, const char *inputQueryFile)
{
	FILE *fpQuery;
	char ch;
	int64_t itemNum, queryLen;

	fpQuery = fopen(inputQueryFile, "r");
	if(fpQuery==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, inputQueryFile);
		return FAILED;
	}

	ch = fgetc(fpQuery);
	while(ch!='\n' && ch!=EOF) ch = fgetc(fpQuery);

	itemNum = 0;
	queryLen = 0;
	while((ch=fgetc(fpQuery))!=EOF)
	{
		if(ch=='>')
		{ // head title
			while(ch!='\n' && ch!=EOF) ch = fgetc(fpQuery);
			itemNum ++;

			sizeArray[itemNum-1] = queryLen;
			queryLen = 0;
		}else
		{ // base sequence
			if(ch!='\n')
				queryLen ++;
		}
	}
	fclose(fpQuery);

	if(queryLen>0)
	{
		itemNum ++;
		sizeArray[itemNum-1] = queryLen;
	}

	if(itemNum!=queryNum)
	{
		printf("line=%d, In %s(), itemNum=%ld, queryNum=%ld, error.\n", __LINE__, __func__, itemNum, queryNum);
		return FAILED;
	}

	//printf("itemNum=%ld\n", itemNum);

	return SUCCESSFUL;
}

/**
 * Determine the query size order.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short determineQuerySizeOrder(int32_t *sizeOrderFlag, int32_t *sizeArray, int64_t queryNum)
{
	int64_t i, j, sumArray[10], itemNumArray[10], sumArraySize, startRow, endRow, rowNum, entryID, queryRow;
	double averEntryNum;

	if(queryNum>=10)
		sumArraySize = 10;
	else
		sumArraySize = queryNum;

	averEntryNum = (double) queryNum / sumArraySize;

	queryRow = 0;
	for(entryID=0; entryID<sumArraySize; entryID++)
	{
		startRow = queryRow;
		endRow = round((entryID + 1) * averEntryNum);

		if(endRow>queryNum-1) endRow = queryNum - 1;
		if(endRow<startRow) endRow = startRow;
		rowNum = endRow - startRow + 1;
		queryRow = endRow + 1;

		//printf("[%ld]: startRow=%ld, endRow=%ld, rowNum=%ld\n", entryID, startRow, endRow, rowNum);

		sumArray[entryID] = 0;
		itemNumArray[entryID] = rowNum;
		for(j=startRow; j<=endRow; j++)
			sumArray[entryID] += sizeArray[j];
	}

//	for(entryID=0; entryID<sumArraySize; entryID++)
//		printf("[%ld]: sum=%ld, itemNum=%ld, aver=%.2f\n", entryID, sumArray[entryID], itemNumArray[entryID], (double)sumArray[entryID]/itemNumArray[entryID]);

	if((double)sumArray[0]/itemNumArray[0]*5<(double)sumArray[entryID-1]/itemNumArray[entryID-1])
		*sizeOrderFlag = 2;
	else
		*sizeOrderFlag = 1;

	return SUCCESSFUL;
}

/**
 * Compute the size of each sub-queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short computeSizeSubQueries(threadPara_t *threadParaArray, int64_t threadNum, int32_t sizeOrderFlag, const char *inputQueryFile)
{
	int64_t sumQueryLen, queryNum;

	// get the sum length, compute the average length
	if(getSumLengthQuery(&sumQueryLen, &queryNum, inputQueryFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the queries, error.\n", __LINE__, __func__);
		return FAILED;
	}

	if(computeQuerySubSum(threadParaArray, threadNum, sumQueryLen, sizeOrderFlag)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the queries, error.\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get the sum length of queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short getSumLengthQuery(int64_t *sumQueryLen, int64_t *queryNum, const char *inputQueryFile)
{
	FILE *fpQuery;
	char ch;

	fpQuery = fopen(inputQueryFile, "r");
	if(fpQuery==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, inputQueryFile);
		return FAILED;
	}

	*queryNum = 0;
	*sumQueryLen = 0;
	while((ch=fgetc(fpQuery))!=EOF)
	{
		if(ch=='>')
		{ // head title
			while(ch!='\n' && ch!=EOF) ch = fgetc(fpQuery);
			(*queryNum) ++;
		}else
		{ // base sequence
			if(ch!='\n')
				(*sumQueryLen) ++;
		}
	}

	fclose(fpQuery);

	if((*sumQueryLen)<=0)
	{
		printf("line=%d, In %s(), cannot get the total query length, error.\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Compute sub-sum size for each sub-query files.
 *  @return:
 *  	if succeed, return SUCCESSFUL; if end of file, return FAILED; else, return ERROR.
 */
short computeQuerySubSum(threadPara_t *threadParaArray, int64_t threadNum, int64_t sumQueryLen, int32_t sizeOrderFlag)
{
	int32_t i, j, averSize, startRow, endRow, rowNum;
	int64_t totalSubQueryLen, tmp;
	double averEntryLen, averRatioSample, sumAdjustRatio;
	double sampleRatio[40]=
	{ 1.4097484277, 1.514527027,  1.4650326797, 1.3584848485, 1.4097484277, 1.3263313609, 1.2452777778, 1.4368589744, 1.4461290323, 1.3342261905,
	  1.4650326797, 1.2452777778, 1.2956647399, 1.1986631016, 1.2452777778, 1.3667682927, 1.3263313609, 1.3185294118, 1.2956647399, 1.3031976744,
	  1.3031976744, 1.3031976744, 1.2663841808, 1.2522346369, 1.1554123711, 1.0673809524, 1.1263819095, 1.1151741294, 0.9300829876, 1.0051569507,
	  0.9745652174, 0.9457805907, 1.1494871795, 0.7547138047, 0.9538297872, 0.8654440154, 0.6157967033, 0.6441091954, 0.4412401575, 0.226185671
	};

	// smooth sample ratio array
	if(smoothSampleRatio(sampleRatio, 40)==FAILED)
	{
		printf("line=%d, In %s(), cannot smooth sample ratio array, error.\n", __LINE__, __func__);
		return FAILED;
	}

	averSize = sumQueryLen / threadNum + 1;
	averEntryLen = 40.0 / threadNum;

	for(i=0; i<threadNum; i++)
	{
		startRow = round(i * averEntryLen);
		endRow = round((i + 1) * averEntryLen) - 1;
		if(startRow<0) startRow = 0;
		if(endRow<0) endRow = 0;
		else if(endRow>39) endRow = 39;

		if(startRow>endRow) endRow = startRow;
		rowNum = endRow - startRow + 1;

		averRatioSample = 0;
		for(j=startRow; j<=endRow; j++)
			averRatioSample += sampleRatio[j];
		averRatioSample /= rowNum;

		threadParaArray[i].initQuerySubSum = averRatioSample * averSize;

		//printf("[%d].initQuerySubSum=%ld\n", i, threadParaArray[i].initQuerySubSum);
	}

	totalSubQueryLen = 0;
	for(i=0; i<threadNum; i++)
	{
		totalSubQueryLen += threadParaArray[i].initQuerySubSum;
	}

	sumAdjustRatio = (double)sumQueryLen / totalSubQueryLen;
	//printf("sumQueryLen=%ld, totalSubQueryLen=%ld, sumAdjustRatio=%f\n", sumQueryLen, totalSubQueryLen, sumAdjustRatio);

	// reverse the array
	if(sizeOrderFlag==2)
	{
		for(i=0; i<threadNum/2; i++)
		{
			tmp = threadParaArray[i].initQuerySubSum;
			threadParaArray[i].initQuerySubSum = threadParaArray[threadNum-i-1].initQuerySubSum;
			threadParaArray[threadNum-i-1].initQuerySubSum = tmp;
		}
	}

	for(i=0; i<threadNum; i++)
	{
		threadParaArray[i].initQuerySubSum *= sumAdjustRatio;
		threadParaArray[i].initQuerySubSum ++;
		//printf("[%d].initQuerySubSum=%ld\n", i, threadParaArray[i].initQuerySubSum);
	}

	return SUCCESSFUL;
}

/**
 * Smooth sample ratio array.
 *  @return:
 *  	if succeed, return SUCCESSFUL; if end of file, return FAILED; else, return ERROR.
 */
short smoothSampleRatio(double *sampleRatioArray, int32_t arraySize)
{
	int32_t i, j, startRow, endRow, rowNum;
	double resultArray[arraySize], sum, aver;

	for(i=0; i<arraySize; i++)
	{
		startRow = i - 1;
		endRow = i + 1;
		if(startRow<0) startRow = 0;
		if(endRow>arraySize-1) endRow = arraySize - 1;

		rowNum = endRow - startRow + 1;
		sum = 0;
		for(j=startRow; j<=endRow; j++)
		{
			sum += sampleRatioArray[j];
		}
		sum /= rowNum;

		resultArray[i] = sum;
	}

	for(i=0; i<arraySize; i++) sampleRatioArray[i] = resultArray[i];

	return SUCCESSFUL;
}

/**
 * Generate sub-query files.
 *  @return:
 *  	if succeed, return SUCCESSFUL; if end of file, return FAILED; else, return ERROR.
 */
short generateSubQueries(threadPara_t *threadParaArray, int64_t threadNum, const char *inputQueryFile)
{
	FILE *fpQuery, *fpSubQuery;
	int64_t queryLen, sumSubQueryLen, maxQueryLen, subFileID, returnFlag;
	char *querySeq, queryHeadTitle[1000];

	fpQuery = fopen(inputQueryFile, "r");
	if(fpQuery==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, inputQueryFile);
		return FAILED;
	}

	if(getMaxQueryLenFromFile(&maxQueryLen, inputQueryFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the maximal query length, error!\n", __LINE__, __func__);
		return FAILED;
	}

	querySeq = (char *) malloc ((maxQueryLen+1)*sizeof(char));
	if(querySeq==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}


	for(subFileID=0; subFileID<threadNum; subFileID++) threadParaArray[subFileID].querySubSum = 0;

	// create the sub-query file
	subFileID = 0;
	fpSubQuery = fopen(threadParaArray[subFileID].queryFileName, "w");
	if(fpSubQuery==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, threadParaArray[0].queryFileName);
		return FAILED;
	}

	//get queries in fasta format
	sumSubQueryLen = 0;
	while((returnFlag=getSingleFastaItemFromFile(fpQuery, queryHeadTitle, querySeq, &queryLen))==SUCCESSFUL)
	{
		sumSubQueryLen += queryLen;
		fprintf(fpSubQuery, ">%s\n%s\n", queryHeadTitle, querySeq);

		if(sumSubQueryLen>threadParaArray[subFileID].initQuerySubSum)
		{
			//printf("sumFile[%ld]: %ld\n", subFileID, sumSubQueryLen);
			fclose(fpSubQuery);
			threadParaArray[subFileID].querySubSum = sumSubQueryLen;
			subFileID ++;
			sumSubQueryLen = 0;

			if(subFileID<threadNum)
			{
				fpSubQuery = fopen(threadParaArray[subFileID].queryFileName, "w");
				if(fpSubQuery==NULL)
				{
					printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, threadParaArray[subFileID].queryFileName);
					return FAILED;
				}
			}
		}
	}
	if(subFileID<threadNum)
	{
		fclose(fpSubQuery);

		if(sumSubQueryLen>0)
		{
			//printf("sumFile[%ld]: %ld\n", subFileID, sumSubQueryLen);
			threadParaArray[subFileID].querySubSum = sumSubQueryLen;
		}else
			remove(threadParaArray[subFileID].queryFileName);
	}else
	{
		printf("line=%d, In %s(), invalid file ID %ld, error!\n", __LINE__, __func__, subFileID);
		return FAILED;
	}

	free(querySeq);
	querySeq = NULL;

	fclose(fpQuery);
	fpQuery = NULL;

	// validate the sub-files
	for(subFileID=0; subFileID<threadNum; subFileID++)
	{
		if(threadParaArray[subFileID].querySubSum>0)
			threadParaArray[subFileID].validFlag = YES;
		else
			threadParaArray[subFileID].validFlag = NO;
	}

	// handle the error situation
	if(returnFlag==ERROR)
	{
		printf("line=%d, In %s(), cannot get query, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get the maximal query length from query file.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxQueryLenFromFile(int64_t *maxQueryLen, const char *queryFileName)
{
	FILE *fpQuery;
	int64_t maxSeqLen, queryLen;
	char ch, headTitle[1000];

	fpQuery = fopen(queryFileName, "r");
	if(fpQuery==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, queryFileName);
		return FAILED;
	}

	*maxQueryLen = 0;
	ch = fgetc(fpQuery);
	if(ch!='>')
	{
		printf("The query file [ %s ] was not in fasta format, error!\n", queryFileName);
		return FAILED;
	}

	// get the head title of a query
	fscanf(fpQuery, "%s\n", headTitle);

	maxSeqLen = 0;
	queryLen = 0;
	while((ch=fgetc(fpQuery))!=EOF)
	{
		if(ch!='>')
		{ // base sequence
			if(ch!='\n')
			{
				queryLen ++;
			}
		}else
		{ // ends of a query
			if(queryLen>maxSeqLen)
			{
				maxSeqLen = queryLen;
			}
			queryLen = 0;
		}
	}

	// process the last query
	if(queryLen>0)
	{
		if(queryLen>maxSeqLen)
		{
			maxSeqLen = queryLen;
		}
		queryLen = 0;
	}else
	{
		printf("line=%d, In %s(), queryLen=%ld, error!\n", __LINE__, __func__, queryLen);
		return FAILED;
	}

	fclose(fpQuery);
	fpQuery = NULL;

	*maxQueryLen = maxSeqLen;
	if(maxQueryLen<=0)
	{
		printf("line=%d, In %s(), maxQueryLen=%ld, error!\n", __LINE__, __func__, *maxQueryLen);
		return FAILED;
	}

	return SUCCESSFUL;
}


/**
 * Get a single fasta item from the fasta file.
 *  @return:
 *  	if succeed, return SUCCESSFUL; if end of file, return FAILED; else, return ERROR.
 */
short getSingleFastaItemFromFile(FILE *fpQuery, char *queryHeadTitle, char *querySeq, int64_t *queryLen)
{
	char ch, *pseq, base;

	*queryLen = 0;
	ch = fgetc(fpQuery);
	if(ch=='>') //start of a new contig
	{
		ch = fgetc(fpQuery);
		pseq = queryHeadTitle;
		while(ch!='\n' && ch!=EOF)
		{
			*pseq = ch;
			pseq ++;

			ch = fgetc(fpQuery);
		}
		*pseq = '\0';

	}else // end of file, return FAILED
	{
		return FAILED;
	}

	pseq = querySeq;
	ch = fgetc(fpQuery);
	while(ch!='>' && ch!=EOF)
	{
		if(ch!='\n')
		{
			switch(ch)
			{
				case 'A':
				case 'C':
				case 'G':
				case 'T':
				case 'N':
					base = ch; break;
				case 'a':
				case 'c':
				case 'g':
				case 't':
				case 'n':
					base = ch - 32; break;
				case '.':
					base = 'N'; break;
				default: base = ch;
			}
			*pseq = base;
			pseq ++;
			(*queryLen) ++;
		}
		ch = fgetc(fpQuery);
	}
	*pseq = '\0';

	if(ch=='>')
	{
		fseek(fpQuery, -1, SEEK_CUR);
	}

	return SUCCESSFUL;
}

/**
 * Create the threads for blastn alignment.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short createThreadsBlastn(pthread_t *threadArray, threadPara_t *threadParaArray, int64_t threadNum)
{
	int32_t i, ret, validNum;

	validNum = 0;
	for(i=0; i<threadNum; i++)
		if(threadParaArray[i].validFlag==YES)
			validNum ++;

	printf("The BLASTN alignment will be performed using %d threads.\n", validNum);

	for(i=0; i<threadNum; i++)
	{
		if(threadParaArray[i].validFlag==YES)
		{
			ret = pthread_create(threadArray+i, NULL, (void  *) generateBlastnResultThread, threadParaArray+i);
			if(ret!=0)
			{
				printf("line=%d, In %s(), cannot create threads for alignment, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Wait the threads for blastn alignment.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short waitThreadsBlastn(pthread_t *threadArray, threadPara_t *threadParaArray, int64_t threadNum)
{
	int32_t i;

	for(i=0; i<threadNum; i++)
	{
		if(threadParaArray[i].validFlag==YES)
			pthread_join(threadArray[i], NULL);
	}

	for(i=0; i<threadNum; i++)
		if(threadParaArray[i].validFlag==YES && threadParaArray[i].successFlag==FAILED)
			return FAILED;

	return SUCCESSFUL;
}

/**
 * Generate the blastn alignment result information.
 *  @return:
 *  	If succeeds, threadPara->successFlag=SUCCESSFUL; otherwise, threadPara->successFlag=FAILED.
 */
void generateBlastnResultThread(threadPara_t *threadPara)
{
	int returnCode;
	char blastnCommand[LINE_CHAR_MAX+1];
	char queryOption[LINE_CHAR_MAX+1];
	char subjectOption[LINE_CHAR_MAX+1];
	char outputOption[LINE_CHAR_MAX+1];
	char ctrlOption[LINE_CHAR_MAX+1];

	strcpy(queryOption, "-query ");
	strcat(queryOption, threadPara->queryFileName);
	strcpy(subjectOption, "-subject ");
	strcat(subjectOption, threadPara->subjectFileName);
	strcpy(outputOption, "-out ");
	strcat(outputOption, threadPara->blastnFileName);
	strcpy(ctrlOption, "-best_hit_overhang 0.1");

	strcpy(blastnCommand, "blastn");
	strcat(blastnCommand, " ");
	strcat(blastnCommand, queryOption);
	strcat(blastnCommand, " ");
	strcat(blastnCommand, subjectOption);
	strcat(blastnCommand, " ");
	strcat(blastnCommand, ctrlOption);
	strcat(blastnCommand, " ");
	strcat(blastnCommand, outputOption);

	printf("BLASTN command: %s\n", blastnCommand);

	returnCode = system(blastnCommand);
	if(returnCode==0)
	{
		threadPara->successFlag = SUCCESSFUL;
	}else
	{
		printf("Please run the correct blastn command or check whether the blastn was correctly installed.\n");
		threadPara->successFlag = FAILED;
	}
}

/**
 * Generate the blastn alignment result information without thread.
 *  @return:
 *  	If succeeds, threadPara->successFlag=SUCCESSFUL; otherwise, threadPara->successFlag=FAILED.
 */
short generateBlastnResultNoThread(const char *inputBlastnFile, const char *inputQueryFile, const char *mergedSegFile)
{
	int returnCode;
	char blastnCommand[LINE_CHAR_MAX+1];
	char queryOption[LINE_CHAR_MAX+1];
	char subjectOption[LINE_CHAR_MAX+1];
	char outputOption[LINE_CHAR_MAX+1];
	char ctrlOption[LINE_CHAR_MAX+1];

	strcpy(queryOption, "-query ");
	strcat(queryOption, inputQueryFile);
	strcpy(subjectOption, "-subject ");
	strcat(subjectOption, mergedSegFile);
	strcpy(outputOption, "-out ");
	strcat(outputOption, inputBlastnFile);
	strcpy(ctrlOption, "-best_hit_overhang 0.1");

	strcpy(blastnCommand, "blastn");
	strcat(blastnCommand, " ");
	strcat(blastnCommand, queryOption);
	strcat(blastnCommand, " ");
	strcat(blastnCommand, subjectOption);
	strcat(blastnCommand, " ");
	strcat(blastnCommand, ctrlOption);
	strcat(blastnCommand, " ");
	strcat(blastnCommand, outputOption);

	printf("BLASTN command: %s\n", blastnCommand);

	returnCode = system(blastnCommand);
	if(returnCode!=0)
	{
		printf("Please run the correct blastn command or check whether the blastn was correctly installed.\n");
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Merge blastn alignment results.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short mergeBlastnResults(const char *inputBlastnFile, threadPara_t *threadParaArray, int64_t threadNum)
{
	FILE *fpBlastn, *fpSubBlastn;
	int32_t i;
	char ch;

	fpBlastn = fopen(inputBlastnFile, "w");
	if(fpBlastn==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, inputBlastnFile);
		return FAILED;
	}

	for(i=0; i<threadNum; i++)
	{
		if(threadParaArray[i].validFlag==YES)
		{
			fpSubBlastn = fopen(threadParaArray[i].blastnFileName, "r");
			if(fpSubBlastn==NULL)
			{
				printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, threadParaArray[i].blastnFileName);
				return FAILED;
			}

			while(!feof(fpSubBlastn))
			{
				ch = fgetc(fpSubBlastn);
				fprintf(fpBlastn, "%c", ch);
			}
			fprintf(fpBlastn, "\n");

			fclose(fpSubBlastn);
		}
	}

	fclose(fpBlastn);

	return SUCCESSFUL;
}
