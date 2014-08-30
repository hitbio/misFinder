/*
 * mf.c
 *
 *  Created on: May 31, 2012
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Compute the mis-assemblies.
 *  @parameters:
 *  	operationFlag: 0 -- all (default); 1 -- only merge subjects; 2 -- only compute metrics.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
int computeGlobalMetrics(int32_t operationMode, char *outputPathName, char *inputQueryFileName, char *subjectsConfigFileName, int minQueryLenThreshold, double matchPercentThreshold, int32_t pairedModePara, char **readFilesPara, int32_t readFileNumPara, int32_t threadNumPara, int32_t indelSizeThresPara)
{
	struct timeval tp_start,tp_end;
	double time_used;
	gettimeofday(&tp_start,NULL);


	// initialize the global parameters
	if(initGlobalParas(operationMode, outputPathName, inputQueryFileName, subjectsConfigFileName, minQueryLenThreshold, matchPercentThreshold, pairedModePara, readFilesPara, readFileNumPara, threadNumPara, indelSizeThresPara)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the global parameters, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// copy queries
	if(copyQueryFile(inputQueryFile, inputQueryFileInit)==FAILED)
	{
		printf("line=%d, In %s(), cannot copy queries, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// merge the subject files
	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_MERGE || operationMode==OPERATION_MODE_METRICS)
	{
		printf("\n========== Begin merging the subjects ...\n");

		if(mergeRefSegmentsFasta(mergedSegFile, mergeSubjectsFile)==FAILED)
		{
			printf("line=%d, In %s(), cannot merge segments in fasta, error!\n", __LINE__, __func__);
			return FAILED;
		}

		printf("========== End merged the subjects.\n");
	}

	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_METRICS)
	{
		printf("\n========== Begin computing the query metrics ...\n");

		// run the blastn command to generate the alignment information
		if(generateBlastnResult(outputPathStr, inputBlastnFile, inputQueryFile, mergedSegFile, threadNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the alignment information, error!\n", __LINE__, __func__);
			return FAILED;
		}


		// parse the input blastn file
		if(parseBlastn(parseResultFile, inputBlastnFile)==FAILED)
		{
			printf("line=%d, In %s(), cannot parse the blastn result file, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// convert the query match information into a binary file
		if(convertMatchInfo(queryMatchInfoFile, parseResultFile)==FAILED)
		{
			printf("line=%d, In %s(), cannot convert the query match information, error!\n", __LINE__, __func__);
			return FAILED;
		}


		// classify the parsed result
		if(classifyQueries(perfectQueryFile, matchedQueryFile, disjunctQueryFile, unmatchedQueryFile, linkedQueryMatchInfoFile, queryMatchInfoFileNew, queryMatchInfoFile, inputQueryFile, mergedSegFile)==FAILED)
		{
			printf("line=%d, In %s(), cannot classify the queries, error!\n", __LINE__, __func__);
			return FAILED;
		}


		// output the statistics of queries
		if(getQueryMetrics(queryStatisticsFile, sortedQueryFile, queryMatchInfoFileNew, inputBlastnFile, mergedSegFile)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the statistics of queries, error!\n", __LINE__, __func__);
			return FAILED;
		}

		printf("========== End computed the query metrics.\n");
	}

	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_MISASS)
	{
		printf("\n========== Begin validating the mis-assemblies ...\n");

		// validate the potential mis-assembled queries
		if(validateMisassQueries(queryMatchInfoFileNew, inputQueryFile, mergedSegFile, readFilesInput, readFileNum, readsFileFormatType)==FAILED)
		{
			printf("line=%d, In %s(), cannot validate potential mis-assembled queries, error!\n", __LINE__, __func__);
			return FAILED;
		}

		printf("========== End validated the mis-assemblies.\n");
	}

	resetGlobalParas();


	gettimeofday(&tp_end,NULL);
	time_used = tp_end.tv_sec-tp_start.tv_sec + (double)(tp_end.tv_usec-tp_start.tv_usec)/1000000;
	printf("\nTotal Used Time: %.2f seconds.\n", time_used);

	return SUCCESSFUL;
}

/**
 * Initialize the global parameters.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short initGlobalParas(int32_t operationMode, char *outputPathName, char *queryFileName, char *mergeSubjectsFileName, int minQueryLenThreshold, double matchPercentThreshold, int32_t pairedModePara, char **readFilesPara, int32_t readFileNumPara, int32_t threadNumPara, int32_t indelSizeThresPara)
{
	int32_t i;

	// global paths and file names
	if(setGlobalPath(outputPathName)==FAILED)
	{
		printf("line=%d, In %s(), cannot set global paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// the subject merge files
	if(mergeSubjectsFileName!=NULL)
	{
		strcpy(mergeSubjectsFile, mergeSubjectsFileName);
	}else
	{
		printf("line=%d, In %s(), please specify the input subjects file, error!\n", __LINE__, __func__);
		return FAILED;
	}
	strcpy(mergedSegFile, outputPathStr);
	strcat(mergedSegFile, "mergedRefSegs");

	if(queryFileName!=NULL)
	{
		strcpy(inputQueryFileInit, queryFileName);

		strcpy(inputQueryFile, outputPathStr);
		strcat(inputQueryFile, "query.fa");
	}else
	{
		printf("line=%d, In %s(), please specify the input queries file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// the blastn files
	strcpy(inputBlastnFile, outputPathStr);
	strcat(inputBlastnFile, "blastnResult");

	// the metrics files
	strcpy(queryStatisticsFile, outputPathStr);
	strcat(queryStatisticsFile, "queryStatistics");

	strcpy(parseResultFile, outputPathStr);
	strcat(parseResultFile, "parseResult");
	strcpy(queryMatchInfoFile, outputPathStr);
	strcat(queryMatchInfoFile, "queryMatchInfo.bin");
	strcpy(queryMatchInfoFileNew, queryMatchInfoFile);
	strcat(queryMatchInfoFileNew, ".new");

	strcpy(linkedQueryMatchInfoFile, outputPathStr);
	strcat(linkedQueryMatchInfoFile, "linkedQueryMatchInfo");
	strcpy(sortedQueryFile, outputPathStr);
	strcat(sortedQueryFile, "sortedQueries");
	strcpy(refDeletionFile, outputPathStr);
	strcat(refDeletionFile, "refDeletion");

	strcpy(perfectQueryFile, outputPathStr);
	strcat(perfectQueryFile, "1_perfectQueries");
	strcpy(matchedQueryFile, outputPathStr);
	strcat(matchedQueryFile, "2_matchedQueries");
	strcpy(disjunctQueryFile, outputPathStr);
	strcat(disjunctQueryFile, "3_disjunctQueries");
	strcpy(unmatchedQueryFile, outputPathStr);
	strcat(unmatchedQueryFile, "4_unmatchedQueries");

	strcpy(errorsFile, outputPathStr);
	strcat(errorsFile, "result_errors");
	strcpy(svFile, outputPathStr);
	strcat(svFile, "result_sv");
	strcpy(misUncertainFile, outputPathStr);
	strcat(misUncertainFile, "result_warning");
	strcpy(gapFile, outputPathStr);
	strcat(gapFile, "result_gap");

	strcpy(newQueryFile, outputPathStr);
	strcat(newQueryFile, "query_cor.fa");

	readFileNum = readFileNumPara;
	for(i=0; i<readFileNum; i++)
	{
		readFilesInput[i] = (char *) calloc (256, sizeof(char));
		if(readFilesInput[i]==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		strcpy(readFilesInput[i], readFilesPara[i]);
	}

	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_MISASS)
	{
		// get reads file format type
		if(getReadsFileFormat(&readsFileFormatType, readFilesInput, readFileNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot get reads file formats, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}
	pairedMode = pairedModePara;
	hashTableSizeReadseq = HASH_TABLE_SIZE;
	reserveHashItemBlocksFlag = RESERVE_HASH_ITEM_READ_SET;

//	kmerSize = KMER_SIZE_DEFAULT;
//	hashTableSize = HASH_TABLE_SIZE;
//	entriesPerKmer = ((kmerSize-1) / 32) + 1;
//	if(kmerSize%32==0)
//	{
//		lastEntryMask = (uint64_t) -1;
//		lastEntryBaseNum = 32;
//	}else
//	{
//		lastEntryMask = (1LLU << ((kmerSize%32)<<1)) - 1;
//		lastEntryBaseNum = kmerSize % 32;
//	}

	// global variables
	if(minQueryLenThreshold>0)
	{
		minQueryLenThres = minQueryLenThreshold;
	}else
	{
		minQueryLenThres = MIN_QUERY_LEN_THRES;
		printf("The minimal contig size will be set to be %d bp by default.\n", minQueryLenThres);
	}

	shortQueryLenThres = SHORT_QUERY_LEN_THRES;

	if(matchPercentThreshold>0)
	{
		matchPercentThres = matchPercentThreshold;
	}else
	{
		matchPercentThres = MATCHED_PERCENT_THRES;
		printf("The matched percent threshold will be set to be %.2f by default.\n", matchPercentThres);
	}

	matchPercentFactor = MATCH_PERCENT_FACTOR;
	varyEndLenThres = VARY_LEN_THRES;
	minDisjunctDistanceThres = MIN_DISJUNCT_DISTANCE_THRES;
	minAlignedSegLenThres = MIN_ALIGNED_SEG_LEN_THRES;

	if(indelSizeThresPara>0)
		indelSizeThres = indelSizeThresPara;
	else
	{
		indelSizeThres = INDEL_SIZE_DEFAULT;
		printf("The minimal indel size will be set to be %d bp by default.\n", INDEL_SIZE_DEFAULT);
	}

	if(threadNumPara>0)
	{
		if(threadNumPara>sysconf(_SC_NPROCESSORS_ONLN))
			threadNum = sysconf(_SC_NPROCESSORS_ONLN);
		else
			threadNum = threadNumPara;
	}else
	{
		threadNum = sysconf(_SC_NPROCESSORS_ONLN);
		if(threadNum<=0)
		{
			printf("line=%d, In %s(), cannot get the thread number from CPU, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}
	//printf("The alignment will be run with %ld threads.\n", threadNum);

	return SUCCESSFUL;
}


/**
 * Set the global input and output path directory.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short setGlobalPath(const char *outPathStr)
{
	int outPathLen;
	struct stat st;

	if(outPathStr!=NULL)
		strcpy(outputPathStr, outPathStr);
	else
		strcpy(outputPathStr, "./");


	outPathLen = strlen(outputPathStr);
	if(outPathLen>=255)
	{
		printf("line=%d, In %s(), output path length=%d, error!\n", __LINE__, __func__, outPathLen);
		return FAILED;
	}

	if(outputPathStr[outPathLen-1]!='/')
	{
		outputPathStr[outPathLen] = '/';
		outputPathStr[outPathLen+1] = '\0';
	}

	if(stat(outputPathStr, &st)==-1)
	{
		if(mkdir(outputPathStr, 0755)==-1)
		{
			printf("line=%d, In %s(), cannot create directory [ %s ], error!\n", __LINE__, __func__, outputPathStr);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

void resetGlobalParas()
{
	int32_t i;

	matchPercentThres = 0;
	shortQueryLenThres = 0;
	minQueryLenThres = 0;

	for(i=0; i<readFileNum; i++)
	{
		free(readFilesInput[i]);
		readFilesInput[i] = NULL;
	}

	//remove(mergedSegFile);
}

/**
 * Get the read file format.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadsFileFormat(int *readsFileFormatType, char **readFilesInput, int32_t readFileNum)
{
	int32_t i, readFormat;
	FILE *fpRead;
	char ch;

	readFormat = -1;
	for(i=0; i<readFileNum; i++)
	{
		fpRead = fopen(readFilesInput[i], "r");
		if(fpRead==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readFilesInput[i]);
			return FAILED;
		}

		ch = fgetc(fpRead);
		if(ch=='>')
		{
			if(readFormat!=-1 && readFormat!=FILE_FORMAT_FASTA)
			{
				printf("Reads files cannot be in multiple format combined together, error!\n");
				return FAILED;
			}

			readFormat = FILE_FORMAT_FASTA;
		}else if(ch=='@')
		{
			if(readFormat!=-1 && readFormat!=FILE_FORMAT_FASTQ)
			{
				printf("Reads files cannot be in multiple format combined together, error!\n");
				return FAILED;
			}

			readFormat = FILE_FORMAT_FASTQ;
		}else
		{
			printf("Reads files must be in fasta or fastq format, error!\n");
			return FAILED;
		}

		fclose(fpRead);
	}

	if(readFormat!=-1)
		*readsFileFormatType = readFormat;
	else
	{
		*readsFileFormatType = -1;
		printf("Reads files must be in fasta or fastq format, error!\n");
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Copy queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short copyQueryFile(char *inputQueryFile, char *inputQueryFileInit)
{
	char *querySeq, queryHeadTitle[1000], *pch;
	FILE *fpQuery, *fpQueryInit;
	int64_t maxQueryLen, secQueryLen, queryID, queryLen, returnFlag;

	fpQueryInit = fopen(inputQueryFileInit, "r");
	if(fpQueryInit==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, inputQueryFileInit);
		return FAILED;
	}

	fpQuery = fopen(inputQueryFile, "w");
	if(fpQuery==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, inputQueryFile);
		return FAILED;
	}

	if(getMaxQueryLenFromFile(&maxQueryLen, inputQueryFileInit)==FAILED)
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

	//get queries in fasta format
	while((returnFlag=getSingleFastaItemFromFile(fpQueryInit, queryHeadTitle, querySeq, &queryLen))==SUCCESSFUL)
	{
		pch = queryHeadTitle;
		while(*pch)
		{
			if((*pch)=='\t')
			{
				*pch = '\0';
				break;
			}
			pch ++;
		}

		fprintf(fpQuery, ">%s\n%s\n", queryHeadTitle, querySeq);
	}

	free(querySeq);
	querySeq = NULL;

	fclose(fpQuery);
	fclose(fpQueryInit);

	return SUCCESSFUL;
}

