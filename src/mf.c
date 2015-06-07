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
int computeGlobalMetrics(int32_t operationMode, char *outputPathName, char *configFilePara, int minQueryLenThreshold, double matchPercentThreshold, int32_t threadNumPara, int32_t indelSizeThresPara)
{
	struct timeval tp_start,tp_end;
	double time_used;
	gettimeofday(&tp_start,NULL);


	// initialize the global parameters
	if(initGlobalParas(operationMode, outputPathName, configFilePara, minQueryLenThreshold, matchPercentThreshold, threadNumPara, indelSizeThresPara)==FAILED)
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

		if(mergeRefSegmentsFasta(mergedSegFile, subjectsFile)==FAILED)
		{
			printf("line=%d, In %s(), cannot merge segments in fasta, error!\n", __LINE__, __func__);
			return FAILED;
		}

		printf("========== End merged the subjects.\n");
	}

	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_METRICS)
	{
		printf("\n========== Begin Blastn alignment ...\n");

		if(generateAlignResult(outputPathStr, queryMatchInfoFile, inputBlastnFile, inputQueryFile, mergedSegFile, threadNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the alignment information, error!\n", __LINE__, __func__);
			return FAILED;
		}

		printf("========== End Blastn alignment.\n");

		if(operationMode==OPERATION_MODE_METRICS)
		{
			// output the statistics of queries
			if(getQueryMetrics(outputPathStr, queryMatchInfoFile, inputBlastnFile, mergedSegFile)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the statistics of queries, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else
		{
			// compute query statistics after validation
			if(computeLenStatisticsFromFastaFile(inputQueryFile, "\nQuery statistics BEFORE mis-assembly validation:")==FAILED)
			{
				printf("line=%d, In %s(), cannot compute the statistics of queries, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_MISASS)
	{
		printf("\n========== Begin mis-assembly validation ...\n");

		// validate the potential mis-assembled queries
		if(validateMisassQueries(outputPathStr, newQueryFile, queryMatchInfoFile, inputQueryFile, mergedSegFile, readFileList)==FAILED)
		{
			printf("line=%d, In %s(), cannot validate potential mis-assembled queries, error!\n", __LINE__, __func__);
			return FAILED;
		}

		printf("========== End mis-assembly validation.\n");

		// compute query statistics after validation
		if(computeLenStatisticsFromFastaFile(newQueryFile, "\nQuery statistics AFTER mis-assembly validation:")==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the statistics of queries, error!\n", __LINE__, __func__);
			return FAILED;
		}
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
short initGlobalParas(int32_t operationMode, char *outputPathName, char *configFilePara, int minQueryLenThreshold, double matchPercentThreshold, int32_t threadNumPara, int32_t indelSizeThresPara)
{
	// global paths and file names
	if(setGlobalPath(outputPathName)==FAILED)
	{
		printf("line=%d, In %s(), cannot set global paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	strcpy(configFile, configFilePara);

	strcpy(subjectsFile, outputPathStr);
	strcat(subjectsFile, "subjects");

	// parse configuration file
	if(parseConfigFile(inputQueryFileInit, subjectsFile, &readFileList, configFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot parse configuration file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// check the files to be validate
	if(checkFilesInConfigfile(operationMode, inputQueryFileInit, subjectsFile, readFileList)==FAILED)
	{
		printf("line=%d, In %s(), cannot check configuration file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// output the configuration information
	if(outputConfigInfo(operationMode, inputQueryFileInit, subjectsFile, readFileList)==FAILED)
	{
		printf("line=%d, In %s(), cannot output configuration information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// set file names
	strcpy(mergedSegFile, outputPathStr);
	strcat(mergedSegFile, "mergedRefSegs");
	strcpy(inputQueryFile, outputPathStr);
	strcat(inputQueryFile, "query.fa");
	strcpy(newQueryFile, outputPathStr);
	strcat(newQueryFile, "query_cor.fa");

	// the blastn files
	strcpy(inputBlastnFile, outputPathStr);
	strcat(inputBlastnFile, "blastnResult");

	strcpy(queryMatchInfoFile, outputPathStr);
	strcat(queryMatchInfoFile, "queryMatchInfo.bin");

	hashTableSizeReadseq = HASH_TABLE_SIZE;

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
	readFile_t *readFileNode, *readFileTmp;

	matchPercentThres = 0;
	shortQueryLenThres = 0;
	minQueryLenThres = 0;

	readFileNode = readFileList;
	while(readFileNode)
	{
		readFileTmp = readFileNode->next;
		for(i=0; i<readFileNode->readFileNum; i++)
			free(readFileNode->readFiles[i]);
		free(readFileNode->readFiles);
		free(readFileNode);
		readFileNode = readFileTmp;
	}
	readFileList = NULL;
}

/**
 * Parse the configuration file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short parseConfigFile(char *inputQueryFileInit, char *subjectSegsFile, readFile_t **readFileList, char *configFile)
{
	FILE *fpConfig, *fpSubjectSegs;
	char line[LINE_CHAR_MAX+1];
	int32_t i, fileStatus, len, startLineRow, despLineFlag, blankLineFlag, querySectionFlag, refSectionFlag, readSectionFlag;

	fpConfig = fopen(configFile, "r");
	if(fpConfig==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, configFile);
		return FAILED;
	}

	fpSubjectSegs = fopen(subjectSegsFile, "w");
	if(fpSubjectSegs==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, subjectSegsFile);
		return FAILED;
	}

	querySectionFlag = NO;
	refSectionFlag = NO;
	readSectionFlag = NO;
	while(1)
	{
		if(readLine(&fileStatus, line, &len, LINE_CHAR_MAX, fpConfig)==FAILED)
		{
			printf("line=%d, In %s(), cannot read a line, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(len==0)
		{
			if(fileStatus==EOF_STATUS)
				break;
			else
				continue;
		}

		despLineFlag = NO;
		blankLineFlag = YES;
		for(i=0; i<len; i++)
		{
			if(line[i]!=' ' && line[i]!='\t')
			{
				if(line[i]=='#')
					despLineFlag = YES;
				blankLineFlag = NO;
				startLineRow = i;
				break;
			}
		}

		if(despLineFlag==NO && blankLineFlag==NO)
		{
			if(querySectionFlag==YES)
			{
				if(strcmp(line+startLineRow, "END")==0)
					querySectionFlag = NO;
				else
					strcpy(inputQueryFileInit, line+startLineRow);
			}else if(refSectionFlag==YES)
			{
				if(strcmp(line+startLineRow, "END")==0)
					refSectionFlag = NO;
				else
					fprintf(fpSubjectSegs, "%s\n", line+startLineRow);
			}else if(readSectionFlag==YES)
			{
				if(strcmp(line+startLineRow, "END")==0)
					readSectionFlag = NO;
				else
				{
					// add the file to read file list
					if(addToReadFileList(readFileList, line+startLineRow)==FAILED)
					{
						printf("line=%d, In %s(), cannot add read file to list, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
			}else if(strcmp(line+startLineRow, "QUERY")==0)
				querySectionFlag = YES;
			else if(strcmp(line+startLineRow, "REF")==0)
				refSectionFlag = YES;
			else if(strcmp(line+startLineRow, "READS")==0)
				readSectionFlag = YES;
		}
	}

	fclose(fpConfig);
	fclose(fpSubjectSegs);

	return SUCCESSFUL;
}

/**
 * Add read files to list.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addToReadFileList(readFile_t **readFileList, char *readFiles)
{
	readFile_t *readFileNode, *readFileNodeTmp;
	char readTypeStr[256], *readFileNames[100];
	int32_t i, readFileNum, len, len_tmp, startRow;
	int32_t readType, pairedMode, readsFileFormatType;

	readFileNum = 0;

	len = strlen(readFiles);
	for(i=0; i<len; i++)
	{
		if(readFiles[i]=='=')
		{
			len_tmp = i;
			strncpy(readTypeStr, readFiles, len_tmp);
			readTypeStr[len_tmp] = '\0';

			break;
		}
	}

	i ++;
	while(i<len)
	{
		for(; i<len; ++i)
			if(readFiles[i]!=' ')
				break;
		startRow = i;
		for(i=startRow; i<len; i++)
		{
			if(i==len-1 || readFiles[i]==' ' || readFiles[i]=='\t' || readFiles[i]=='\n')
			{
				readFileNames[readFileNum] = (char*) calloc (256, sizeof(char));
				if(readFileNames[readFileNum]==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(i==len-1)
					len_tmp = i - startRow + 1;
				else
					len_tmp = i - startRow;
				strncpy(readFileNames[readFileNum], readFiles+startRow, len_tmp);
				readFileNames[readFileNum][len_tmp] = '\0';
				readFileNum ++;

				i ++;
				break;
			}
		}
	}

	if(strcmp(readTypeStr, "PE")==0)
		readType = PE_READ_TYPE;
	else if(strcmp(readTypeStr, "MP")==0)
		readType = MP_READ_TYPE;
	else if(strcmp(readTypeStr, "PE_IN")==0)
		readType = MP_IN_READ_TYPE;
	else if(strcmp(readTypeStr, "LONG_PE")==0)
		readType = LONG_PE_READ_TYPE;
	else if(strcmp(readTypeStr, "LONG_SE")==0)
		readType = LONG_SE_READ_TYPE;

	// get reads file format type
	if(getReadsFileFormat(&readsFileFormatType, readFileNames, readFileNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot get reads file formats, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(readFileNum==1)
	{
		if(readType==PE_READ_TYPE || readType==MP_READ_TYPE || readType==MP_IN_READ_TYPE || readType==LONG_PE_READ_TYPE)
			pairedMode = 2;
		else if(readType==LONG_SE_READ_TYPE)
			pairedMode = 0;

	}else if(readFileNum==2)
	{
		if(readType==PE_READ_TYPE || readType==MP_READ_TYPE || readType==MP_IN_READ_TYPE || readType==LONG_PE_READ_TYPE)
			pairedMode = 1;
		else if(readType==LONG_SE_READ_TYPE)
		{
			printf("Please specify the correct read files.\n");
			return FAILED;
		}
	}else
	{
		printf("Please specify the correct read files.\n");
		return FAILED;
	}

	// add new node to list
	readFileNode = (readFile_t*) calloc (1, sizeof(readFile_t));
	if(readFileNode==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	readFileNode->readsFileFormatType = readsFileFormatType;
	readFileNode->readType = readType;
	readFileNode->pairedMode = pairedMode;
	readFileNode->readFileNum = readFileNum;
	readFileNode->next = NULL;
	readFileNode->readFiles = (char**) calloc (readFileNum, sizeof(char*));
	if(readFileNode->readFiles==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	for(i=0; i<readFileNum; i++)
	{
		readFileNode->readFiles[i] = (char*) calloc (256, sizeof(char));
		if(readFileNode->readFiles[i]==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		strcpy(readFileNode->readFiles[i], readFileNames[i]);
	}

	if((*readFileList)==NULL)
		*readFileList = readFileNode;
	else
	{
		readFileNodeTmp = *readFileList;
		while(readFileNodeTmp)
		{
			if(readFileNodeTmp->next==NULL)
			{
				readFileNodeTmp->next = readFileNode;
				break;
			}
			readFileNodeTmp = readFileNodeTmp->next;
		}
	}

	// free memory
	for(i=0; i<readFileNum; i++)
		free(readFileNames[i]);

	return SUCCESSFUL;
}

/**
 * Get the read file format.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadsFileFormat(int32_t *readsFileFormatType, char **readFilesInput, int32_t readFileNum)
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
 * Check the files in Config file, including the query file, subject file, and read file list.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short checkFilesInConfigfile(int32_t operationMode, char *inputQueryFileInit, char *subjectsFile, readFile_t *readFileList)
{
	int32_t i, j;
	struct stat st;
	FILE *subjectFileTmp;
	int32_t len, fileStatus, fileNumTmp;
	char fileNameTmp[LINE_CHAR_MAX+1];
	readFile_t *readFileNode;

	// check the query and subjects file names
	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_MERGE || operationMode==OPERATION_MODE_METRICS)
	{
		// check the query file name
		if(inputQueryFileInit==NULL || inputQueryFileInit[0]=='\0')
		{
			printf("line=%d, In %s(), please specify the query in configuration file.\n", __LINE__, __func__);
			return FAILED;
		}else if(stat(inputQueryFileInit, &st)==-1)
		{
			printf("line=%d, In %s(), query file [ %s ] does not exist, error! Please specify the correct query file in configuration file.\n", __LINE__, __func__, inputQueryFileInit);
			return FAILED;
		}

		// check the subjects file names
		if(subjectsFile==NULL || subjectsFile[0]=='\0')
		{
			printf("line=%d, In %s(), please specify the subjects in configuration file.\n", __LINE__, __func__);
			return FAILED;
		}else if(stat(subjectsFile, &st)==-1)
		{
			printf("line=%d, In %s(), please specify the correct subjects in configuration file.\n", __LINE__, __func__);
			return FAILED;
		}else
		{
			subjectFileTmp = fopen(subjectsFile, "r");
			if(subjectFileTmp==NULL)
			{
				printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, subjectsFile);
				return FAILED;
			}

			fileNumTmp = 0;
			while(1)
			{
				if(readLine(&fileStatus, fileNameTmp, &len, LINE_CHAR_MAX, subjectFileTmp)==FAILED)
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

				if(stat(fileNameTmp, &st)==-1)
				{
					printf("line=%d, In %s(), please specify the correct subject [ %s ] in configuration file.\n", __LINE__, __func__, fileNameTmp);
					return FAILED;
				}

				fileNumTmp ++;
			}

			if(fileNumTmp==0)
			{
				printf("line=%d, In %s(), please specify the subjects in configuration file.\n", __LINE__, __func__);
				return FAILED;
			}

			fclose(subjectFileTmp);
		}
	}

	// check the read file names
	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_MISASS)
	{
		if(readFileList==NULL)
		{
			printf("line=%d, In %s(), please specify the read files in configuration file.\n", __LINE__, __func__);
			return FAILED;
		}

		fileNumTmp = 0;
		readFileNode = readFileList;
		while(readFileNode)
		{
			i = 0;
			while(i<readFileNode->readFileNum)
			{
				if(stat(readFileNode->readFiles[i], &st)==-1)
				{
					printf("line=%d, In %s(), please specify the correct read file [ %s ] in configuration file.\n", __LINE__, __func__, readFileNode->readFiles[i]);
					return FAILED;
				}

				if(i==1)
				{
					if(strcmp(readFileNode->readFiles[0], readFileNode->readFiles[1])==0)
					{
						printf("line=%d, In %s(), please specify the correct read file [ %s ] in configuration file.\n", __LINE__, __func__, readFileNode->readFiles[1]);
						return FAILED;
					}
				}

				i ++;
			}

			fileNumTmp ++;

			readFileNode = readFileNode->next;
		}

		if(fileNumTmp==0)
		{
			printf("line=%d, In %s(), please specify the read files in configuration file.\n", __LINE__, __func__);
			return FAILED;
		}
		else if(fileNumTmp>=2)
		{
			printf("line=%d, In %s(), one paired-end library is supported in current version, and please specify the correct read files in configuration file.\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Output the configuration information.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputConfigInfo(int32_t operationMode, char *inputQueryFileInit, char *subjectsFile, readFile_t *readFileList)
{
	int32_t i, j;
	struct stat st;
	FILE *subjectFileTmp;
	int32_t len, fileStatus, fileNumTmp;
	char fileNameTmp[LINE_CHAR_MAX+1];
	readFile_t *readFileNode;

	printf("\nConfiguration information:\n");

	// check the query and subjects file names
	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_MERGE || operationMode==OPERATION_MODE_METRICS)
	{
		printf("  Query file: %s\n", inputQueryFileInit);

		printf("  Subject files:\n");

		subjectFileTmp = fopen(subjectsFile, "r");
		if(subjectFileTmp==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, subjectsFile);
			return FAILED;
		}

		fileNumTmp = 0;
		while(1)
		{
			if(readLine(&fileStatus, fileNameTmp, &len, LINE_CHAR_MAX, subjectFileTmp)==FAILED)
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

			printf("    [%d]: %s\n", fileNumTmp, fileNameTmp);

			fileNumTmp ++;
		}

		fclose(subjectFileTmp);
	}

	// check the read file names
	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_MISASS)
	{
		printf("  Read files:\n");

		fileNumTmp = 0;
		readFileNode = readFileList;
		while(readFileNode)
		{
			printf("    [%d]: ", fileNumTmp);

			i = 0;
			while(i<readFileNode->readFileNum)
			{
				if(i==0)
					printf("%s", readFileNode->readFiles[i]);
				else
					printf(", %s", readFileNode->readFiles[i]);

				i ++;
			}
			printf("\n");

			fileNumTmp ++;

			readFileNode = readFileNode->next;
		}
	}

	printf("\n");

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
	int64_t maxQueryLen, queryLen, returnFlag;

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

