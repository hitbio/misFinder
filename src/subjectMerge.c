/*
 * subjectMerge.c
 *
 *  Created on: May 29, 2012
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Merge reference segments into a single file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short mergeRefSegmentsFasta(const char *mergedSegFile, const char *mergeSubjectsFile)
{
	// initialize the memory for reference segments mergence
	if(initMemRefSegmentsFasta(mergeSubjectsFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory for reference segments mergence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// merge reference segments
	if(mergeRefSegsFasta(mergedSegFile, segFileArr, segFileNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot merge segments in fasta, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// free memory for reference segments merging
	freeMemRefSegmentsFasta();

	return SUCCESSFUL;
}

/**
 * Initialize the memory for reference segments mergence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short initMemRefSegmentsFasta(const char *mergeSubjectsFile)
{
	FILE *fpSubjects;
	char line[LINE_CHAR_MAX+1];
	int fileStatus, len;
	int64_t tmpFileNum;

	// get the total segment number
	if(getSegmentFileNum(&segFileNum, mergeSubjectsFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the segment file number, error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(segFileNum<=0)
	{
		printf("There are no segment files in [ %s ], please configure the segment files first.\n", mergeSubjectsFile);
		return FAILED;
	}

	// allocate the segment file memory
	segFileArr = (char **) malloc (segFileNum * sizeof(char*));
	if(segFileArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	for(tmpFileNum=0; tmpFileNum<segFileNum; tmpFileNum++)
	{
		segFileArr[tmpFileNum] = (char *) calloc (LINE_CHAR_MAX+1, sizeof(char));
		if(segFileArr[tmpFileNum]==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	fpSubjects = fopen(mergeSubjectsFile, "r");
	if(fpSubjects==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, mergeSubjectsFile);
		return FAILED;
	}

	tmpFileNum = 0;
	while(1)
	{
		if(readLine(&fileStatus, line, &len, LINE_CHAR_MAX, fpSubjects)==FAILED)
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

		strcpy(segFileArr[tmpFileNum], line);
		tmpFileNum ++;
	}

	fclose(fpSubjects);
	fpSubjects = NULL;

	// ####################### Debug information #######################
	if(tmpFileNum!=segFileNum)
	{
		printf("line=%d, In %s(), tmpFileNum=%ld !=segFileNum=%ld, error!\n", __LINE__, __func__, tmpFileNum, segFileNum);
		freeMemRefSegmentsFasta();
		return FAILED;
	}
	// ####################### Debug information #######################

	return SUCCESSFUL;
}

/**
 * Release the memory of global parameters.
 */
void freeMemRefSegmentsFasta()
{
	int64_t i;

	for(i=0; i<segFileNum; i++)
	{
		free(segFileArr[i]);
		segFileArr[i] = NULL;
	}
	free(segFileArr);
	segFileArr = NULL;

	segFileNum = 0;
}

/**
 * Get the subject segments file number in mergeSubjectsFile.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getSegmentFileNum(int64_t *segmentFileNum, const char *mergeSubjectsFileName)
{
	FILE *fpSubjects;
	char line[LINE_CHAR_MAX+1];
	int fileStatus, len;

	fpSubjects = fopen(mergeSubjectsFileName, "r");
	if(fpSubjects==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, mergeSubjectsFileName);
		return FAILED;
	}

	*segmentFileNum = 0;
	while(1)
	{
		if(readLine(&fileStatus, line, &len, LINE_CHAR_MAX, fpSubjects)==FAILED)
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

		(*segmentFileNum) ++;
	}

	fclose(fpSubjects);
	fpSubjects = NULL;

	return SUCCESSFUL;
}

/**
 * Merge the segments in fasta.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mergeRefSegsFasta(const char *mergedSegFile, char **segFileArr, int64_t segFileNum)
{
	FILE *fpMergedSeg;
	int64_t i, maxFileSizeByte, bufLen;
	char *segBuf;

	fpMergedSeg = fopen(mergedSegFile, "w");
	if(fpMergedSeg==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, mergedSegFile);
		return FAILED;
	}

	if(getMaxFileSizeByte(&maxFileSizeByte, segFileArr, segFileNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the maximal file size, error!\n", __LINE__, __func__);
		return FAILED;
	}

	segBuf = (char *) malloc (maxFileSizeByte * sizeof(char));
	if(segBuf==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<segFileNum; i++)
	{
		// fill the segment buffer
		if(fillSegBuf(segBuf, &bufLen, segFileArr[i])==FAILED)
		{
			printf("line=%d, In %s(), cannot fill segment buffer, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(bufLen>0)
		{
			// trim the tail non-base characters
			if(filterTailNonBases(segBuf, &bufLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot filter the tail non-base characters of segment buffer, error!\n", __LINE__, __func__);
				free(segBuf); segBuf = NULL;
				return FAILED;
			}

			fprintf(fpMergedSeg, "%s\n", segBuf);
		}
	}

	free(segBuf);
	segBuf = NULL;
	fclose(fpMergedSeg);
	fpMergedSeg = NULL;

	return SUCCESSFUL;
}

/**
 * Get the maximal file size.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxFileSizeByte(int64_t *maxFileSizeByte, char **segFileArray, int64_t fileNums)
{
	int64_t i;
	struct stat st;

	*maxFileSizeByte = 0;
	for(i=0; i<fileNums; i++)
	{
		if(stat(segFileArray[i], &st)==-1)
		{
			printf("The file [ %s ] does not exist, please conform whether the file is correctly spelled!\n", segFileArray[i]);
			return FAILED;
		}

		if(st.st_size>(*maxFileSizeByte))
			*maxFileSizeByte = st.st_size;
	}

	if((*maxFileSizeByte)<=0)
	{
		printf("All the segments contain no sequence data, error!\n");
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Filter the tail non-base characters in the segment buffer.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short filterTailNonBases(char *segBuf, int64_t *bufLen)
{
	int64_t i, tmpBufLen, validStartPos, validEndPos;

	tmpBufLen = *bufLen;

	// get the valid start character position
	i = 0;
	while(i<tmpBufLen)
	{
		if(segBuf[i]=='\n')
		{
			validStartPos = i + 1;
			break;
		}
		i++;
	}

	validEndPos = -1;
	for(i=tmpBufLen-1; i>=validStartPos; i--)
	{
		switch(segBuf[i])
		{
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case 'a':
			case 'c':
			case 'g':
			case 't':
			case 'N':
			case 'n':
			case '.':
				validEndPos = i;
				break;
		}

		if(validEndPos>=0)
			break;
	}

	if(validEndPos<validStartPos)
	{
		printf("line=%d, In %s(), validEndPos=%ld < validStartPos=%ld, error!\n", __LINE__, __func__, validEndPos, validStartPos);
		return FAILED;
	}

	segBuf[validEndPos+1] = '\0';
	*bufLen = validEndPos + 1;

	return SUCCESSFUL;
}

/**
 * fill the segment buffer.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillSegBuf(char *segBuf, int64_t *bufLen, char *segmentFile)
{
	FILE *fpSegment;
	char ch;

	fpSegment = fopen(segmentFile, "r");
	if(fpSegment==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, segmentFile);
		return FAILED;
	}

	*bufLen = 0;
	while((ch=fgetc(fpSegment))!=EOF) segBuf[(*bufLen)++] = ch;

	fclose(fpSegment);
	fpSegment = NULL;

	return SUCCESSFUL;
}

/**
 * Fill subject sequences.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillSubjectSeqs(queryMatchInfo_t *queryMatchInfoSet, char *mergedSegFile)
{
	int64_t i, maxSubjectLen, subjectID, subjectLen, returnFlag;
	subject_t *subjectArray;
	char *subjectSeq, subjectHeadTitle[1000];
	FILE *fpSubject;

	fpSubject = fopen(mergedSegFile, "r");
	if(fpSubject==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, mergedSegFile);
		return FAILED;
	}

	subjectArray = queryMatchInfoSet->subjectArray;

	maxSubjectLen = 0;
	for(i=0; i<queryMatchInfoSet->itemNumSubjectArray; i++)
		if(maxSubjectLen<subjectArray[i].subjectLen)
			maxSubjectLen = subjectArray[i].subjectLen;

	subjectSeq = (char *) malloc ((maxSubjectLen+1)*sizeof(char));
	if(subjectSeq==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	subjectID = 1;
	while((returnFlag=getSingleFastaItemFromFile(fpSubject, subjectHeadTitle, subjectSeq, &subjectLen))==SUCCESSFUL)
	{
		if(strcmp(subjectHeadTitle, subjectArray[subjectID-1].subjectTitle)==0)
		{
			subjectArray[subjectID-1].subjectSeq = (char*) calloc (subjectArray[subjectID-1].subjectLen+1, sizeof(char));
			if(subjectArray[subjectID-1].subjectSeq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			strcpy(subjectArray[subjectID-1].subjectSeq, subjectSeq);
		}else
		{
			printf("line=%d, In %s(), cannot get the subject base sequences, error!\n", __LINE__, __func__);
			return FAILED;
		}

		subjectID ++;
	}


	free(subjectSeq);
	subjectSeq = NULL;

	fclose(fpSubject);

	return SUCCESSFUL;
}
