/*
 * accuracyStastics.c
 *
 *  Created on: May 19, 2012
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Get the metrics for queries, including the length metrics and the accuracy metrics.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short getQueryMetrics(char *queryStatisticFile, char *sortedQueryFile, char *queryMatchInfoFile, char *blastnResultFile, char *mergedSegmentsFile)
{
	// initialize the metrics memory
	if(initMemQueryMetrics(queryMatchInfoFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory for metrics of queries, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the length statistics for queries
	if(queryLenStatistics(queryMetrics, sortedQueryFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the length statistics of queries, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the accuracy metrics for queries
	if(queryAccuracyStatistics(queryMetrics, queryMatchInfoSet->queryArray, queryMatchInfoSet->itemNumQueryArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the accuracy statistics of queries, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// compute the reference covered ratio
	if(computeReferenceCovRatio(queryMetrics, queryMatchInfoSet->subjectArray, queryMatchInfoSet->itemNumSubjectArray, queryMatchInfoSet->queryArray, queryMatchInfoSet->itemNumQueryArray, blastnResultFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the genome covered ratio from queries, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the GC content ratio
	if(computeGCRatio(queryMetrics, mergedSegmentsFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the statistics of queries, error!\n", __LINE__, __func__);
		return FAILED;
	}

/*
	// ################################ Debug information ############################
	if(checkRefCoveredRatio(queryMetrics, queryMatchInfoSet->subjectArray, queryMatchInfoSet->itemNumSubjectArray, queryMatchInfoSet->queryArray, queryMatchInfoSet->itemNumQueryArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot check the genome covered ratio, error!\n", __LINE__, __func__);
		return FAILED;
	}
	// ################################ Debug information ############################
*/
	// save the query metrics to file
	if(saveQueryStatisticsToFile(queryStatisticFile, queryMetrics)==FAILED)
	{
		printf("line=%d, In %s(), cannot save the statistics of queries to file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// free the memory of metrics
	freeMemQueryMetrics();

	return SUCCESSFUL;
}

/**
 * Initialize the memory for query metrics.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short initMemQueryMetrics(char *queryMatchInfoFile)
{
	// load the query match information from the binary file
	if(loadQueryMatchInfoFromFile(&queryMatchInfoSet, queryMatchInfoFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the query match information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate memory for the query metrics
	if(allocateQueryMetrics(&queryMetrics, queryMatchInfoSet->subjectArray, queryMatchInfoSet->itemNumSubjectArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Allocate the memory for query metrics.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short allocateQueryMetrics(metrics_t **queryMetrics, subject_t *subjectArray, int64_t itemNumSubjectArray)
{
	int i;

	*queryMetrics = (metrics_t *) calloc (1, sizeof(metrics_t));
	if((*queryMetrics)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	(*queryMetrics)->baseBitArrays = (baseBitArray_t*) calloc(itemNumSubjectArray, sizeof(baseBitArray_t));
	if((*queryMetrics)->baseBitArrays==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	(*queryMetrics)->subjectNum = itemNumSubjectArray;

	for(i=0; i<(*queryMetrics)->subjectNum; i++)
	{
		(*queryMetrics)->baseBitArrays[i].subjectBitArray = (uint8_t *) calloc(subjectArray[i].subjectLen, sizeof(uint8_t));
		if((*queryMetrics)->baseBitArrays[i].subjectBitArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		(*queryMetrics)->baseBitArrays[i].totalRefBaseNum = subjectArray[i].subjectLen;
		(*queryMetrics)->baseBitArrays[i].coveredBaseNum = 0;
		(*queryMetrics)->baseBitArrays[i].coveredRatio = 0;
	}

	return SUCCESSFUL;
}

/**
 * Release query metrics.
 */
void releaseQueryMetrics(metrics_t **queryMetrics)
{
	int i;

	for(i=0; i<(*queryMetrics)->subjectNum; i++)
	{
		free((*queryMetrics)->baseBitArrays[i].subjectBitArray);
		(*queryMetrics)->baseBitArrays[i].subjectBitArray = NULL;
		(*queryMetrics)->baseBitArrays[i].totalRefBaseNum = 0;
		(*queryMetrics)->baseBitArrays[i].coveredBaseNum = 0;
		(*queryMetrics)->baseBitArrays[i].coveredRatio = 0;
	}

	free((*queryMetrics)->baseBitArrays);
	(*queryMetrics)->baseBitArrays = NULL;
	(*queryMetrics)->subjectNum = 0;

	free(*queryMetrics);
	*queryMetrics = NULL;
}

/**
 * Free the memory of query metrics.
 */
void freeMemQueryMetrics()
{
	// free the memory of query match information
	releaseQueryMatchInfo(&queryMatchInfoSet);

	// free the memory of query metrics
	releaseQueryMetrics(&queryMetrics);
}

/**
 * Save the query metrics to file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short saveQueryStatisticsToFile(char *queryStatisticFile, metrics_t *queryMetrics)
{
	fpStatistics = fopen(queryStatisticFile, "w");
	if(fpStatistics==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, queryStatisticFile);
		return FAILED;
	}

/*
	// print the result on screen
	printf("\nStatistics of queries:\n");
	printf("Queries number   : %d\n", queryMetrics->lengthMetrics.totalNum);
	printf("Total length     : %lu\n", queryMetrics->lengthMetrics.totalLen);
	printf("Total ref. len   : %lu\n", queryMetrics->lengthMetrics.totalRefBaseNum);
	printf("Assem Len Ratio  : %.2f %%\n", queryMetrics->lengthMetrics.lengthRatio * 100);
	printf("Maximal length   : %d\n", queryMetrics->lengthMetrics.maxSize);
	printf("N50 size         : %d\n", queryMetrics->lengthMetrics.N50);
	printf("Mean size        : %.2f\n", queryMetrics->lengthMetrics.meanSize);
	printf("Median size      : %d\n", queryMetrics->lengthMetrics.medianSize);
	printf("Ref. cov ratio   : %.2f %%\n", queryMetrics->lengthMetrics.coveredRatio * 100);
	printf("Organism GC ratio: %.2f %%\n", queryMetrics->lengthMetrics.GCRatio * 100);

	printf("\nPerfect queries  : lenRatio=%.2f %%, totalLen=%ld, Num=%d, meanSize=%.2f\n", queryMetrics->accuracyMetrics[0].lengthRatio * 100, queryMetrics->accuracyMetrics[0].totalLen, queryMetrics->accuracyMetrics[0].totalNum, queryMetrics->accuracyMetrics[0].meanSize);
	printf("Matched queries  : lenRatio=%.2f %%, totalLen=%ld, Num=%d, meanSize=%.2f\n", queryMetrics->accuracyMetrics[1].lengthRatio * 100, queryMetrics->accuracyMetrics[1].totalLen, queryMetrics->accuracyMetrics[1].totalNum, queryMetrics->accuracyMetrics[1].meanSize);
	printf("Disjunct queries : lenRatio=%.2f %%, totalLen=%ld, Num=%d, meanSize=%.2f\n", queryMetrics->accuracyMetrics[2].lengthRatio * 100, queryMetrics->accuracyMetrics[2].totalLen, queryMetrics->accuracyMetrics[2].totalNum, queryMetrics->accuracyMetrics[2].meanSize);
	printf("Unmatched queries: lenRatio=%.2f %%, totalLen=%ld, Num=%d, meanSize=%.2f\n\n", queryMetrics->accuracyMetrics[3].lengthRatio * 100, queryMetrics->accuracyMetrics[3].totalLen, queryMetrics->accuracyMetrics[3].totalNum, queryMetrics->accuracyMetrics[3].meanSize);
*/

	// output the result to file
	fprintf(fpStatistics, "Statistics of queries:\n");
	fprintf(fpStatistics, "Queries number   : %d\n", queryMetrics->lengthMetrics.totalNum);
	fprintf(fpStatistics, "Total length     : %lu\n", queryMetrics->lengthMetrics.totalLen);
	fprintf(fpStatistics, "Total ref. len   : %lu\n", queryMetrics->lengthMetrics.totalRefBaseNum);
	fprintf(fpStatistics, "Assem Len Ratio  : %.2f %%\n", queryMetrics->lengthMetrics.lengthRatio * 100);
	fprintf(fpStatistics, "Maximal length   : %d\n", queryMetrics->lengthMetrics.maxSize);
	fprintf(fpStatistics, "N50 size         : %d\n", queryMetrics->lengthMetrics.N50);
	fprintf(fpStatistics, "Mean size        : %.2f\n", queryMetrics->lengthMetrics.meanSize);
	fprintf(fpStatistics, "Median size      : %d\n", queryMetrics->lengthMetrics.medianSize);
	fprintf(fpStatistics, "Ref. cov ratio   : %.2f %%\n", queryMetrics->lengthMetrics.coveredRatio * 100);
	fprintf(fpStatistics, "Organism GC ratio: %.2f %%\n", queryMetrics->lengthMetrics.GCRatio * 100);

	fprintf(fpStatistics, "\nPerfect queries  : lenRatio=%.2f %%, totalLen=%ld, Num=%d, meanSize=%.2f\n", queryMetrics->accuracyMetrics[0].lengthRatio * 100, queryMetrics->accuracyMetrics[0].totalLen, queryMetrics->accuracyMetrics[0].totalNum, queryMetrics->accuracyMetrics[0].meanSize);
	fprintf(fpStatistics, "Matched queries  : lenRatio=%.2f %%, totalLen=%ld, Num=%d, meanSize=%.2f\n", queryMetrics->accuracyMetrics[1].lengthRatio * 100, queryMetrics->accuracyMetrics[1].totalLen, queryMetrics->accuracyMetrics[1].totalNum, queryMetrics->accuracyMetrics[1].meanSize);
	fprintf(fpStatistics, "Disjunct queries : lenRatio=%.2f %%, totalLen=%ld, Num=%d, meanSize=%.2f\n", queryMetrics->accuracyMetrics[2].lengthRatio * 100, queryMetrics->accuracyMetrics[2].totalLen, queryMetrics->accuracyMetrics[2].totalNum, queryMetrics->accuracyMetrics[2].meanSize);
	fprintf(fpStatistics, "Unmatched queries: lenRatio=%.2f %%, totalLen=%ld, Num=%d, meanSize=%.2f\n", queryMetrics->accuracyMetrics[3].lengthRatio * 100, queryMetrics->accuracyMetrics[3].totalLen, queryMetrics->accuracyMetrics[3].totalNum, queryMetrics->accuracyMetrics[3].meanSize);

	fclose(fpStatistics);
	fpStatistics = NULL;

	return SUCCESSFUL;
}


/**
 * Get the accuracy statistics for queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short queryAccuracyStatistics(metrics_t *queryMetrics, query_t *queryArray, int64_t itemNumQueryArray)
{
	int64_t i;

	for(i=0; i<4; i++) queryMetrics->accuracyMetrics[i].totalNum = queryMetrics->accuracyMetrics[i].totalLen = queryMetrics->accuracyMetrics[i].meanSize = queryMetrics->accuracyMetrics[i].lengthRatio = 0;

	for(i=0; i<itemNumQueryArray; i++)
	{
		if(queryArray[i].queryLen>=minQueryLenThres)
		{
			switch(queryArray[i].globalMatchKind)
			{
				case PERFECT_MATCH_KIND:
					queryMetrics->accuracyMetrics[0].totalNum ++;
					queryMetrics->accuracyMetrics[0].totalLen += queryArray[i].queryLen;
					break;
				case MATCHED_KIND:
					queryMetrics->accuracyMetrics[1].totalNum ++;
					queryMetrics->accuracyMetrics[1].totalLen += queryArray[i].queryLen;
					break;
				case DISJUNCT_MATCH_KIND:
					queryMetrics->accuracyMetrics[2].totalNum ++;
					queryMetrics->accuracyMetrics[2].totalLen += queryArray[i].queryLen;
					break;
				case UNMATCHED_KIND:
					queryMetrics->accuracyMetrics[3].totalNum ++;
					queryMetrics->accuracyMetrics[3].totalLen += queryArray[i].queryLen;
					break;
				default:
					printf("In %s(), unknown match kind %d, error!\n", __func__, queryArray[i].globalMatchKind);
					return FAILED;
			}
		}
	}

	for(i=0; i<4; i++)
	{
		if(queryMetrics->accuracyMetrics[i].totalNum>0)
		{
			queryMetrics->accuracyMetrics[i].meanSize = (double)queryMetrics->accuracyMetrics[i].totalLen / queryMetrics->accuracyMetrics[i].totalNum;
			if(queryMetrics->lengthMetrics.totalRefBaseNum>0)
				queryMetrics->accuracyMetrics[i].lengthRatio = (double)queryMetrics->accuracyMetrics[i].totalLen / queryMetrics->lengthMetrics.totalRefBaseNum;
			else
			{
				printf("line=%d, In %s(), the total reference base count is %ld, error!\n", __LINE__, __func__, queryMetrics->lengthMetrics.totalRefBaseNum);
				return FAILED;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the GC ratio of reference.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short computeGCRatio(metrics_t *queryMetrics, char *mergedSegmentsFile)
{
	int64_t GC_BaseNum, AT_BaseNum;
	FILE *fpSegs;
	char ch;

	fpSegs = fopen(mergedSegmentsFile, "r");
	if(fpSegs==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, mergedSegmentsFile);
		return FAILED;
	}

	GC_BaseNum = 0;
	AT_BaseNum = 0;

	// skip the head descriptions
	while(1)
	{
		ch = fgetc(fpSegs);
		if(ch=='>')
		{
			fseek(fpSegs, -1, SEEK_CUR);
			break;
		}
	}

	while(1)
	{
		ch = fgetc(fpSegs);
		if(ch==EOF)
		{ // end of file
			break;
		}else if(ch=='>')
		{ // header line
			while(ch!='\n') ch = fgetc(fpSegs);
		}else
		{ // base part
			switch(ch)
			{
				case 'G':
				case 'C':
				case 'g':
				case 'c':
					GC_BaseNum ++;
					break;
				case 'A':
				case 'T':
				case 'a':
				case 't':
					AT_BaseNum ++;
					break;
			}
		}
	}

	if(GC_BaseNum+AT_BaseNum>0)
		queryMetrics->lengthMetrics.GCRatio = (double)GC_BaseNum / (GC_BaseNum + AT_BaseNum);
	else
	{
		printf("Please input the correct genome reference when calculate the GC content ratio, error!\n");
		return FAILED;
	}

	fclose(fpSegs);
	fpSegs = NULL;

	return SUCCESSFUL;
}
