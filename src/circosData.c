/*
 * circosData.c
 *
 *  Created on: Aug 8, 2014
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Generate the data for Circos.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateCircosData(queryMatchInfo_t *queryMatchInfoSet, readSet_t *readSet)
{
	int32_t i, queryLen;
	query_t *queryItem;
	baseCov_t *baseCovArray;
	char newQueryLabel[256];
	char queryFileCircos[256], covFileCircos[256], disagreeFileCircos[256];
	char zeroCovFileCircos[256], discorFileCircos[256], orphanedCovFileCircos[256];
	char SPRatioFileCircos[256], SMinusRatioFileCircos[256], SPlusRatioFileCircos[256];
	char resultFileCircos[256];
	FILE *fpQueryCovData, *fpQueryDisNum, *fpZeroCov, *fpDiscorCov, *fpOrphanedCov;
	FILE *fpSPRatio, *fpSMinusRatio, *fpSPlusRatio, *fpResult;

	strcpy(queryFileCircos, outputPathStr);
	strcat(queryFileCircos, "circos_queries.txt");

	strcpy(covFileCircos, outputPathStr);
	strcat(covFileCircos, "circos_coverage.txt");

	strcpy(zeroCovFileCircos, outputPathStr);
	strcat(zeroCovFileCircos, "circos_zeroCov.txt");

	strcpy(discorFileCircos, outputPathStr);
	strcat(discorFileCircos, "circos_discorCov.txt");

	strcpy(orphanedCovFileCircos, outputPathStr);
	strcat(orphanedCovFileCircos, "circos_orphanedCov.txt");

	strcpy(SPRatioFileCircos, outputPathStr);
	strcat(SPRatioFileCircos, "circos_SPRatio.txt");

	strcpy(SMinusRatioFileCircos, outputPathStr);
	strcat(SMinusRatioFileCircos, "circos_SMinusRatio.txt");

	strcpy(SPlusRatioFileCircos, outputPathStr);
	strcat(SPlusRatioFileCircos, "circos_SPlusRatio.txt");

	strcpy(disagreeFileCircos, outputPathStr);
	strcat(disagreeFileCircos, "circos_disagreeNum.txt");

	strcpy(resultFileCircos, outputPathStr);
	strcat(resultFileCircos, "circos_result.txt");


	fpQueryCovData = fopen (covFileCircos, "w");
	if(fpQueryCovData==NULL)
	{
		printf("line=%d, In %s(), cannot create circos query file [%s], error!\n", __LINE__, __func__, covFileCircos);
		return FAILED;
	}

	fpZeroCov = fopen (zeroCovFileCircos, "w");
	if(fpZeroCov==NULL)
	{
		printf("line=%d, In %s(), cannot create circos query file [%s], error!\n", __LINE__, __func__, zeroCovFileCircos);
		return FAILED;
	}

	fpDiscorCov = fopen (discorFileCircos, "w");
	if(fpDiscorCov==NULL)
	{
		printf("line=%d, In %s(), cannot create circos query file [%s], error!\n", __LINE__, __func__, discorFileCircos);
		return FAILED;
	}

	fpOrphanedCov = fopen (orphanedCovFileCircos, "w");
	if(fpOrphanedCov==NULL)
	{
		printf("line=%d, In %s(), cannot create circos query file [%s], error!\n", __LINE__, __func__, orphanedCovFileCircos);
		return FAILED;
	}

	fpSPRatio = fopen (SPRatioFileCircos, "w");
	if(fpSPRatio==NULL)
	{
		printf("line=%d, In %s(), cannot create circos query file [%s], error!\n", __LINE__, __func__, SPRatioFileCircos);
		return FAILED;
	}

	fpSMinusRatio = fopen (SMinusRatioFileCircos, "w");
	if(fpSMinusRatio==NULL)
	{
		printf("line=%d, In %s(), cannot create circos query file [%s], error!\n", __LINE__, __func__, SMinusRatioFileCircos);
		return FAILED;
	}

	fpSPlusRatio = fopen (SPlusRatioFileCircos, "w");
	if(fpSPlusRatio==NULL)
	{
		printf("line=%d, In %s(), cannot create circos query file [%s], error!\n", __LINE__, __func__, SPlusRatioFileCircos);
		return FAILED;
	}

	fpQueryDisNum = fopen (disagreeFileCircos, "w");
	if(fpQueryDisNum==NULL)
	{
		printf("line=%d, In %s(), cannot create circos query file [%s], error!\n", __LINE__, __func__, disagreeFileCircos);
		return FAILED;
	}

	fpResult = fopen (resultFileCircos, "w");
	if(fpResult==NULL)
	{
		printf("line=%d, In %s(), cannot create circos query file [%s], error!\n", __LINE__, __func__, resultFileCircos);
		return FAILED;
	}

	// generate circos query file
	if(generateCircosQueryData(queryFileCircos, queryMatchInfoSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate circos base coverage data for single query, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		queryItem = queryMatchInfoSet->queryArray + i;
		queryLen = queryItem->queryLen;

//		if(queryItem->queryID==16 || queryItem->queryID==44 || queryItem->queryID==59
//			|| queryItem->queryID==73 || queryItem->queryID==231 || queryItem->queryID==262
//			|| queryItem->queryID==27 || queryItem->queryID==147 || queryItem->queryID==62 || queryItem->queryID==156 || queryItem->queryID==309)
		{

			// ########################### Debug information ##############################
			//if(queryItem->queryID==4 || strcmp(queryItem->queryTitle, "ctg7180000002351")==0)
			{
				printf("^^^^^^ queryID=%d, queryTitle=%s, queryLen=%d, subjectNum=%d\n", queryItem->queryID, queryItem->queryTitle, queryItem->queryLen, queryItem->querySubjectNum);
			}
			// ########################### Debug information ##############################

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

			strcpy(newQueryLabel, "Q");
			sprintf(newQueryLabel+strlen(newQueryLabel), "%d", i+1);

			// generate base coverage file
			if(fillBaseCovCircosDataSingleQuery(fpQueryCovData, fpZeroCov, newQueryLabel, baseCovArray, queryLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate circos base coverage data for single query, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// generate disagreeNum data
			if(generateCircosDisNumSingleQuery(fpQueryDisNum, newQueryLabel, baseCovArray, queryLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate circos disagreeNum data for single query, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// generate coverage of discordant reads data
			if(generateCircosHeatmapAndRatioDataSingleQuery(fpDiscorCov, fpOrphanedCov, fpSPRatio, fpSMinusRatio, fpSPlusRatio, newQueryLabel, queryItem, readSet)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate circos heat map for discordant pairs for single query, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// generate result data
			if(generateResultCircosDataSingleQuery(fpResult, newQueryLabel, queryItem)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate circos result data for single query, error!\n", __LINE__, __func__);
				return FAILED;
			}

			free(baseCovArray);
		}
	}

	fclose(fpQueryCovData);
	fclose(fpZeroCov);
	fclose(fpDiscorCov);
	fclose(fpOrphanedCov);
	fclose(fpSPRatio);
	fclose(fpSMinusRatio);
	fclose(fpSPlusRatio);
	fclose(fpQueryDisNum);
	fclose(fpResult);

	return SUCCESSFUL;
}

/**
 * Generate circos query data.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateCircosQueryData(char *queryFileCircos, queryMatchInfo_t *queryMatchInfoSet)
{
	int32_t i, j, queryLen, unknownBaseNum, startRow, endRow, gapID;
	query_t *queryItem;
	char newQueryLabel[256], gapLabel[256], *querySeq;
	FILE *fpQueryData;

	fpQueryData = fopen (queryFileCircos, "w");
	if(fpQueryData==NULL)
	{
		printf("line=%d, In %s(), cannot create circos query file [%s], error!\n", __LINE__, __func__, queryFileCircos);
		return FAILED;
	}

	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		queryItem = queryMatchInfoSet->queryArray + i;
		queryLen = queryItem->queryLen;

//		if(queryItem->queryID==16 || queryItem->queryID==44 || queryItem->queryID==59
//			|| queryItem->queryID==73 || queryItem->queryID==231 || queryItem->queryID==262
//			|| queryItem->queryID==27 || queryItem->queryID==147 || queryItem->queryID==62 || queryItem->queryID==156 || queryItem->queryID==309)
		{

			strcpy(newQueryLabel, "Q");
			sprintf(newQueryLabel+strlen(newQueryLabel), "%d", i+1);

			fprintf(fpQueryData, "chr - %s\t%d\t%d\t%d\t%s\n", newQueryLabel, i+1, 0, queryLen-1, "green");
		}
	}

	// generate circos query band data
	gapID = 1;
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		queryItem = queryMatchInfoSet->queryArray + i;
		queryLen = queryItem->queryLen;

//		if(queryItem->queryID==16 || queryItem->queryID==44 || queryItem->queryID==59
//			|| queryItem->queryID==73 || queryItem->queryID==231 || queryItem->queryID==262
//			|| queryItem->queryID==27 || queryItem->queryID==147 || queryItem->queryID==62 || queryItem->queryID==156 || queryItem->queryID==309)
		{

			strcpy(newQueryLabel, "Q");
			sprintf(newQueryLabel+strlen(newQueryLabel), "%d", i+1);

			querySeq = queryItem->querySeq;
			unknownBaseNum = 0;
			for(j=0; j<queryLen; j++)
			{
				if(querySeq[j]=='N' || querySeq[j]=='n' || querySeq[j]=='.')
				{
					unknownBaseNum ++;
				}else if(unknownBaseNum>0)
				{
					startRow = j - unknownBaseNum;
					endRow = j - 1;

					strcpy(gapLabel, "gap");
					sprintf(gapLabel+strlen(gapLabel), "%d", gapID);

					fprintf(fpQueryData, "band\t%s\t%s\t%s\t%d\t%d\t%s\n", newQueryLabel, gapLabel, gapLabel, startRow, endRow, "gpos75");

					unknownBaseNum = 0;
					gapID ++;
				}
			}
		}
	}

	fclose(fpQueryData);

	return SUCCESSFUL;
}

/**
 * Fill the base coverage data for Circos from single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillBaseCovCircosDataSingleQuery(FILE *fpCovData, FILE *fpZeroCov, char *queryLabel, baseCov_t *baseCovArray, int32_t arraySize)
{
	int32_t i, subRegSize, startRow, endRow, itemNum, sumCov, covNum, maxValue;
	double averCov;

	subRegSize = 50;
	maxValue = 200;
	sumCov = itemNum = 0;
	for(i=0; i<arraySize; i++)
	{
		covNum = baseCovArray[i].baseNumArray[5];

		if(covNum==0)
		{
			fprintf(fpZeroCov, "%s\t%d\t%d\t%d\n", queryLabel, i, i, 1);
		}

		sumCov += covNum;
		itemNum ++;
		if(itemNum==subRegSize || i==arraySize-1)
		{
			if(itemNum>0)
			{
				averCov = (double)sumCov / itemNum;
				if(averCov>maxValue)
					averCov = maxValue;
				endRow = i;
				startRow = endRow - itemNum + 1;
				fprintf(fpCovData, "%s\t%d\t%d\t%.4f\n", queryLabel, startRow, endRow, averCov);
			}

			itemNum = 0;
			sumCov = 0;
		}
	}

	return SUCCESSFUL;
}

/**
 * Generate Circos disagreeNum for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateCircosDisNumSingleQuery(FILE *fpQueryDisNum, char *queryLabel, baseCov_t *baseCovArray, int32_t arraySize)
{
	int32_t i, j, maxID, maxValue, total;
	double maxRatio;

	for(i=0; i<arraySize; i++)
	{
		if(baseCovArray[i].baseNumArray[5]>0)
		{
			maxID = -1;
			maxValue = 0;
			total = baseCovArray[i].baseNumArray[5];
			for(j=0; j<5; j++)
			{
				if(baseCovArray[i].baseNumArray[j]>maxValue)
				{
					maxValue = baseCovArray[i].baseNumArray[j];
					maxID = j;
				}
			}

			maxRatio = (double)maxValue / total;
			if(maxRatio<DISAGREE_RATIO_THRES)
			{
				fprintf(fpQueryDisNum, "%s\t%d\t%d\t%d\n", queryLabel, i, i, 1);
			}
		}else
		{
			fprintf(fpQueryDisNum, "%s\t%d\t%d\t%d\n", queryLabel, i, i, 1);
		}
	}

	return SUCCESSFUL;
}

/**
 * Generate Circos heat map for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateCircosHeatmapAndRatioDataSingleQuery(FILE *fpDiscorCov, FILE *fpOrphanedCov, FILE *fpSPRatio, FILE *fpSMinusRatio, FILE *fpSPlusRatio, char *queryLabel, query_t *queryItem, readSet_t *readSet)
{
	ratioRegion_t *ratioRegionArray;
	int32_t ratioRegionNum, subRegSize;

	subRegSize = 500;

	// initialize the ratio region array
	if(initRatioRegArrayCircos(&ratioRegionArray, &ratioRegionNum, subRegSize, queryItem)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the ratio array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the ratioRegion array
	if(fillRatioRegionArray(ratioRegionArray, ratioRegionNum, queryItem, readSet, insertSize, standDev)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the ratioRegion array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// compute their values
	if(computeRatios(ratioRegionArray, ratioRegionNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the ratios, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// output the discordant, orphaned data
	if(outputCircosHeatmapSingleQuery(fpDiscorCov, fpOrphanedCov, queryLabel, ratioRegionArray, ratioRegionNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot output the circos heat map data, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// output ratios
	if(outputCircosRatiosSingleQuery(fpSPRatio, fpSMinusRatio, fpSPlusRatio, queryLabel, ratioRegionArray, ratioRegionNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot output the circos ratio data, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// free memory
	if(ratioRegionNum>0)
		free(ratioRegionArray);

	return SUCCESSFUL;
}

/**
 * Initialize ratio region array for Circos data for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initRatioRegArrayCircos(ratioRegion_t **ratioRegionArray, int32_t *ratioRegionNum, int32_t subRegSize, query_t *queryItem)
{
	int32_t i, itemNum, startRegPos, endRegPos, midPos;

	*ratioRegionNum = (queryItem->queryLen - 1) / subRegSize + 1;
	if((*ratioRegionNum)<=0)
	{
		return SUCCESSFUL;
	}

	*ratioRegionArray = (ratioRegion_t*) calloc (*ratioRegionNum, sizeof(ratioRegion_t));
	if((*ratioRegionArray)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	itemNum = 0;
	for(i=1; i<queryItem->queryLen; i+=subRegSize)
	{
		startRegPos = i;
		endRegPos = startRegPos + subRegSize - 1;
		if(endRegPos>queryItem->queryLen)
			endRegPos = queryItem->queryLen;

		midPos = (startRegPos + endRegPos) / 2;
		(*ratioRegionArray)[itemNum].midQPos = midPos;
		(*ratioRegionArray)[itemNum].startQPosLHalf = startRegPos;
		(*ratioRegionArray)[itemNum].endQPosLHalf = midPos;
		(*ratioRegionArray)[itemNum].startQPosRHalf = midPos + 1;
		(*ratioRegionArray)[itemNum].endQPosRHalf = endRegPos;

		(*ratioRegionArray)[itemNum].disagreeNum = 0;
		(*ratioRegionArray)[itemNum].zeroCovNum = 0;
		(*ratioRegionArray)[itemNum].discorNum = 0;

		(*ratioRegionArray)[itemNum].SPRatio = -1;
		(*ratioRegionArray)[itemNum].singleMinusRatio = -1;
		(*ratioRegionArray)[itemNum].singlePlusRatio = -1;
		(*ratioRegionArray)[itemNum].discorRatio = -1;
		itemNum ++;
	}

	if(itemNum>(*ratioRegionNum))
	{
		printf("line=%d, In %s(), itemNum=%d, ratioRegionNum=%d, error!\n", __LINE__, __func__, itemNum, *ratioRegionNum);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Output Circos heat map data for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputCircosHeatmapSingleQuery(FILE *fpDiscorCov, FILE *fpOrphanedCov, char *queryLabel, ratioRegion_t *ratioRegionArray, int32_t ratioRegionNum)
{
	int32_t i, startRow, endRow, discorNum, singleNum;

	for(i=0; i<ratioRegionNum; i++)
	{
		startRow = ratioRegionArray[i].startQPosLHalf - 1;
		endRow = ratioRegionArray[i].endQPosRHalf - 1;
		discorNum = ratioRegionArray[i].discorNum;
		singleNum = ratioRegionArray[i].singleNum;

		if(discorNum>0)
			fprintf(fpDiscorCov, "%s\t%d\t%d\t%d\n", queryLabel, startRow, endRow, discorNum);
		if(singleNum>0)
			fprintf(fpOrphanedCov, "%s\t%d\t%d\t%d\n", queryLabel, startRow, endRow, singleNum);
	}

	return SUCCESSFUL;
}

/**
 * Output Circos ratio data for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputCircosRatiosSingleQuery(FILE *fpSPRatio, FILE *fpSMinusRatio, FILE *fpSPlusRatio, char *queryLabel, ratioRegion_t *ratioRegionArray, int32_t ratioRegionNum)
{
	int32_t i, startRow, endRow;
	double SPRatio, singleMinusRatio, singlePlusRatio;

	for(i=0; i<ratioRegionNum; i++)
	{
		startRow = ratioRegionArray[i].startQPosLHalf - 1;
		endRow = ratioRegionArray[i].endQPosRHalf - 1;
		SPRatio = ratioRegionArray[i].SPRatio;
		singleMinusRatio = ratioRegionArray[i].singleMinusRatio;
		singlePlusRatio = ratioRegionArray[i].singlePlusRatio;

		if(SPRatio<0)
			SPRatio = 0;
		if(singleMinusRatio<0)
			singleMinusRatio = 0;
		if(singlePlusRatio<0)
			singlePlusRatio = 0;

		fprintf(fpSPRatio, "%s\t%d\t%d\t%.8f\n", queryLabel, startRow, endRow, SPRatio);
		fprintf(fpSMinusRatio, "%s\t%d\t%d\t%.8f\n", queryLabel, startRow, endRow, singleMinusRatio);
		fprintf(fpSPlusRatio, "%s\t%d\t%d\t%.8f\n", queryLabel, startRow, endRow, singlePlusRatio);
	}

	return SUCCESSFUL;
}

/**
 * Generate Circos result data for single query.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateResultCircosDataSingleQuery(FILE *fpResult, char *queryLabel, query_t *queryItem)
{
	misInfo_t *misInfo;
	misassSeq_t *misassSeq;
	char colorStr[256];
	int32_t startRow, endRow;

	misInfo = queryItem->misInfoList;
	while(misInfo)
	{
		if(misInfo->misassFlag==TRUE_MISASS)
			strcpy(colorStr, "color=red");
		else if(misInfo->misassFlag==STRUCTURE_VARIATION)
			strcpy(colorStr, "color=blue");
		else if(misInfo->misassFlag==UNCERTAIN_MISASS)
			strcpy(colorStr, "color=yellow");

		misassSeq = misInfo->misassSeqList;
		while(misassSeq)
		{
			startRow = misassSeq->startQueryPos - 1;
			endRow = misassSeq->endQueryPos - 1;
			fprintf(fpResult, "%s\t%d\t%d\t%d\t%s\n", queryLabel, startRow, endRow, 1, colorStr);

			misassSeq = misassSeq->next;
		}

		misInfo = misInfo->next;
	}

	return SUCCESSFUL;
}
