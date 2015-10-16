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
	char zeroCovFileCircos[256], discorFileCircos[256];
	char multiRatioFileCircos[256], resultFileCircos[256];
	FILE *fpQueryCovData, *fpQueryDisNum, *fpZeroCov;
	FILE *fpMultiRatio, *fpDiscorRatio, *fpResult;

	strcpy(queryFileCircos, outputPathStr);
	strcat(queryFileCircos, "circos_queries.txt");

	strcpy(covFileCircos, outputPathStr);
	strcat(covFileCircos, "circos_coverage.txt");

	strcpy(zeroCovFileCircos, outputPathStr);
	strcat(zeroCovFileCircos, "circos_zeroCov.txt");

	strcpy(discorFileCircos, outputPathStr);
	strcat(discorFileCircos, "circos_discorRatio.txt");

	strcpy(multiRatioFileCircos, outputPathStr);
	strcat(multiRatioFileCircos, "circos_multiRatio.txt");

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
		printf("line=%d, In %s(), cannot create file [%s], error!\n", __LINE__, __func__, zeroCovFileCircos);
		return FAILED;
	}

	fpDiscorRatio = fopen (discorFileCircos, "w");
	if(fpDiscorRatio==NULL)
	{
		printf("line=%d, In %s(), cannot create file [%s], error!\n", __LINE__, __func__, discorFileCircos);
		return FAILED;
	}

	fpMultiRatio = fopen (multiRatioFileCircos, "w");
	if(fpMultiRatio==NULL)
	{
		printf("line=%d, In %s(), cannot create file [%s], error!\n", __LINE__, __func__, multiRatioFileCircos);
		return FAILED;
	}

	fpQueryDisNum = fopen (disagreeFileCircos, "w");
	if(fpQueryDisNum==NULL)
	{
		printf("line=%d, In %s(), cannot create file [%s], error!\n", __LINE__, __func__, disagreeFileCircos);
		return FAILED;
	}

	fpResult = fopen (resultFileCircos, "w");
	if(fpResult==NULL)
	{
		printf("line=%d, In %s(), cannot create file [%s], error!\n", __LINE__, __func__, resultFileCircos);
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

		// data for ecoli
//		if(queryItem->queryID==1 || queryItem->queryID==3 || queryItem->queryID==8 || queryItem->queryID==22
//			|| queryItem->queryID==28 || queryItem->queryID==36 || queryItem->queryID==38
//			|| queryItem->queryID==41 || queryItem->queryID==42 || queryItem->queryID==43 || queryItem->queryID==45 || queryItem->queryID==46
//			|| queryItem->queryID==47 || queryItem->queryID==54 || queryItem->queryID==56 || queryItem->queryID==60
//			|| queryItem->queryID==61 || queryItem->queryID==64 || queryItem->queryID==67 || queryItem->queryID==71)

		// data for spombe
		if(strcmp(queryItem->queryTitle, "scf7180000013910")==0 || strcmp(queryItem->queryTitle, "scf7180000013931")==0
			|| strcmp(queryItem->queryTitle, "scf7180000013975")==0 || strcmp(queryItem->queryTitle, "scf7180000014158")==0
			|| strcmp(queryItem->queryTitle, "scf7180000014214")==0 || strcmp(queryItem->queryTitle, "scf7180000014218")==0
			|| strcmp(queryItem->queryTitle, "scf7180000013964")==0 || strcmp(queryItem->queryTitle, "scf7180000014209")==0
			|| strcmp(queryItem->queryTitle, "scf7180000014217")==0 || strcmp(queryItem->queryTitle, "scf7180000014226")==0)
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
			if(generateCircosHeatmapAndRatioDataSingleQuery(fpDiscorRatio, fpMultiRatio, newQueryLabel, queryItem, readSet)==FAILED)
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
	fclose(fpDiscorRatio);
	fclose(fpMultiRatio);
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

		// data for ecoli
//		if(queryItem->queryID==1 || queryItem->queryID==3 || queryItem->queryID==8 || queryItem->queryID==22
//			|| queryItem->queryID==28 || queryItem->queryID==36 || queryItem->queryID==38
//			|| queryItem->queryID==41 || queryItem->queryID==42 || queryItem->queryID==43 || queryItem->queryID==45 || queryItem->queryID==46
//			|| queryItem->queryID==47 || queryItem->queryID==54 || queryItem->queryID==56 || queryItem->queryID==60
//			|| queryItem->queryID==61 || queryItem->queryID==64 || queryItem->queryID==67 || queryItem->queryID==71)

		// data for spombe
		if(strcmp(queryItem->queryTitle, "scf7180000013910")==0 || strcmp(queryItem->queryTitle, "scf7180000013931")==0
			|| strcmp(queryItem->queryTitle, "scf7180000013975")==0 || strcmp(queryItem->queryTitle, "scf7180000014158")==0
			|| strcmp(queryItem->queryTitle, "scf7180000014214")==0 || strcmp(queryItem->queryTitle, "scf7180000014218")==0
			|| strcmp(queryItem->queryTitle, "scf7180000013964")==0 || strcmp(queryItem->queryTitle, "scf7180000014209")==0
			|| strcmp(queryItem->queryTitle, "scf7180000014217")==0 || strcmp(queryItem->queryTitle, "scf7180000014226")==0)
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

		// data for ecoli
//		if(queryItem->queryID==1 || queryItem->queryID==3 || queryItem->queryID==8 || queryItem->queryID==22
//			|| queryItem->queryID==28 || queryItem->queryID==36 || queryItem->queryID==38
//			|| queryItem->queryID==41 || queryItem->queryID==42 || queryItem->queryID==43 || queryItem->queryID==45 || queryItem->queryID==46
//			|| queryItem->queryID==47 || queryItem->queryID==54 || queryItem->queryID==56 || queryItem->queryID==60
//			|| queryItem->queryID==61 || queryItem->queryID==64 || queryItem->queryID==67 || queryItem->queryID==71)

		// data for spombe
		if(strcmp(queryItem->queryTitle, "scf7180000013910")==0 || strcmp(queryItem->queryTitle, "scf7180000013931")==0
			|| strcmp(queryItem->queryTitle, "scf7180000013975")==0 || strcmp(queryItem->queryTitle, "scf7180000014158")==0
			|| strcmp(queryItem->queryTitle, "scf7180000014214")==0 || strcmp(queryItem->queryTitle, "scf7180000014218")==0
			|| strcmp(queryItem->queryTitle, "scf7180000013964")==0 || strcmp(queryItem->queryTitle, "scf7180000014209")==0
			|| strcmp(queryItem->queryTitle, "scf7180000014217")==0 || strcmp(queryItem->queryTitle, "scf7180000014226")==0)
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
short generateCircosHeatmapAndRatioDataSingleQuery(FILE *fpDiscorRatio, FILE *fpMultiRatio, char *queryLabel, query_t *queryItem, readSet_t *readSet)
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
	if(fillRatioRegionArray(ratioRegionArray, ratioRegionNum, queryItem, readSet, NO)==FAILED)
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

	// output the discordant, multi-align ratio data
	if(outputCircosHeatmapSingleQuery(fpDiscorRatio, fpMultiRatio, queryLabel, ratioRegionArray, ratioRegionNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot output the circos heat map data, error!\n", __LINE__, __func__);
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
short outputCircosHeatmapSingleQuery(FILE *fpDiscorRatio, FILE *fpMultiRatio, char *queryLabel, ratioRegion_t *ratioRegionArray, int32_t ratioRegionNum)
{
	int32_t i, startRow, endRow, discorNum;
	int32_t endHeadRow, startTailRow;
	double multiRatio, discorRatio;

	for(i=0; i<ratioRegionNum; i++)
	{
		startRow = ratioRegionArray[i].startQPosLHalf - 1;
		endRow = ratioRegionArray[i].endQPosRHalf - 1;
		discorNum = ratioRegionArray[i].discorNum;

		if(discorNum>0)
			discorRatio = ratioRegionArray[i].discorRatio;
		else
			discorRatio = 0;

		if(discorRatio>0)
		{
			fprintf(fpDiscorRatio, "%s\t%d\t%d\t%.8f\n", queryLabel, startRow, endRow, discorRatio);
		}else
			fprintf(fpDiscorRatio, "%s\t%d\t%d\t%d\n", queryLabel, startRow, endRow, 0);
	}

	endHeadRow = startTailRow = -1;
	for(i=0; i<ratioRegionNum-1; i++)
	{
		if(ratioRegionArray[i].multiReadsRatio>0 && ratioRegionArray[i+1].multiReadsRatio==0)
		{
			endHeadRow = i;
			break;
		}
	}
	for(i=ratioRegionNum-1; i>=0; i--)
	{
		if(ratioRegionArray[i].multiReadsRatio>0 && ratioRegionArray[i-1].multiReadsRatio==0)
		{
			startTailRow = i;
			break;
		}
	}

	for(i=0; i<ratioRegionNum; i++)
	{
		startRow = ratioRegionArray[i].startQPosLHalf - 1;
		endRow = ratioRegionArray[i].endQPosRHalf - 1;

		if(i<=endHeadRow || i>=startTailRow)
			multiRatio = 0;
		else
			multiRatio = ratioRegionArray[i].multiReadsRatio;

		if(multiRatio>0)
			fprintf(fpMultiRatio, "%s\t%d\t%d\t%.8f\n", queryLabel, startRow, endRow, multiRatio);
		else
			fprintf(fpMultiRatio, "%s\t%d\t%d\t%d\n", queryLabel, startRow, endRow, 0);
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
		{
			if(misInfo->misType==QUERY_MISJOIN_KIND)
				strcpy(colorStr, "color=vdred");
			else if(misInfo->misType==QUERY_INDEL_KIND)
				strcpy(colorStr, "color=lorange");
			else
			{
				printf("line=%d, In %s(), error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else if(misInfo->misassFlag==STRUCTURE_VARIATION)
		{
			if(misInfo->misType==QUERY_MISJOIN_KIND)
				strcpy(colorStr, "color=dgreen");
			else if(misInfo->misType==QUERY_INDEL_KIND)
				strcpy(colorStr, "color=dblue");
			else
			{
				printf("line=%d, In %s(), error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else if(misInfo->misassFlag==UNCERTAIN_MISASS)
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
