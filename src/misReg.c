/*
 * misReg.c
 *
 *  Created on: Aug 8, 2014
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Extract the mis-assembly regions.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short extractMisReg(queryMatchInfo_t *queryMatchInfoSet)
{
	int32_t i, processedNum;
	query_t *queryArray;

	printf("Begin extracting regions for assembly errors and correct assemblies due to structural variations ...\n");

	processedNum = 0;
	queryArray = queryMatchInfoSet->queryArray;
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		// ########################### Debug information ##############################
		//if(queryArray[i].queryID==8 || strcmp(queryArray[i].queryTitle, "scf7180000013826")==0)
		//{
		//	printf("****** queryID=%d, queryTitle=%s, queryLen=%d, subjectNum=%d\n", queryArray[i].queryID, queryArray[i].queryTitle, queryArray[i].queryLen, queryArray[i].querySubjectNum);
		//}
		// ########################### Debug information ##############################

		// compute the structural variations in single query
		if(extractMisRegSingleQuery(queryArray+i, queryMatchInfoSet->subjectArray)==FAILED)
		{
			printf("line=%d, In %s(), cannot determine structural variations in single query, error!\n", __LINE__, __func__);
			return FAILED;
		}

		processedNum ++;
//		if(processedNum%100==0)
//			printf("Queries processed: %d\n", processedNum);
	}

//	if(processedNum%100!=0)
//		printf("Queries processed: %d\n", processedNum);

	return SUCCESSFUL;
}

/**
 * Extract the mis-assembly regions for single query.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short extractMisRegSingleQuery(query_t *queryItem, subject_t *subjectArray)
{
	misInfo_t *misInfo;
	queryMargin_t *queryMargin;
	queryIndel_t *queryIndel;
	int32_t startQueryPos, endQueryPos, startSubjectPos, endSubjectPos, subjectID, misassKind;
	int32_t leftSegRow, rightSegRow, tmp;
	globalValidSeg_t *globalSegArray;

	globalSegArray = queryItem->globalValidSegArray;

	misInfo = queryItem->misInfoList;
	while(misInfo)
	{
		misInfo->misassSeqList = misInfo->tailMisassSeq = NULL;
		misInfo->misassSeqNum = 0;

		if(misInfo->misType==QUERY_MISJOIN_KIND)
		{ // misjoin
			queryMargin = misInfo->queryMargin;

			subjectID = -1;
			if(queryMargin->leftMargin<queryMargin->rightMargin)
			{
				startQueryPos = queryMargin->leftMargin;
				endQueryPos = queryMargin->rightMargin;
			}else
			{
				endQueryPos = queryMargin->leftMargin;
				startQueryPos = queryMargin->rightMargin;
			}

			if(misInfo->misassFlag==TRUE_MISASS)
				misassKind = ERR_MISJOIN;
			else if(misInfo->misassFlag==STRUCTURE_VARIATION)
				misassKind = SV_MISJOIN;
			else if(misInfo->misassFlag==UNCERTAIN_MISASS)
				misassKind = UNCER_MISJOIN;
			else
			{
				printf("line=%d, In %s(), invalid misassFlag=%d, error!\n", __LINE__, __func__, misInfo->misassFlag);
				return FAILED;
			}

			if(addMisassSeqNodeToMisInfoNode(misInfo, misassKind, subjectID, startQueryPos, endQueryPos, -1, -1)==FAILED)
			{
				printf("line=%d, In %s(), cannot add misseqNode to misInfo, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else
		{ // indel

			queryIndel = misInfo->queryIndel;
			leftSegRow = misInfo->leftSegRow;
			rightSegRow = misInfo->rightSegRow;

			if(leftSegRow==-1 && rightSegRow>=0)
			{
				startQueryPos = 1;
				endQueryPos = queryIndel->leftMargin;
				subjectID = -1;
				startSubjectPos = endSubjectPos = -1;
			}else if(leftSegRow>=0 && rightSegRow==-1)
			{
				startQueryPos = queryIndel->rightMargin;
				endQueryPos = queryItem->queryLen;
				subjectID = -1;
				startSubjectPos = endSubjectPos = -1;
			}else if((leftSegRow>=0 && rightSegRow>=0) || misInfo->innerFlag==YES)
			{
				if(queryIndel->leftMargin<queryIndel->rightMargin)
				{
					startQueryPos = queryIndel->leftMargin;
					endQueryPos = queryIndel->rightMargin;
				}else
				{
					endQueryPos = queryIndel->leftMargin;
					startQueryPos = queryIndel->rightMargin;
				}

				subjectID = queryIndel->subjectID;
				startSubjectPos = queryIndel->leftMarginSubject;
				endSubjectPos = queryIndel->rightMarginSubject;
				if(startSubjectPos>endSubjectPos)
				{
					tmp = startSubjectPos;
					startSubjectPos = endSubjectPos;
					endSubjectPos = tmp;
				}
			}else
			{
				startQueryPos = queryIndel->leftMargin;
				endQueryPos = queryIndel->rightMargin;

				subjectID = -1;
				startSubjectPos = endSubjectPos = -1;
			}

			if(misInfo->misassFlag==TRUE_MISASS)
			{
				if(misInfo->gapFlag==NO)
				{
					if(queryIndel->queryIndelKind==QUERY_INSERT)
						misassKind = ERR_INSERT;
					else if(queryIndel->queryIndelKind==QUERY_DEL)
						misassKind = ERR_DEL;
					else if(queryIndel->queryIndelKind==QUERY_GAP)
					{
						if(queryIndel->difQuery>queryIndel->difSubject)
							misassKind = ERR_INSERT;
						else
							misassKind = ERR_DEL;
					}else
					{
						printf("line=%d, In %s(), invalid queryIndelKind=%d, error!\n", __LINE__, __func__, queryIndel->queryIndelKind);
						return FAILED;
					}
				}else
				{
					if(queryIndel->queryIndelKind==QUERY_INSERT)
						misassKind = GAP_INSERT;
					else if(queryIndel->queryIndelKind==QUERY_DEL)
						misassKind = GAP_DEL;
					else if(queryIndel->queryIndelKind==QUERY_GAP)
					{
						if(queryIndel->difQuery>queryIndel->difSubject)
							misassKind = GAP_INSERT;
						else
							misassKind = GAP_DEL;
					}else
					{
						printf("line=%d, In %s(), invalid queryIndelKind=%d, error!\n", __LINE__, __func__, queryIndel->queryIndelKind);
						return FAILED;
					}
				}
			}else if(misInfo->misassFlag==STRUCTURE_VARIATION)
			{
				if(queryIndel->queryIndelKind==QUERY_INSERT)
					misassKind = SV_INSERT;
				else if(queryIndel->queryIndelKind==QUERY_DEL)
					misassKind = SV_DEL;
				else if(queryIndel->queryIndelKind==QUERY_GAP)
				{
					if(queryIndel->difQuery>queryIndel->difSubject)
						misassKind = SV_INSERT;
					else
						misassKind = SV_DEL;
				}else
				{
					printf("line=%d, In %s(), invalid queryIndelKind=%d, error!\n", __LINE__, __func__, queryIndel->queryIndelKind);
					return FAILED;
				}
			}else if(misInfo->misassFlag==UNCERTAIN_MISASS)
			{
				if(queryIndel->queryIndelKind==QUERY_INSERT)
					misassKind = UNCER_INSERT;
				else if(queryIndel->queryIndelKind==QUERY_DEL)
					misassKind = UNCER_DEL;
				else if(queryIndel->queryIndelKind==QUERY_GAP)
				{
					if(queryIndel->difQuery>queryIndel->difSubject)
						misassKind = UNCER_INSERT;
					else
						misassKind = UNCER_DEL;
				}else
				{
					printf("line=%d, In %s(), invalid queryIndelKind=%d, error!\n", __LINE__, __func__, queryIndel->queryIndelKind);
					return FAILED;
				}
			}

			if(addMisassSeqNodeToMisInfoNode(misInfo, misassKind, subjectID, startQueryPos, endQueryPos, startSubjectPos, endSubjectPos)==FAILED)
			{
				printf("line=%d, In %s(), cannot add misseqNode to misInfo, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		misInfo = misInfo->next;
	}

	return SUCCESSFUL;
}

/**
 * Add misassSeq node to the list.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short addMisassSeqNodeToMisInfoNode(misInfo_t *misInfo, int32_t misassKind, int32_t subjectID, int32_t startQueryPos, int32_t endQueryPos, int32_t startSubjectPos, int32_t endSubjectPos)
{
	misassSeq_t *misassSeq;

	misassSeq = (misassSeq_t*) calloc (1, sizeof(misassSeq_t));
	if(misassSeq==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	misassSeq->misassKind = misassKind;
	misassSeq->startQueryPos = startQueryPos;
	misassSeq->endQueryPos = endQueryPos;
	misassSeq->subjectID = subjectID;
	misassSeq->startSubjectPos = startSubjectPos;
	misassSeq->endSubjectPos = endSubjectPos;
	misassSeq->next = NULL;

	if(misInfo->misassSeqList)
		misInfo->tailMisassSeq->next = misassSeq;
	else
		misInfo->misassSeqList = misassSeq;
	misInfo->tailMisassSeq = misassSeq;
	misInfo->misassSeqNum ++;

	return SUCCESSFUL;
}

