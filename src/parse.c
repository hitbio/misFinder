/*
 * parse.c
 *
 *  Created on: Nov 21, 2011
 *      Author: xiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Parse blastn result of queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short parseBlastn(char *parseResultFile, const char *blastnResultFile)
{
	int j, len, fileStatus, matchFlag;
	char line[LINE_CHAR_MAX+1], *pch;
	FILE *fpBlastnRe, *fpParseResult;
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


	fpBlastnRe = fopen(blastnResultFile, "r");
	if(fpBlastnRe==NULL)
	{
		printf("In %s(), cannot open file [ %s ], error!\n", __func__, blastnResultFile);
		return FAILED;
	}

	fpParseResult = fopen(parseResultFile, "w");
	if(fpParseResult==NULL)
	{
		printf("In %s(), cannot open file [ %s ], error!\n", __func__, parseResultFile);
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

	stage = MATCH_START_STAGE;
	step = MATCH_INFO_START_STEP;
	stepMatchFlag = MATCH_INFO_START_STEP;
	while(1)
	{
		if(readLine(&fileStatus, line, &len, LINE_CHAR_MAX, fpBlastnRe)==FAILED)
		{
			printf("In %s(), cannot read a line, error!\n", __func__);
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
					strcat(queryTitle, " ");
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

					fprintf(fpParseResult, ">%s\t%d\t%s\t%d\n", queryTitle, queryLen, subjectTitle, subjectLen);

					stage = SUBJECT_LEN_STAGE;
					continue;
				}
				else
				{
					strcat(subjectTitle, " ");
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
							queryPos_ch[j++] = *pch;
							pch ++;
						}
						queryPos_ch[j] = '\0';
						startQueryPos = atoi(queryPos_ch);
						//printf("startQueryPos=%d\n", startQueryPos);

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

						// get the end subject pos
						pch = line + len - 1;
						while(*pch!=' ')
						{
							pch --;
						}
						pch ++;

						endSubjectPos = atoi(pch);
						//printf("endSubjectPos=%d\n", endSubjectPos);

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
						pch = line + len - 1;
						while(*pch!=' ')
						{
							pch --;
						}
						pch ++;
						endSubjectPos = atoi(pch);
						//printf("endSubjectPos=%d\n", endSubjectPos);

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

					// print the result
					fprintf(fpParseResult, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\n", startSubjectPos, endSubjectPos, startQueryPos, endQueryPos, strand, matchLen, totalMatchLen, gapNum, matchPercent);


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
					// print the result
					fprintf(fpParseResult, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\n", startSubjectPos, endSubjectPos, startQueryPos, endQueryPos, strand, matchLen, totalMatchLen, gapNum, matchPercent);

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

	fclose(fpBlastnRe);
	fclose(fpParseResult);

	return SUCCESSFUL;
}

/**
 * Read a line in a file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short readLine(int *fileStatus, char *str, int *len, int maxChNum, FILE *fp)
{
	int i;
	char ch;

	i = 0;
	ch = fgetc(fp);
	while(ch!=EOF && ch!='\n')
	{
		str[i++] = ch;
		ch = fgetc(fp);
	}
	str[i] = '\0';

	if(i>maxChNum)
	{
		printf("Too small character numbers, please set a new one that larger than %d\n", i);
		return FAILED;
	}else
	{
		*len = i;
	}

	if(feof(fp))
		*fileStatus = EOF_STATUS;
	else
		*fileStatus = READING_STATUS;

	return SUCCESSFUL;
}

