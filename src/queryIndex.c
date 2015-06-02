/*
 * queryIndex.c
 *
 *  Created on: May 29, 2014
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Build the queryIndex.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short buildQueryIndex(queryIndex_t **queryIndex, queryMatchInfo_t *queryMatchInfoSet)
{
	// initialize the queryIndex
	if(initQueryIndex(queryIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the queryIndex, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// count the queryKmer occurrences
	if(countQueryKmerOccs(*queryIndex, queryMatchInfoSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot count the kmers from queries, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// add the query pos information
	if(addQueryKmerQuerypos(*queryIndex, queryMatchInfoSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot count the kmers from queries, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Initialize queryIndex.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initQueryIndex(queryIndex_t **queryIndex)
{
	*queryIndex = (queryIndex_t *) calloc(1, sizeof(queryIndex_t));
	if((*queryIndex)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize k-mer block
	if(initQueryKmerBlock(*queryIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the queryKmer blocks in k-mer hash table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize kmerseq block
	if(initQueryKmerseqBlock(*queryIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the kmerseq blocks in k-mer hash table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Release the query index.
 */
void releaseQueryIndex(queryIndex_t **queryIndex)
{
	int32_t i;

	if((*queryIndex)==NULL)
		return;

	// free k-mer blocks
	if((*queryIndex)->kmerBlockArr)
	{
		for(i=0; i<(*queryIndex)->blocksNumKmer; i++)
			free((*queryIndex)->kmerBlockArr[i].kmerArr);
		(*queryIndex)->blocksNumKmer = 0;
		free((*queryIndex)->kmerBlockArr);
		(*queryIndex)->kmerBlockArr = NULL;
	}
	if((*queryIndex)->kmerHashtable)
	{
		free((*queryIndex)->kmerHashtable);
		(*queryIndex)->kmerHashtable = NULL;
	}

	// free kmerseq blocks
	if((*queryIndex)->kmerSeqBlockArr)
	{
		for(i=0; i<(*queryIndex)->blocksNumKmerSeq; i++)
			free((*queryIndex)->kmerSeqBlockArr[i].kmerSeqArr);
		(*queryIndex)->blocksNumKmerSeq = 0;
		free((*queryIndex)->kmerSeqBlockArr);
		(*queryIndex)->kmerSeqBlockArr = NULL;
	}

	// free query pos blocks
	if((*queryIndex)->queryPosBlockArr)
	{
		for(i=0; i<(*queryIndex)->blocksNumQuerypos; i++)
			free((*queryIndex)->queryPosBlockArr[i].queryPosArr);
		(*queryIndex)->blocksNumQuerypos = 0;
		free((*queryIndex)->queryPosBlockArr);
		(*queryIndex)->queryPosBlockArr = NULL;
	}

	// free index node
	free(*queryIndex);
	*queryIndex = NULL;
}

/**
 * Initialize queryKmer blocks.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initQueryKmerBlock(queryIndex_t *queryIndex)
{
	queryIndex->bytesPerKmer = sizeof(queryKmer_t);
	queryIndex->totalItemNumKmer = 0;
	queryIndex->maxBlocksNumKmer = MAX_BLOCKS_NUM_KMER;
	queryIndex->maxItemNumPerKmerBlock = BLOCK_SIZE_PER_KMER / queryIndex->bytesPerKmer;

	queryIndex->hashTableSize = HASH_TABLE_SIZE;
	queryIndex->kmerSize = KMER_SIZE_DEFAULT;

	queryIndex->kmerHashtable = (kmerHashBucket_t *) calloc (queryIndex->hashTableSize, sizeof(kmerHashBucket_t));
	if( queryIndex->kmerHashtable == NULL )
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate queryKmer blocks
	queryIndex->kmerBlockArr = (queryKmerBlock_t *) calloc(queryIndex->maxBlocksNumKmer, sizeof(queryKmerBlock_t));
	if( queryIndex->kmerBlockArr == NULL )
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	queryIndex->blocksNumKmer = 0;

	// add new queryKmer block
	if(addNewBlockQueryKmer(queryIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Add new queryKmer block.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewBlockQueryKmer(queryIndex_t *queryIndex)
{
	if(queryIndex->blocksNumKmer>=queryIndex->maxBlocksNumKmer)
	{
		queryIndex->kmerBlockArr = (queryKmerBlock_t *) realloc (queryIndex->kmerBlockArr, queryIndex->maxBlocksNumKmer*2*sizeof(queryKmerBlock_t));
		if(queryIndex->kmerBlockArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		queryIndex->maxBlocksNumKmer *= 2;
	}

	queryIndex->kmerBlockArr[queryIndex->blocksNumKmer].kmerArr = (queryKmer_t *) calloc (queryIndex->maxItemNumPerKmerBlock, queryIndex->bytesPerKmer);
	if(queryIndex->kmerBlockArr[queryIndex->blocksNumKmer].kmerArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	queryIndex->kmerBlockArr[queryIndex->blocksNumKmer].blockID = queryIndex->blocksNumKmer + 1;
	queryIndex->kmerBlockArr[queryIndex->blocksNumKmer].itemNum = 0;

	queryIndex->blocksNumKmer ++;

	return SUCCESSFUL;
}

/**
 * Initialize kmerseq blocks in query index.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initQueryKmerseqBlock(queryIndex_t *queryIndex)
{
	queryIndex->entriesPerKmer = ((queryIndex->kmerSize-1) / 32) + 1;
	if(queryIndex->kmerSize%32==0)
	{
		queryIndex->lastEntryMask = (uint64_t) -1;
		queryIndex->lastEntryBaseNum = 32;
	}else
	{
		queryIndex->lastEntryMask = (1LLU << ((queryIndex->kmerSize%32)<<1)) - 1;
		queryIndex->lastEntryBaseNum = queryIndex->kmerSize % 32;
	}
	queryIndex->bytesPerKmerseq = queryIndex->entriesPerKmer * sizeof(uint64_t);
	queryIndex->totalItemNumKmerSeq = 0;
	queryIndex->maxBlocksNumKmerSeq = MAX_BLOCKS_NUM_KMER_SEQ;
	queryIndex->maxItemNumPerKmerSeqBlock = BLOCK_SIZE_PER_KMER_SEQ / queryIndex->bytesPerKmerseq;

	// allocate kmerseq blocks
	queryIndex->kmerSeqBlockArr = (kmerseqBlock_t *) malloc(queryIndex->maxBlocksNumKmerSeq * sizeof(kmerseqBlock_t));
	if( queryIndex->kmerSeqBlockArr == NULL )
	{
		printf("line=%d, In %s(), cannot allocate the kmerSeqBlockArr for the k-mer hash table, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	queryIndex->blocksNumKmerSeq = 0;

	// add new kmerseq block
	if(addNewBlockQueryKmerSeq(queryIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Add new kmerseq block in query index.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewBlockQueryKmerSeq(queryIndex_t *queryIndex)
{
	if(queryIndex->blocksNumKmerSeq>=queryIndex->maxBlocksNumKmerSeq)
	{
		queryIndex->kmerSeqBlockArr = (kmerseqBlock_t *) realloc (queryIndex->kmerSeqBlockArr, queryIndex->maxBlocksNumKmerSeq*2*sizeof(kmerseqBlock_t));
		if(queryIndex->kmerSeqBlockArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		queryIndex->maxBlocksNumKmerSeq *= 2;
	}

	queryIndex->kmerSeqBlockArr[queryIndex->blocksNumKmerSeq].kmerSeqArr = (uint64_t *) calloc (queryIndex->maxItemNumPerKmerSeqBlock, queryIndex->bytesPerKmerseq);
	if(queryIndex->kmerSeqBlockArr[queryIndex->blocksNumKmerSeq].kmerSeqArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	queryIndex->kmerSeqBlockArr[queryIndex->blocksNumKmerSeq].blockID = queryIndex->blocksNumKmerSeq + 1;
	queryIndex->kmerSeqBlockArr[queryIndex->blocksNumKmerSeq].itemNum = 0;

	queryIndex->blocksNumKmerSeq ++;

	return SUCCESSFUL;
}

/**
 * Count the queryKmer occurrences.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short countQueryKmerOccs(queryIndex_t *queryIndex, queryMatchInfo_t *queryMatchInfoSet)
{
	int32_t i, j, queryLen, startBasePos, endBasePos, basePos, baseInt, newBasePos, validFlag;
	int32_t kmerSize, hashTableSize, entriesPerKmer, lastEntryBaseNum;
	char *querySeq;
	uint64_t *pKmerSeqIntDone, *pKmerSeqIntDoing, hashcode, lastEntryMask;
	kmerseqBlock_t *pKmerseqBlock;
	query_t *queryArray;

	pKmerseqBlock = queryIndex->kmerSeqBlockArr + queryIndex->blocksNumKmerSeq - 1;
	queryIndex->pKmerSeqAdding = pKmerseqBlock->kmerSeqArr;

	kmerSize = queryIndex->kmerSize;
	entriesPerKmer = queryIndex->entriesPerKmer;
	lastEntryBaseNum = queryIndex->lastEntryBaseNum;
	lastEntryMask = queryIndex->lastEntryMask;
	hashTableSize = queryIndex->hashTableSize;

	queryArray = queryMatchInfoSet->queryArray;
	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		querySeq = queryArray[i].querySeq;
		queryLen = queryArray[i].queryLen;

		// get start and end query positions
		startBasePos = endBasePos = -1;
		if(getStartEndQueryPos(&startBasePos, &endBasePos, querySeq, queryLen, kmerSize)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the start and end query positions, error!\n", __LINE__, __func__);
			return FAILED;
		}

		while(startBasePos>=0 && endBasePos>=0)
		{
			pKmerSeqIntDoing = queryIndex->pKmerSeqAdding;
			// generate the first kmerseqInt
			if(generateReadseqInt(pKmerSeqIntDoing, querySeq+startBasePos, kmerSize, entriesPerKmer)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate the integer sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}
			pKmerSeqIntDone = pKmerSeqIntDoing;

			// count the kmer
			hashcode = kmerhashInt(pKmerSeqIntDoing, entriesPerKmer, lastEntryBaseNum, hashTableSize);
			if(countQueryKmer(hashcode, pKmerSeqIntDoing, queryIndex)==FAILED)
			{
				printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// update startBasePos
			startBasePos += kmerSize;
			if(startBasePos>endBasePos)
			{
				// recalculate the start and end query positions
				if(getStartEndQueryPos(&startBasePos, &endBasePos, querySeq, queryLen, kmerSize)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the start and end query positions, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}else if(startBasePos+kmerSize-1>endBasePos)
				startBasePos = endBasePos - kmerSize + 1;
		}

/*
		// get start and end query positions
		if(getStartEndQueryPos(&startBasePos, &endBasePos, querySeq, queryLen, kmerSize)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the start and end query positions, error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(startBasePos<0 || endBasePos<0)
		{
			continue;
		}

		pKmerSeqIntDoing = queryIndex->pKmerSeqAdding;
		// generate the first kmerseqInt
		if(generateReadseqInt(pKmerSeqIntDoing, querySeq+startBasePos, kmerSize, entriesPerKmer)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the integer sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}
		pKmerSeqIntDone = pKmerSeqIntDoing;

		// count the first scafKmer
		hashcode = kmerhashInt(pKmerSeqIntDoing, entriesPerKmer, lastEntryBaseNum, hashTableSize);
		if(countQueryKmer(hashcode, pKmerSeqIntDoing, queryIndex)==FAILED)
		{
			printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// process other scafKmers
		for(basePos=startBasePos+kmerSize; basePos<=endBasePos; basePos++)
		{
			validFlag = YES;
			switch(querySeq[basePos])
			{
				case 'A':
				case 'a': baseInt = 0; break;
				case 'C':
				case 'c': baseInt = 1; break;
				case 'G':
				case 'g': baseInt = 2; break;
				case 'T':
				case 't': baseInt = 3; break;
				default: validFlag = NO;
			}

			if(validFlag==YES)
			{
				pKmerSeqIntDoing = queryIndex->pKmerSeqAdding;

				// generate the kmer integer sequence
				if(entriesPerKmer>=2)
				{
					for(j=0; j<entriesPerKmer-2; j++)
					{
						pKmerSeqIntDoing[j] = (pKmerSeqIntDone[j] << 2) | (pKmerSeqIntDone[j+1] >> 62);
					}
					pKmerSeqIntDoing[entriesPerKmer-2] = (pKmerSeqIntDone[entriesPerKmer-2] << 2) | (pKmerSeqIntDone[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
				}
				pKmerSeqIntDoing[entriesPerKmer-1] = ((pKmerSeqIntDone[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

				pKmerSeqIntDone = pKmerSeqIntDoing;

				hashcode = kmerhashInt(pKmerSeqIntDoing, entriesPerKmer, lastEntryBaseNum, hashTableSize);
				if(countQueryKmer(hashcode, pKmerSeqIntDoing, queryIndex)==FAILED)
				{
					printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}else
			{
				// get next valid query base position
				if(getNextValidQueryBasePos(&newBasePos, querySeq, basePos, endBasePos, kmerSize)==FAILED)
				{
					printf("line=%d, In %s(), cannot get next valid query base position, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(newBasePos<0)
					break;

				pKmerSeqIntDoing = queryIndex->pKmerSeqAdding;
				// generate the first kmerseqInt
				if(generateReadseqInt(pKmerSeqIntDoing, querySeq+newBasePos, kmerSize, entriesPerKmer)==FAILED)
				{
					printf("line=%d, In %s(), cannot generate the integer sequence, error!\n", __LINE__, __func__);
					return FAILED;
				}
				pKmerSeqIntDone = pKmerSeqIntDoing;

				// count the first scafKmer
				hashcode = kmerhashInt(pKmerSeqIntDoing, entriesPerKmer, lastEntryBaseNum, hashTableSize);
				if(countQueryKmer(hashcode, pKmerSeqIntDoing, queryIndex)==FAILED)
				{
					printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
					return FAILED;
				}

				basePos += kmerSize - 1;
			}
		}
*/
	}

	// the last k-mer block is empty, remove it
	if(queryIndex->kmerBlockArr[queryIndex->blocksNumKmer-1].itemNum==0)
	{
		free(queryIndex->kmerBlockArr[queryIndex->blocksNumKmer-1].kmerArr);
		queryIndex->kmerBlockArr[queryIndex->blocksNumKmer-1].kmerArr = NULL;
		queryIndex->blocksNumKmer --;
	}

	// the last kmerseq block is empty, remove it
	if(queryIndex->kmerSeqBlockArr[queryIndex->blocksNumKmerSeq-1].itemNum==0)
	{
		free(queryIndex->kmerSeqBlockArr[queryIndex->blocksNumKmerSeq-1].kmerSeqArr);
		queryIndex->kmerSeqBlockArr[queryIndex->blocksNumKmerSeq-1].kmerSeqArr = NULL;
		queryIndex->blocksNumKmerSeq --;
	}

	return SUCCESSFUL;
}

/**
 * Get the start and end query positions (start from 0).
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getStartEndQueryPos(int32_t *startBasePos, int32_t *endBasePos, char *querySeq, int32_t queryLen, int32_t kmerSize)
{

	int32_t i, j, baseNum, startRow, endRow, validStartRow, validBaseFlag;

	// compute the start base position
	if((*endBasePos)==-1)
		startRow = 0;
	else
		startRow = (*endBasePos) + 1;

	baseNum = 0;
	validStartRow = NO;
	for(i=startRow; i<=queryLen-kmerSize; i++)
	{
		if(checkValidQueryBase(&validBaseFlag, querySeq[i])==FAILED)
		{
			printf("line=%d, In %s(), cannot check the valid query base, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(validBaseFlag==YES)
		{
			if(baseNum==0)
				startRow = i;

			baseNum ++;
			if(baseNum>=kmerSize)
			{
				validStartRow = YES;
				break;
			}
		}else
		{
			baseNum = 0;
		}
	}

	if(validStartRow==YES)
		*startBasePos = startRow;
	else
		*startBasePos = *endBasePos = -1;

	// compute the end base position
	if((*startBasePos)>=0)
	{
		endRow = -1;
		startRow = (*startBasePos) + kmerSize;
		for(i=startRow; i<queryLen; i++)
		{
			if(checkValidQueryBase(&validBaseFlag, querySeq[i])==FAILED)
			{
				printf("line=%d, In %s(), cannot check the valid query base, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(validBaseFlag==NO)
			{
				endRow = i - 1;
				break;
			}
		}
		if(i>=queryLen)
			endRow = queryLen - 1;

		if(endRow>=0)
			*endBasePos = endRow;
		else
			*endBasePos = -1;
	}else
	{
		*endBasePos = -1;
	}

/*
	int32_t i, j, validFlag;
	char base1, base2;

	*startBasePos = *endBasePos = -1;
	for(i=0; i<=queryLen-kmerSize; i++)
	{
		base1 = querySeq[i];
		if(checkValidQueryBase(&validFlag, base1)==FAILED)
		{
			printf("line=%d, In %s(), cannot check the valid query base, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(validFlag==YES)
		{
			for(j=1; j<kmerSize; j++)
			{
				base2 = querySeq[i+j];
				if(checkValidQueryBase(&validFlag, base2)==FAILED)
				{
					printf("line=%d, In %s(), cannot check the valid query base, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(validFlag==NO)
				{
					validFlag = NO;
					break;
				}
			}
			if(validFlag==YES)
			{
				*startBasePos = i;
				break;
			}else
				i += j;
		}
	}

	for(i=queryLen-1; i>=kmerSize-1; i--)
	{
		base1 = querySeq[i];
		if(checkValidQueryBase(&validFlag, base1)==FAILED)
		{
			printf("line=%d, In %s(), cannot check the valid query base, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(validFlag==YES)
		{
			*endBasePos = i;
			break;
		}
	}
*/

	return SUCCESSFUL;
}

/**
 * Check the valid query base.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short checkValidQueryBase(int32_t *validFlag, char base)
{
	if(base=='A' || base=='a' || base=='C' || base=='c' || base=='G' || base=='g' || base=='T' || base=='t')
		*validFlag = YES;
	else
		*validFlag = NO;

	return SUCCESSFUL;
}

/**
 * Get next valid query base position.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNextValidQueryBasePos(int32_t *newBasePos, char *querySeq, int32_t basePos, int32_t endBasePos, int32_t kmerSize)
{
	int32_t i, j, validFlag;
	char base1, base2;

	*newBasePos = -1;
	for(i=basePos+1; i<=endBasePos-kmerSize; i++)
	{
		base1 = querySeq[i];
		if(checkValidQueryBase(&validFlag, base1)==FAILED)
		{
			printf("line=%d, In %s(), cannot check the valid query base, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(validFlag==YES)
		{
			for(j=1; j<kmerSize; j++)
			{
				base2 = querySeq[i+j];
				if(checkValidQueryBase(&validFlag, base2)==FAILED)
				{
					printf("line=%d, In %s(), cannot check the valid query base, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(validFlag==NO)
				{
					validFlag = NO;
					break;
				}
			}
			if(validFlag==YES)
			{
				*newBasePos = i;
				break;
			}else
				i += j;
		}
	}

	return SUCCESSFUL;
}

/**
 * Compute the kmer hash code.
 *  @return:
 *    The hash code.
 */
uint64_t kmerhashInt(uint64_t *seqInt, int32_t entriesPerKmer, int32_t lastEntryBaseNum, int64_t hashArraySize)
{
	uint64_t hashcode;
	int i, j;

	hashcode = 5381;
	for(i=0; i<entriesPerKmer-1; i++)
	{
		for(j=0; j<32; j++)
		{
			hashcode += (hashcode << 5) | ((seqInt[i] >> (62-2*j)) & 3);
		}
	}

	for(j=0; j<lastEntryBaseNum; j++)
	{
		hashcode += (hashcode << 5) | ((seqInt[entriesPerKmer-1] >> (2*(lastEntryBaseNum-j-1))) & 3);
	}

	//return (hashcode & 0x7FFFFFFF) % hashTableSize;
	return hashcode % hashArraySize;
}

/**
 * Count the k-mer occurrences.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short countQueryKmer(uint64_t hashcode, uint64_t *kmerSeqInt, queryIndex_t *queryIndex)
{
	queryKmer_t *kmer;
	queryKmerBlock_t *pKmerBlock;
	kmerseqBlock_t *pKmerseqBlock;

	kmer = getQueryKmerByHash(hashcode, kmerSeqInt, queryIndex);
	if(!kmer)
	{
		// process k-mer blocks
		pKmerBlock = queryIndex->kmerBlockArr + queryIndex->blocksNumKmer - 1;
		kmer = pKmerBlock->kmerArr + pKmerBlock->itemNum;

		//kmer->kmerseqBlockID = pKmerBlock->blockID;
		//kmer->itemRowKmerseqBlock = pKmerBlock->itemNum;
		kmer->arraysize = 1;
		kmer->nextKmerBlockID = queryIndex->kmerHashtable[hashcode].kmerBlockID;
		kmer->nextItemRowKmerBlock = queryIndex->kmerHashtable[hashcode].itemRowKmerBlock;
		queryIndex->kmerHashtable[hashcode].kmerBlockID = pKmerBlock->blockID;
		queryIndex->kmerHashtable[hashcode].itemRowKmerBlock = pKmerBlock->itemNum;

		pKmerBlock->itemNum ++;
		queryIndex->totalItemNumKmer ++;

		if(pKmerBlock->itemNum >= queryIndex->maxItemNumPerKmerBlock)
		{
			// add new kmer block
			if(addNewBlockQueryKmer(queryIndex)==FAILED)
			{
				printf("line=%d, In %s(), cannot add new k-mer block, Error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		// process kmerseq blocks
		pKmerseqBlock = queryIndex->kmerSeqBlockArr + queryIndex->blocksNumKmerSeq - 1;
		kmer->kmerseqBlockID = pKmerseqBlock->blockID;
		kmer->itemRowKmerseqBlock = pKmerseqBlock->itemNum;
		pKmerseqBlock->itemNum ++;
		queryIndex->totalItemNumKmerSeq ++;
		queryIndex->pKmerSeqAdding += queryIndex->entriesPerKmer;

		if(pKmerseqBlock->itemNum >= queryIndex->maxItemNumPerKmerSeqBlock)
		{
			// add new kmerseq block
			if(addNewBlockQueryKmerSeq(queryIndex)==FAILED)
			{
				printf("line=%d, In %s(), cannot add new kmerseq block, Error!\n", __LINE__, __func__);
				return FAILED;
			}
			queryIndex->pKmerSeqAdding = queryIndex->kmerSeqBlockArr[queryIndex->blocksNumKmerSeq-1].kmerSeqArr;
		}
	}else
	{
		kmer->arraysize ++;
	}

	return SUCCESSFUL;
}

queryKmer_t *getQueryKmerByHash(uint64_t hashvalue, uint64_t *kmerSeqInt, queryIndex_t *queryIndex)
{
	queryKmer_t *kmer;
	uint64_t *kmerseq;
	kmerHashBucket_t *pKmerBucket;

	pKmerBucket = queryIndex->kmerHashtable + hashvalue;
	if(pKmerBucket->kmerBlockID>0)
	{
		kmer = queryIndex->kmerBlockArr[pKmerBucket->kmerBlockID-1].kmerArr + pKmerBucket->itemRowKmerBlock;
		while(kmer)
		{
			kmerseq = queryIndex->kmerSeqBlockArr[kmer->kmerseqBlockID-1].kmerSeqArr + kmer->itemRowKmerseqBlock * queryIndex->entriesPerKmer;
			if(identicalKmerSeq(kmerSeqInt, kmerseq, queryIndex->entriesPerKmer)==YES)
				break;

			if(kmer->nextKmerBlockID>0)
				kmer = queryIndex->kmerBlockArr[kmer->nextKmerBlockID-1].kmerArr + kmer->nextItemRowKmerBlock;
			else
				kmer = NULL;
		}
	}else
	{
		kmer = NULL;
	}

	return kmer;
}

/**
 * Check whether the two sequence is identical.
 *  @return:
 *  	If identical, return YES; otherwise, return NO.
 */
short identicalKmerSeq(uint64_t *kmerSeqInt1, uint64_t *kmerSeqInt2, int32_t entriesNum)
{
	int32_t i;
	for(i=0; i<entriesNum; i++)
	{
		if(kmerSeqInt1[i] != kmerSeqInt2[i])
			return NO;
	}

	return YES;
}

/**
 * Add the k-mer query pos information.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addQueryKmerQuerypos(queryIndex_t *queryIndex, queryMatchInfo_t *queryMatchInfoSet)
{
	int32_t i, j, queryID, queryLen, newBasePos, startBasePos, queryPos, endBasePos, basePos, baseInt, validFlag;
	int32_t kmerSize, hashTableSize, entriesPerKmer, lastEntryBaseNum;
	char *querySeq;
	uint64_t pKmerSeqInt[queryIndex->entriesPerKmer], hashcode, lastEntryMask;
	query_t *queryArray;

	queryArray = queryMatchInfoSet->queryArray;
	kmerSize = queryIndex->kmerSize;
	entriesPerKmer = queryIndex->entriesPerKmer;
	lastEntryBaseNum = queryIndex->lastEntryBaseNum;
	lastEntryMask = queryIndex->lastEntryMask;
	hashTableSize = queryIndex->hashTableSize;


	// initialize contigpos blocks and the ridpos regions in that blocks
	if(initQueryposBlocksInQueryIndex(queryIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<queryMatchInfoSet->itemNumQueryArray; i++)
	{
		queryID = queryArray[i].queryID;
		querySeq = queryArray[i].querySeq;
		queryLen = queryArray[i].queryLen;

		// get start and end query positions
		startBasePos = endBasePos = -1;
		if(getStartEndQueryPos(&startBasePos, &endBasePos, querySeq, queryLen, kmerSize)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the start and end query positions, error!\n", __LINE__, __func__);
			return FAILED;
		}

		queryPos = startBasePos + 1;
		while(startBasePos>=0 && endBasePos>=0)
		{
			// generate the first kmerseqInt
			if(generateReadseqInt(pKmerSeqInt, querySeq+startBasePos, kmerSize, entriesPerKmer)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate the integer sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// count the kmer
			hashcode = kmerhashInt(pKmerSeqInt, entriesPerKmer, lastEntryBaseNum, hashTableSize);
			if(addQueryKmer(hashcode, pKmerSeqInt, queryID, queryPos, queryIndex)==FAILED)
			{
				printf("line=%d, In %s(), cannot add kmer occurrence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// update startBasePos
			startBasePos += kmerSize;
			queryPos = startBasePos + 1;
			if(startBasePos>endBasePos)
			{
				// recalculate the start and end query positions
				if(getStartEndQueryPos(&startBasePos, &endBasePos, querySeq, queryLen, kmerSize)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the start and end query positions, error!\n", __LINE__, __func__);
					return FAILED;
				}
				queryPos = startBasePos + 1;
			}else if(startBasePos+kmerSize-1>endBasePos)
			{
				startBasePos = endBasePos - kmerSize + 1;
				queryPos = startBasePos + 1;
			}
		}

/*
		// get start and end query positions
		if(getStartEndQueryPos(&startBasePos, &endBasePos, querySeq, queryLen, kmerSize)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the start and end query positions, error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(startBasePos<0 || endBasePos<0)
		{
			continue;
		}

		queryPos = startBasePos + 1;

		// generate the first kmerseqInt
		if(generateReadseqInt(pKmerSeqInt, querySeq+startBasePos, kmerSize, entriesPerKmer)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the integer sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// count the first scafKmer
		hashcode = kmerhashInt(pKmerSeqInt, entriesPerKmer, lastEntryBaseNum, hashTableSize);
		if(addQueryKmer(hashcode, pKmerSeqInt, queryID, queryPos, queryIndex)==FAILED)
		{
			printf("line=%d, In %s(), cannot add kmer occurrence, error!\n", __LINE__, __func__);
			return FAILED;
		}
		queryPos ++;

		// process other scafKmers
		for(basePos=startBasePos+kmerSize; basePos<=endBasePos; basePos++, queryPos++)
		{
			validFlag = YES;
			switch(querySeq[basePos])
			{
				case 'A':
				case 'a': baseInt = 0; break;
				case 'C':
				case 'c': baseInt = 1; break;
				case 'G':
				case 'g': baseInt = 2; break;
				case 'T':
				case 't': baseInt = 3; break;
				default: validFlag = NO;
			}

			if(validFlag==YES)
			{
				// generate the kmer integer sequence
				if(entriesPerKmer>=2)
				{
					for(j=0; j<entriesPerKmer-2; j++)
					{
						pKmerSeqInt[j] = (pKmerSeqInt[j] << 2) | (pKmerSeqInt[j+1] >> 62);
					}
					pKmerSeqInt[entriesPerKmer-2] = (pKmerSeqInt[entriesPerKmer-2] << 2) | (pKmerSeqInt[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
				}
				pKmerSeqInt[entriesPerKmer-1] = ((pKmerSeqInt[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

				hashcode = kmerhashInt(pKmerSeqInt, entriesPerKmer, lastEntryBaseNum, hashTableSize);
				if(addQueryKmer(hashcode, pKmerSeqInt, queryID, queryPos, queryIndex)==FAILED)
				{
					printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}else
			{
				// get next valid query base position
				if(getNextValidQueryBasePos(&newBasePos, querySeq, basePos, endBasePos, kmerSize)==FAILED)
				{
					printf("line=%d, In %s(), cannot get next valid query base position, error!\n", __LINE__, __func__);
					return FAILED;
				}
				if(newBasePos<0)
					break;

				queryPos = newBasePos + 1;

				// generate the first kmerseqInt
				if(generateReadseqInt(pKmerSeqInt, querySeq+newBasePos, kmerSize, entriesPerKmer)==FAILED)
				{
					printf("line=%d, In %s(), cannot generate the integer sequence, error!\n", __LINE__, __func__);
					return FAILED;
				}

				hashcode = kmerhashInt(pKmerSeqInt, entriesPerKmer, lastEntryBaseNum, hashTableSize);
				if(addQueryKmer(hashcode, pKmerSeqInt, queryID, queryPos, queryIndex)==FAILED)
				{
					printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
					return FAILED;
				}

				basePos += kmerSize - 1;
			}
		}
*/
	}

	return SUCCESSFUL;
}

/**
 * Initialize queryPos block in query index.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initQueryposBlocksInQueryIndex(queryIndex_t *queryIndex)
{
	int64_t i, j, totalNum, sum;
	queryKmer_t *kmer;
	uint32_t blocksNumKmer, itemNumKmerBlock;

	queryIndex->bytesPerQuerypos = sizeof(queryPosBlock_t);
	queryIndex->blocksNumQuerypos = 0;
	queryIndex->totalItemNumQuerypos = 0;
	queryIndex->maxBlocksNumQuerypos = MAX_BLOCKS_NUM_RIDPOS;
	queryIndex->maxItemNumPerQueryposBlock = BLOCK_SIZE_PER_RIDPOS / queryIndex->bytesPerQuerypos;

	// add querypos block head nodes
	queryIndex->queryPosBlockArr = (queryPosBlock_t *) malloc (queryIndex->maxBlocksNumQuerypos * sizeof(queryPosBlock_t));
	if(queryIndex->queryPosBlockArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// add new querypos block
	if(addNewBlockQuerypos(queryIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// set start querypos regions and allocate querypos blocks
	totalNum = sum = 0;
	blocksNumKmer = queryIndex->blocksNumKmer;
	for(i=0; i<blocksNumKmer; i++)
	{
		kmer = queryIndex->kmerBlockArr[i].kmerArr;
		itemNumKmerBlock = queryIndex->kmerBlockArr[i].itemNum;
		for(j=0; j<itemNumKmerBlock; j++)
		{
			kmer->ppos = queryIndex->queryPosBlockArr[queryIndex->blocksNumQuerypos-1].queryPosArr + sum;

			sum += kmer->arraysize;
			if(sum >= queryIndex->maxItemNumPerQueryposBlock)
			{
				if(sum > queryIndex->maxItemNumPerQueryposBlock)
				{
					sum -= kmer->arraysize;
					queryIndex->queryPosBlockArr[queryIndex->blocksNumQuerypos-1].itemNum = sum;
					totalNum += sum;

					// add new querypos block
					if(addNewBlockQuerypos(queryIndex)==FAILED)
					{
						printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
						return FAILED;
					}

					kmer->ppos = queryIndex->queryPosBlockArr[queryIndex->blocksNumQuerypos-1].queryPosArr;

					sum = kmer->arraysize;
				}else
				{
					queryIndex->queryPosBlockArr[queryIndex->blocksNumQuerypos-1].itemNum = sum;
					totalNum += sum;

					// add new querypos block
					if(addNewBlockQuerypos(queryIndex)==FAILED)
					{
						printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
						return FAILED;
					}
					sum = 0;
				}
			}
			kmer ++;
		}
	}

	// process the last querypos block
	if(sum==0)
	{ // remove the last querypos block if it is empty
		free(queryIndex->queryPosBlockArr[queryIndex->blocksNumQuerypos-1].queryPosArr);
		queryIndex->queryPosBlockArr[queryIndex->blocksNumQuerypos-1].queryPosArr = NULL;
		queryIndex->queryPosBlockArr[queryIndex->blocksNumQuerypos-1].blockID = 0;
		queryIndex->queryPosBlockArr[queryIndex->blocksNumQuerypos-1].itemNum = 0;

		queryIndex->blocksNumQuerypos --;
	}else
	{
		queryIndex->queryPosBlockArr[queryIndex->blocksNumQuerypos-1].itemNum = sum;

		totalNum += sum;
	}

	queryIndex->totalItemNumQuerypos = totalNum;

	return SUCCESSFUL;
}

/**
 * Add new querypos block.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewBlockQuerypos(queryIndex_t *queryIndex)
{
	if(queryIndex->blocksNumQuerypos>=queryIndex->maxBlocksNumQuerypos)
	{
		queryIndex->queryPosBlockArr = (queryPosBlock_t *) realloc (queryIndex->queryPosBlockArr, queryIndex->maxBlocksNumQuerypos*2*sizeof(queryPosBlock_t));
		if(queryIndex->queryPosBlockArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		queryIndex->maxBlocksNumQuerypos *= 2;
	}

	queryIndex->queryPosBlockArr[queryIndex->blocksNumQuerypos].queryPosArr = (queryPos_t *) calloc (queryIndex->maxItemNumPerQueryposBlock, queryIndex->bytesPerQuerypos);
	if(queryIndex->queryPosBlockArr[queryIndex->blocksNumQuerypos].queryPosArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	queryIndex->queryPosBlockArr[queryIndex->blocksNumQuerypos].blockID = queryIndex->blocksNumQuerypos + 1;
	queryIndex->queryPosBlockArr[queryIndex->blocksNumQuerypos].itemNum = 0;

	queryIndex->blocksNumQuerypos ++;

	return SUCCESSFUL;
}

/**
 * Add a kmer to query index.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addQueryKmer(uint64_t hashcode, uint64_t *kmerSeqInt, int32_t queryID, int32_t queryPos, queryIndex_t *queryIndex)
{
	queryKmer_t *kmer;
	queryPos_t *querypos;

	kmer = getQueryKmerByHash(hashcode, kmerSeqInt, queryIndex);
	if(kmer)
	{
		querypos = kmer->ppos + kmer->multiplicity;
		querypos->queryID = queryID;
		querypos->queryPos = queryPos;

		kmer->multiplicity ++;
	}else
	{
		printf("line=%d, In %s(), cannot get k-mer, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get the kmer bases from integer.
 */
short getKmerBaseByInt(char *baseSeq, uint64_t *kmerSeqInt, int32_t entriesPerKmer, int32_t lastEntryBaseNum)
{
	int32_t i, j, k, baseInt;

	k = 0;
	for(i=0; i<entriesPerKmer-1; i++)
	{
		for(j=0; j<32; j++)
		{
			baseInt = (kmerSeqInt[i] >> (62-2*j)) & 3;
			switch(baseInt)
			{
				case 0: baseSeq[k]='A'; break;
				case 1: baseSeq[k]='C'; break;
				case 2: baseSeq[k]='G'; break;
				case 3: baseSeq[k]='T'; break;
			}
			k ++;
		}
	}

	for(j=0; j<lastEntryBaseNum; j++)
	{
		baseInt = (kmerSeqInt[entriesPerKmer-1] >> (2*(lastEntryBaseNum-j-1))) & 3;
		switch(baseInt)
		{
			case 0: baseSeq[k]='A'; break;
			case 1: baseSeq[k]='C'; break;
			case 2: baseSeq[k]='G'; break;
			case 3: baseSeq[k]='T'; break;
		}
		k ++;
	}

	baseSeq[k] = '\0';

	return SUCCESSFUL;
}

/**
 * Fill queries.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillQueries(queryMatchInfo_t *queryMatchInfoSet, char *inputQueryFile)
{
	int64_t i, maxQueryLen, queryID, queryLen, returnFlag;
	query_t *queryArray;
	char *querySeq, queryHeadTitle[1000], *pch;
	FILE *fpQuery;

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

	queryArray = queryMatchInfoSet->queryArray;

	queryID = 1;
	while((returnFlag=getSingleFastaItemFromFile(fpQuery, queryHeadTitle, querySeq, &queryLen))==SUCCESSFUL)
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

		if(strcmp(queryHeadTitle, queryArray[queryID-1].queryTitle)==0)
		{
			queryArray[queryID-1].querySeq = (char*) calloc (queryArray[queryID-1].queryLen+1, sizeof(char));
			if(queryArray[queryID-1].querySeq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			strcpy(queryArray[queryID-1].querySeq, querySeq);
		}else
		{
			printf("line=%d, In %s(), cannot get the potential mis-assembled queries, error!\n", __LINE__, __func__);
			return FAILED;
		}

		queryID ++;
	}

	free(querySeq);
	querySeq = NULL;

	fclose(fpQuery);

	return SUCCESSFUL;
}
