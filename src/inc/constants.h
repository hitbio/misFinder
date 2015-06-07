/*
 * constants.h
 *
 *  Created on: May 23, 2012
 *      Author: zhuxiao
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_ 1

#define __USE_FILE_OFFSET64
#define __USE_LARGEFILE64

#ifndef NULL
#define NULL ((void *)0)
#endif

#define MISFINDER_VERSION_STR			("v0.4.05.02")
#define MISFINDER_RELEASE_DATE_STR		("Jun 7, 2015")


#define DEBUG_SCAF_OVERLAP_FLAG		(NO)


#define SUCCESSFUL					0
#define FAILED						1
#define ERROR						-1

#define YES							1
#define NO							0


#define OPERATION_MODE_ALL				0
#define OPERATION_MODE_MERGE			1
#define OPERATION_MODE_METRICS			2
#define OPERATION_MODE_MISASS			3


#define ONLY_MATCHED_FOR_COV		(YES)

#define LINE_CHAR_MAX				1000

#define PLUS_STRAND					0
#define MINUS_STRAND				1

#define ORIENTATION_PLUS			0
#define ORIENTATION_MINUS			1

// file status
#define READING_STATUS				0
#define EOF_STATUS					1

#define COL_NUM						9


// match head information
// 0--match start, 1--queryTitle, 2--queryLen, 4--subjectHead, 4--subjectLen, 5--matchInfo, 6--unmatchInfo, 7--lambda, 8--finished
#define MATCH_START_STAGE			0
#define QUERY_TITLE_STAGE			1
#define QUERY_LEN_STAGE				2
#define SUBJECT_HEAD_STAGE			3
#define SUBJECT_LEN_STAGE			4
#define MATCH_INFO_STAGE			5
#define UNMATCH_INFO_STAGE			6
#define LAMBDA_STAGE				7
#define MATCH_FINISHED_STAGE		8

// step head information
// 0--start, 1--Score, 2--Identities/Gaps, 3--Strand, 4--step match info, 5--finish
#define MATCH_INFO_START_STEP		0  /////////
#define SCORE_STEP					1
#define IDENTITY_GAP_STEP			2
#define STRAND_STEP					3
#define STEP_INFO_STEP				4
#define MATCH_INFO_FINSISH_STEP		5

// step match information
// 0--start, 1--query, 2--middle, 3--subject, 4--finish
#define STEP_START_FLAG		0
#define STEP_QUERY_FLAG		1
#define STEP_MIDDLE_FLAG	2
#define STEP_SUBJECT_FLAG	3
#define STEP_FINISH_FLAG	4

// skips of characters for the first line
#define QUERY_ID_SKIP_NUM			6
#define QUERY_LEN_SKIP_NUM			7
#define SUBJECT_HEAD_SKIP_NUM		8
#define SUBJECT_LEN_SKIP_NUM		7
#define SCORE_SKIP_NUM				8
#define IDENTITY_SKIP_NUM			13
//#define GAPS_SKIP_NUM				6
#define STRAND_SKIP_NUM				8
#define STEP_QUERY_SKIP_NUM			5
#define STEP_SUBJECT_SKIP_NUM		5
#define LAMBDA_SKIP_NUM				6

#define MAX_MATCH_ITEM_NUM		100000
#define MAX_ARR_SIZE_SEG_LINK	100000

// Match kinds: 0-- perfect match, 1--matched, 2--disjunct match, 3--unmatched
#define PERFECT_MATCH_KIND			0
#define MATCHED_KIND				1
#define DISJUNCT_MATCH_KIND			2
#define UNMATCHED_KIND				3

#define UNCIRCULAR_FLAG				0
#define CIRCULAR_FLAG				1

#define MIN_QUERY_LEN_THRES				100
#define SHORT_QUERY_LEN_THRES			200
#define MATCHED_PERCENT_THRES			0.95f
#define MATCH_PERCENT_FACTOR			0.8f
//#define MIN_TOTAL_MATCH_LEN_THRES		30
#define VARY_LEN_THRES					1000
//#define END_IGNORE_LEN					30	// deleted 2013-04-23
#define END_IGNORE_LEN					300		// added 2013-04-23
#define MIN_DISJUNCT_DISTANCE_THRES		5000
//#define MIN_DISJUNCT_DISTANCE_THRES		10000
#define MIN_ALIGNED_SEG_LEN_THRES		100

//#define FRAGMENT_SIZE		200
#define VARY_GAPSIZE		100

#define MIN_BASE_FLAG_DELETION		27


//=============== misass ================
#define MAX_READ_BUF_SIZE				10000
#define MAX_READ_LEN_IN_BUF				5000

#define FILE_FORMAT_FASTA				1
#define FILE_FORMAT_FASTQ				2

#define PE_READ_TYPE					1
#define MP_READ_TYPE					2
#define MP_IN_READ_TYPE					3
#define LONG_PE_READ_TYPE				4
#define LONG_SE_READ_TYPE				5

#define BLOCK_SIZE_PER_READ				(1 << 26)	// 64 MB
#define MAX_BLOCKS_NUM_READ				512

#define BLOCK_SIZE_PER_READ_SEQ			(1 << 26)	// 64 MB
#define MAX_BLOCKS_NUM_READ_SEQ			512

#define BLOCK_SIZE_PER_READ_SEQ_HASH	(1 << 26)	// 64 MB
#define MAX_BLOCKS_NUM_READ_SEQ_HASH	512

#define BLOCK_SIZE_PER_KMER_SEQ			(1 << 26)	// 64 MB
#define MAX_BLOCKS_NUM_KMER_SEQ			512

#define BLOCK_SIZE_PER_KMER				(1 << 26)	// 64 MB
#define MAX_BLOCKS_NUM_KMER				512

#define BLOCK_SIZE_PER_RIDPOS			(1 << 27)	// 128 MB
#define MAX_BLOCKS_NUM_RIDPOS			1024

#define READS_NUM_PER_FILE_SAMPLE		10000
#define UNKNOWN_BASE_REPLACE_CHAR		'C'
#define MAX_UNKNOWN_BASE_NUM			0

#define ARTIFACTS_BASE_A_THRESHOLD 		0.9f

#define HASH_TABLE_SIZE					(15728681LLU)
#define KMER_SIZE_DEFAULT				25

//=============== misass ===============
#define MATCH_SCORE						1
#define MISMATCH_SCORE					-2
#define GAP_SCORE						-4

#define MISMATCH_WIN_SIZE				100
#define MIS_DENSITY_THRES				0.1f

#define DISAGREE_RATIO_THRES			0.9f

#define UNUSED_MISASS					0
#define POTENTIAL_MISASS				1
#define TRUE_MISASS						2
#define STRUCTURE_VARIATION				3
#define UNCERTAIN_MISASS				4

#define QUERY_MISJOIN_KIND				1
#define QUERY_INDEL_KIND				2

#define QUERY_INDEL_UNCERTAIN			0
#define QUERY_INSERT					1
#define QUERY_DEL						2
#define QUERY_GAP						3

#define INDEL_SIZE_DEFAULT				5
//#define INDEL_SIZE_DEFAULT				1
#define END_BACK_CHECK_SIZE				100

#define ERR_MISJOIN						1
#define ERR_INSERT						2
#define ERR_DEL							3
#define SV_MISJOIN						4
#define SV_INSERT						5
#define SV_DEL							6
#define GAP_INSERT						7
#define GAP_DEL							8
#define UNCER_MISJOIN					9
#define UNCER_INSERT					10
#define UNCER_DEL						11


#endif /* CONSTANTS_H_ */
