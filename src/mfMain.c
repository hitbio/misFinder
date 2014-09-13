/*
 * mfMain.c
 *
 *  Created on: Nov 21, 2011
 *      Author: xiao
 */

#include "inc/stdinc.h"
#include "inc/global.h"


int main(int argc, char **argv)
{
	if(parseCommandParasAndExe(argc, argv)==ERROR)
	{
		printf("line=%d, In %s(), cannot parse the command parameters, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Parse the input parameters.
 *  @parameters:
 *  	(1) -q <file>
 *  			Query file in fasta format.
 *  	(2) -s <file>
 *  			Subjects configuration file. It is necessary for -op=1 and op=1.
 *  	(3) -o(-out) <directory>
 *  			Output directory for the results. If it does not assigned, the default current command path will be assigned instead.
 *  	(4) -m <int>
 *  			The minimal query length to be considered in globalMetrics. Default 100 bp.
 *  	(5) -p <value>
 *  			The minimal identity percentage for matched queries and matched segments.
 *  	(6) -t <int>
 *  			The thread number to run the alignment between queries and subjects.
 *  	(6) -h(--help)
 *  			Show help information.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL;
 *  	if the job is executed and failed, return ERROR;
 *  	otherwise, return FAILED.
 */
short parseCommandParasAndExe(int argc, char **argv)
{
	int32_t i, operationMode;
	char queryFilePara[256];
	char subjectsFilePara[256];
	char outputDirPara[256];
	char *readFilesPara[256];
	char readFilesBuf[256][256];
	int64_t minQueryLenPara, threadNumPara;
	double minIdentityPercentPara;
	int32_t readFileNumPara, pairedModePara, indelSizeThresPara;

	if(argc==1)
	{
		if(showUsageInfo()==FAILED)
		{
			printf("line=%d, In %s(), cannot show the usage information for globalMetrics, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else if(argc==2)
	{
		i = 1;
		while(i<argc)
		{
			if(strcmp(argv[i], "all")==0 || strcmp(argv[i], "merge")==0 || strcmp(argv[i], "metrics")==0)
			{
				if(showUsageInfo()==FAILED)
				{
					printf("line=%d, In %s(), cannot show the usage information for globalMetrics, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}else if(strcmp(argv[i], "-h")==0 || strcmp(argv[i], "-help")==0)
			{
				if(showUsageInfo()==FAILED)
				{
					printf("line=%d, In %s(), cannot show the usage information for globalMetrics, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}else
			{
				if(argv[i][0]=='-')
				{
					printf("%s : unknown argument\n", argv[i]);
					printf("Please use -h or -help for more information.\n");
				}else
				{
					printf("%s : invalid argument\n", argv[i]);
					printf("Please use -h or -help for more information.\n");
				}

				return FAILED;
			}

			i ++;
		}

		return FAILED;
	}else
	{
		// reset the parameters
		operationMode = -1;
		queryFilePara[0] = '\0';
		subjectsFilePara[0] = '\0';
		outputDirPara[0] = '\0';
		minQueryLenPara = 0;
		minIdentityPercentPara = 0;
		threadNumPara = 0;
		readFileNumPara = 0;
		pairedModePara = 0;
		indelSizeThresPara = 0;

		if(strcmp(argv[1], "all")==0)
		{
			operationMode = OPERATION_MODE_ALL;
		}else if(strcmp(argv[1], "merge")==0)
		{
			operationMode = OPERATION_MODE_MERGE;
		}else if(strcmp(argv[1], "metrics")==0)
		{
			operationMode = OPERATION_MODE_METRICS;
		}else if(strcmp(argv[1], "misass")==0)
		{
			operationMode = OPERATION_MODE_MISASS;
		}else
		{
			if(strcmp(argv[1], "-h")==0 || strcmp(argv[1], "-help")==0)
			{
				if(showUsageInfo()==FAILED)
				{
					printf("line=%d, In %s(), cannot show the usage information for globalMetrics, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}else
			{
				if(argv[1][0]=='-')
				{
					printf("%s : unknown argument\n", argv[1]);
					printf("Please use -h or -help for more information.\n");
				}else
				{
					printf("%s : invalid argument\n", argv[1]);
					printf("Please use -h or -help for more information.\n");
				}
			}

			return FAILED;
		}

		i = 2;
		while(i<argc)
		{

			if(strcmp(argv[i], "-q")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						strcpy(queryFilePara, argv[i+1]);
					}else
					{
						printf("Exception: please specify the correct query file.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the query file.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-s")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						strcpy(subjectsFilePara, argv[i+1]);
					}else
					{
						printf("Exception: please specify the correct subject configuration file.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the subject configuration file.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-o")==0 || strcmp(argv[i], "-out")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						strcpy(outputDirPara, argv[i+1]);
					}else
					{
						printf("Exception: please specify the correct output directory.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the output directory.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-m")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						minQueryLenPara = atol(argv[i+1]);
						if(minQueryLenPara<=0)
						{
							printf("Exception: please specify the correct minimal query length.\n");
							return FAILED;
						}
					}else
					{
						printf("Exception: please specify the correct minimal query length.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the minimal query length.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-pt")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						minIdentityPercentPara = atof(argv[i+1]);
						if(minIdentityPercentPara<=0 || minIdentityPercentPara>1.0)
						{
							printf("Exception: please specify the correct minimal aligned identity percent value, error and exit!\n");
							return FAILED;
						}
					}else
					{
						printf("Exception: please specify the correct minimal aligned identity percent value.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the minimal aligned identity percent value.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-p")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						pairedModePara = atol(argv[i+1]);
						if(pairedModePara<1 || pairedModePara>2)
						{
							printf("Exception: please specify the correct paired mode.\n");
							return FAILED;
						}
					}else
					{
						printf("Exception: please specify the paired mode.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the paired mode.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-f")==0)
			{ // read files
				readFileNumPara = 0;

				i++;
				while(i<argc)
				{
					if(argv[i][0]!='-')
					{
						strcpy(readFilesBuf[readFileNumPara], argv[i]);
						readFilesPara[readFileNumPara] = readFilesBuf[readFileNumPara];
						readFileNumPara ++;
					}else
					{ // next option
						break;
					}

					i ++;
				}

				if(readFileNumPara==0)
				{
					printf("Exception: please specify the read files.\n");
					return FAILED;
				}
			}else if(strcmp(argv[i], "-t")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						threadNumPara = atol(argv[i+1]);
						if(threadNumPara<=0)
						{
							printf("Exception: please specify the correct number of threads.\n");
							return FAILED;
						}
					}else
					{
						printf("Exception: please specify the correct number of threads.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the number of threads.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-i")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						indelSizeThresPara = atol(argv[i+1]);
						if(indelSizeThresPara<=0)
						{
							printf("Exception: please specify the minimal indel size.\n");
							return FAILED;
						}
					}else
					{
						printf("Exception: please specify the correct minimal indel size.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the minimal indel size.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-h")==0 || strcmp(argv[i], "-help")==0)
			{
				if(showUsageInfo()==FAILED)
				{
					printf("line=%d, In %s(), cannot show the usage information for globalMetrics, error!\n", __LINE__, __func__);
					return FAILED;
				}

				return FAILED;
			}else
			{
				if(argv[i][0]=='-')
				{
					printf("%s : unknown argument\n", argv[i]);
					printf("Please use -h or -help for more information.\n");
				}else
				{
					printf("%s : invalid argument\n", argv[i]);
					printf("Please use -h or -help for more information.\n");
				}

				return FAILED;
			}
		}

		// check whether the parameter is valid
		if(operationMode!=OPERATION_MODE_ALL && operationMode!=OPERATION_MODE_MERGE && operationMode!=OPERATION_MODE_METRICS && operationMode!=OPERATION_MODE_MISASS)
		{
			printf("Exception: please specify the correct operations: all, merge, metrics or misass\n");
			return FAILED;
		}else if(operationMode==OPERATION_MODE_ALL)
		{
			if(strlen(queryFilePara)==0)
			{
				printf("Exception: please specify the query file\n");
				return FAILED;
			}

			if(strlen(subjectsFilePara)==0)
			{
				printf("Exception: please specify the subject configuration file\n");
				return FAILED;
			}
			if(readFileNumPara<=0)
			{
				printf("Exception: please specify the correct reads file\n");
				return FAILED;
			}
		}else if(operationMode==OPERATION_MODE_MERGE)
		{
			if(strlen(subjectsFilePara)==0)
			{
				printf("Exception: please specify the subject configuration file\n");
				return FAILED;
			}
		}else if(operationMode==OPERATION_MODE_METRICS)
		{
			if(strlen(queryFilePara)==0)
			{
				printf("Exception: please specify the query file\n");
				return FAILED;
			}
			if(strlen(subjectsFilePara)==0)
			{
				printf("Exception: please specify the subject configuration file\n");
				return FAILED;
			}
		}else if(operationMode==OPERATION_MODE_MISASS)
		{
			if(strlen(queryFilePara)==0)
			{
				printf("Exception: please specify the query file\n");
				return FAILED;
			}
			if(strlen(subjectsFilePara)==0)
			{
				printf("Exception: please specify the subject configuration file\n");
				return FAILED;
			}
			if(readFileNumPara<=0)
			{
				printf("Exception: please specify the correct reads file\n");
				return FAILED;
			}
		}

		// begin to do the job
		if(computeGlobalMetrics(operationMode, outputDirPara, queryFilePara, subjectsFilePara, minQueryLenPara, minIdentityPercentPara, pairedModePara, readFilesPara, readFileNumPara, threadNumPara, indelSizeThresPara)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the metrics, error!\n", __LINE__, __func__);
			return ERROR;
		}
	}

	return SUCCESSFUL;
}

/**
 * Show the usage information.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short showUsageInfo()
{
	// version information
	printf("misFinder: %s\n", MISFINDER_VERSION_STR);
	printf("Released : %s\n", MISFINDER_RELEASE_DATE_STR);

	printf("\nUsage: mf <command> [option]\n");
	printf("    merge       merge multiple subjects into a multi-fasta format file\n");
	printf("    metrics     compute the metrics\n");
	printf("    misass      compute mis-assemblies\n");
	printf("    all         do all the above in turn\n");

	printf("\nPROGRAM OPTIONS:\n");
	printf("  1) merge -- merge subjects:\n");
	printf("    -s <FILE>          Subjects configuration file. It is required for \n"
		   "                       subjects merge.\n");
	printf("    -o <STR>\n");
	printf("    -out <STR>         Output directory for the output files. Default is \"./\"\n");
	printf("    -h\n");
	printf("    -help              Show help information.\n");

	printf("  2) metrics -- compute the metrics:\n");
	printf("    -q <FILE>          Query file in fasta format.\n");
	printf("    -s <FILE>          Subjects configuration file. It is required for \n"
		   "                       subjects merge.\n");
	printf("    -m <INT>           The minimal query length. Default is 100.\n");
	printf("    -pt <FLOAT>        The minimal identity percentage for matched queries and \n"
		   "                       matched segments. Default is 0.95.\n");
	printf("    -t <INT>           The number of threads for the alignment between queries \n"
		   "                       and subjects. Default is the number of CPU cores.\n");
	printf("    -o <STR>\n");
	printf("    -out <STR>         Output directory for the output files. Default is \"./\"\n");
	printf("    -h\n");
	printf("    -help              Show help information.\n");

	printf("  3) misass -- compute mis-assemblies:\n");
	printf("    -q <FILE>          Query file in fasta format.\n");
	printf("    -s <FILE>          Subjects configuration file. It is required for \n"
		   "                       subjects merge.\n");
	printf("    -p <INT>           Paired end mode for read files.\n"
		   "                       0 - do not treat reads as paired (default);\n"
		   "                       1 - reads are paired with the first read in the first\n"
		   "                       file, and the second read in the second file;\n"
		   "                       2 - reads are paired and interleaved within a single\n"
		   "                       file.\n");
	printf("    -f <FILES>         Read files. It is necessary to be specified for commands\n"
		   "                       'all' and 'misass'.\n");
	printf("    -i <INT>           Minimal indel size. Default is 5 bp.\n");
	printf("    -o <STR>\n");
	printf("    -out <STR>         Output directory for the output files. Default is \"./\"\n");
	printf("    -h\n");
	printf("    -help              Show help information.\n");

	printf("\nREPORT BUGS:\n");
	printf("    zhuxiao.hit@gmail.com\n");

	return SUCCESSFUL;
}
