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
 *  @return:
 *  	If succeeds, return SUCCESSFUL;
 *  	if the job is executed and failed, return ERROR;
 *  	otherwise, return FAILED.
 */
short parseCommandParasAndExe(int argc, char **argv)
{
	int32_t i, operationMode;
	char configFilePara[256];
	char outputDirPara[256];
	int32_t minQueryLenPara, threadNumPara;
	double minIdentityPercentPara;
	int32_t indelSizeThresPara;

#if __WORDSIZE!=64
	printf("Please compile and run the tool on x86_64 Linux system.\n");
	return FAILED;
#endif

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
		configFilePara[0] = '\0';
		outputDirPara[0] = '\0';
		minQueryLenPara = 0;
		minIdentityPercentPara = 0;
		threadNumPara = 0;
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

			if(strcmp(argv[i], "-o")==0 || strcmp(argv[i], "-out")==0)
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
			}
			else if(strcmp(argv[i], "-conf")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						strcpy(configFilePara, argv[i+1]);
					}else
					{
						printf("Exception: please specify the correct configuration file.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the configuration file.\n");
					return FAILED;
				}
				i += 2;
			}
			else if(strcmp(argv[i], "-t")==0)
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
		}else if(strlen(configFilePara)==0)
		{
			printf("Exception: please specify the correct configuration file.\n");
			return FAILED;
		}

		// begin to do the job
		if(computeGlobalMetrics(operationMode, outputDirPara, configFilePara, minQueryLenPara, minIdentityPercentPara, threadNumPara, indelSizeThresPara)==FAILED)
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
	printf("    metrics     Compute the assembly metrics\n");
	printf("    misass      Compute mis-assemblies\n");
	printf("    all         Do all the above in turn\n");

	printf("\nPROGRAM OPTIONS:\n");
	printf("  1) metrics -- compute the metrics:\n");
	printf("    -conf <FILE>       Configuration file. It is required.\n");
	printf("    -m <INT>           The minimal query length. Default is 100.\n");
	printf("    -pt <FLOAT>        The minimal identity percentage for matched queries and \n"
		   "                       matched segments. Default is 0.95.\n");
	printf("    -t <INT>           The number of threads for the alignment between queries \n"
		   "                       and subjects. Default is the number of CPU cores.\n");
	printf("    -o <STR>\n");
	printf("    -out <STR>         Output directory for the output files. Default is \"./\"\n");
	printf("    -h\n");
	printf("    -help              Show help information.\n");

	printf("  2) misass -- compute mis-assemblies:\n");
	printf("    -conf <FILE>       Configuration file. It is required.\n");
	printf("    -i <INT>           Minimal indel size. Default is 5 bp.\n");
	printf("    -t <INT>           The number of threads for reads alignment. Default is\n"
		   "                       the number of CPU cores.\n");
	printf("    -o <STR>\n");
	printf("    -out <STR>         Output directory for the output files. Default is \"./\"\n");
	printf("    -h\n");
	printf("    -help              Show help information.\n");

	printf("\nREPORT BUGS:\n");
	printf("    zhuxiao.hit@gmail.com\n");

	return SUCCESSFUL;
}
