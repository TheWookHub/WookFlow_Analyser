#!/usr/bin/env Rscript

## Libraries required to run
suppressWarnings(
    suppressPackageStartupMessages(
        {
            library(optparse)
            library(tidyverse)
            # library(data.table)
            # library(checkmate)
        }
    )
)

# Load all the supporting functions to be used
# in main function
source("WookFlow_Functions.R")

###################################
# Define Global Variables for use #
###################################

tool_name = 'WookFlow Analyser'
tool_ver = 'v0.1'

# THis is a list of options so that we can figure out
# later if we're missing
option_names = c(
    "sample_table","peptide_table","hits_table",
    "counts_table", "min_pos_agreement_count",
    "min_pos_agreement_percent"
)

# path to write files to
DEFAULT_OUTPUTDIR = getwd()

# default_prefix is an empty string
DEFAULT_PREFIX = NULL

# Define info message just to state some info for this tool
info_msg = paste0(
    "***********************************************************************************\n",
    "HI! HOW ARE YOU?\n",
    "I'M FINE THANK YOU, AND YOU?\n",
    "Welcome and thank you for using ", tool_name, " - ",tool_ver, "\n",
    "This tool:\n",
    "  - Takes Wookflow output and conducts post pipeline report about the data.\n",
    "  - Tool mainly focuses on VirScan. HuScan to be implemented over time.\n",
    "  - Lenient filtering to be done when Ki Wook can be in office for 10 days straight.\n",
    "------------------------------------------\n",
    "  - Tool constructed by Preston L.\n",
    "***********************************************************************************\n"
)

##################################
# Defining user input parameters #
##################################

# Options for the code to run
option_list = list(
    make_option(
        c("--sample_table"),
        action = "store",
        help = "Sample Table file (.csv) containing information on samples of PhIPSeq Output.",
        type = "character"
    ),
    make_option(
        c("--peptide_table"),
        action = "store",
        help = "Peptide table file (.csv) containing information peptide info such as species, protein sequence etc.",
        type = "character"
    ),
    make_option(
        c("--hits_table"),
        action = "store",
        help = "Hits file (.csv) containing information of peptide hits for each peptide for all samples.",
        type = "character"
    ),
    make_option(
        c("--counts_table"),
        action = "store",
        help = "Counts table (.csv) containing information of peptide counts for each peptide for all samples.",
        type = "character"
    ),
    make_option(
        c("--min_pos_agreement_count"),
        action = "store",
        default = 15,
        help = paste(
            "Minimum positive agreement counts for a species to be acknowledged as identified.",
            "[Default: 15]"
        ),
        type = "integer"
    ),
    make_option(
        c("--min_pos_agreement_percent"),
        action = "store",
        default = 0.03,
        help = paste(
            "Minimum positive agreement percent (out of total agreements) for a",
            "species to be acknowledged as identified. [Default: 0.03]"
        ),
        type = "double"
    ),
    make_option(
        c("--file_prefix"),
        action = "store",
        help = paste(
            "An optional user defined prefix character string to add to output",
            "file names."
        ),
        type = 'character'
    ),
    make_option(
        c("--outpath"),
        action = "store",
        help = paste(
            "An optional user defined output path to write output",
            "file to. [Default = Current Work Directory]"
        ),
        type = 'character'
    ),
    make_option(
        c("--comprehensive"),
        action = "store_true",
        help = "If flagged, the program will generate additional outputs. E.g. Mid and Low priority findings."
    )
)

# Here we parse the params from option list
option_object = OptionParser(
    option_list = option_list,
    description = info_msg
)

usr_params = parse_args(option_object)

####################################
# Main Process where stuff happens #
####################################

# check required params, otherwise print help msg
missing = setdiff(option_names, names(usr_params))

if(length(missing) > 0){
    print_help(option_object)
    message("-------------------------------------------")
    message("-------------------------------------------")
    message(paste("Missing required parameter(s):",paste(missing, collapse = ', ')))
    message("Program will terminate.")
}else{
    # Print welcome msg
    message(info_msg)
    
    # Check and change DEFAULT_OUTPUTDIR and/or DEFAULT_PREFIX
    # if provided and valid
    if(!is.null(usr_params$outpath)){
        if(dir.exists(usr_params$outpath)){
            paste0("Modfied file name DEFAULT_OUTPUTDIR: ", DEFAULT_OUTPUTDIR)
            DEFAULT_OUTPUTDIR = usr_params$outpath
        }else{
            message("-------------------------------------------")
            message(
                paste0("Outpath not found: ", usr_params$outpath)
            )
            message(
                paste0("Using DEFAULT_OUTPUTDIR: ", DEFAULT_OUTPUTDIR)
            )
            message("-------------------------------------------")
        }
    }
    
    if(!is.null(usr_params$file_prefix)){
        DEFAULT_PREFIX = usr_params$file_prefix
        message(
            paste0("Modfied file name DEFAULT_PREFIX: ", DEFAULT_PREFIX)
        )
        message("-------------------------------------------")
    }
    
    # Code main process here
    main_function(
        sample_table_path = usr_params$sample_table,
        peptide_table_path = usr_params$peptide_table,
        hits_table_path = usr_params$hits_table,
        counts_table_path = usr_params$counts_table,
        min_pos_percent = usr_params$min_pos_agreement_percent,
        min_pos_count = usr_params$min_pos_agreement_count,
        comprehensive_output = usr_params$comprehensive
    )
    
    message(
        paste0("Done!")
    )
}

###############
# END OF CODE #
###############

