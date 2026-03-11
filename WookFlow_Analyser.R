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

#------------------------------------*
# The Main Function to do everything *
#------------------------------------*
main_function = function(
        sample_table_path,
        peptide_table_path,
        hits_table_path,
        counts_table_path,
        min_pos_percent,
        min_pos_count,
        comprehensive_output
){
    #-----------------------------------------------------------#
    # Stage 0 - Reading files and generate required data tables #
    #-----------------------------------------------------------#
    
    # If check doesn't complain then the file exists
    # and we can continue to extract sample table data
    sample_table = check_and_read_file(sample_table_path)
    
    sample_prefix = sample_table %>% 
        filter(control_status == 'empirical') %>% 
        select(sample_source) %>% distinct()
    
    sample_repeat_list = repeat_finder(sample_table,'empirical')
    
    # check if there are any samples that have < 2 replicates
    insufficient_replicates = sample_repeat_list %>% filter(n < 2)
    if(nrow(insufficient_replicates) > 0){
        message(paste("These samples did not have at least 2 repeats and have been excluded:"))
        message(paste(insufficient_replicates$sample_source,"\n"))
        message(paste('----------------------------------------',"\n"))
    }
    sample_repeat_list = sample_repeat_list %>% filter(n >= 2)
    beads_only_repeat = repeat_finder(sample_table, 'beads_only')
    
    # extract peptide table data
    peptide_table = check_and_read_file(peptide_table_path)
    
    # extract hits table data
    hits_table = check_and_read_file(hits_table_path)
    
    # extract counts table data
    counts_table = check_and_read_file(counts_table_path)
    
    
    if(!is.integer(min_pos_count) | min_pos_count < 0){
        message(paste("Non-integer provided or value is less than 0 in min_pos_count:",min_pos_count))
        message(paste("Reverting to default value. min_pos_count = 15"))
        min_pos_count = as.integer(15)
    }
    
    if(!is.double(min_pos_percent) | min_pos_percent < 0){
        message(paste("Non-double provided or value is less than 0 in min_pos_percent:",min_pos_percent))
        message(paste("Reverting to default value. min_pos_percent = 0.03"))
        min_pos_percent = as.double(0.03)
    }
    
    # we calculate which samples have multi-hits for the peptides
    multi_hits = generate_multi_hit_table(sample_prefix, hits_table)
    
    #------------------------------------------------------------#
    # Stage 1A - Strict Filtering & generating strict agreements #
    #------------------------------------------------------------#
    # prefix s_ denotes strict
    # suffix _pep dataframe contains per peptide information
    # suffix _spe dataframe contains per species information
    
    s_list = hits_filtering(multi_hits,sample_repeat_list,mode = 'strict')
    
    # this gives us the hit peptides that passed strict filtering conditions
    s_hits_pep = s_list[[1]]
    # grouped by species and then summed up number of hits per species
    s_pos_agreements_spe = s_list[[2]]
    
    
    # now we want to retreive all types of agreements
    s_agreement_list = calculate_strict_agreements(
        multi_hits,
        sample_repeat_list,
        include_disagreement = FALSE # Switch to TRUE if you want disagreement table.
    )
    # This table gives peptides that had agreement for positives and negatives.
    # NOTE: this is different to s_hits_pep.
    #   1)  It's about agreements, not hits.
    #   2)  it consider replicates agreeing a peptide has found to NOT BE a hit too.
    s_agreements_pep = s_agreement_list[[1]]
    
    # Only runs if disagreement parameter was TRUE in
    # calculate_strict_agreement function
    if(length(s_agreement_list) == 2){
        # the opposite of s_agreements_pep
        s_disagreements_pep = s_agreement_list[[2]]
    }
    # create a per species tables of how many agreements were there
    s_agreements_spe = s_agreements_pep %>%
        select(Species, sample_prefix$sample_source) %>%
        group_by(Species) %>%
        summarise(across(everything(), sum))
    
    # calculate how many total peptides represent each species
    max_peptide_per_species = s_agreements_pep %>% select(Species) %>% count(Species)
    colnames(max_peptide_per_species)[2] = "Max_Possible_Agreement"
    
    # create a table so we can normalise the agreements by max possible agreement
    # easier to see with % agreed in running the experiment
    div_matrix = matrix(
        max_peptide_per_species$Max_Possible_Agreement,
        nrow = s_agreements_spe %>% select(sample_prefix$sample_source) %>% nrow(),
        ncol = s_agreements_spe %>% select(sample_prefix$sample_source) %>% ncol()
    )
    
    # adding the max possible peptide next to s_agreement_spe
    s_agreements_spe = inner_join(max_peptide_per_species,s_agreements_spe, by = "Species")
    
    # normalising (omitting the non-numeric columns first) and call it a s_agreements_spe_pct
    s_agreements_spe_pct = as_tibble(
        s_agreements_spe %>% select(sample_prefix$sample_source) / div_matrix
    )
    
    # we'll re-add the species and max possible agreement column
    s_agreements_spe_pct = add_column(
        s_agreements_spe_pct,
        s_agreements_spe[,1:2],
        .before = 1
    )
    
    remove(s_list,s_agreement_list)
    
    #--------------------------------------------------------------#
    # Stage 1B - Lenient Filtering & generating lenient agreements #
    #--------------------------------------------------------------#
    
    # TODO Space to later implement Lenient Filtering if desired
    
    # l_list = hits_filtering(multi_hits,sample_repeat_list,mode = 'lenient')
    # l_hits_pep = l_list[[1]]
    # l_combined_pos_agreements_spe = l_list[[2]]
    
    #-----------------------------------------------------------------#
    # Stage 2A - Selecting peptides and species using agreement ratio # 
    #            Strict Filtering                                     #
    #-----------------------------------------------------------------#
    
    # Calculate the percentage of positive agreements out of all agreements
    # The reason why we used "all-agreements" rather than total number of peptide
    # representatives of each species was because disagreements between replicates yields
    # unreliable information under Strict Filtering conditions. A peptide that's found to
    # be a hit in a subset of the replicates is neither a confirmation of a "found hit" or
    # a "found non-hit". We exclude this gray-zone of a finding.
    s_pos_agreements_pct = as_tibble(
        s_pos_agreements_spe %>% select(all_of(sample_repeat_list$sample_source)) / 
            s_agreements_spe %>% select(all_of(sample_repeat_list$sample_source))
    ) 
    
    # re-adding species names to the dataframe
    s_pos_agreements_pct = add_column(
        s_pos_agreements_pct, 
        s_pos_agreements_spe %>% select(Species), 
        .before = 1
    )
    
    #-----------------------------------------------------------------#
    # Stage 2B - Selecting peptides and species using agreement ratio # 
    #            Lenient Filtering                                     #
    #-----------------------------------------------------------------#
    # TODO Space to later implement Lenient Filtering if desired
    
    #----------------------------------------------------#
    # Stage 3A - Cut off Reporting for Strict Filtering  #
    #----------------------------------------------------#
    
    # Pivot the both pct and spe tables so we can filter then join them
    # into one df for further filtering and plotting.
    # Note: ct -> count, pct -> percent
    s_ct_vs_pct = inner_join(
        make_species_pivot_tables(s_pos_agreements_spe, 'Pos_Ct_Agreement'),
        make_species_pivot_tables(s_pos_agreements_pct, 'Pos_Pct_Agreement'),
        by = c("Species","Samples"))
    
    cutoff_done = prioritise_species(s_ct_vs_pct,min_pos_percent,min_pos_count)
    h_priority = cutoff_done[[1]] # Good probability that these species were detected
    m_priority = cutoff_done[[2]] # a bit dodgy, investigate with caution
    l_priority = cutoff_done[[3]] # we think these are absolutely dodgy
    nan_entries = cutoff_done[[4]] # ratio was NaN, indicating dividing by 0.
    non_high = bind_rows(m_priority,l_priority)
    
    # We filter viruses that have been found in multiple samples (at least 2). extract which samples
    # to which species retaining both sample name and species name
    h_priority_info = find_common_species(h_priority,min_shared = 2)  # increase min_shared to make it more stringent
    h_priority_shared = h_priority_info[[1]]
    h_priority_samples_species = h_priority_info[[2]]
    
    # detailing which species as well as which samples with what species
    # did not make it as a high priority
    exclusion_info_list = get_excluded_species_info(non_high,h_priority_shared,h_priority_samples_species)
    excluded_species = exclusion_info_list[[1]]
    excluded_species_samples = exclusion_info_list[[2]]
    num_notFound_samples = exclusion_info_list[[3]]
    
    remove(cutoff_done,exclusion_info_list)
    
    # Generate the plots here if needed. We choose everything above 0 so
    # we can actually see the min_pos_percent and and min_pos_count cut offs
    myGraph1 = make_PosAgreementPlot(s_ct_vs_pct, min_pos_count, min_pos_percent, h_priority_shared)
    
    
    #----------------------------------------------------#
    # Stage 3B - Cut off Reporting for Lenient Filtering #
    #----------------------------------------------------#
    # TODO Space to later implement Lenient Filtering if desired
    
    
    #----------------------------------------------------#
    # Stage 4A - Strict Filtering peptide level analysis #
    #----------------------------------------------------#
    
    # Here we're retrieving all the peptides from species that have passed the
    # cutoffs in each sample. In other words, all the peptides that were:
    #   1)  strict positive agreements
    #   2)  greater than min_pos_count & min_pos_percent
    s_pos_peptide_dt = s_hits_pep %>% 
        pivot_longer(
            -c(peptide_id,original_id,Species),
            names_to = "Samples", 
            values_to = 'Strict_Pos_Agreement'
        ) %>% 
        filter(
            Species %in% unique(h_priority_samples_species$Species) & 
                Samples %in% unique(h_priority_samples_species$Samples) & 
                Strict_Pos_Agreement > 0
        )
    
    # Retrieving peptide information of peptides from s_pos_peptide_dt
    s_pos_peptide_info = merge(
        peptide_table %>% select(
            peptide_id,
            original_id,
            Prot,
            Prot_Start,
            Prot_End,
            UniProtEntry,
            taxon_ids,
            common_taxonomies,
            Species
        ),
        s_pos_peptide_dt %>% select(-Strict_Pos_Agreement,-Species),
        by = c("peptide_id","original_id")
    )
    
    # generating the number of times a peptide has been detected in the samples
    s_pos_peptide_occurrence = as_tibble(
        merge(
            peptide_table %>% select(
                peptide_id,
                original_id,
                Prot,
                Prot_Start,
                Prot_End,
                UniProtEntry,
                taxon_ids,
                common_taxonomies,
                Species
            ),
            s_pos_peptide_dt %>% select(peptide_id) %>% count(peptide_id)
        )
    )
    
    #-----------------------------------------------------#
    # Stage 4B - Lenient Filtering peptide level analysis #
    #-----------------------------------------------------#
    # TODO Space to later implement Lenient Filtering if desired
    
    
    #---------------------------------------------#
    # Stage 5A - Strict Filtering counts analysis #
    #---------------------------------------------#
    
    normalised_counts = rpm_normalise(counts_table,4)
    hits_names = names(hits_table %>% select(-peptide_id,-original_id,-Species))
    
    # Now to find the ones that were hits and retain their read count
    normalised_hits_counts = (normalised_counts %>% select(all_of(hits_names))) * (hits_table %>%  select(all_of(hits_names)))
    normalised_hits_counts = attach_pep_info(counts_table,normalised_hits_counts)
    
    # 1)    For each sample, we'll look at all replicates for every individual peptide and 
    #       find peptides that are > 0 in `normalised_hits_counts`. 
    # 2)    If all replicates are > 0  we'll take the mean of hits counts of that peptide 
    #       using values from all replicates.
    # 3)    If any of the replicate is 0, we'll return a 0 because not all replicates have 
    #       agreed they found sufficient reads to be a hit.
    # 4)    This will be our `strict_pos_mean_hits_counts` table.
    s_pos_mean_hits_counts = bind_cols(
        lapply(
            sample_prefix$sample_source,
            function(i){
                b = normalised_hits_counts %>% 
                    select(starts_with(i)) %>% 
                    mutate(
                        value = ifelse( 
                            # look for any numeric columns which are 0 and sum them. If at least
                            # one column is a 0, then rowSums would > 0.
                            rowSums(across(where(is.numeric)) == 0) > 0, 
                            0, # If rowSums > 0, we'll put it as zero
                            rowMeans(across(where(is.numeric)),na.rm = TRUE) # else get rowMean from the 3 triplicates
                        )
                    ) %>% 
                    select(value)
                colnames(b) = i
                return(b)
            }
        )
    )
    # adding back some peptide info
    s_pos_mean_hits_counts = attach_pep_info(counts_table,s_pos_mean_hits_counts)
    
    #-------------------------------#
    # FINAL Stage - Writing outputs #
    #-------------------------------#
    
    # s_pos results (Species level info)
    output_csv(h_priority_shared,"HighPriority_Species_Shared.csv",DEFAULT_PREFIX,DEFAULT_OUTPUTDIR)
    output_csv(h_priority_samples_species,'HighPriority_Sample_to_Species.csv',DEFAULT_PREFIX,DEFAULT_OUTPUTDIR)
    
    # s_pos results (Peptide Level info)
    output_csv(s_pos_peptide_info, "HighPriority_Species_Peptide_Info.csv",DEFAULT_PREFIX,DEFAULT_OUTPUTDIR)
    
    # s_pos_hits_counts results
    output_csv(s_pos_mean_hits_counts,"Strict_Positive_Mean_Hits_Counts.csv",DEFAULT_PREFIX,DEFAULT_OUTPUTDIR)
    
    if(!is.null(comprehensive_output)){
        # raw l_priority
        output_csv(h_priority,"RawHighPriority_List.csv",DEFAULT_PREFIX,DEFAULT_OUTPUTDIR)
        # raw m_priority
        output_csv(m_priority,"RawMediumPriority_List.csv",DEFAULT_PREFIX,DEFAULT_OUTPUTDIR)
        # raw l_priority
        output_csv(l_priority,"RawLowPriority_List.csv",DEFAULT_PREFIX,DEFAULT_OUTPUTDIR)
        # Nan agreements
        output_csv(nan_entries,"Nan_Agreements.csv.csv",DEFAULT_PREFIX,DEFAULT_OUTPUTDIR)
        # s_pos results for counting peptide identification per sample
        output_csv(s_pos_peptide_occurrence,"HighPriority_Species_Peptide_Count_Per_Sample.csv",DEFAULT_PREFIX,DEFAULT_OUTPUTDIR)
    }
    
    # Pos Agreement Plot
    ggsave(
        file.path(
            DEFAULT_OUTPUTDIR,
            paste0(DEFAULT_PREFIX,"_","PostiveAgreement_Plot.pdf")
        ),
        plot = myGraph1 + theme(aspect.ratio = 0.7),
        device = cairo_pdf,
        width = 1920,
        height = 1200,
        units = "px",
        scale = 3
    ) 
}

##################################
# Defining user input parameters #
##################################

# Options for the code to run
option_list = list(
    # Will probably use this option later. But First make the params work
    # make_option(c("-p", "--param_file"), action="store", default=NULL, type='character',
    #             help="parameter file to configure SKAT-O run. Must be in comma separated format."),  
    # make_option(
    #     c("--cores"),
    #     action = "store",
    #     # action = "store_true",
    #     default = 4,
    #     help = "Number of cores to use. [Default = 4].",
    #     type = "integer"
    # ),
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

