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
        # raw h_priority
        output_csv(h_priority,"RawHighPriority_List.csv",DEFAULT_PREFIX,DEFAULT_OUTPUTDIR)
        # raw m_priority
        output_csv(m_priority,"RawMediumPriority_List.csv",DEFAULT_PREFIX,DEFAULT_OUTPUTDIR)
        # raw l_priority
        output_csv(l_priority,"RawLowPriority_List.csv",DEFAULT_PREFIX,DEFAULT_OUTPUTDIR)
        # Nan agreements
        output_csv(nan_entries,"Nan_Agreements.csv",DEFAULT_PREFIX,DEFAULT_OUTPUTDIR)
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

#################################
# Define supporting function(s) #
#################################

# if file doesn't exist, we stop and produce error message
# if it does exist, then we read in the file

check_and_read_file = function(filePath){
    if(!file.exists(filePath)){
        stop(paste("File does not exist:",filePath))
    }
    infile = read_csv(filePath, progress = FALSE, show_col_types = FALSE)
    return(infile)
}

# lazy function to add peptide info to dataframes
attach_pep_info = function(info_df,to_attach_df){
    if(nrow(info_df) != nrow(to_attach_df)){
        message("Dataframes do not have matching rows for adding column!")   
    }else{
        donaldTrump = add_column(
            as_tibble(to_attach_df),
            .before = 1, 
            peptide_id = info_df$peptide_id,
            original_id = info_df$original_id,
            Species = info_df$Species
        )
        return(donaldTrump)
    }
}

# performs reads per million normalisation
rpm_normalise = function(dt, start_col){
    scaling_factor <- apply(dt[,start_col:ncol(dt)], 2, sum)/1000000
    dt_RPM = as_tibble(mapply('/', dt[,start_col:ncol(dt)], scaling_factor))
    dt_RPM = add_column(dt_RPM,dt[,1:(start_col-1)],.before = 1)
    return(dt_RPM)
}

# how to figure how many repeats there are
# assumes the column sample_source exists
repeat_finder = function(sample_table,type){
    rpt_list = sample_table %>% 
        filter(control_status %in% c(type)) %>%
        select(sample_source) %>% 
        count(sample_source)
    return(rpt_list)
}

# function to generate a table that records how many hits among
# repeats per peptide
generate_multi_hit_table = function(prefix_dt, hits_dt){
    results = bind_cols(
        lapply(
            prefix_dt$sample_source,
            function(i){
                sample_hits = hits_dt %>% select(starts_with(i))
                x = data.frame(value = rowSums(sample_hits))
                colnames(x) = i
                return(x)
            }
        )
    )
    # adding back some peptide info
    results = attach_pep_info(hits_dt,results)
    return(results)
}

# Calculate strict and lenient filtering
hits_filtering = function(multi_hits_df, replicate_info,mode){
    # Pre-defined as NULL because nothing exists instead of NA
    # which signifies something should be here but is missing.
    filtered_hits = NULL 
    filtered_species = NULL
    
    sample_names = replicate_info$sample_source
    relevant_data = multi_hits_df %>% select(all_of(sample_names))
    mode = tolower(mode)
    
    if(!mode %in% c('strict','lenient')){
        stop(paste("No such mode for hits_filtering function:",mode))
    }
    
    if(mode == 'strict'){
        filtered_hits = bind_cols(
            lapply(
                sample_names,
                function(i){
                    count_condition = replicate_info %>% 
                        filter(sample_source == i) %>% 
                        select(n) %>% 
                        as.integer()
                    x_hits = ifelse(multi_hits_df[,c(i)] >= count_condition, 1,0)
                    return(x_hits)
                }
            )
        )
    }else if(mode == 'lenient'){
        filtered_hits = bind_cols(
            lapply(
                sample_names,
                function(i){
                    count_condition = replicate_info %>% 
                        filter(sample_source == i) %>% 
                        select(n) %>% 
                        as.integer()
                    count_condition = ceiling(count_condition / 2)
                    x_hits = ifelse(multi_hits_df[,c(i)] >= count_condition, 1,0)
                    return(x_hits)
                }
            )
        )
    }
    filtered_hits = attach_pep_info(multi_hits_df, filtered_hits)
    filtered_species = filtered_hits %>%
        select(Species,all_of(sample_names)) %>%
        group_by(Species) %>%
        summarise(across(everything(),sum))
    return(list(filtered_hits,filtered_species))
}

# Calculate agreements and disagreements for strict filtering 
#   -   Agreement is defined as all replicates agree in either 
#       yes to a peptide as a hit or no it is not a hit.
#   -   Disagreement is defined as subset of replicates agree it is a hit
#       and the remainder replicates do not.
# Note: lenient partial agreements is the same as strict disagreements (but we don't really need this here).
calculate_strict_agreements = function(multi_hits_df,replicate_info, include_disagreement = FALSE){
    sample_names = replicate_info$sample_source
    # Calculate agreements
    s_agreements = bind_cols(
        lapply(
            sample_names,
            function(i){
                count_condition = replicate_info %>% 
                    filter(sample_source == i) %>% 
                    select(n) %>% 
                    as.integer()
                x_agreements = ifelse(
                    multi_hits_df[,c(i)] >= count_condition | multi_hits_df[,c(i)] == 0, 1,0
                )
                return(x_agreements)
            }
        )
    )
    
    if(include_disagreement){
        # Quick maffs to make it into disagreements
        s_disagreements = s_agreements * -1 + 1
        s_disagreements = attach_pep_info(multi_hits_df, s_disagreements)
        s_agreements = attach_pep_info(multi_hits_df, s_agreements)
        return(list(s_agreements, s_disagreements))
    }else{
        s_agreements = attach_pep_info(multi_hits_df, s_agreements)
        return(list(s_agreements))
    }
}

# Turn any species table into a long table for plotting purposes
make_species_pivot_tables = function(speciesPivot_df, value_names){
    pivot_table = speciesPivot_df %>% 
        pivot_longer(
            -Species, 
            names_to = "Samples", 
            values_to = value_names
        )
    return(pivot_table)
}

# takes a count vs pct dataframe, based on min count and min ratio cutoff
# return subsets of the dataframe that shows high/mid/low priority species.
prioritise_species = function(ct_v_pct_df, min_percent, min_count){
    # record any nan entires
    nan_entries = ct_v_pct_df %>% filter(is.nan(Pos_Pct_Agreement))
    
    # removes any rows where pct results in a NaN due to divide by 0
    # i.e. no agreements were found at all, so 0 / 0 gives NaN.
    # Usually happens with viruses with very low number of representative peptides
    clean_df = drop_na(ct_v_pct_df)
    
    h_priority = clean_df %>% 
        filter(Pos_Pct_Agreement > min_percent & Pos_Ct_Agreement > min_count)
    m_priority = bind_rows(
        clean_df %>% 
            filter(Pos_Pct_Agreement > min_percent & Pos_Ct_Agreement <= min_count),
        clean_df %>% 
            filter(Pos_Pct_Agreement <= min_percent & Pos_Ct_Agreement > min_count),
    )
    l_priority = clean_df %>% 
        filter(Pos_Pct_Agreement <= min_percent & Pos_Ct_Agreement <= min_count)
    return(list(h_priority,m_priority,l_priority,nan_entries))
}

# returns a dataframe for species that (after cutoff) have been shared by at least
# min_shared number of samples [Default = 2]. Also returns a dataframe that shows
# which viruses were identified with which sample.
find_common_species = function(ct_vs_pct,min_shared = 2){
    # Then we filter viruses that have been found in multiple samples.
    shared_species = ct_vs_pct %>% 
        count(Species) %>% 
        arrange(desc(n)) %>% 
        filter(n >= min_shared)
    
    samples_species_info = unique(
        ct_vs_pct %>% 
            filter(Species %in% shared_species$Species) %>% 
            select(Species,Samples)
    )
    return(list(shared_species, samples_species_info))
}

# retrieve extra information for species and samples that didn't make it
# to the high priority (for comprehensiveness of results)
get_excluded_species_info = function(excluded_dt,included_shared, included_samples_species){
    # extracting species that were not part of the h_priority_shared list    
    excluded_species = excluded_dt %>% 
        filter(!Species %in% included_shared$Species) %>%
        count(Species) %>% 
        arrange(n)
    
    # In this dataframe, the "n" column isn't too important
    # some species will have less than total samples because they
    # had NaN and were excluded due to divide by 0.    
    excluded_species_samples = bind_rows(
        lapply(
            unique(excluded_dt$Samples),
            function(i){
                x = excluded_dt %>% filter(Samples %in% i)
                y = included_samples_species %>% filter(Samples %in% i)
                return(x %>% filter(!Species %in% y$Species))
            }
        )
    )
    # This table shows how many species was not accepted per sample
    num_notFound_samples = excluded_species_samples %>% 
        group_by(Samples) %>% 
        count() %>% 
        arrange(n)
    return(list(excluded_species,excluded_species_samples,num_notFound_samples))    
}

# generate the Postive Agreement Plot. Produces a ggplot object with:
#    -  y-axis being pos. agreement count
#           *   defined as the total counts of peptides of a species collectively agreed
#               to be detected across replicates. 
#    -  x-axis being pos. agreement percentage
#           *   defined as pos. agreement count divided by the total agreements (which 
#               include negatives agreements).
#    -  top 5 [default] species of high_priority highlighted in colour for the plot
make_PosAgreementPlot = function(ct_vs_pct_dt, min_count,min_percent,high_priority_info, top = 5){
    # some basic filtering 
    toPlot_nonZero = ct_vs_pct_dt %>% filter(Pos_Pct_Agreement > 0.00 & Pos_Ct_Agreement > 0)
    # determining a valid top number in case there are less than 5 [default].
    top = ifelse(length(high_priority_info$Species) < top, length(high_priority_info$Species), top)
    topList = high_priority_info$Species[1:top]
    # Here we choose top 5 (arbitrary value) those in h_priority_shared to highlight in the plot.
    # Since we will never know how many h_priority_shared will have, to prevent highlighting
    # too many (loses the purpose), we'll just settle with top 5. 
    toPlot_nonZero_2 = toPlot_nonZero %>% 
        mutate(
            group = if_else(
                Species %in% topList, 
                Species, 
                "Other"
            )
        )
    
    # Give it some fancy colours
    highlight_cols <- setNames(scales::hue_pal()(top), topList)
    highlight_cols <- c(highlight_cols, Other = "lightgray")
    
    # Do the ploties
    myPlot = ggplot(toPlot_nonZero_2, aes(x = Pos_Pct_Agreement, y = Pos_Ct_Agreement, colour = group)) +
        geom_point(
            size = 2
        )+
        scale_color_manual(
            name = "Species",
            values = highlight_cols, 
            breaks = c(topList, "Other")
        ) +
        labs(
            title = "Strict Positive Agreement Plot: Count vs Percent",
            subtitle = paste("Highlighting Top 5 out of",nrow(high_priority_info),"high priority species."),
            x = 'Pos. Agreement Percent', 
            y = 'Pos. Agreement Count'
        ) +
        geom_vline(xintercept = min_percent,colour = "red") +
        annotate(
            "text",
            x = min_percent, 
            y = Inf, 
            label = paste0((min_percent) * 100,"% Pos. Agreement"),
            vjust = -1, 
            hjust = 3, 
            angle = 90
        ) +
        geom_hline(yintercept = min_count,colour = "red",) +
        annotate(
            "text",
            x = Inf, 
            y = min_count, 
            label = paste(min_count, "counts of Pos. Agreement"),
            vjust = -0.5, 
            hjust = 3
        ) +
        theme(legend.position = "bottom",legend.key.size = unit(1,'line'))
    return(myPlot)
} 

# Write output function
output_csv = function(dt,filename,prefix,write_path){
    if(is.null(prefix)){
        write_csv( # raw h_priority
            dt,
            file.path(write_path,filename)
        )    
    }else{
        write_csv( # raw h_priority
            dt,
            file.path(
                write_path,
                paste0(prefix,"_",filename)
            )
        )
    }
}


