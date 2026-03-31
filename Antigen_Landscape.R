# Trying out Bea / Legana code for antigen landscape
# 

library(tidyverse)
library(ampir)

X

calculate_mean_rpk_difference <- function(data, sample_id_col, condition_col, pep_id_col, abundance_col, blastp_data) {
    # Calculate RPK and mean RPK per peptide
    rpk_mean_peptide_hits <- data %>%
        group_by({{ sample_id_col }}) %>% 
        mutate(rpk = {{ abundance_col }} / sum({{ abundance_col }}) * 100000) %>%
        ungroup() %>% 
        group_by({{ condition_col }}, {{ pep_id_col }}) %>% 
        mutate(mean_rpk_per_peptide = mean(rpk)) %>% 
        ungroup()
    
    # Separate cases and controls to calculate mean_rpk per group
    mean_rpk_cc <- rpk_mean_peptide_hits %>%
        pivot_wider(names_from = {{ condition_col }}, values_from = mean_rpk_per_peptide, names_prefix = "mean_rpk_per_pep", values_fill = 0)  %>%
        select({{ pep_id_col }}, mean_rpk_per_pepCase, mean_rpk_per_pepControl) %>%
        distinct() %>%
        filter(mean_rpk_per_pepCase != 0 | mean_rpk_per_pepControl != 0) #remove rows where both case and control are 0

    # Nifty code trick to get cases and controls in the same row per peptide
    mean_rpk_cc_reduced <- mean_rpk_cc %>%
        group_by({{ pep_id_col}}) %>%
        mutate(
            mean_rpk_per_pepControl = if_else(mean_rpk_per_pepControl == 0 & n() > 1, NA_real_, mean_rpk_per_pepControl),
            mean_rpk_per_pepCase = if_else(mean_rpk_per_pepCase == 0 & n() > 1, NA_real_, mean_rpk_per_pepCase)) %>%
        fill(mean_rpk_per_pepControl, mean_rpk_per_pepCase, .direction = "downup") %>%
        distinct() %>%
        ungroup()

    # Join BLAST results to mean RPK calculations for case and control and calculate mean_rpk_difference
    blastp_data %>%
        left_join(mean_rpk_cc_reduced, by = join_by(qaccver == {{ pep_id_col }})) %>%
        group_by(saccver) %>% #protein ID
        mutate(mean_rpk_difference = mean_rpk_per_pepCase - mean_rpk_per_pepControl) %>%
        drop_na() %>%
        ungroup() %>%
        select("seqid" = qaccver, "start" = sstart, "end" = send, mean_rpk_per_pepCase, mean_rpk_per_pepControl, mean_rpk_difference, saccver)
}



test_window_size <- function(start, end, win_size, step_size) {
    if ((end - start + 1) >= win_size) {
        seq(start, end - win_size + 1, by = step_size)
    } else {
        NA
    }
}

calculate_moving_sum <- function(df, value_column, win_size = NULL, step_size = NULL) {
    # define default values
    default_win_size <- 4
    default_step_size <- 1
    
    # use defaults if parameters are NULL
    win_size <- win_size %||% default_win_size
    step_size <- step_size %||% default_step_size
    
    if(missing(value_column)) {
        stop("The `value_column` must be specified!")
    }
    if (win_size == default_win_size && step_size == default_step_size) {
        message("Default values used: win_size = ", default_win_size, ", step_size = ", default_step_size)
    } else if (win_size == default_win_size) {
        message("Default value used for win_size: ", default_win_size)
    } else if (step_size == default_step_size) {
        message("Default value used for step_size: ", default_step_size)
    }
    
    df %>%
        mutate(windows = pmap(list(start, end), test_window_size, win_size = win_size, step_size = step_size)) %>%
        unnest(windows) %>%
        filter(!is.na(windows)) %>%
        mutate(window_start = windows,
               window_end = windows + win_size - 1) %>%
        rowwise() %>%
        mutate(moving_sum = sum(df %>% filter(start <= window_start, 
                                              end >= window_end) %>% 
                                    pull({{ value_column }}))) %>%
        ungroup() %>%
        select(-windows)
}


# Plot the moving sum in cases and controls 
ms_plot_clean <- function(moving_sum_dataframe){
    #plot without x axis label or legend position for easier plot stacking
    moving_sum_dataframe %>% 
        distinct(window_start, .keep_all = TRUE) %>% #remove duplicate windows to avoid overstacking
        mutate(Condition = if_else(moving_sum > 0, "Case", "Control")) %>% 
        ggplot(aes(x = (window_start + window_end) / 2, y = moving_sum, fill = Condition)) +
        geom_bar(stat = "identity") +
        labs(x = "", fill = "", y = "") +
        theme_minimal(base_size = 14) +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank()) +
        scale_fill_manual(values = c("Case" = "#d73027", "Control" = "#4575b4"), labels = c("Case", "Control")) +
        theme(legend.position = "none")
}



# Function to read BLAST output
read_blast <- function(input_blast_result) {
    read_tsv(input_blast_result, col_names = c("qaccver", "saccver", "pident", "nident", "length", "evalue", "bitscore", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "qseq", "sseq", "ppos", "stitle", "frames"))
}

# function to read enterovirus polyprotein metadata obtained from UniProt
# note this function is optimised for coxsackievirus B1. Protein names and positions in other ev species MAY differ
read_ev_polyprotein_uniprot_metadata <- function(path_to_enterovirus_metadata.tsv) {
    `%notin%` <- Negate(`%in%`)
    
    overlapping_ev_proteins <- c("P1",
                                 "Genome polyprotein",
                                 "Capsid protein VP0",
                                 "P2",
                                 "P3",
                                 "Protein 3A",
                                 "Viral protein genome-linked",
                                 "Protein 3CD")
    
    ev_proteins <- c("VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3AB", "3C", "3D")
    
    read_tsv(path_to_enterovirus_metadata.tsv) %>%
        separate_longer_delim(Chain, delim = "; CHAIN") %>% 
        select(-c(Entry, Reviewed)) %>% 
        mutate(
            start = as.numeric(str_extract(Chain, "(?<=^|\\s)\\d+")),
            end = as.numeric(str_extract(Chain, "(?<=\\.\\.)\\d+")),
            note = str_extract(Chain, '(?<=/note=\").+?(?=\")'),
            id = str_extract(Chain, '(?<=/id=\").+?(?=\")')
        ) %>% 
        filter(note %notin% overlapping_ev_proteins) %>% # filter out overlapping proteins
        mutate(ev_proteins = str_extract(note, paste(ev_proteins, collapse = "|")),
               ev_proteins = if_else(is.na(ev_proteins), "3D", ev_proteins), #tidy up names in bulk
               start = ifelse(start == 2, 1, start), #start start position from 1 to match epitopes starting at 1
               protein_aa = str_sub(Sequence, start, end)) # add column containing the sequence of each protein
}


# This function uses the output of the `read_ev_polyprotein_uniprot_metadata()` function
# to plot the enterovirus polyprotein labeled with proteins

plot_ev_polyprotein <- function(ev_polyprotein_uniprot_metadata) {
    
    ev_proteins <- c("VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3AB", "3C", "3D")
    
    protein_sections <- c(ev_polyprotein_uniprot_metadata$start, ev_polyprotein_uniprot_metadata$end)
    
    protein_colours <- c(
        "VP4" = "#428984",
        "VP2" = "#6FC0EE",
        "VP3" = "#26DED8E6",
        "VP1" = "#C578E6",
        "2A" = "#F6F4D6",
        "2B" = "#D9E8E5",
        "2C" = "#EBF5D8",
        "3AB" = "#EDD9BA",
        "3C" = "#EBD2D0",
        "3D" = "#FFB19A")
    
    ggplot(ev_polyprotein_uniprot_metadata) +
        geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 0.1, fill = ev_proteins), alpha = 1) +
        geom_text(aes(x = (start + end) / 2, y = 0.05, label = ev_proteins), color = "black", size = 3, fontface = "bold", family = "Verdana") +
        geom_segment(data = data.frame(x = protein_sections), aes(x = x, xend = x, y = 0, yend = 0.1), linetype = "solid", color = "black", linewidth = 0.2) +
        geom_segment(aes(x = start, xend = end, y = 0, yend = 0), color = "black", linewidth = 0.5) +  # Bottom line
        geom_segment(aes(x = start, xend = end, y = 0.1, yend = 0.1), color = "black", linewidth = 0.5) +  # Top line
        theme_minimal() +
        labs(title = "", x = "", y = "") +
        theme_bw(base_size = 12) +
        scale_y_continuous(expand = expansion(mult = c(0, 0))) +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              panel.border = element_blank(),
              legend.position = "none",
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank()) +
        annotate("text", x = min(ev_polyprotein_uniprot_metadata$start) - 5, y = 0.05, label = "5'", size = 3, hjust = 1) +
        annotate("text", x = max(ev_polyprotein_uniprot_metadata$end) + 5, y = 0.05, label = "3'", size = 3, hjust = 0) +
        scale_fill_manual(values = protein_colours) +
        geom_vline(xintercept = min(ev_polyprotein_uniprot_metadata$start), color = "black", linewidth = 0.2) + # Left vertical line
        geom_vline(xintercept = max(ev_polyprotein_uniprot_metadata$start), color = "black", linewidth = 0.2) # Right vertical line
}



# What bea's data looks like'
endia_virscan_hits_counts_annot1_3 <- read_tsv("CDI_VirScan_Outputs/run2/phip11_plate1-3_v_kiwook_2/CDIVirScan_000/phip11_plate1-3_v_kiwook_2_CDIVirScan_000_Hits_counts_annotated.tsv") %>% 
    select(-contains("Beads_Only")) %>% 
    pivot_longer(cols = starts_with("KiWook"), names_to = "participant_id", values_to = "abundance") 

endia_virscan_hits_counts_annot4_6 <- read_tsv("CDI_VirScan_Outputs/run2/phip11_plate4-6_v_kiwook_2/CDIVirScan_000/phip11_plate4-6_v_kiwook_2_CDIVirScan_000_Hits_counts_annotated.tsv") %>% 
    select(-contains("Beads_Only"),
           -KiWook_PhIP11_Plate6_81_G9_POWH_SAMPLE) %>% #remove positive control 
    pivot_longer(cols = starts_with("KiWook"), names_to = "participant_id", values_to = "abundance")

endia_virscan_long <- bind_rows(endia_virscan_hits_counts_annot1_3, endia_virscan_hits_counts_annot4_6) %>% 
    group_by(pep_id) %>%
    filter(any(abundance > 0)) %>% # remove peptides that are 0 in all samples
    ungroup() %>% 
    group_by(participant_id) %>% 
    filter(sum(abundance) != 0) %>% # remove samples with 0 peptides detected 
    ungroup() %>% 
    mutate(sample_id = str_remove(participant_id, "^(?:[^_]+_){5}"),   # keep everything after the 5th _
           timepoint = str_extract(sample_id, "[^_]+$"),           # extract last segment from sample_id
           participant_id = str_remove(sample_id, "_[^_]+$"),        # remove timepoint from sample_id
           .keep = "unused")


# Extract peptide info and export to fa file
h_priority_species_peptide_counts = read_csv("MANGO_HighPriority_Species_Peptide_Count_Per_Sample.csv")
h_priority_species_peptide_info = read_csv("MANGO_HighPriority_Species_Peptide_Info.csv")
h_priority_species_peptide_info %>% 
    select(peptide_id,Prot) %>% 
    rename(peptide_seq = Prot) %>%
    distinct() %>%
    as.data.frame() %>%
    df_to_faa(file = "MANGO_High_Priority_Peptides.fasta")


h_priority_species_peptide_info %>% 
    filter(grepl("^Enterovirus",Species)) %>%
    select(peptide_id,Prot) %>% 
    rename(peptide_seq = Prot) %>%
    distinct() %>%
    as.data.frame() %>%
    df_to_faa(file = "MANGO_High_Priority_Peptides_Enterovirus.fasta")


# have to do the blast outside R
ENDIA_blastp_evB1 <- read_blast("blast_results/blastp_endia_evB1_all_virscan_peps.blast")

# Trying to make my data look the same
sample_table = read_csv("../../../NextFlow_Storage/KIM17394/full_sample_table.csv")
sample_prefix = sample_table %>% 
    filter(control_status == 'empirical') %>% 
    select(sample_source) %>% distinct()
sample_prefix_2 = sample_prefix %>% filter(!grepl("^CMV",sample_source))
sample_status = sample_table %>% 
    select(sample_source, sample_status) %>% 
    filter(!sample_source %in% c("BG","CMV_1","CMV_2")) %>% 
    distinct() %>% 
    mutate(sample_status = str_to_title(sample_status))


s_pos_mean_hits_counts = read_csv("MANGO_Strict_Positive_Mean_Hits_Counts.csv")
peptide_table = read_csv("../../../NextFlow_Storage/KIM17394/Part1/peptide_table.csv")


s_pos_mean_hits_counts_2 = merge(
    s_pos_mean_hits_counts,
    peptide_table %>% select(peptide_id,UniProtEntry,Prot_Start,Prot_End,Peptide_Length),
    by = "peptide_id"
) %>%
    pivot_longer(cols = sample_prefix_2$sample_source, names_to = "participant_id", values_to = "abundance") %>%    
    group_by(peptide_id) %>%
    filter(any(abundance > 0)) %>%
    ungroup() %>%
    group_by(participant_id) %>%
    filter(sum(abundance) != 0) %>%
    ungroup() %>%
    mutate(
        sample_id = participant_id,
        timepoint = str_extract(sample_id, "[^_]+$"),
        participant_id = str_remove(participant_id,"_[^_]+$"),
        .keep = "unused"
    )

s_pos_mean_hits_counts_2 = merge(
    s_pos_mean_hits_counts_2,
    sample_status, 
    by.x = "sample_id", 
    by.y = 'sample_source'
)

s_pos_fc = calculate_mean_rpk_difference(s_pos_mean_hits_counts_2, sample_id, sample_status, peptide_id, abundance, ENDIA_blastp_evB1)


s_pos_ms = calculate_moving_sum(s_pos_fc, value_column = mean_rpk_difference, win_size = 32, step_size = 4)
s_pos_ms_plot <- s_pos_ms %>% 
    ms_plot_clean() +
    ggtitle("ENDIA")

# curl cmd to download coxsackievirus B1 genome metadata
# curl -X 'GET' 'https://rest.uniprot.org/uniprotkb/stream?query=accession%3AP08291%20AND%20reviewed%3Atrue&fields=id%2Csequence%2Caccession%2Cprotein_name%2Cft_binding%2Cft_chain%2Corganism_name%2Clength%2Creviewed&sort=accession%20desc' -H 'accept: text/plain;format=tsv'


coxsackievirusB1_P08291  = read_ev_polyprotein_uniprot_metadata("coxsackievirusB1_P08291.tsv")
EV_B1_plot <- plot_ev_polyprotein(coxsackievirusB1_P08291)

combined_ms_plots <- wrap_plots(EV_B1_plot,
                                s_pos_ms_plot + labs(x = "Position in sequence (amino acids)") + theme(legend.position = "bottom"),
                                ncol = 1, heights = c(0.3, 3, 3)) &
    theme(plot.margin = margin(5.5, 5.5, 5.5, 0),
          plot.title = element_text(size = 10)) 


figure_1 <- wrap_elements(combined_ms_plots) +
    geom_vline(xintercept = 848.5, linetype = "dashed", linewidth = 0.2) +
    labs(tag = expression(sum((bar(X)[rpk_cases] - bar(X)[rpk_controls])))) +
    theme(
        plot.tag = element_text(size = rel(1), angle = 90),
        plot.tag.position = "left"
    )
