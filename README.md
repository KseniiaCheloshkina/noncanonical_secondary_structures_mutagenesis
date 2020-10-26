# noncanonical_secondary_structures_mutagenesis
Reproduce results of paper Georgakopoulos-Soares et al. "Noncanonical secondary structures arising from non-B DNA motifs are determinants of mutagenesis" (2018)

Structure of repository:
 - run: all scripts used in paper. Main scripts:
     - filter_bins.R: script for getting the list of all genome windows for specified window length excluding first and last windows, centromeres, telomeres and blacklisted regions
     - target.R: script for processing breakpoints data
 - cbp_data: common tools for R bioinformatics projects
