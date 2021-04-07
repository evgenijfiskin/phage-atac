# Custom code repository for PHAGE ATAC manuscript

## About

Each folder corresponds to a particular dataset or analysis mode that is presented in this work. Rscripts can be found in the `code` folder that can be executed sequentially to reproduce analyses in part or in whole. Intermediate output files are also hosted for ease of use. 

## Setup

To best use this resource, we recommend pairing with large data files (that are not compatible with github as they exceed 100Mb). These files are available from the [Open Science Framework].

Once one downloads the zip archieve from OSF (xx Gb), place the extracted folder in the same directory as this repository named `phage_atac_large_data_files` (as shown below). This will enable running custom code to reproduce items in the output folders. 


```
.
├── phage_atac_large_data_files
│   ├── input
│   └── output
└── phage-atac
    ├── README.md
    ├── cd8_hashing
    ├...
```

## Figure to directory mapping

**Note** this is for the updated version under review and _not_ the bioRxiv verison.
```
Figure 1 | species_mix_asapseq
Figure 2 | phage_asap_comparison
Figure 3 | cd8_hashing
Figure 4 | full_phage_nb_sequencing,spike_mix
```
The `intracellular` analysis is part of ED Figure 6. 

## Contact

For questions related to this work, please raise an issue on this repository. 

<br><br>

