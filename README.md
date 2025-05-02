# Tropics to the poles: A Snapshot of Coastal Eukaryotic Marine Microalgae Diversity Across Five Ecoregions.  
### Contents
[Pipeline Overview](#pipeline-overview)   |   [Abstract](#abstract)   |   [Script files](#script-files)   |   [Metadata examples](#metadata-examples)   |   [R requirements](#r-requirements)

## Pipeline Overview
_Authors:_  
Jacqui Stuart (1,2), Ken Ryan (1), Natalie Robinson (3), Svenja Halfter (3), John K. Pearman (1), 
Jacob Thomson-Laing (2), Kirsty F. Smith (2)

_Affiliations:_   
1.  School of Biological Sciences, Victoria University of Wellington, New Zealand.  
2.	Cawthron Institute, Private Bag 2, Nelson 7042, New Zealand.  
3.	National Institute of Water and Atmospheric Research (NIWA),  Wellington

**This analysis pipeline includes:**   
- Processing Raw Metabarcoding files
- Rarefaction Curves
- Alpha diversity analysis
- Beta diversity analysis
- Venn diagram construction
- Community composition and
- Dendrogram construction
- Indicator Species analysis 

All analysis are to compare the community composition of coastal marine microalgae from five ecoregions; Tropical, subtropical, temperate, sub-polar and polar. An additional analysis of assembly processes between and within ecoregions was complete, the script for this is in a separate file as production of the map kept causing the markdown file to crash.


## Abstract
The objectives of this study are to analyse and compare the community composition and diversity of planktonic marine microalgae in a latitudinal gradient spanning the tropical to polar ecoregions in the South Pacific. This was achieved using eDNA samples and multi-region small subunit rDNA metabarcoding. Our focus was to identify similarities and differences at the functional group level across a large latitudinal gradient and assess diversity.

Samples were collected from sites in tropical, sub-tropical, and temperate ecoregions (Figure 1a). All sites were coastal or near shore, without significant rainfall 2-3 days before sampling. All samples from the tropical, sub-tropical and temperate sites were collected in spring, with tropical collected in the spring of 2022 between 28th November – 30th November, sub-tropical and temperate samples collected between 20th October – 29th November 2021, and polar samples were between 10th October to 24th November, respectively. Sampling in the sub-polar regions was undertaken on 14th February 2023, and polar sites were sampled via holes in the sea ice (10th October - 24th November 2022), with an additional five open water polar site samples during the National Institute of Water and Atmospheric Research (NIWA) Tangaroa voyage to the Ross Sea  6th – 11th February 2023.

<table align="center" border="0">
  <tr>
    <td width="10%"></td>
    <td width="40%">
      <img src="images/figure3.jpg"/>
      <br></br>
      <sub>
        <strong>Figure 1.</strong> a) Map depicting sampling sites in tropical, sub-tropical to temperate, sub-polar and polar ecoregions. b) The taxonomic-functional groups of eukaryotic marine microalgae used in this study. The other categories include golden, yellow, and brown classes. Divisions are based on PR2 database classifications (Vaulot et al., 2022), as proposed by Adl et al. (2012), and the World Register of Marine Species (WoRMS 2023) database. Image from Stuart et al., 2024 DOI: [10.1038/s41598-024-56993-4](https://www.nature.com/articles/s41598-024-56993-4)
      </sub>
    </td>
    <td width="10%"></td>
  </tr>
</table>
The community composition of the EMC was divided into taxon-based functional groups: dinoflagellates, pennate diatoms, radial centric diatoms, polar centric diatoms, chlorophytes, haptophytes, and ‘Other’ eukaryotic microalgae (Figure 1b). Taxonomic ranks were identified using the revised classification of eukaryotic groups proposed by Adl et al. (2012), which is also used for taxonomic assignment in the PR2 database. This resolution was based on taxon-based functional divisions used to assess the microalgal community structure (Figure 1b; Kruk, 2002; Litchman et al., 2007; Litchman and Klausmeier, 2008; Edwards et al., 2013; Wentzky et al., 2020).  


## Script Files  

- **Scripts:** R scripts used to process raw metabarcoding files, assign taxonomy, assess microalgal community composition, assembly processes and indicator species across five ecoregions.
  - `Step1_PhD_Chp3_ProcessingMetabarcodingData`: This R script delineates a comprehensive metabarcoding data processing workflow, encompassing library installation, data organization, primer analysis, Cutadapt-based trimming, *DADA2* quality checks, filtering, and chimera removal. Utilising checkpoints ensures adaptability in handling interruptions or errors. Subsequently, the script navigates through pivotal steps: "Merging seqtab files," which consolidates tables from various sequencing runs based on gene region (v9), and saves them as .rds files. The "Taxonomy" section processes merged tables, metadata, and taxonomy files, crafting phyloseq objects. Addressing contamination, the "Subtracting negative controls" segment deducts sequences from negative controls, adjusting the main data accordingly. The script then eliminates zero samples, creating a new phyloseq object without negative controls. Concluding with a comparison of count sums, the script saves the cleaned dataset as .rds and tracks removed sequences in a CSV. It culminates by emphasizing the preparedness of the cleaned data for community and diversity analysis.
  - `Step2_PhD_Chp3_DataAnalysis_CommunityComp_and_Diversity`: 
  - `Step3_PhD_Chp3_DataAnalysis_AssemblyProcess_forGitHub`: 


## Metadata Examples



## R Requirements




  
