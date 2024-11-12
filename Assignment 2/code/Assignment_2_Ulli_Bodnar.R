# Setup --------------------------------------------------------

# Import the libraries so that their functions can be used
# Just a heads up, I find R kinda unpredictible with the cache and holding libraries loaded so if this doesn't work by just the ones i've uncommented, please just uncomment all of these below and run it again. It will for sure work with all loaded, but I am pretty sure I only need the ones that are labelled 'in use'
# MIGHT NEED LIBRARIES:
# library(stats)
# library(sf)
# library(randomForest)
# library(viridis)
# library(seqinr)

# Libraries in use
library(tidyverse)
library(rentrez)
library(Biostrings)
library(fmsb)
library(ape)
library(dendextend)
library(phytools)
library(DECIPHER)
library(muscle)

#The method to load libraries above is standard, however with number of libraries used for this project it will be best to use lapply or purr. I am gonna use lapply as I am familiar with that
# packages <- c('tidyverse', 'rentrez','Biostrings', 'fmsb',' ape', 'dendextend', 'phytools', 'DECIPHER', 'muscle')
# lapply(packages, library, character.only = T)
# rm(packages)

#The code above will make your loading of packages in 2 clicks rather than 9.
#If packages mentioned above aren't installed you can use the if function below to install them
#  installed_packages <- packages %in% rownames(installed.packages())
# if (any(!installed_packages)) {
#   install.packages(packages[!installed_packages])
# }
#Now its 3 :)

# Set working directory to utilize data from datasets
setwd("/Users/ullibodnar/Documents/School/Guelph Masters/Bioinformatics Software Tools/Assignment 2/code")

#Good practice to remind user to set working directory, however it would be best to just remind the user rather than giving path to your working directory such as

#Set working directory according to your system
# setwd('./Assignment_2')

#Not sure you were there but we were told by Dr. Cottenie that we aren't suppose to throw in setwd() in our code, but its still a good reminder for the user



# Global variables ---------------------------------------------------------------

missingData <- 0.01 #can potentially use a better variable name-M
lengthVar <- 50
chosenModel <- "K80" # K2P
clusteringMethod <- "ML"

#excellent use of global variables.


# Global functions --------------------------------------------------------

# Use the species name from BOLD dataset to obtain corresponding mass (in grams) from Pantheria dataset
getAverageMass <- function (df, speciesName) {
  averageMass <- df |>
    filter(MSW05_Binomial == speciesName) |>
    pull("5-1_AdultBodyMass_g")
  
  # Substitute -999 masses
  ifelse(averageMass |>
           length() == 0
         , 0, averageMass)
}



# Importing data ----------------------------------------------------------

# BOLD DB
# API CALL TO SHOW I KNOW HOW TO DO IT
# Lepus <- read_tsv("http://v3.boldsystems.org/index.php/API_Public/combined?taxon=Lepus&format=tsv")
# write_tsv(Lepus, "lepus_bold_data.txt")

# Import lepus data from BOLD database to have larger selection of COI-5P sequences
rawBoldLepus <- read_tsv(file = "../data/lepus_bold_data.txt")

# Import Pantheria DB for body mass in grams data
pantheriaData <- read_tsv(file = "../data/Pantheria.tsv")
#where did you get Pantheria dataset, also you mispelled it when writing it, so it wouldnt load :(.
#pantheriaData <- read_tsv(file = "../data/PanTHERIA.tsv")
# NCBI's nucleotide ---
# Determine possible database search locations
entrez_dbs()

# Sequence length range chosen becuase NCBI returned 671bp length for this sequence of interest when inspecting the database. Range to catch ones with missing data
# Query NCBI's nucleotide database to test the waters of possible data
lepusSearch <- entrez_search(db = "nucleotide", term = "Lepus[ORGN] AND COI AND 600:800[SLEN]", retmax = 100)

# Determine result count and change entrez_search to provide unique IDs for each possible returned result
lepusRetmax <- lepusSearch$count
lepusSearch <- entrez_search(db = "nucleotide", term = "Lepus[ORGN] AND COI AND 600:800[SLEN]", retmax = lepusRetmax)

# Fetch the data as fasta
lepusFetch <- entrez_fetch(db = "nucleotide", id = lepusSearch$ids, rettype = "fasta")

# Write the data to fasta file; commented out because the data have already been imported
# write(lepusFetch, "lepus_fetch.fasta", sep = "\n") --- Done on Oct. 18

#good job so far, moving forward from here it will be good practice to rm variables you are no longer using
rm(lepusFetch,lepusRetmax,lepusSearch, packages)
# Import already created NCBI data file
nucleotideLepusStringSet <- readDNAStringSet("../data/lepus_fetch.fasta")

# Convert to a dataframe object for easy manipulation
nucleotideLepus <- data.frame(title = names(nucleotideLepusStringSet), nucleotides = paste(nucleotideLepusStringSet))



# Formatting data ---------------------------------------------------------
# Subset bold lepus for only sections needed to reduce cognitive load when viewing dataframe
# remove ones that don't contain nucleotide data and ones that arent COI
boldLepus <- rawBoldLepus[, c("processid", "species_name", "markercode", "nucleotides")] |>
  filter(!is.na(nucleotides)) |>
  filter(markercode == "COI-5P")
#you are subsetting then using tidyverse, you can tidyverse the whole thing in one shot
# boldLepus <- rawBoldLepus %>%
#   select(processid,species_name,markercode,nucleotides) %>%
#   filter(!is.na(nucleotides)) %>%
#   filter(markercode == 'COI-5P')

# Change BOLD's processid column name to "id" to match NCBI database
names(boldLepus)[names(boldLepus) == "processid"] <- "id"
  
# Manipulate the NCBI dataframe to have correct column names for later merge with boldLepus
nucleotideLepus$id <- word(nucleotideLepus$title, 1L)
nucleotideLepus$species_name <- word(nucleotideLepus$title, 2L, 3L)

# filter out any non cytochrome sequences
nucleotideLepus <- filter(nucleotideLepus, grepl("cytochrome oxidase subunit", title))

# Add marker code column and rearrange columns to allow for clean merge with bold
nucleotideLepus$markercode <- "COI-5P"
nucleotideLepus <- nucleotideLepus[ , c("id", "species_name", "markercode", "nucleotides")]


# MERGE DATAFRAMES
lepusSeq <- merge(nucleotideLepus, boldLepus, all = T)

#Here would be a good place to remove all variables that are no longer in use moving forward to clean up the workspace
#rm(boldLepus,nucleotideLepus, nucleotideLepusStringSet, rawBoldLepus)

# Map body masses from Pantheria to the subsetted column for downstream analysis of body masses
lepusSeq$mass_g <- purrr::map(lepusSeq$species_name, function (x) {getAverageMass(pantheriaData, x)}) |>
  as.numeric()

#very nice and fancy, professional right here-M

# Remove duplicates, trim Ns from the ends and remove gaps, remove entries with NA species_name, and remove entries with no mass data
lepusSeq <- lepusSeq[!duplicated(lepusSeq$nucleotides), ]
#you can use distinct(nucleotides, keep_all = T) below to remove duplicates-M
lepusSeq <- lepusSeq |>
  filter(!is.na(species_name)) |>
  filter(mass_g > 0) |>
  mutate(nucleotides2 = str_remove_all(nucleotides, "^N+|N+$|-")) |>
  filter(str_count(nucleotides2, "N") <= (missingData * str_count(nucleotides))) |>
  filter(str_count(nucleotides2) >= median(str_count(nucleotides2)) - lengthVar & str_count(nucleotides2) <= median(str_count(nucleotides2)) + lengthVar)

# Create a subset dataframe containing a random sequence from each species; this allows for comparison between just one sample per species
set.seed(1234) # so that we get the same result

lepusSeqSubset <- lepusSeq |>
  group_by(species_name) |>
  slice_sample(n = 1) |> # Randomly selects one row per species, because Karl told me to do it :)
  ungroup() |>
  as.data.frame()
#LOL @ Karl's comment would probably help as to why its needed-M


# View data in radar chart ------------------------------------------------
# The following chart idea taken from https://r-graph-gallery.com/142-basic-radar-chart.html

# Map the data as names and mass to visualize in a radar chart
massAndNames <- as.data.frame(matrix(lepusSeqSubset$mass_g , ncol=12))
colnames(massAndNames) <- paste("L.", word(lepusSeqSubset$species_name, 2L), sep = " ")

# Add in the upper and lower limits to comply to the formatting requirements for fmsb library
massAndNames <- rbind(rep(5000,12), rep(1000,12), massAndNames)

# check data formatted properly
# head(massAndNames)

# Remove styling from other graph so it doesn't display weird when running
par(mar=c(1,1,1,1))

# Render the radar chart
radarchart( massAndNames, axistype=1, 
            #customize the polygon, grid, and labels
            pcol=rgb(1,0.6,0,0.9) , pfcol=rgb(1,0.6,0,0.5) , plwd=4 , 
            cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(1000,5000, 1000), cglwd=0.8,
            vlcex=0.8 
)



# Aligning sequences ------------------------------------------------------
# put entire lepus subset into DNAStringSet to work with the library
lepusSeqSubset$nucleotides2 <- DNAStringSet(lepusSeqSubset$nucleotides2)

# Map the names, substituting L. for Lepus for readability. 
names(lepusSeqSubset$nucleotides2) <- paste("L.", word(lepusSeqSubset$species_name, 2L), sep = " ")
# Check that it worked
names(lepusSeqSubset$nucleotides2)

# Conduct alignment with muscle
lepusSeqSubsetAlignment <- DNAStringSet(muscle::muscle(lepusSeqSubset$nucleotides2))
# Check it out in the browser to see if anything is out of place --> originally, I saw that one of the names was NA and I had to go back to remove NA species_name entries
#BrowseSeqs(lepusSeqSubsetAlignment)

#nice check! -M



# Clustering --------------------------------------------------------------
# Convert to a dataclass used by Ape for distance clustering
lepusSeqBin <- as.DNAbin(lepusSeqSubsetAlignment)
# Check conversion is correct
class(lepusSeqBin)

# Create distance matrix
distanceMatrix <- dist.dna(lepusSeqBin, model = chosenModel, as.matrix = TRUE, pairwise.deletion = TRUE)
# Check out distance matrix worked
head(distanceMatrix)

# PHYLOGENETIC TREE
clustersLepusCOI <- DECIPHER::TreeLine(lepusSeqSubsetAlignment,
                                       myDistMatrix = distanceMatrix,
                                       method = clusteringMethod,
                                       model = chosenModel,
                                       reconstruct = TRUE,
                                       maxTime = 0.01)

#Try to add in a metric to measure the model efficiency-M

# Bottom, left, top, right margins
par(mar=c(5,5,1,10))

# Create a vector of masses scaled relative to the max mass 
maxMass <- max(lepusSeqSubset$mass_g)
scaledMass <- lepusSeqSubset$mass_g / maxMass * 2

# Plot the phylogenetic tree
clustersLepusCOI |>
  set("leaves_pch", 20)  |>
  set("leaves_cex", scaledMass) |>
  set("nodes_col", "orange") |>
  set("labels_col", "black") |>
  plot(horiz = TRUE, xlab = "Distance", ylab = "Species")


# Test for phylogenetic conservatism --------------------------------------

# Export tree as a format used by ape package
tree_phylo <- as.phylo(clustersLepusCOI)

# Lambda estimation will test the extent to which body mass shows phylogenetic signal. 
# If it is close to 1, the trait follows the phylogeny
# If it is close to 0, the trait does not really follow the phylogenetic signal and does not exhibit phylogenetic conservatism
lambda_estimation <- phylosig(tree_phylo, lepusSeqSubset$mass_g, method = "lambda")

print(lambda_estimation)

#Excellent

