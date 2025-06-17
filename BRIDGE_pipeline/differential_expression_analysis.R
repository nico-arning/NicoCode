#.libPaths("./datapond/Dependencies/R_libs")
suppressPackageStartupMessages(library("findpython"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("BiocParallel"))

parser = ArgumentParser()

# specify our desired options
parser$add_argument(
  "-pl",
  "--ReadPools",
  type = "character",
  help = "List seperated by commas denoting by which regular start of the file to group the reads by"
)

parser$add_argument(
  "-acc",
  "--accession",
  type = "character",
  help = "Accession of the experiment"
)

parser$add_argument(
  "-st",
  "--SheetTreatments",
  type = "character",
  help = "List seperated by commas denoting which treatment sheet to use for each plate"
)

parser$add_argument(
  "-sc",
  "--SheetConcentrations",
  type = "character",
  help = "List seperated by commas denoting which concentration sheet to use for each plate"
)

parser$add_argument(
  "-bl",
  "--BridgeLibraryR",
  default = "./datapond/Scripts/differential_expression_library.R",
  type = "character",
  help = "Location of the differential expression library"
)

parser$add_argument(
  "-md",
  "--MetaData",
  type = "character",
  help = "Location of the metadata"
)

parser$add_argument(
  "-c",
  "--cores",
  type = "integer",
  default = 2,
  help = "Number of cores"
)

parser$add_argument(
  "-odb",
  "--orgdb",
  default = "org.Sfrugiperda.eg.db",
  type = "character",
  help = "orgdb made with eggnogg"
)

parser$add_argument(
  "-gtf",
  "--gtf",
    default = "/workdir/RefGenomes/SpodopteraFrugiperda/GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.gtf",
  type = "character",
  help = "GTF file"
)

parser$add_argument(
  "-en",
  "--eggnog",
  default = "/workdir/RefGenomes/SpodopteraFrugiperda/SpodopteraFrugiperda.emapper.annotations",
  type = "character",
  help = "Eggnnog .annotations file"
)


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args = parser$parse_args()

############# This is just for testing the analysis
#args$ReadPools="placeholder"
#args$SheetTreatments="placeholder"
#args$SheetConcentrations = "placeholder"
#args$MetaData="./MetaData/Lyon_bulk_raw.csv"
#args$BridgeLibraryR="/home/ubuntu/Scripts/differential_expression_library.R"
#args$c = 20
#args$gtf = "/workdir_nextflow/RefGenomes/BotrytisCinerea/Botrytis_cinerea.ASM83294v1.59.gtf"
#args$eggnog = "/workdir_nextflow/RefGenomes/BotrytisCinerea/BotrytisCinerea.emapper.annotations"
#args$orgdb = "/workdir_nextflow/RefGenomes/BotrytisCinerea/org.Bcinerea.eg.db"
#############go_enrichment
print(args)
# Source functions from our library
source(args$BridgeLibraryR)

# Make directories to dump data in
system("mkdir -p ./DeSeq/Tables")
system("mkdir -p ./DeSeq/Tables/DoseResponses")
system("mkdir -p ./DeSeq/Plots")
system("pwd")

# Register cores for parallelisation
register(MulticoreParam(workers = as.numeric(args$c)))

# Get the different plates
read_pools = unlist(strsplit(args$ReadPools, ","))

# Get the number of plates
plate_number = length(read_pools)

# Get which treatment sheets to use for each read pool
sheet_treatments = unlist(strsplit(args$SheetTreatments, ","))
sheet_concentrations = unlist(strsplit(args$SheetConcentrations, ","))

if (grepl("\\.csv$", args$MetaData)) {

  data <- read.csv(args$MetaData)
  data = data[, 2:ncol(data)]
    plate_matrix <- as.data.frame(data[, grepl("^Metadata", colnames(data))])
    matrix <- as.data.frame(t(data[, !grepl("^Metadata", colnames(data))]))
    gene_names = rownames(matrix)
    matrix <- as.data.frame(apply(matrix, 2, as.integer))
    rownames(matrix) = gene_names

    # Make the substitutions in the column names for plate_matrix
    colnames(plate_matrix) <- colnames(plate_matrix) %>%
      gsub("Metadata_treatments", "Treatments", .) %>%
      gsub("Metadata_treatment", "Treatments", .) %>%
      gsub("Metadata_plate", "Plate", .) %>%
      gsub("Metadata_Plate", "Plate", .) %>%
      gsub("Metadata_concentration", "Concentrations", .) %>%
      gsub("Metadata_Concentrations", "Concentrations", .) %>%
      #gsub("Metadata_timepoint", "Timepoints", .) %>%
      gsub("Metadata_", "", .)
    rownames(plate_matrix) = plate_matrix$Metadata_Wells
    colnames(matrix) = rownames(plate_matrix)
    plate_matrix$Concentrations = as.numeric(plate_matrix$Concentrations)
    concentrations_array = as.array(plate_matrix$Concentrations)
    concentrations_array = concentrations_array[concentrations_array != 0]
    concentrations_array = concentrations_array[!is.na(concentrations_array)]

    plate_matrix$Concentrations[is.na(plate_matrix$Concentrations)] = 1
    plate_matrix$Concentrations[plate_matrix$Concentrations == 0] = 1
    plate_matrix$Concentrations = as.numeric(plate_matrix$Concentrations)
    plate_matrix$logdose = log10(plate_matrix$Concentrations)
    plate_matrix$ReadSumsPerWell = colSums(matrix)
    plate_matrix$Wellnames = plate_matrix$Wells
    # =========================================================================
    #plate_matrix = plate_matrix[plate_matrix$Timepoints == 24, ]
    # =========================================================================

    if ("Timepoints" %in% colnames(plate_matrix) == 1) {
      plate_matrix$Treatments = paste0(plate_matrix$Treatments, " at ", plate_matrix$Timepoints)
      # Use gsub to replace strings that start with "DMSO" with "DMSO"
      plate_matrix$Treatments <- gsub("^DMSO.*", "DMSO", plate_matrix$Treatments)
    }
    if ("Wells" %in% colnames(plate_matrix) == 1) {
      plate_matrix$Wellnames = plate_matrix$Wells
    } else {
      plate_matrix$Wellnames = paste(plate_matrix$Treatments, plate_matrix$Concentrations, sep = "_")
    }
  
} else {
  data = read_in_meta(read_pools, sheet_treatments, sheet_concentrations, args$MetaData, args$accession)
  matrix = data$matrix
  plate_matrix = data$plate_matrix
  plate_matrix = plate_matrix[rowSums(is.na(plate_matrix)) == 0, ]
}

print(plate_matrix)
write.csv(plate_matrix, "./DeSeq/Tables/plate_matrix.csv")

# Check for treatments with only one unique logdose
# single_logdose_treatments <- plate_matrix %>%
#   group_by(Treatments) %>%
#   summarise(unique_logdose_count = n_distinct(logdose)) %>%
#   filter(unique_logdose_count == 1) %>%
#   pull(Treatments)



#print("The following treatments only have one unique Concentration")
#print(Treatments)

# Remove samples with treatments that have only one unique logdose
#plate_matrix <- plate_matrix[!plate_matrix$Treatments %in% single_logdose_treatments, ]

# Ensure plate_matrix and matrix have matching columns
matrix = matrix[, rownames(plate_matrix)]

plate_matrix_copy = plate_matrix
colnames(plate_matrix_copy) = paste("Metadata_", colnames(plate_matrix_copy))
colnames(plate_matrix_copy) <- colnames(plate_matrix_copy) %>%
  gsub("Metadata_ ", "Metadata_", .)
complete_matrix = cbind(t(as.matrix(matrix)), plate_matrix_copy)
write.csv(complete_matrix, "full_raw_matrix_and_meta.csv")

system("pwd")

# # Plot the distrubtion of Read counts per gene
# sparse_plot(matrix, plate_matrix, "Raw_counts")

# Run DeSeq for analysing differential expression
de_results = differential_expression(matrix, plate_matrix, length(read_pools))
dds = de_results$dds
vsd = de_results$vsd

# Parse the wellnames
vsd$Wellnames = sapply(strsplit(as.character(vsd$Wellnames), "_"), `[`, 1)
corrected_counts = counts(dds, normalized = TRUE)

complete_matrix = cbind(t(as.matrix(corrected_counts)), plate_matrix_copy)
write.csv(complete_matrix, "full_corrected_matrix_and_meta.csv")

# Initialise the outline columns
vsd$Outlier = ""

# Make PCA before batch correction
make_pca(vsd, plate_matrix, "./DeSeq/Plots/Plate_uncorrected_pca.png", vsd$Plate)

print("Removing batch effects")
# Correct for batch effects
if (plate_number > 1) {
  assay(vsd) = limma::removeBatchEffect(
    assay(vsd),
    batch = vsd$Plate,
  )
}
print("Done removing batch effects")

# Get correct counts of reads per well
plate_matrix$CorrectedReadsPerWell = colSums(assay(vsd))

# Plot the sums of reads per well
reads_per_well_plot(plate_matrix)

# # Perform automated outlier removal by GridPCA
outliers = outlier_removal_pca(vsd)
vsd$Outlier[colnames(vsd) %in% outliers] = "Outlier"

# Write the corrected counts and variance stabilised counts to file
# For further processing with python
write.csv(corrected_counts, "./DeSeq/Tables/corrected_count_matrix.csv")
tmp = assay(vsd)
tmp$Metadata_treatments = vsd$Treatments
tmp$Metadata_concentration = vsd$Concentrations
write.csv(tmp, "./DeSeq/Tables/variance_stabilised_counts.csv")

# Make PCA after batch correction
make_pca(vsd, plate_matrix, "./DeSeq/Plots/Plate_corrected_treatments_pca.png", vsd$Treatments)
make_pca(vsd, plate_matrix, "./DeSeq/Plots/Plate_corrected_readsums_pca.png", vsd$ReadSumsPerWell)

# # Remove outliers

# plate_matrix = plate_matrix[!(rownames(plate_matrix) %in% outliers), ]

# # Check for treatments with only one unique logdose
# single_logdose_treatments <- plate_matrix %>%
#   group_by(Treatments) %>%
#   summarise(unique_logdose_count = n_distinct(logdose)) %>%
#   filter(unique_logdose_count == 1) %>%
#   pull(Treatments)

# # Print the treatments that have only one unique logdose
# print("The following treatments only have one unique logdose:")


# #Remove samples with treatments that have only one unique logdose
# plate_matrix <- plate_matrix[!plate_matrix$Treatments %in% single_logdose_treatments, ]


# matrix = matrix[, rownames(plate_matrix)]
# dds = dds[, !(colnames(dds) %in% outliers)]
# vsd = vsd[, !(colnames(vsd) %in% outliers)]

# de_results = differential_expression(matrix, plate_matrix, length(read_pools))
# dds = de_results$dds
# vsd = de_results$vsd
# vsd$Outlier=""

# Make another PCA without outliers
# make_pca(vsd, plate_matrix, "./DeSeq/Plots/Outlier_removed_pca.png", vsd$Treatments)

# Write to file for downstream processing
corrected_counts = counts(dds, normalized = TRUE)

plate_matrix = plate_matrix[colnames(corrected_counts),]

corrected_counts = as.data.frame(t(corrected_counts))

corrected_counts$Treatments = plate_matrix$Treatments
corrected_counts$Concentration = plate_matrix$Concentrations
write.csv(corrected_counts, "./DeSeq/Tables/corrected_count_matrix.csv")

# tmp = assay(vsd)
# tmp$Metadata_treatments = vsd$Treatments
# tmp$Metadata_concentration = vsd$Concentrations
# write.csv(tmp, "./DeSeq/Tables/variance_stabilised_counts.csv")

if (args$orgdb == "org.Hsapiens.eg.db") { # This is human
    print("Using human annotation")
    setwd("./DeSeq/Plots")

    for (treatment in unique(dds$Treatments)) {
    if (treatment != "DMSO") {
        print(treatment)
        volcano_plot(dds, treatment, "Treatments", "DMSO")
        functional_analysis_human(dds)
    }
    }
    
} else { # Everything else
    if (basename(args$orgdb) == "org.Bcinerea.eg.db") {
        gene2geneID =  gsub("g", "p", rownames(matrix))
        names(gene2geneID) = rownames(matrix)
    } else {
        gene2geneID = GetGene2ID(args$gtf)
    }

    kegg = PrepEggnog(args$eggnog)
    GetOrgdb(args$orgdb)

    setwd("./DeSeq/Plots")

    for (treatment in unique(dds$Treatments)) {
    if (treatment != "DMSO") {
        print(treatment)
        volcano_plot(dds, treatment, "Treatments", "DMSO")
        functional_analysis(dds, args$orgdb, gene2geneID, treatment, kegg)
    }
}

}
