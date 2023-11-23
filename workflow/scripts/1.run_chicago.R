#!/usr/bin/env Rscript

# Load environment and libraries ----
library(Chicago)
library(optparse)


# Define parser ----
parser <- optparse::OptionParser()

parser <- optparse::add_option(
    parser, c("--inputfolder"),
    type = "character",
    help = "Path to folder containing input files",
    metavar = "character"
)

parser <- optparse::add_option(
    parser, c("--designfolder"),
    type = "character",
    help = "Path to folder containing Chicago design files",
    metavar = "character"
)

parser <- optparse::add_option(
    parser, c("-o", "--out"),
    type = "character",
    help = "Path to output folder",
    metavar = "character"
)

parser <- optparse::add_option(
    parser, c("--samples"),
    type = "character",
    help = "Samples to use",
    metavar = "character"
)

parser <- optparse::add_option(
    parser, c("-t", "--tissue"),
    type = "character",
    help = "Tissue to test",
    metavar = "character"
)


# Get arguments. ----
args <- optparse::parse_args(parser)


# Retrieve the correct files. ----
regex_samples <- paste0(base::strsplit(args$samples, split = " ")[[1]],  collapse = "|")
inputs <- list.files(path = args$inputfolder, recursive = TRUE, full.names = TRUE, pattern = regex_samples)
inputs <- inputs[grepl("chinput$", inputs)]

# Run Chicago ----
chicago_data <- Chicago::setExperiment(designDir = args$designfolder)
chicago_data <- Chicago::readAndMerge(files = inputs, cd = chicago_data)
chicago_data <- Chicago::chicagoPipeline(chicago_data, outprefix = paste0(args$out, args$tissue))

# Export results and RDS.
Chicago::exportResults(chicago_data, sprintf("%s/%s_chicago", args$out, args$tissue))
saveRDS(chicago_data, file = sprintf("%s/%s_chicago.rds", args$out, args$tissue))


# Save full Chicago table without cutoff score ----
utils::write.table(
    x = chicago_data@x,
    file = sprintf("%s/%s_chicago_full_table.txt", args$out, args$tissue),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# Plot QC figures. ----
for (i in c(1:6)) {
    pngName <- paste0(args$out, "_", i, ".pdf")
    pdf(pngName)
    plotBaits(chicago_data, n = 1, plotBprof = T)
    dev.off()
}

# args <- optparse::parse_args(
#     parser,
#     args = c(
#         "--inputfolder=~/odomLab/capturehic/data/workflow/chicago/input/canis_familiaris/",
#         "--designfolder=~/odomLab/capturehic/data/workflow/chicago/design/canis_familiaris/",
#         "--samples=do26735 do26772",
#         "--tissue=brain",
#         "--out=~/odomLab/capturehic/data/workflow/chicago/results/canis_familiaris/"
#     )
# )
