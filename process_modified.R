args <- commandArgs(trailingOnly = TRUE)

mapping_stats <- data.frame(sample = character(), reads_in_input = integer(), reads_after_processing = integer(), reads_aligned = integer(), aligned_percentage = numeric(), modified_by_indel = numeric(), modified_by_3mer_indel = numeric(), modified_by_non3mer_indel = numeric())
i <- 1
for (j in args) {
	print(paste0("processing samples in ",j))
	setwd(j)
	htmls <- list.files(pattern = "CRISPResso_on.*\\.html")
	samples <- unlist(lapply(htmls, function(x) {
		tmp <- substr(x, 15, nchar(x) - 5)
	}))
	for (s in samples) {
		mapping_stats[i, 1] <- s
		tryCatch({
			data0 <- read.csv(paste0("CRISPResso_on_", s, "/CRISPResso_mapping_statistics.txt"), sep = "\t", header = T, stringsAsFactors = F)

			data1 <- as.data.frame(t(read.csv(paste0("CRISPResso_on_", s, "/Quantification_window_nucleotide_percentage_table.txt"), sep = "\t", stringsAsFactors = F, header = F, row.names = 1)))
			ind <- which(data1$V1 == "A")
			AG <- 100 * as.numeric(data1$G[data1$V1 == "A"])

			file <- list.files(paste0("CRISPResso_on_", s), pattern = "Alleles_frequency_table_around_sgRNA.*.txt")
			data <- read.csv(paste0("CRISPResso_on_", s, "/", file), sep = "\t", header = T, stringsAsFactors = F)
			modified_by_indel <- sum(data$X.Reads.1[unlist(apply(data[, c("n_deleted", "n_inserted")], 1, function(x) {
				return(any(x != 0))
			}))])
			modified_by_3mer_indel <- sum(data$X.Reads.1[unlist(apply(data[, c("n_deleted", "n_inserted")], 1, function(x) {
				return(any(x != 0) & (x[1] - x[2]) %% 3 == 0)
			}))])
			modified_by_non3mer_indel <- sum(data$X.Reads.1[unlist(apply(data[, c("n_deleted", "n_inserted")], 1, function(x) {
				return(any(x != 0) & (x[1] - x[2]) %% 3 != 0)
			}))])
			mapping_stats[i, 2:4] <- data0[, c("READS.IN.INPUTS", "READS.AFTER.PREPROCESSING", "READS.ALIGNED")]
			mapping_stats[i, 5] <- 100 * data0[1, "READS.ALIGNED"] / data0[1, "READS.IN.INPUTS"]
			mapping_stats[i, 6] <- modified_by_indel
			mapping_stats[i, 7] <- modified_by_3mer_indel
			mapping_stats[i, 8] <- modified_by_non3mer_indel
			if (length(AG) >= 1) {
				mapping_stats[i, 9:(8 + length(AG))] <- AG
			}
		},error=function(e){
			cat("发生错误: ", e$message, "\n")
		})
		
		i <- i + 1
	}
}

names(mapping_stats)[1:8] <- c("sample", "reads_in_input", "reads_after_processing", "reads_aligned", "aligned_percentage (%)", "modified_by_indel (%)", "modified_by_3mer_indel (%)", "modified_by_non3mer_indel (%)")
if (ncol(mapping_stats) >= 9) {
	names(mapping_stats)[9:ncol(mapping_stats)] <- paste("A", 1:(ncol(mapping_stats) - 8), " to G (%)", sep = "")
}
setwd(args[1])
write.table(mapping_stats, "sample_summary.tsv", sep = "\t", row.names = F, quote = F, na = "")
