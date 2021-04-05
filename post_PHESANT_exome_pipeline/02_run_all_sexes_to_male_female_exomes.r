source("functions_for_run_all_sexes_male_female.r")

# Need the phenotype summary file as well.
# Use 01_run_summarise_phenotypes_initial_PHESANT_run.r if these have not yet been created.

# Note that this currently needs to be performed locally - awk on the cluster segfaults out, I think due to differences 
# between linux and mac versions of awk and how many columns are dealt with.

# Ensure that the cts_single and cts_multi tsv files are created locally and copied up to the cloud.
{

	if (!file.exists(paste0(QCed_io_name, ".tsv"))) {
		system(paste0("gsutil cp ", QCed_io_name_cloud, " ../../"))
	}

	# Pull out the make and female userIds for our application. Do this using the initially create file (paste0(QCed_io_name, "tsv")).
	sex_column <- grep("x31_", strsplit(system(paste0("head -1 ", QCed_io_name, ".tsv"), intern=TRUE), split='\t')[[1]])
	userID_column <- grep("userId", strsplit(system(paste0("head -1 ", QCed_io_name, ".tsv"), intern=TRUE), split='\t')[[1]])

	if (!file.exists(sex_df_path)) {
		cat('Creating data table containing the sexes of samples in the entire data set to restrict to...\n')
		sex_df <- data.table(userID = system(paste0("awk 'BEGIN { FS=\"\t\" } { print $", userID_column, " }' ", QCed_io_name, ".tsv"), intern=TRUE)[-1],
							 sex = system(paste0("awk 'BEGIN { FS=\"\t\" } { print $", sex_column, " }' ", QCed_io_name, ".tsv"), intern=TRUE)[-1])
		cat('Created data table containing the sexes of samples in the entire data set to restrict to...\n')
		fwrite(sex_df, file = sex_df_path, sep='\t')
	} else {
		sex_df <- fread(sex_df_path)
	}

	males <- (sex_df %>% filter(sex == "1") %>% select('userID'))$userID
	females <- (sex_df %>% filter(sex == "0") %>% select('userID'))$userID

	# Let's read in the phenotype information file
	phenofile <- fread(variable_info, sep='\t', header=TRUE)

	# Want to pull out variables that end up as cts variables in the both_sex PHESANT file.
	# First, let's get the header.
	file <- paste0(QCed_io_name , ".tsv")
	header <- strsplit(system(paste0("head -n 1 ", file), intern=TRUE), split='\t')[[1]]

	single_cts_columns <- c()
	multi_cts_columns <- c()

	for(i in 1:n_chunks) {
		pheno_summary <- paste0(filename_root, i, "_phenosummary.tsv")
		
		if(!file.exists(pheno_summary)) {
			system(paste0("gsutil cp ", intermediate_output_location, "*", i, "_phenosummary.tsv", " ../../"))
		}

		cts_variables <- system(paste("grep IRNT", pheno_summary, "| cut -f1 -d'\t'"), intern=TRUE)

		# These columns can be written to a new file as is (no need to perform any averaging)
		single_cts_columns <- c(single_cts_columns, unlist(lapply(cts_variables, get_raw_cts, header, single=TRUE)))
		multi_cts_columns <- c(multi_cts_columns, unlist(lapply(cts_variables, get_raw_cts, header, single=FALSE)))
	}

	# awk out the single cts columns and write them to a file.
	# Include the userId
	single_cts_columns <- c(1, single_cts_columns)
	multi_cts_columns <- c(1, multi_cts_columns)
	outfile_single <- paste0(QCed_io_name , "_cts_single.tsv")

	# This is the part that can't be run on the cloud due to awk differences between OS X and linux
	if (!file.exists(outfile_single)) {
		if (local == TRUE) {
			system(paste0("awk -F $'\t' -v OFS=$'\t' '{print ", paste0("$", single_cts_columns, collapse=","), "}' ", file, " > ", outfile_single))
			system(paste0("gsutil cp ", outfile_single, " ", intermediate_output_location))
		} else {
			cat("Need to run locally to use awk script, before copying result to the cloud for the remainder of this script.\n")
			cat(paste0("Copy across using:\n", "gsutil cp ", intermediate_output_location, "*_cts_single.tsv ../../"))
			system(paste0("gsutil cp ", intermediate_output_location, "*_cts_single.tsv ", "../../"))
			stop()
		}

	}

	outfile_multi <- paste0(QCed_io_name , "_cts_multi.tsv")

	# and this too
	if (!file.exists(outfile_multi)) { 
		if (local == TRUE) {
			system(paste0("awk -F $'\t' -v OFS=$'\t' '{print ", paste0("$", multi_cts_columns, collapse=","), "}' ", file, " > ", outfile_multi))
			system(paste0("gsutil cp ", outfile_multi, " ", intermediate_output_location))
		} else {
			cat("Need to run locally to use awk script, before copying result to the cloud for the remainder of this script.\n")
			cat(paste0("Copy across using:\n", "gsutil cp ", intermediate_output_location, "*_cts_multi.tsv ../../"))
			system(paste0("gsutil cp ", intermediate_output_location, "*_cts_multi.tsv ", "../../"))
			stop()
		}
	}

	# Now, create the average columns for the other cts variables...
	# Get the names of the variables, and then grep for them and take the average.
	cts_multi <- fread(outfile_multi, header=TRUE, sep='\t', data.table=FALSE)
	cts_multi_to_average_over <- unique(gsub("_.*", "_", names(cts_multi)[-1]))
	codings <- reassign_all_values(substr(cts_multi_to_average_over, 2, nchar(cts_multi_to_average_over)-1), phenofile)

	# Need extra step for cts_multi
	if(length(codings) > 0) {
		new_codings <- list()
		for(i in 1:length(codings)) {
			for(j in 1:length(codings[[i]])) {
				if(j == 1) {
					new_codings[[i]] <- names(cts_multi)[grep(paste0("x", codings[[i]][j]), names(cts_multi))]
				} else {
					new_codings[[i]] <- c(new_codings[[i]], names(cts_multi)[grep(paste0("x", codings[[i]][j]), names(cts_multi))])
				}
			}
		}

		names(new_codings) <- names(codings)
		codings <- new_codings
		cts_multi <- change_values(codings, cts_multi, coding_info_file)
		cts_multi <- data.frame(userId=cts_multi$userId, sapply(cts_multi_to_average_over, average_over_cts_multi, cts_multi))
	}

	names(cts_multi) <- c("userId", gsub("x(.*)_.*", "\\1", cts_multi_to_average_over))

	# Now, write these, along with userID, to disk.
	cts_columns <- fread(outfile_single, header=TRUE, sep='\t', data.table=FALSE)
	names(cts_columns) <- gsub("x([^_]*)_.*", "\\1", names(cts_columns))

	# Make the alterations before the transform
	codings <- reassign_all_values(names(cts_columns), phenofile)
	if(length(codings) > 0)
		cts_columns <- change_values(codings, cts_columns, coding_info_file)

	cts_columns <- merge(cts_multi, cts_columns, by='userId')

	if(	!file.exists(paste0(QCed_io_name , "_cts_irnt.tsv")) |
		!file.exists(paste0(QCed_io_name , "_cts_irnt_males.tsv")) |
		!file.exists(paste0(QCed_io_name , "_cts_irnt_females.tsv")) )
	{
		cat('Continuous IRNT variables...\n')

		# As a check, I want to IRNT the cts variables in all sexes, this should be the same as before.
		# Not anymore - as we've made some changes to the encodings and PHESANT!
		# The male and female IRNT should be different.

		cts_output <- data.frame(sapply(cts_columns, irnt))
		cts_output_males <- data.frame(sapply(cts_columns[cts_columns$userId %in% males,], irnt))
		cts_output_females <- data.frame(sapply(cts_columns[cts_columns$userId %in% females,], irnt))

		names(cts_output) <- names(cts_columns)
		names(cts_output_males) <- names(cts_columns)
		names(cts_output_females) <- names(cts_columns)

		# Finally, need to remove all the cts variables that have less than 5000 instances of non-NAs.
		cat(paste("Started with", ncol(cts_output), "cts.\n"))
		to_keep_all_sexes_cts <- apply(cts_output, 2, keep_cat_ordered, n_cts_min)
		cat(paste(sum(to_keep_all_sexes_cts), "cts remain for both sexes.\n"))
		to_keep_all_sexes_cts[1] <- TRUE
		to_keep_males_cts <- apply(cts_output_males, 2, keep_cat_ordered, n_cts_min)
		cat(paste(sum(to_keep_males_cts), "cts remain for males.\n"))
		to_keep_males_cts[1] <- TRUE
		to_keep_females_cts <- apply(cts_output_females, 2, keep_cat_ordered, n_cts_min)
		cat(paste(sum(to_keep_females_cts), "cts remain for females.\n"))
		to_keep_females_cts[1] <- TRUE

		fwrite(cts_output[,to_keep_all_sexes_cts], file=paste0(QCed_io_name , "_cts_irnt.tsv") , sep='\t')
		fwrite(cts_output_males[, to_keep_males_cts], file=paste0(QCed_io_name , "_cts_irnt_males.tsv") , sep='\t')
		fwrite(cts_output_females[, to_keep_females_cts], file=paste0(QCed_io_name , "_cts_irnt_females.tsv") , sep='\t')

		# Perform the check
		for(i in 1:n_chunks) {
			all_sexes <- fread(paste0(filename_root, i, ".tsv"), header=TRUE, sep='\t', data.table=FALSE)
			to_look <- names(all_sexes)[names(all_sexes) %in% names(cts_output)]
			all_sexes_check <- all_sexes[,names(all_sexes) %in% names(cts_output)]
			if(length(to_look) > 1) {
				for(j in 2:ncol(all_sexes_check)){
					damn <- max(abs(cts_output[,to_look[j]] - all_sexes_check[,j]),na.rm=TRUE)
					if(damn > 1e-14) {
						print(paste(to_look[j], names(all_sexes_check)[j]))
						print(damn)
					}
				}
			}
		}

		rm("cts_output", "cts_output_males", "cts_output_females")
	}

	if(	!file.exists(paste0(QCed_io_name , "_cts_raw.tsv")) |
		!file.exists(paste0(QCed_io_name , "_cts_raw_males.tsv")) |
		!file.exists(paste0(QCed_io_name , "_cts_raw_females.tsv")) ) {
		
		cat('Continuous raw variables...\n')

		# Now create the raw cts files.

		cts_output_raw <- cts_columns
		cts_output_males_raw <- cts_columns[cts_columns$userId %in% males,]
		cts_output_females_raw <- cts_columns[cts_columns$userId %in% females,]

		# Finally, need to remove all the cts variables that have less than 5000 instances of non-NAs.
		cat(paste("Started with", ncol(cts_output_raw), "cts.\n"))
		to_keep_all_sexes_cts_raw <- apply(cts_output_raw, 2, keep_cat_ordered, n_cts_min)
		cat(paste(sum(to_keep_all_sexes_cts_raw), "cts remain for both sexes.\n"))
		to_keep_all_sexes_cts_raw[1] <- TRUE
		to_keep_males_cts_raw <- apply(cts_output_males_raw, 2, keep_cat_ordered, n_cts_min)
		cat(paste(sum(to_keep_males_cts_raw), "cts remain for males.\n"))
		to_keep_males_cts_raw[1] <- TRUE
		to_keep_females_cts_raw <- apply(cts_output_females_raw, 2, keep_cat_ordered, n_cts_min)
		cat(paste(sum(to_keep_females_cts_raw), "cts remain for females.\n"))
		to_keep_females_cts_raw[1] <- TRUE

		fwrite(cts_output_raw[,to_keep_all_sexes_cts_raw], file=paste0(QCed_io_name , "_cts_raw.tsv"), sep='\t')
		fwrite(cts_output_males_raw[, to_keep_males_cts_raw], file=paste0(QCed_io_name , "_cts_raw_males.tsv") , sep='\t')
		fwrite(cts_output_females_raw[, to_keep_females_cts_raw], file=paste0(QCed_io_name , "_cts_raw_females.tsv") , sep='\t')

		rm("cts_output_raw", "cts_output_males_raw", "cts_output_females_raw")

	}

	for(i in 1:n_chunks)
	{
		cat(paste0('\nExamining chunk ', i, '\n'))
		all_sexes <- fread(paste0(filename_root, i, ".tsv"), header=TRUE, sep='\t', data.table=FALSE)

		# If there's just the age, sex and userId, skip to the next iteration of the loop.
		if(ncol(all_sexes) == 3) {
			next
		}

		log_file_both_sex <- paste0(filename_root, i, ".log")
		pheno_summary <- paste0(filename_root, i, "_phenosummary.tsv")
		cts_variables <- system(paste("grep IRNT", pheno_summary, "| cut -f1 -d'\t'"), intern=TRUE)

		# Now, look in the PHESANT file for IRNT - and exclude those from the male and females 
		all_sexes <- all_sexes[,-which(names(all_sexes) %in% c(cts_variables, "age", "sex")), drop=FALSE]
		if (ncol(all_sexes) == 1) {
			print("No categorical variables")
			break
		}

		PHESANT_males <- all_sexes[all_sexes$userId %in% males,]
		PHESANT_females <- all_sexes[all_sexes$userId %in% females,]

		# Check to make sure that all remaining variables are accounted for - they should be mapped to either Binary, or ordered categorical.
		cat_ordered <- system(paste("grep CAT-ORD", pheno_summary, "| cut -f1 -d'\t'"), intern=TRUE)
		cat_unordered <- system(paste("grep CAT-SINGLE-BINARY", pheno_summary, "| cut -f1 -d'\t'"), intern=TRUE)
		cat_mul_unordered <- system(paste("grep CAT-MUL-BINARY-VAR", pheno_summary, "| cut -f1 -d'\t'"), intern=TRUE)
		edge_cases_cat_unordered <- system(paste("grep 'Combine .* two bins and treat as binary'", pheno_summary, "| cut -f1 -d'\t'"), intern=TRUE)
		# Extra edge case - Integer phenotypes where there are only 2 available phenotypes are mapped to case-control, and we name them INT-BINARY-VAR.
		int_unordered <- system(paste("grep INT-BINARY-VAR", pheno_summary, "| cut -f1 -d'\t'"), intern=TRUE)

		# Check that the grepping for cat-unordered, cat-ordered, and cat-mul-unordered are disjoint.
		if(length(c(cat_unordered, cat_ordered, cat_mul_unordered, int_unordered))!= length(unique(c(cat_unordered, cat_ordered, cat_mul_unordered, int_unordered))))
			cat("Error: non-unique greps for CAT-ORD, CAT-SINGLE-BINARY and CAT-MUL-BINARY-VAR.\n")

		cat_unordered <- c(cat_unordered, cat_mul_unordered, edge_cases_cat_unordered, int_unordered)

		cat_unordered_all_sexes <- all_sexes[,which(names(all_sexes) %in% cat_unordered), drop=FALSE]
		cat_unordered_males <- PHESANT_males[,which(names(PHESANT_males) %in% cat_unordered), drop=FALSE]
		cat_unordered_females <- PHESANT_females[,which(names(PHESANT_females) %in% cat_unordered), drop=FALSE]

		cat_ordered_all_sexes <- all_sexes[,which(names(all_sexes) %in% cat_ordered), drop=FALSE]
		cat_ordered_males <- PHESANT_males[,which(names(PHESANT_males) %in% cat_ordered), drop=FALSE]
		cat_ordered_females <- PHESANT_females[,which(names(PHESANT_females) %in% cat_ordered), drop=FALSE]

		cat(paste("Started with", ncol(cat_unordered_all_sexes), "categorical unordered.\n"))
		to_keep_all_sexes <- apply(cat_unordered_all_sexes, 2, minimum_bin, n_cat_min)
		cat(paste(sum(to_keep_all_sexes), "categorical unordered remain for both sexes.\n"))
		to_keep_male <-  apply(cat_unordered_males, 2, minimum_bin, n_cat_min)
		cat(paste(sum(to_keep_male), "categorical unordered remain for males.\n"))
		to_keep_female <- apply(cat_unordered_females, 2, minimum_bin, n_cat_min)
		cat(paste(sum(to_keep_female), "categorical unordered remain for females.\n"))

		# Finally, need to remove ordered variables that are categorical ordered and have less that 5000 (the default in PHESANT) non-NA values.
		cat(paste("Started with", ncol(cat_ordered_all_sexes), "categorical ordered.\n"))
		to_keep_all_sexes_ordered <- apply(cat_ordered_all_sexes, 2, keep_cat_ordered, n_cat_ordered)
		cat(paste(sum(to_keep_all_sexes_ordered), "categorical ordered remain for both sexes.\n"))
		to_keep_males_ordered <- apply(cat_ordered_males, 2, keep_cat_ordered, n_cat_ordered)
		cat(paste(sum(to_keep_males_ordered), "categorical ordered remain for males.\n"))
		to_keep_females_ordered <- apply(cat_ordered_females, 2, keep_cat_ordered, n_cat_ordered)
		cat(paste(sum(to_keep_females_ordered), "categorical ordered remain for females.\n"))

		fwrite(cbind(all_sexes$userId, cat_ordered_all_sexes[, to_keep_all_sexes_ordered, drop=FALSE], cat_unordered_all_sexes[, to_keep_all_sexes, drop=FALSE]),
			sep='\t', file=paste0(QCed_io_name, "_cat_variables_both_sexes.", i, '.tsv'))
		fwrite(cbind(PHESANT_males$userId, cat_ordered_males[, to_keep_males_ordered, drop=FALSE], cat_unordered_males[, to_keep_male, drop=FALSE]),
			sep='\t', file=paste0(QCed_io_name, "_cat_variables_males.", i, '.tsv'))
		fwrite(cbind(PHESANT_females$userId, cat_ordered_females[, to_keep_females_ordered, drop=FALSE], cat_unordered_females[, to_keep_female, drop=FALSE]),
			sep='\t', file=paste0(QCed_io_name, "_cat_variables_females.", i, '.tsv'))

		all_sexes_check <- all_sexes[,-c(1,which(names(all_sexes) %in% c(cat_unordered, cat_ordered, cat_mul_unordered)))]
		if(ncol(all_sexes_check) != 0) {
			cat("Error: there are variables that have not been accounted for!\n")
			cat(names(all_sexes_check), "\n")
		}

		rm("all_sexes", "PHESANT_males", "PHESANT_females")

	}
}
