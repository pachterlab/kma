## ------------------------------------------------------------------------
system.file("pre-process", package="kma")

## ----, eval=FALSE--------------------------------------------------------
#  devtools::install_github("http://github.com/pachterlab/kma")

## ----, eval=FALSE--------------------------------------------------------
#  system.file("pre-process", package = "kma")
#  

## ------------------------------------------------------------------------
base_dir <- system.file("example", package="kma")
xprs_fnames <- Sys.glob(file.path(base_dir, "experiment/*/*/xprs_out/results.xprs"))
xprs_fnames

## ------------------------------------------------------------------------
sample_names <- sub(file.path(base_dir, "experiment/[a-z]+/"), "", xprs_fnames) %>%
    sub("xprs_out/results.xprs", "", .) %>%
    gsub("/", "", .)
sample_names

## ------------------------------------------------------------------------
condition_names <- sub("[0-9]+", "", sample_names)
condition_names

## ------------------------------------------------------------------------
xprs <- read_express(xprs_fnames, sample_names, condition_names)

## ------------------------------------------------------------------------
names(xprs)

## ------------------------------------------------------------------------
intron_to_trans <- data.table::fread(file.path(base_dir, "kma_pre-process_out",
    "intron_to_transcripts.txt"), data.table = FALSE)
head(intron_to_trans)

## ------------------------------------------------------------------------
ir <- newIntronRetention(xprs$tpm, intron_to_trans, xprs$condition,
    xprs$uniq_counts)

## ------------------------------------------------------------------------
print(ir)

## ------------------------------------------------------------------------
ir <- ir %>%
    filter_low_tpm(1) %>%
    filter_perfect_psi() %>%
    filter_low_frags(3)
colnames(ir$flat)

## ------------------------------------------------------------------------
zc_fnames <- Sys.glob(file.path(base_dir, "experiment/*/*/zero_coverage.txt"))
zc_samples <- sub(file.path(base_dir, "experiment/[a-z]+/"), "", zc_fnames) %>%
    sub("zero_coverage.txt", "", .) %>%
    gsub("/", "", .)
zc_conditions <- sub("[0-9]+", "", zc_samples)
all_zc <- get_batch_intron_zc(zc_fnames, zc_samples, zc_conditions)
head(all_zc)

## ------------------------------------------------------------------------
ir <- summarize_zero_coverage(ir, all_zc)

## ------------------------------------------------------------------------
colnames(ir$flat)

## ------------------------------------------------------------------------
ir_test <- retention_test(ir)
head(ir_test)

## ------------------------------------------------------------------------
ir_test %>%
    filter(qvalue <= 0.10) %>%
    select(-c(pvalue))

