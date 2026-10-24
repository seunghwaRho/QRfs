# adrian.R: builds data/adrian.rda; run from the package root (QRfs/).
# The numbers are those embedded in qrfs_adrian_data() of the reviewed library (bit-identical to
# DataVulnerability.xls of the authors' replication package, openICPSR 113169, CC BY 4.0; checked
# 2026-09-21). The library file is sourced from the review repository.
lib <- "/Users/shr/Dropbox/QuantileFixedb/codeReviewFinal/library/qrfs_functions_reduced.R"
env <- new.env(); source(lib, local = env)
adrian <- env$qrfs_adrian_data()
stopifnot(is.data.frame(adrian), nrow(adrian) == 171L,
          identical(names(adrian), c("quarter", "ggdp_h1", "ggdp", "NFCI")))
save(adrian, file = "data/adrian.rda", compress = "xz", version = 2)
