# Deploy without using renv.lock
# This will use the currently installed packages instead

# First, temporarily rename renv.lock to bypass it
if (file.exists("renv.lock")) {
  file.rename("renv.lock", "renv.lock.backup")
}

# Deploy the app
rsconnect::deployApp(
  appFiles = c(
    "app.R",
    "DESCRIPTION",
    "R/",
    "all_domains.csv",
    "all_yEvo_vcf.csv",
    "chromosome_info.csv",
    "gene_info.csv",
    "motifs.csv",
    "pathways.csv",
    "pfam.csv",
    "prosite.csv",
    "img/",
    "molstar-custom/",
    "static/"
  ),
  forceUpdate = TRUE
)

# Restore renv.lock after deployment
if (file.exists("renv.lock.backup")) {
  file.rename("renv.lock.backup", "renv.lock")
}