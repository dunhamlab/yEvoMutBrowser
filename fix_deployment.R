# Fix deployment issues by updating renv and removing deprecated packages

# Install remotes if missing
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Restore renv environment
renv::restore()

# Update all packages
renv::update()

# Remove plogr from the lockfile if it exists
lockfile <- renv::lockfile_read()
if ("plogr" %in% names(lockfile$Packages)) {
  lockfile$Packages$plogr <- NULL
  renv::lockfile_write(lockfile)
}

# Snapshot the current state without deprecated packages
renv::snapshot()

# Try deployment again
rsconnect::deployApp()