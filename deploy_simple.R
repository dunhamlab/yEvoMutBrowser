# Simple deployment without renv complications

# Move renv.lock temporarily
if (file.exists("renv.lock")) {
  message("Moving renv.lock temporarily...")
  file.rename("renv.lock", "renv.lock.temp")
}

# Deploy
tryCatch({
  rsconnect::deployApp()
}, finally = {
  # Restore renv.lock
  if (file.exists("renv.lock.temp")) {
    message("Restoring renv.lock...")
    file.rename("renv.lock.temp", "renv.lock")
  }
})