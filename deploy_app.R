library(rsconnect)

rsconnect::deployApp(
  appDir = ".",                  # root folder of your app
  appName = "yEvoMutBrowser2",    # or whatever name you want
  account = "leahmarieanderson"  # your shinyapps.io account name
)
