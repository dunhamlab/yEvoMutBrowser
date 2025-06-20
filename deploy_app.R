library(rsconnect)

rsconnect::deployApp(
  appDir = "/Users/leahanderson/Documents/GitHub/yEvoMutBrowser",                  # root folder of your app
  appName = "yEvoMutBrowser2",    # or whatever name you want
  account = "leahmarieanderson"  # your shinyapps.io account name
)
