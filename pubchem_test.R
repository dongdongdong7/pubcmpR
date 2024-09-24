# devtools::install_github("selcukorkmaz/PubChemR")
library(magrittr)
compound_info <- PubChemR::get_compounds(c(2244, 305))
page <- xml2::read_html("https://pubchem.ncbi.nlm.nih.gov/compound/305")
r <- httr::GET("https://pubchem.ncbi.nlm.nih.gov/compound/305")
page <- rawToChar(r$content)
page <- rvest::read_html(page)
page %>% rvest::html_element(css = "body") %>%
  rvest::html_elements(css = "div#root") %>%
  rvest::html_text()

install.packages("RSelenium")

xml2::xml_text(html)
# 静态html没有想要的信息
library(RSelenium)
rD <- rsDriver()
remDr <- rD[["client"]]
remDr$navigate("http://www.google.com/ncr")
remDr$navigate("http://www.bbc.com")
remDr$close()
# stop the selenium server
rD[["server"]]$stop()


sdf <- ChemmineR::smiles2sdf("CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O")
ChemmineR::plot(sdf[1])

