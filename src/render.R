library(rmarkdown)
args = commandArgs(trailingOnly=TRUE)
render("LRR_structure.Rmd", output_format = "html_document", output_file = "LRR_structure.html", params = list(wd=args[1]))
