library(rmarkdown)
args = commandArgs(trailingOnly=TRUE)

#report_path <- tempfile(fileext = ".Rmd")
#file.copy("LRR_structure.Rmd", report_path, overwrite = TRUE)

render("LRR_structure.Rmd", output_format = "html_document", output_file = "LRR_structure.html", params = list(wd=args[1]))
