library(rmarkdown)
render("LRR_structure.Rmd", output_format = "html_document", output_file = "$(pwd)/LRR_structure_${NAME}.html", params = list(wd="$(pwd)"))
