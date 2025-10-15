suppressWarnings(suppressMessages({
    library(argparser)
    library(tidyverse)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--df", help = "df")
argv <- add_argument(argv, "--spname", help = "spname, default: spname")
argv <- parse_args(argv)

df <- argv$df
spname <- argv$spname

if(!dir.exists(spname)){
    dir.create(spname, recursive = T)
}

# read
df <- read.table(df, sep = ",", header = T)

data <- data[ data$type != "intergenic",]
data$type <- factor(data$type, levels = data$type)

options(repr.plot.height = 5, repr.plot.width = 4)
ggplot(data, aes(x = spname, y= reads, fill=type)) +
    geom_col(position = "fill", color="black", width = 0.3) + 
    theme_classic() +
    ylab("Proportion(% intergenic)") +
    theme(axis.title.x = element_blank(), legend.title = element_blank())

ggsave(str_glue("{spname}/{spname}_intergenic_proportion.png"), height = 5, width = 4)