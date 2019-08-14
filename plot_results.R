# load libraries
library(ggplot2)
library(readr)
library(cowplot)

# set working directory
setwd('~/Desktop/FGM/1nFGM4conflict/')

# get files
# first file is conflict and standard results
# second file is abiotic change
# mutation size and amount of conflict are used to name figure
args = commandArgs(trailingOnly=TRUE)
df1 <- read_csv(args[1])
df2 <- read_csv(args[2])
mut_size <- args[3]
conflict <- args[4]

# merge dataframes
df <- rbind(df1[, -7], df2)
# plot position and fitness over time and save figure to directory
# Party 2 is exluded from first plot because z is the same as that of party 1
a <- ggplot(df[!df$Population == 'Party 2', ], aes(Iteration,z,color = Population)) + geom_line()
b <- ggplot(df, aes(Iteration,Fitness,color = Population)) + geom_line()
plot_grid(a,b,nrow = 2, ncol = 1)
ggsave(paste(mut_size, conflict, 'plot.png', sep = '_'), width = 5, height = 5)


