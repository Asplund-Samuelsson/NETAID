#!/usr/bin/Rscript

# Load command line arguments
args = commandArgs(trailingOnly=T)
infiles = head(args, length(args)-1)
outfile = tail(args, 1)

# Check that the outfile is pdf
if (!endsWith(outfile, "pdf")) {
  message("Error: Outfile must end with 'pdf'.")
  quit("no")
}

# Read data
for (infile in infiles) {
  cd_new = read.table(infile, stringsAsFactors=F, quote="", sep="\t", header=T, fill=T)
  if (exists("cd")) {
    cd = rbind(cd, cd_new)
  } else {
    cd = cd_new
  }
}

# Remove the ID column
cd = cd[,grep("^ID$", colnames(cd), invert=T)]

# Melt data
library(reshape2)
cd = melt(cd)

# Change from mM to µM
cd$value = cd$value / 1000

# Create metabolite tag
cd$Metabolite = paste(
  cd$Name, rep(" [", nrow(cd)), cd$Compartment, rep("]", nrow(cd)), sep=""
)

# Create measured vs NET tag
cd$Source = ifelse(cd$variable %in% c("LowIn", "HighIn"), "Bounds", "NET-opt.")
cd$variable = ifelse(cd$variable %in% c("LowIn", "LowOut"), "Lower", "Upper")

# Cast data
cd_lower = subset(cd, variable == "Lower")[c("Metabolite", "Label", "Source", "value")]
colnames(cd_lower)[length(colnames(cd_lower))] = "Lower"

cd_upper = subset(cd, variable == "Upper")[c("Metabolite", "Label", "Source", "value")]
colnames(cd_upper)[length(colnames(cd_upper))] = "Upper"

cd = merge(cd_lower, cd_upper)

# Plot the data
library(ggplot2)
library(RColorBrewer)

gp = ggplot(cd, aes(x=Label, ymin=Lower, ymax=Upper, color=Label, linetype=Source))
gp = gp + geom_errorbar(position = "dodge", width = 0.6)
gp = gp + geom_hline(yintercept=1e-7, linetype="dashed", color="grey")
gp = gp + geom_hline(yintercept=0.1, linetype="dashed", color="grey")
gp = gp + facet_wrap("Metabolite", ncol=3)
gp = gp + theme_bw()
gp = gp + scale_y_log10()
gp = gp + scale_color_manual(values=brewer.pal(8, 'Paired'), guide=guide_legend(reverse=TRUE))
gp = gp + ylab("Concentration range (M)")
gp = gp + theme(axis.title.x=element_blank())
gp = gp + theme(strip.text.x = element_text(size=7))
gp = gp + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gp = gp + coord_flip()

# Calculate optimal height of plot
height_mm = 20 + (10 + 5 * length(unique(cd$Label))) * ceiling(length(unique(cd$Metabolite)) / 3)

ggsave(outfile, gp, width=210/25.4, height=height_mm/25.4, limitsize=F)
