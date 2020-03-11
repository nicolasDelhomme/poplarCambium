

# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000087","mitotic M phase",0.004,22.0000,0.962,0.000,"mitotic M phase"),
c("GO:0002252","immune effector process",0.185,5.0000,0.864,0.000,"immune effector process"),
c("GO:0009141","nucleoside triphosphate metabolic process",1.605,300.0000,0.689,0.000,"nucleoside triphosphate metabolism"),
c("GO:0006523","alanine biosynthetic process",0.050,26.0000,0.791,0.273,"nucleoside triphosphate metabolism"),
c("GO:0006021","inositol biosynthetic process",0.016,10.0000,0.817,0.331,"nucleoside triphosphate metabolism"),
c("GO:0030388","fructose 1,6-bisphosphate metabolic process",0.014,25.0000,0.836,0.402,"nucleoside triphosphate metabolism"),
c("GO:0031326","regulation of cellular biosynthetic process",10.816,113.0000,0.808,0.551,"nucleoside triphosphate metabolism"),
c("GO:0006102","isocitrate metabolic process",0.022,10.0000,0.839,0.276,"nucleoside triphosphate metabolism"),
c("GO:0006417","regulation of translation",0.692,119.0000,0.786,0.550,"nucleoside triphosphate metabolism"),
c("GO:0010608","posttranscriptional regulation of gene expression",0.719,124.0000,0.887,0.267,"nucleoside triphosphate metabolism"),
c("GO:0006414","translational elongation",0.777,300.0000,0.836,0.157,"nucleoside triphosphate metabolism"),
c("GO:0006183","GTP biosynthetic process",0.053,144.0000,0.670,0.520,"nucleoside triphosphate metabolism"),
c("GO:0006662","glycerol ether metabolic process",0.122,142.0000,0.847,0.295,"nucleoside triphosphate metabolism"),
c("GO:0022402","cell cycle process",1.053,22.0000,0.798,0.215,"nucleoside triphosphate metabolism"),
c("GO:0009142","nucleoside triphosphate biosynthetic process",0.651,300.0000,0.646,0.662,"nucleoside triphosphate metabolism"),
c("GO:0009635","response to herbicide",0.003,31.0000,0.892,0.000,"response to herbicide"),
c("GO:0010447","response to acidic pH",0.006,14.0000,0.863,0.546,"response to herbicide"),
c("GO:0009268","response to pH",0.011,14.0000,0.860,0.188,"response to herbicide"),
c("GO:0009581","detection of external stimulus",0.057,5.0000,0.850,0.595,"response to herbicide"),
c("GO:0002526","acute inflammatory response",0.016,6.0000,0.852,0.662,"response to herbicide"),
c("GO:0009595","detection of biotic stimulus",0.005,8.0000,0.861,0.676,"response to herbicide"),
c("GO:0009597","detection of virus",0.001,8.0000,0.884,0.632,"response to herbicide"),
c("GO:0009611","response to wounding",0.127,5.0000,0.872,0.437,"response to herbicide"),
c("GO:0009615","response to virus",0.117,8.0000,0.857,0.232,"response to herbicide"),
c("GO:0009416","response to light stimulus",0.157,5.0000,0.836,0.658,"response to herbicide"),
c("GO:0006954","inflammatory response",0.110,6.0000,0.851,0.341,"response to herbicide"),
c("GO:0032505","reproduction of a single-celled organism",0.181,14.0000,0.955,0.000,"reproduction of a single-celled organism"),
c("GO:0034220","ion transmembrane transport",3.528,300.0000,0.912,0.000,"ion transmembrane transport"),
c("GO:0055085","transmembrane transport",8.916,300.0000,0.932,0.497,"ion transmembrane transport"),
c("GO:0006818","hydrogen transport",1.149,300.0000,0.866,0.366,"ion transmembrane transport"),
c("GO:0015797","mannitol transport",0.005,11.0000,0.903,0.377,"ion transmembrane transport"),
c("GO:0030255","protein secretion by the type IV secretion system",0.013,18.0000,0.848,0.404,"ion transmembrane transport"),
c("GO:0042273","ribosomal large subunit biogenesis",0.151,297.0000,0.887,0.000,"ribosomal large subunit biogenesis"),
c("GO:0034728","nucleosome organization",0.129,31.0000,0.859,0.342,"ribosomal large subunit biogenesis"),
c("GO:0042255","ribosome assembly",0.164,297.0000,0.833,0.649,"ribosomal large subunit biogenesis"),
c("GO:0051262","protein tetramerization",0.044,84.0000,0.864,0.657,"ribosomal large subunit biogenesis"),
c("GO:0051259","protein oligomerization",0.188,83.0000,0.852,0.652,"ribosomal large subunit biogenesis"),
c("GO:0052126","movement in host environment",0.022,8.0000,0.891,0.000,"movement in host environment"),
c("GO:0042026","protein refolding",0.069,204.0000,0.957,0.029,"protein refolding"),
c("GO:0009853","photorespiration",0.005,11.0000,0.927,0.042,"photorespiration"),
c("GO:0019684","photosynthesis, light reaction",0.069,86.0000,0.852,0.052,"photosynthesis, light reaction"),
c("GO:0043467","regulation of generation of precursor metabolites and energy",0.030,9.0000,0.859,0.502,"photosynthesis, light reaction"),
c("GO:0019646","aerobic electron transport chain",0.042,14.0000,0.829,0.513,"photosynthesis, light reaction"),
c("GO:0006808","regulation of nitrogen utilization",0.079,88.0000,0.890,0.072,"regulation of nitrogen utilization"),
c("GO:0006879","cellular iron ion homeostasis",0.110,84.0000,0.808,0.180,"regulation of nitrogen utilization"),
c("GO:0050821","protein stabilization",0.045,7.0000,0.919,0.450,"regulation of nitrogen utilization"),
c("GO:0031647","regulation of protein stability",0.070,7.0000,0.917,0.464,"regulation of nitrogen utilization"),
c("GO:0019740","nitrogen utilization",0.085,88.0000,0.934,0.073,"nitrogen utilization"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "REVIGO Gene Ontology treemap",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
