#defining some variables which will be used throughout the script
#this is the working directory.
# VARIABLES ####
path = "C:\\Users\\2873826G\\Documents\\My Masters Degree UofG\\Data Exploration\\Assignment"
p_threshold = 0.001
fold_threshold = 2
sample_colours = c("pink4","navyblue","limegreen")  #these are the colours which will be used to color the three sample groups
colour1 = "pink4"
colour2 = "navyblue"
colour3 = "limegreen"
default_plot_width = (763/96)
default_plot_height = (370/96) #using standard export height and width, these are in pixels in R so dividing by 96 to convert to inches, which is what the functions use.
#for consistency, we'll think of proliferating to be sample group I, senescent to be sample group II, and senescent-mitochondria-depleted to be sample group III.

# LIBRARIES ####
#loading all the required libraries
library(ggplot2)
library(ggrepel)
library(reshape2)
library(amap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db) #even though the dataset is human, some functions are written to accept an "organism" argument and run using either human or mouse as specified.
library(STRINGdb)
library(devEMF)
library(grDevices)
library(ggeasy)

# THEME ####
#defining a theme for consistency in plots
customized_theme = theme(plot.title = element_text(size = 20),axis.text = element_text(size=10),axis.title = element_text(size = 15),legend.text = element_text(size = 10),legend.position = "bottom")

# FUNCTIONS ####
#creating all required functions
DE_table_creator = function(de_filepath,em_annotated_tablename,separator,samples,p_threshold,fold_threshold)
{
  #loading it
  if(separator == "tab") {DE_table = read.table(de_filepath,header=TRUE, row.names=1, sep= "\t")}
  if(separator == "comma") {DE_table = read.table(de_filepath,header=TRUE, row.names=1)}
  
  #merging it
  DE_table = merge(em_annotated_tablename,DE_table,by.x=0,by.y=0)
  
  #renaming first column as ENSEMBL
  names(DE_table)[1] = "ENSEMBL"
  
  #adding a sig column
  DE_table$sig = as.factor(DE_table$p.adj < p_threshold & abs(DE_table$log2fold) > fold_threshold)
  
  #adding mlog10p column
  DE_table$mlog10p = -log10(DE_table$p.adj)
  
  #adding a column containing mean expression values for each gene across all replicates
  DE_table$mean_exp = rowMeans(DE_table[,samples])
  
  #setting row names to be gene symbols
  row.names(DE_table) = DE_table$symbol
  
  #removing NAs
  DE_table = na.omit(DE_table)
  
  #ordering the table by p.adj value
  sorted_order = order(DE_table[,"p.adj"], decreasing=FALSE)
  DE_table = DE_table[sorted_order,]
  
  #removing the duplicate ENSEMBL column 
  DE_table = DE_table[,-2]
  
  #returning the table
  return(DE_table)
}

#creating a function to return a vector containing the significant genes from parsed DE table
sig_genes_lister = function(parsed_DE_table,p_threshold,fold_threshold) 
{
  return(row.names(subset(parsed_DE_table,p.adj < p_threshold & abs(log2fold) > fold_threshold)))
}

#creating volcano plot function that just takes the DE table as actual data and the rest of the arguments are parameters.
volcano_plotter = function(DE_table,p_threshold = 0.05,fold_threshold = 1,non_sig_colour = "black",up_colour = "red",down_colour = "blue",plot_title="",x_axis_label="Log2 of Fold Change",y_axis_label="Minus Log10 of p Adjusted Value",threshold_linetype = "dashed",threshold_linecolour = "gray",threshold_linesize = 0.5,xlimits = c(-20,20), ylimits = c(0,100),down_labels = 5, up_labels = 5) 
{
  #creating a subset containing significant genes only
  DE_sig = subset(DE_table, p.adj < p_threshold & abs(log2fold) > fold_threshold)
  
  #creating two separate tables for up-regulated and down-regulated genes (significant only)
  DE_up = subset(DE_sig, log2fold > 0)
  DE_down = subset(DE_sig, log2fold < 0)
  
  #sorting these tables by p value
  DE_up = DE_up[order(DE_up$p.adj,decreasing = FALSE),]
  DE_down = DE_down[order(DE_down$p.adj,decreasing = FALSE),]
  
  #selecting the top genes for each
  DE_up_top = DE_up[1:(up_labels),]
  DE_down_top = DE_down[1:(down_labels),]
  
  #loading the necessary libraries
  library(ggplot2)
  library(ggrepel)
  
  #creating the plot object
  ggp = ggplot(DE_table,aes(x=log2fold,y=mlog10p)) + 
    
    #plotting the base layer
    geom_point(aes(colour="a")) +
    
    #plotting the colored up-regulated significant genes
    geom_point(aes(colour="b"),data=DE_up) +
    
    #plotting the colored down-regulated significant genes
    geom_point(aes(colour="c"),data=DE_down) +
    
    #labelling
    labs(title=plot_title,x=x_axis_label,y=y_axis_label) +
    
    #adding themes
    theme_classic() +
    customized_theme + 
    
    #drawing thresholds
    geom_vline(xintercept= -(fold_threshold),linetype=threshold_linetype, color = threshold_linecolour, size=threshold_linesize) + 
    geom_vline(xintercept= fold_threshold,linetype=threshold_linetype, color = threshold_linecolour, size=threshold_linesize) + 
    geom_hline(yintercept=-log10(p_threshold),linetype=threshold_linetype, color = threshold_linecolour, size=threshold_linesize) +
    
    #labelling the top 10 most significant genes
    geom_text_repel(data=DE_up_top, aes(label=symbol),colour=up_colour,show.legend = FALSE) +
    geom_text_repel(data=DE_down_top, aes(label=symbol),colour=down_colour,show.legend = FALSE) +
    
    #adding colours to the layers and labels to the legend appropriately
    scale_colour_manual(values = c(non_sig_colour, up_colour, down_colour), labels=c("Non-Significant","Up-Regulated","Down-Regulated")) + 
    
    #adding limits on axes
    xlim(xlimits) +
    ylim(ylimits)
  
  #the plot is ready, returning this
  return(ggp)
}

#creating function for MA plot
MA_plotter = function(DE_table,p_threshold=0.05,fold_threshold=1,non_sig_colour="black",sig_colour="violetred",threshold_linetype = "dashed",threshold_linecolour = "gray",threshold_linesize = 0.5,plot_title="",x_axis_label="Log10 Mean Expression",y_axis_label="Log2 Fold Change",y_axis_limits=c(-50,50))
{
  #creating a subset containing significant genes only
  DE_sig = subset(DE_table, p.adj < p_threshold & abs(log2fold) > fold_threshold)
  
  #loading necessary library
  library(ggplot2)
  
  #tables are ready, plotting
  ggp = ggplot(DE_table,aes(x=log10(mean_exp),y=log2fold)) +
    
    #applying the first layer
    geom_point(aes(colour="a")) +
    
    #coloring dots by significance
    geom_point(data = DE_sig,aes(colour="b")) +
    
    #adding fold change thresholds
    geom_hline(yintercept= fold_threshold,linetype=threshold_linetype, color = threshold_linecolour, size=threshold_linesize) +
    geom_hline(yintercept= -(fold_threshold),linetype=threshold_linetype, color = threshold_linecolour, size=threshold_linesize) +
    
    #labelling
    labs(title=plot_title,x=x_axis_label,y=y_axis_label) +
    
    #applying themes
    theme_classic() + customized_theme +
    
    #adding limits
    ylim(y_axis_limits) +
    
    #adding appropriate labels and labelling the legend appropriately.
    scale_colour_manual(values = c(non_sig_colour, sig_colour), labels=c("Non-Significant","Significant"))
  
  #our plot is now ready.
  return(ggp)
}

#making a PCA plotter function
PCA_plotter = function(em_symbols_table,ss_table, colour1,colour2,colour3,point_size=2,plot_title="") 
{
  #first we do PCA 
  pca =  prcomp(t(as.matrix(sapply(em_symbols_table, as.numeric)))) #first em_symbols is converted into a numeric dataframe, then it is converted to a matrix, then it is transposed, then the PCA is done.
  #extracting the co-ordinates
  pca_coordinates = data.frame(pca$x)
  #adding variance percentage to the axes
  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  x__label = paste("PC1 ", " (",prop_x, "%)",sep="")
  y__label = paste("PC2 ", " (",prop_y, "%)",sep="")
  
  #plotting 
  library(ggplot2)
  ggp = ggplot(pca_coordinates,aes(x=PC1, y=PC2)) +
    #plotting the points
    geom_point(aes(colour = ss_table$SAMPLE_GROUP),size=point_size) +
    #defining color groups
    scale_color_manual(values = c(colour1,colour2,colour3)) +
    #adding text labelling each sample
    geom_text_repel(aes(label=ss_table$SAMPLE)) +
    #labelling
    labs(title = plot_title,x=x__label,y=y__label) +
    #applying themes
    theme_classic() + customized_theme
  
  #the plot is ready.
  return(ggp)
}

#creating a function for expression matrices
expression_matrices_plotter = function(em_table, colour1,colour2 ,colour3 ,alpha=0.5,plot_title = "", x_axis_label = "Log10 of Expression")
{
  library(reshape2)
  #melting em table
  em.m = melt(em_table)
  
  #creating the density plot
  ggp = ggplot(em.m,aes(x=log10(value+0.01))) +
    #making the density plot
    geom_density(alpha=alpha,aes(fill = variable)) +
    #making facets
    facet_wrap(~variable, ncol=ncol(em_table)) +
    #adding themes
    theme(panel.grid = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5),strip.text.x = element_text(size = 14, family="Arial", face="bold", vjust=1),panel.spacing = unit(1, "lines"),legend.position = "none") + customized_theme +
    #labelling axes
    labs(title=plot_title,x=x_axis_label) + 
    #adding colors 
    scale_fill_manual(values =  rep(c(colour1,colour2,colour3),each = 3))
  
  #plot is ready
  return(ggp)
}

#making boxplot functions

#for an individual boxplot
single_boxplot_plotter = function(em_table,gene_of_interest,sample_group_vector,factor_levels,alpha=0.5,boxplot_size = 1,colour1 ,colour2,colour3,median_line_size=1.5,jitter_colour = "black",plot_title="",x_axis_label="Sample Groups",y_axis_label="Expression Value")
{
  #creating the needed table for the boxplot
  gene_table = em_table[gene_of_interest,]
  #transposing it
  gene_table = data.frame(t(gene_table))
  #appending sample group information
  gene_table$sample_group = sample_group_vector
  #renaming the first column to 'expression_value'
  names(gene_table)[1] = "expression_value"
  #reordering levels
  gene_table$sample_group = factor(gene_table$sample_group,levels = factor_levels)
  #table is ready
  
  #plotting
  library(ggplot2)
  ggp = ggplot(gene_table,aes(x=sample_group,y=expression_value,colour=sample_group),alpha = alpha) +
    #creating boxplot
    geom_boxplot(size = boxplot_size,fatten=median_line_size) +
    #applying themes
    theme_classic() + customized_theme +
    #defining colors
    scale_colour_manual(values = c(colour1,colour2,colour3)) +
    #adding the jitter plot
    geom_jitter(colour = jitter_colour) +
    #adding labels
    labs(title = plot_title,x=x_axis_label,y=y_axis_label)
  
  #plot is ready
  return(ggp)
}

#for faceted boxplots
faceted_boxplot_plotter = function(em_scaled_table,candidate_genes_list,sample_group_vector,factor_levels = "default",alpha=0.5,boxplot_size=1,median_line_size=1,colour1 ,colour2,colour3,plot_title="",x_axis_label="Sample Groups",y_axis_label="Expression Value") 
{
  #creating a table for boxplots 
  candidate_genes_table = em_scaled_table[candidate_genes_list,]
  #transposing and adding sample information
  candidate_genes_table = data.frame(t(candidate_genes_table))
  #adding sample group data
  candidate_genes_table$sample_group = sample_group_vector
  #melting the table
  library(reshape2)
  candidate_genes_table.m = melt(candidate_genes_table)
  #setting levels of the melted table, if default is overrided
  if(factor_levels != "default") {candidate_genes_table.m$sample_group = factor(candidate_genes_table.m$sample_group,levels = factor_levels)}
  
  #plotting
  #the logic is very similar to the individual box-plotting logic
  library(ggplot2)
  ggp = ggplot(candidate_genes_table.m,aes(x=sample_group,y=value,facet = variable,colour=sample_group),alpha = alpha) + 
    geom_boxplot(size = boxplot_size,fatten=median_line_size) + 
    theme_classic() + customized_theme +
    geom_jitter(show.legend = FALSE) +
    scale_colour_manual(values = c(colour1,colour2,colour3)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    facet_wrap(~variable, ncol=length(candidate_genes_list)) +
    #adding labels
    labs(title = plot_title,x=x_axis_label,y=y_axis_label)
  
  return(ggp)
}

#making a heatmap function
heatmap_plotter = function(em_scaled_table,sig_genes_vector,y_cluster = TRUE,x_cluster = FALSE,correlation_method="spearman",clustering_method="average",down_colour="royalblue3",up_colour="red",plot_title="",x_axis_label="Samples",y_axis_label="Genes")
{
  #making a table to create a heatmap
  library(reshape2)
  library(amap)
  #creating a matrix out of the expression table containing significant genes. omitting NAs. 
  hm.matrix = as.matrix(na.omit(em_scaled_table[sig_genes_vector,]))
  
  if(y_cluster == TRUE)
  { 
    #clustering by Y-axis
    #calculating correlation by spearman method
    y.dist = Dist(hm.matrix, method=correlation_method)
    #performing clustering/building the dendrogram
    y.cluster = hclust(y.dist, method=clustering_method)
    #extracting the dendrogram
    y.dd = as.dendrogram(y.cluster)
    #reordering the dendrogram
    y.dd.reorder = reorder(y.dd,0,FUN="average")
    #getting the clustered row indexes
    y.order = order.dendrogram(y.dd.reorder)
    #re-ordering the heatmap matrix according to the calculated clustered indexes
    hm.matrix = hm.matrix[y.order,] #by Y-axis 
  }
  
  if(x_cluster == TRUE)
  { 
    #clustering by X-axis #same logic
    x.dist = Dist(t(na.omit(em_scaled_table[sig_genes_vector,])), method=clustering_method)
    x.cluster = hclust(x.dist, method=clustering_method)
    x.dd = as.dendrogram(x.cluster)
    x.dd.reorder = reorder(x.dd,0,FUN="average")
    x.order = order.dendrogram(x.dd.reorder)
    #re-ordering the heatmap matrix according to the calculated clustered indexes
    hm.matrix = hm.matrix[,x.order] #by X-axis
  }
  
  #melting the table
  hm.matrix = melt(hm.matrix)
  #creating a color palette from blue to red (i.e. down to up)
  colours = c(down_colour,up_colour)
  
  #creating the plot
  #the x axis are the samples, the y axis is the gene, and the tile colour corresponds to the expression value. lower is bluer, higher is redder, as is the convention.
  ggp = ggplot(hm.matrix,aes(x=Var2,y=Var1,fill=value)) + 
    #the heatmap is plotted
    geom_tile() + 
    #the colors are filled appropriately using a gradient created from pre-defined colours  
    scale_fill_gradientn(colours = colorRampPalette(colours)(500)) +
    #automatic Y axis label is removed for the sake of neatness.
    ylab("") + 
    #themes are added  
    theme_classic() + customized_theme + 
    #tilting the x_axis labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #theme is further modified for neatness
    theme(axis.text.y = element_blank(), axis.ticks=element_blank(), legend.title = element_blank(), legend.spacing.x = unit(0.25, 'cm')) +
    #adding labels
    labs(title = plot_title,x=x_axis_label,y=y_axis_label)
  
  #returning the plot
  return(ggp)
}

#making a rug plotting function
rug_plotter = function(sample_group_vector,rug_colours,factor_levels = "default")
{
  library(reshape2)
  library(ggplot2)
  library(ggeasy)
  groups_data = as.matrix(as.numeric(as.factor(sample_group_vector)))
  groups_data = melt(groups_data)
  # heatmap
  #resetting levels to match plots, if factor_levels default is overrided.
  if(factor_levels != "default") {groups_data$value = factor(groups_data$value,levels = factor_levels)}
  #creating plot
  ggp = ggplot(groups_data, aes(y=Var2,x=Var1,fill=value)) + geom_tile(linetype="blank") + scale_fill_gradientn(colours = rug_colours) + labs(x = "", y = "") + theme_classic() + theme(legend.position="none", legend.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks=element_blank()) + easy_remove_axes() #there is nothing in this plot except the actual rug itself, i.e. no labels etc.
  return(ggp)
}

#making function to create an ORA object. 
ORA_calculator = function(master_table,p_threshold = 0.05,fold_threshold = 1,gene_name_type = "SYMBOL",ont="BP",q_value_cutoff = 0.1,analysis="all",organism = "human")
{
  #subsetting the table
  master_sig = subset(master_table,p.adj < p_threshold & abs(log2fold) > fold_threshold)
  
  #loading necessary libraries
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  
  #using the correct dictionary
  if(organism == "mouse") {dictionary = org.Mm.eg.db}
  if(organism == "human") {dictionary = org.Hs.eg.db}
  
  #if we want to analyse all sig genes - 
  if(analysis=="all") 
  {
    #creating a vector containing ENTREZ IDs of sig genes
    sig_genes_entrez = bitr(row.names(master_sig), fromType = gene_name_type, toType = "ENTREZID", OrgDb = dictionary)
    #running over-representation analysis
    ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = dictionary, readable = T, ont = ont, pvalueCutoff = p_threshold, qvalueCutoff = q_value_cutoff)
  }
  
  #if we want to analyse up-regulated sig genes only
  if(analysis=="up")
  {
    sig_genes = row.names(subset(master_sig,log2fold>0))
    #creating a vector containing ENTREZ IDs of sig genes
    sig_genes_entrez = bitr(sig_genes, fromType = gene_name_type, toType = "ENTREZID", OrgDb = dictionary)
    #running over-representation analysis
    ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = dictionary, readable = T, ont = ont, pvalueCutoff = p_threshold, qvalueCutoff = q_value_cutoff)
  }
  
  #if we want to analyse down-regulated sig genes only
  if(analysis=="down")
  {
    sig_genes = row.names(subset(master_sig,log2fold<0))
    #creating a vector containing ENTREZ IDs of sig genes
    sig_genes_entrez = bitr(sig_genes, fromType = gene_name_type, toType = "ENTREZID", OrgDb = dictionary)
    #running over-representation analysis
    ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = dictionary, readable = T, ont = ont, pvalueCutoff = p_threshold, qvalueCutoff = q_value_cutoff)
  }
  
  
  #returning the ora results object
  return(ora_results)
  
}

#creating an ora barplot plotter, using barplot only, the information depicted is redundant for the rest so we don't need all the plots in any case.
ORA_barplot_plotter = function(ora_results_object,how_many_genesets_to_plot=10)
{
  #plotting
  ggp = barplot(ora_results_object, showCategory=how_many_genesets_to_plot)
  return(ggp)
}

#creating function to organize ORA results
ora_results_table_creator = function(ora_results_object)
{
  gene_sets = ora_results_object$geneID
  description = ora_results_object$Description
  p.adj = ora_results_object$p.adjust
  #creating a combined table
  ora_results_table = data.frame(cbind(gene_sets,p.adj))
  row.names(ora_results_table) = description
  
  #returning this
  return(ora_results_table)
}

#writing a function to return candidate genes from an ora results object
ora_candidate_genes_lister = function(ora_results_object,which_ontology = 1) #which_ontology parameter set to 1 will select the most enriched ontology, set to 2 would select the second most enriched ontology, and so on.
{
  table = ora_results_table_creator(ora_results_object)
  enriched_gene_set = as.character(table [(which_ontology),1])
  candidate_genes = unlist(strsplit(enriched_gene_set, "/"))
  
  return(candidate_genes)
}

#writing a function to select most significantly differentially expressed genes in an ORA gene set for its respective DE comparison
candidate_genes_selector = function(DE_table,genes_of_interest,number_to_select = 10)
{
  DE_sig = DE_table[genes_of_interest,]
  #ordering this
  sorted_order = order(DE_sig[,"p.adj"], decreasing=FALSE)
  DE_sig = DE_sig[sorted_order,]
  #selecting the top genes
  DE_top = DE_sig[1:number_to_select,]
  
  #returning the names of genes selected
  return(row.names(DE_top))
}

#writing a function to create a string network
string_network_plotter = function(list_of_candidate_genes,organism = "mouse")
{
  #loading stringdb library
  library(STRINGdb)
  #creating data frame containing a column for the genes and calling that column "gene" (this is compulsory)
  candidate_genes_table = data.frame(list_of_candidate_genes)
  names(candidate_genes_table)[1] = "gene"
  
  #creating variable to store species information
  if(organism == "mouse") {species = 10090}
  if(organism == "human") {species = 9606}
  
  # load the database
  string_db = STRINGdb$new( version="11.5", species=species, score_threshold=200, network_type="full", input_directory="")
  # map our genes to the database
  string_mapped = string_db$map(candidate_genes_table, "gene", removeUnmappedRows = TRUE )
  # now we plot
  return(string_db$plot_network(string_mapped))
}

#creating a gene to plot metagene boxplot
metagene_boxplot_plotter = function(em_scaled_table,genes_of_interest,ss_table,colours,alpha=0.5,boxplot_size = 1,median_line_size=1.5,jitter_colour = "black",plot_title="Metagene Plot",x_axis_label="Sample Groups",y_axis_label="Expression Value")
{
  #creating the required table
  em_scaled_table = na.omit(em_scaled_table)
  em_scaled_table = na.omit(em_scaled_table[genes_of_interest,])
  metagene = data.frame(colMeans(em_scaled_table))
  metagene$sample_group = ss_table$SAMPLE_GROUP
  #renaming column 1 to something sensible
  names(metagene)[1] = "expression_value"
  
  #creating the boxplot
  library(ggplot2)
  ggp = ggplot(metagene,aes(x=sample_group,y=expression_value,colour=sample_group),alpha = alpha) +
    #creating boxplot
    geom_boxplot(size = boxplot_size,fatten=median_line_size) +
    #applying themes
    theme_classic() + customized_theme +
    #defining colors
    scale_colour_manual(values = colours) +
    #adding the jitter plot
    geom_jitter(colour = jitter_colour,width = 0.1, height = 0.1) +
    #adding labels
    labs(title = plot_title,x=x_axis_label,y=y_axis_label)
  
  #returning the plot
  return(ggp)
}

#creating a function to plot fold vs fold plot
fold_vs_fold_scatterplot_plotter = function(master,x_axis_label,y_axis_label,plot_title = "Fold Vs Fold Plot (Cor = ", dot_colour = "black")
{
  #testing the correlation
  cor_value = cor(master$log2fold.smdvs,master$log2fold.svp,method = "spearman")

  #adding it to plot title
  plot_title = paste(plot_title,cor_value,")",sep="")
  
  #creating the plot
  ggp = ggplot(master,aes(x=log2fold.svp,y=log2fold.smdvs),colour = dot_colour) +
  #creating the scatterplot
  geom_point() +
  #adding a theme
  theme_classic() + customized_theme +
  #adding labels
  labs(title = plot_title,x=x_axis_label,y=y_axis_label)
  
  #returning the plot
  return(ggp)
}

#creating a function to save plots
plot_saver = function(path,desired_filename,plot_height,plot_width,plot_object,colour_type)
{
  library(devEMF)
  library(ggplot2)
  library(grDevices)
  #creating the filename
  path = paste(path,"\\",desired_filename,".emf",sep="")
  #different functions create emf files with different efficiency, e.g. the emf() function can plot many different types of colours but cannot create gradients with colors near the edges of the palette (so it has trouble creating heatmaps), on the other hand, savePlot() can handle gradients with edge colours (creating high contrast heatmaps), but not many different types of colors (it was writing white-colored expression density plots for me). for this reason, I have included a "colortype" metric, which calls the desired function.
  if(colour_type == "block"){
  emf(file = path, width = plot_width, height = plot_height)
  print(plot_object)
  # closing the device
  dev.off()
  }
  
  if(colour_type == "gradient"){
    # creating new device so we can modify dimensions
    dev.new(width=plot_width, height=plot_height)
    # Redrawing the plot 
    print(plot_object)
    # saving
    savePlot(path, type = "emf", device = dev.cur())
    # closing the device
    dev.off()
  }
}

# TABLES ####
#loading all the tables from the directory
em = read.table((paste(path,"\\EM.csv",sep = "")), header=TRUE, row.names=1, sep= "\t") #this is the expressions matrix
annotations = read.table((paste(path,"\\Human_Background_GRCh38.p13.csv",sep = "")),header=TRUE, row.names=1, sep= "\t")
#changing annotations table column names
new_names = c("symbol","chromosome","start","stop","type")
names(annotations) = new_names
rm(new_names) #removing variables not needed anymore to keep environment clean 
ss = read.table((paste(path,"\\sample_sheet.csv",sep = "")), header=TRUE, sep= "\t") #this is the sample sheet.

#merging em and annotations
em_annotated = merge(em, annotations, by.x=0, by.y=0)
row.names(em_annotated) = em_annotated$Row.names
names(em_annotated)[1] = "ENSEMBL"

de_senes_v_prolif = DE_table_creator(de_filepath = paste(path,"\\DE_Senes_vs_Prolif.csv",sep = ""),em_annotated_tablename = em_annotated,separator = "tab",p_threshold = p_threshold,fold_threshold = fold_threshold,samples = names(em))
de_senes_mtd_v_senes = DE_table_creator(de_filepath = paste(path,"\\DE_Senes_MtD_vs_Senes.csv",sep = ""),em_annotated_tablename = em_annotated,separator = "tab",p_threshold = p_threshold,fold_threshold = fold_threshold,samples = names(em))
de_senes_mtd_v_prolif = DE_table_creator(de_filepath = paste(path,"\\DE_Senes_MtD_vs_Prolif.csv",sep = ""),em_annotated_tablename = em_annotated,separator = "tab",p_threshold = p_threshold,fold_threshold = fold_threshold,samples = names(em))

#making em symbols and em scaled tables
#all DE tables contain expression values so any one can be used for the subsets
em_symbols = de_senes_v_prolif[,2:10]
#scaling em_symbols by row
em_scaled = data.frame(t(scale(data.frame(t(em_symbols)))))

#for signature analysis in part V, we will need a master table combining senes_v_prolif and senes_mtd_v_senes DE analyses. creating this.
master = merge(de_senes_v_prolif, de_senes_mtd_v_senes, by.x=0, by.y=0, suffixes=c(".svp",".smdvs"))
#setting gene symbols to the gene symbols and removing unneeded Row.names column
row.names(master) = master$symbol.svp
master = master[,-1]

#creating a vector containing sig genes from all the DE comparisons (including genes that are in only one or two of the three)
sig_genes = c((row.names(subset(de_senes_v_prolif,sig == TRUE))),(row.names(subset(de_senes_mtd_v_prolif,de_senes_mtd_v_prolif$sig == TRUE))),row.names(subset(de_senes_mtd_v_senes, sig == TRUE))) #this will be used to plot the global heatmap showing significant differences between the three samples groups
#removing duplicates
sig_genes = unique(sig_genes)

# PART I - QUALITY CHECK PLOTS ####
#for quality checking we have MA and expression density plots

#plotting MA first
ggp = MA_plotter(de_senes_v_prolif,p_threshold,fold_threshold,plot_title = "MA Plot Senes v Prolif",y_axis_label = "Log2 Fold Change Senes v Prolif",y_axis_limits = c(-10,10))
plot_saver(path,"MA1",(450/96),(600/96),ggp,"block") 

ggp = MA_plotter(de_senes_mtd_v_senes,p_threshold,fold_threshold,plot_title = "MA Plot SenesMtD v Senes",y_axis_label = "Log2 Fold Change SenesMtD v Senes",y_axis_limits = c(-10,10))
plot_saver(path,"MA2",(450/96),(600/96),ggp,"block")

ggp = MA_plotter(de_senes_mtd_v_prolif,p_threshold,fold_threshold,plot_title = "MA Plot SenesMtD v Prolif",y_axis_label = "Log2 Fold Change SenesMtD v Prolif",y_axis_limits = c(-10,10))
plot_saver(path,"MA3",(450/96),(600/96),ggp,"block")

#expression density plots now
ggp = expression_matrices_plotter(em,colour1,colour2,colour3,plot_title = "Expression Density Plots")
plot_saver(path,"EM Density",default_plot_height ,(1200/96),ggp,"block")

# PART II - GLOBAL ANALYSIS ####
#this will have the global heatmap and the PCA plot.

#heatmap

#this first one will plot all the genes, across all samples (including non-sig)
ggp = heatmap_plotter(em_scaled,row.names(em_symbols),plot_title = "Global Heatmap (All Genes)",y_axis_label = "Genes")
plot_saver(path,"heatmap_global_all_genes",default_plot_height ,(550/96),ggp,"gradient")

#this second one will plot only genes significant in at least one DE comparison
ggp = heatmap_plotter(em_scaled,sig_genes,plot_title = "Global Heatmap (Sig Genes)",y_axis_label = "Significant Genes")
plot_saver(path,"heatmap_global_sig_genes",default_plot_height ,(550/96),ggp,"gradient")

#now the rug
ggp = rug_plotter(ss$SAMPLE_GROUP,rug_colours = sample_colours) 
plot_saver(path,"rug",(100/96) ,default_plot_width,ggp,"block")

#the PCA plot
ggp = PCA_plotter(em_scaled,ss,colour1,colour2,colour3,plot_title = "PCA Plot")
plot_saver(path,"PCA",(280/96),(600/96),ggp,"block")

# PART III - DE SENES V PROLIF ####

#for a global look, creating a heatmap of all genes significant in this DE analysis
genes_of_interest = row.names(subset(de_senes_v_prolif,sig == TRUE))
ggp = heatmap_plotter(em_scaled[,1:6],genes_of_interest,plot_title = "DE Senes v Prolif Heatmap (Sig Genes)")
plot_saver(path,"heatmap_svp_sig_genes",default_plot_height ,(550/96),ggp,"gradient")

#creating volcano plot
ggp = volcano_plotter(de_senes_v_prolif,p_threshold,fold_threshold,plot_title = "Volcano Plot Senes v Prolif",ylimits = c(0,300),xlimits = c(-10,10),down_labels = 13, up_labels = 10)
plot_saver(path,"volcano svp",default_plot_height,(660/96),ggp,"block")

#plotting individual boxplots of 10 most significant genes
genes_of_interest = row.names(de_senes_v_prolif[1:10,])
ggp = faceted_boxplot_plotter(em_scaled,genes_of_interest,ss$SAMPLE_GROUP,plot_title = "Most Significantly Differentially Expressed Genes",colour1 = colour1,colour2 = colour2, colour3 = colour3)
plot_saver(path,"boxplot sig svp",default_plot_height,default_plot_width,ggp,"block")

#doing pathway analysis using GO (it is better than GSEA because GSEA is blind to gene sets within which some genes go up and some go down, which is a common biological phenomenon.)
#the expectation here is that genes related to stress and inflammation would go up, genes related to proliferation would go down.
ora_results = ORA_calculator(de_senes_v_prolif,p_threshold,fold_threshold)

#plotting the ORA barplot
ggp = ORA_barplot_plotter(ora_results) + customized_theme + theme(legend.position = "right")
plot_saver(path,"ora barplot svp",(450/96),(725/96),ggp,"gradient")

#plotting a string plot
genes_of_interest = ora_candidate_genes_lister(ora_results)
plot_saver(path, "string network svp",(650/96),(950/96),(string_network_plotter(genes_of_interest,organism = "human")),"block")

#heatmap of "nuclear division" genes (i.e. most enriched gene set)
ggp = heatmap_plotter(em_scaled[,1:6],genes_of_interest,plot_title = "Senes v Prolif (Most Enriched Gene-Set)")
plot_saver(path,"heatmap_svp_go_genes",(330/96) ,(630/96),ggp,"gradient")

#creating a boxplot of most significantly different expressed genes present in the most enriched gene-set
genes_of_interest = candidate_genes_selector(de_senes_v_prolif,genes_of_interest = genes_of_interest)
ggp = faceted_boxplot_plotter(em_scaled,genes_of_interest,ss$SAMPLE_GROUP,plot_title = "Most Significantly DE ORA Enriched Genes",colour1 = colour1,colour2 = colour2, colour3 = colour3)
plot_saver(path,"boxplot go svp",(325/96),default_plot_width,ggp,"block")

#repeating this workflow for the other 2 DE comparisons

# PART IV - DE SENESMTD V SENES ####

#for a global look, creating a heatmap of all genes significant in this DE analysis
genes_of_interest = row.names(subset(de_senes_mtd_v_senes,sig == TRUE))
ggp = heatmap_plotter(em_scaled[,4:9],genes_of_interest,plot_title = "DE SenesMtD v Senes Heatmap (Sig Genes)")
plot_saver(path,"heatmap_smdvs_sig_genes",default_plot_height ,(550/96),ggp,"gradient")

#creating volcano plot
ggp = volcano_plotter(de_senes_mtd_v_senes,p_threshold,fold_threshold,plot_title = "Volcano Plot SenesMtD v Senes",ylimits = c(0,325),xlimits = c(-11,11),down_labels = 34, up_labels = 22)
plot_saver(path,"volcano smdvs",default_plot_height,(575/96),ggp,"block")

#plotting individual boxplots of 10 most significant genes
genes_of_interest = row.names(de_senes_mtd_v_senes[1:10,])
ggp = faceted_boxplot_plotter(em_scaled,genes_of_interest,ss$SAMPLE_GROUP,plot_title = "Most Significantly Differentially Expressed Genes",colour1 = colour1,colour2 = colour2, colour3 = colour3)
plot_saver(path,"boxplot sig smdvs",default_plot_height,default_plot_width,ggp,"block")

#doing pathway analysis using GO (it is better than GSEA because GSEA is blind to gene sets within which some genes go up and some go down, which is a common biological phenomenon.)
#the expectation here is that genes related to stress and inflammation would go down, genes related to proliferation would stay down, genes related to aging go down.
ora_results = ORA_calculator(de_senes_mtd_v_senes,p_threshold,fold_threshold)

#plotting the ORA barplot
ggp = ORA_barplot_plotter(ora_results) + customized_theme + theme(legend.position = "right")
plot_saver(path,"ora barplot smdvs",(450/96),(725/96),ggp,"gradient")

#plotting a string plot
genes_of_interest = ora_candidate_genes_lister(ora_results)
plot_saver(path, "string network smdvs",(650/96),(950/96),(string_network_plotter(genes_of_interest,organism = "human")),"block")

#heatmap of "extra-cellular matrix organization" genes (i.e. most enriched gene set)
ggp = heatmap_plotter(em_scaled[,4:9],genes_of_interest,plot_title = "SenesMtD v Senes (Most Enriched Gene-Set)")
plot_saver(path,"heatmap_smdvs_go_genes",(330/96) ,(630/96),ggp,"gradient")

#creating a boxplot of most significantly different expressed genes present in the most enriched gene-set
genes_of_interest = candidate_genes_selector(de_senes_mtd_v_senes,genes_of_interest = genes_of_interest)
ggp = faceted_boxplot_plotter(em_scaled,genes_of_interest,ss$SAMPLE_GROUP,plot_title = "Most Significantly DE ORA Enriched Genes",colour1 = colour1,colour2 = colour2, colour3 = colour3)
plot_saver(path,"boxplot go smdvs",(325/96),(825/325),ggp,"block")

# PART V - DE SENESMTD V PROLIF ####

#for a global look, creating a heatmap of all genes significant in this DE analysis
genes_of_interest = row.names(subset(de_senes_mtd_v_prolif,sig == TRUE))
index1 = c(1:3) #creating the right indexes to select appropriate columns from em table
index2 = c(7:9)
indexes = c(index1,index2)
#removing unnecessary objects from environment
rm(index1)
rm(index2)
ggp = heatmap_plotter(em_scaled[,indexes],genes_of_interest,plot_title = "DE SenesMtD v Prolif Heatmap (Sig Genes)")
plot_saver(path,"heatmap_smdvp_sig_genes",default_plot_height ,(550/96),ggp,"gradient")

#creating volcano plot
ggp = volcano_plotter(de_senes_mtd_v_prolif,p_threshold,fold_threshold,plot_title = "Volcano Plot SenesMtD v Prolif",ylimits = c(0,310),xlimits = c(-10,12),down_labels = 31, up_labels = 24)
plot_saver(path,"volcano smdvp",default_plot_height,(575/96),ggp,"block")

#plotting individual boxplots of 10 most significant genes
genes_of_interest = row.names(de_senes_mtd_v_prolif[1:10,])
ggp = faceted_boxplot_plotter(em_scaled,genes_of_interest,ss$SAMPLE_GROUP,plot_title = "Most Significantly Differentially Expressed Genes",colour1 = colour1,colour2 = colour2, colour3 = colour3)
plot_saver(path,"boxplot sig smdvp",default_plot_height,default_plot_width,ggp,"block")

#doing pathway analysis using GO (it is better than GSEA because GSEA is blind to gene sets within which some genes go up and some go down, which is a common biological phenomenon.)
#the expectation here is that genes related to proliferation would stay down.
ora_results = ORA_calculator(de_senes_mtd_v_prolif,p_threshold,fold_threshold)

#plotting the ORA barplot
ggp = ORA_barplot_plotter(ora_results) + customized_theme + theme(legend.position = "right")
plot_saver(path,"ora barplot smdvp",(450/96),(725/96),ggp,"gradient")

#plotting a string plot
genes_of_interest = ora_candidate_genes_lister(ora_results)
plot_saver(path, "string network smdvp",(650/96),(950/96),(string_network_plotter(genes_of_interest,organism = "human")),"block")

#heatmap of "nuclear division" genes (i.e. most enriched gene set)
ggp = heatmap_plotter(em_scaled[,indexes],genes_of_interest,plot_title = "SenesMtD v Prolif (Most Enriched Gene-Set)")
plot_saver(path,"heatmap_smdvp_go_genes",(330/96) ,(630/96),ggp,"gradient")

#creating a boxplot of most significantly different expressed genes present in the most enriched gene-set
genes_of_interest = candidate_genes_selector(de_senes_mtd_v_prolif,genes_of_interest = genes_of_interest)
ggp = faceted_boxplot_plotter(em_scaled,genes_of_interest,ss$SAMPLE_GROUP,plot_title = "Most Significantly DE ORA Enriched Genes",colour1 = colour1,colour2 = colour2, colour3 = colour3)
plot_saver(path,"boxplot go smdvp",(325/96),default_plot_width,ggp,"block")

# PART VI - SIGNATURES ####
#this analysis will be done between senes v prolif and senesmd v senes, as our research question questions what happens in the transcriptome when cells become senescent, and when senescent cells have their mitochondria depleted.

# first creating a fold vs fold plot
ggp = fold_vs_fold_scatterplot_plotter(master,"Log2 Fold Change Senes v Prolif","Log2 Fold Change SenesMtD v Senes")
plot_saver(path,"foldvfold",(400/96),(620/96),ggp,"block")

#we will take the significant genes and extract 4 signatures, as follows
#signature 1 will be down in senes v prolif and down in senesmd v senes (down, down)
#signature 2 will be up in senes v prolif and up in senesmd v senes (up, up)
#signature 3 will be down in senes v prolif and up in senesmd v senes (down, up)
#signature 4 will be up in senes v prolif and down in senesmd v senes (up, down)

# PART VIa SIGNATURE 1 ####
#signature 1 will be down in senes v prolif and down in senesmd v senes (down, down)
signature1 = c(row.names(subset(master,(master$sig.svp == TRUE & master$log2fold.svp < 0) & (master$sig.smdvs == TRUE & master$log2fold.smdvs < 0))))
#removing duplicates
signature1 = unique(signature1)

#creating a heatmap of the 4 signature genes
ggp = heatmap_plotter(em_scaled,signature1,plot_title = "Signature 1 Heatmap")
plot_saver(path,"sign1 heatmap",(300/96),(350/96),ggp,"gradient")

#creating metagene boxplot
ggp = metagene_boxplot_plotter(em_scaled,signature1,ss,sample_colours)
plot_saver(path,"sign1 metagene",(300/96),(500/96),ggp,"block")

#creating a boxplot
ggp = faceted_boxplot_plotter(em_scaled,signature1,ss$SAMPLE_GROUP,plot_title = "Signature 1 Genes",colour1 = colour1, colour2 = colour2, colour3 = colour3)
plot_saver(path,"sign1 boxplot",(300/96),(450/96),ggp,"block")

#doing ORA, not using function as it is not adaptable here
sig_genes_entrez = bitr(signature1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = "all", pvalueCutoff = p_threshold, qvalueCutoff = 0.1)

#upon inspection, it was found that there are no results that cross the p threshold, and so there is no gene ontology to plot.
#hence GO cannot be shown.

# PART VIb SIGNATURE 2 ####
#signature 2 will be up in senes v prolif and up in senesmd v senes (up, up)
signature2 = c(row.names(subset(master,(master$sig.svp == TRUE & master$log2fold.svp > 0) & (master$sig.smdvs == TRUE & master$log2fold.smdvs > 0))))
#removing duplicates
signature2 = unique(signature2)

#creating a heatmap of the signature genes
ggp = heatmap_plotter(em_scaled,signature2,plot_title = "Signature 2 Heatmap")
plot_saver(path,"sign2 heatmap",default_plot_height,(550/96),ggp,"gradient")

#creating metagene boxplot
ggp = metagene_boxplot_plotter(em_scaled,signature2,ss,sample_colours)
plot_saver(path,"sign2 metagene",(300/96),(500/96),ggp,"block")

#creating a boxplot
ggp = faceted_boxplot_plotter(em_scaled,signature2,ss$SAMPLE_GROUP,plot_title = "Signature 2 Genes",colour1 = colour1, colour2 = colour2, colour3 = colour3)
plot_saver(path,"sign2 boxplot",(350/96),(1200/96),ggp,"block")

#doing ORA, not using function as it is not adaptable here
sig_genes_entrez = bitr(signature2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = "all", pvalueCutoff = p_threshold, qvalueCutoff = 0.1)

#upon inspection, it was found that there are no results that cross the p threshold, and so there is no gene ontology to plot.
#hence GO cannot be shown.

# PART VIc SIGNATURE 3 ####
#this is signature 1 in the report 
#signature 3 will be down in senes v prolif and up in senesmd v senes (down, up)
#we expect that this since depleting mitochondria lowers SASP (senescence-associated secretory phenotype), health markers like bone density which would go down with senescence would increase after depleting mitochondria.
signature3 = c(row.names(subset(master,(master$sig.svp == TRUE & master$log2fold.svp < 0) & (master$sig.smdvs == TRUE & master$log2fold.smdvs > 0))))
#removing duplicates
signature3 = unique(signature3)

#creating a heatmap of the signature genes
ggp = heatmap_plotter(em_scaled,signature3,plot_title = "Signature 3 Heatmap")
plot_saver(path,"sign3 heatmap",default_plot_height,(550/96),ggp,"gradient")

#creating metagene boxplot
ggp = metagene_boxplot_plotter(em_scaled,signature3,ss,sample_colours,plot_title = "Signature 1 Metagene Plot")
plot_saver(path,"sign3 metagene",(300/96),(500/96),ggp,"block")

#using a random subset of genes to create boxplot
genes_of_interest = sample(signature3,10,replace=FALSE)
ggp = faceted_boxplot_plotter(em_scaled,genes_of_interest,ss$SAMPLE_GROUP,plot_title = "Signature 3 Genes",colour1 = colour1, colour2 = colour2, colour3 = colour3)
plot_saver(path,"sign3 boxplot",(300/96),(700/96),ggp,"block")

#doing ORA, not using function as it is not adaptable here
sig_genes_entrez = bitr(signature3, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = "all", pvalueCutoff = p_threshold, qvalueCutoff = 0.1)

#upon inspection, it was found that there are no results that cross the p threshold, and so there is no gene ontology to plot.
#hence GO cannot be shown.

# PART VId SIGNATURE 4 ####
#this is signature 2 in the report.
#signature 4 will be up in senes v prolif and down in senesmd v senes (up, down)
#from literature, we can clearly expect to find that stress and inflammation related genes go up during senescence, and then go down after mitochondria depletion.
#we expect that this since depleting mitochondria lowers SASP (senescence-associated secretory phenotype), health markers like bone density which would go down with senescence would increase after depleting mitochondria.
signature4 = c(row.names(subset(master,(master$sig.svp == TRUE & master$log2fold.svp > 0) & (master$sig.smdvs == TRUE & master$log2fold.smdvs < 0))))
#removing duplicates
signature4 = unique(signature4)

#creating a heatmap of the signature genes
ggp = heatmap_plotter(em_scaled,signature4,plot_title = "Signature 4 Heatmap")
plot_saver(path,"sign4 heatmap",default_plot_height,(550/96),ggp,"gradient")

#creating metagene boxplot
ggp = metagene_boxplot_plotter(em_scaled,signature4,ss,sample_colours,plot_title = "Signature 2 Metagene Plot")
plot_saver(path,"sign4 metagene",(300/96),(500/96),ggp,"block")

#doing ORA, not using function as it is not adaptable here
sig_genes_entrez = bitr(signature4, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = "all", pvalueCutoff = p_threshold, qvalueCutoff = 0.1)

#plotting GO results in a barplot
ggp = ORA_barplot_plotter(ora_results)
plot_saver(path,"sign4 ora barplot",default_plot_height,(600/96),ggp,"gradient")

#plotting a string plot
ggp = string_network_plotter(ora_candidate_genes_lister(ora_results),organism = "human")
plot_saver(path,"sign4 ora string network",(650/96),(950/96),(string_network_plotter(ora_candidate_genes_lister(ora_results),organism = "human")),"block")

#heatmap of "granulocyte chemotaxis" genes (i.e. most enriched gene set)
genes_of_interest = ora_candidate_genes_lister(ora_results)
ggp = heatmap_plotter(em_scaled,genes_of_interest,plot_title = "Signature 4 (Most Enriched Gene-Set)")
plot_saver(path,"sign4 ora heatmap",(330/96) ,(630/96),ggp,"gradient")

#creating a boxplot of genes present in the most enriched gene-set
ggp = faceted_boxplot_plotter(em_scaled,genes_of_interest,ss$SAMPLE_GROUP,plot_title = "Signature 4 Most Enriched Genes",colour1 = colour1,colour2 = colour2, colour3 = colour3)
plot_saver(path,"sign4 ora boxplot",(325/96),(850/96),ggp,"block")

# and done :) ####