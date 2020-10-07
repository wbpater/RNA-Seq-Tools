PCAplotR = function (x, PCX=1, PCY=2, ntop = 500, returnData = F, returnLoadings = F,                     # define function with default arguments
                     palette = NULL, colour = NULL, label = NULL, shape = NULL,                           # define function with default arguments
                     file = NULL, plotsize = c(20,25))                                                    # define function with default arguments
{
  ## PCAplotR is a heavily modified plotPCA function that is designed to plug in to DESeq2. 
  ## The function takes vst normalised Large DESeqTransform objects as input (x) and can output a plot
  ## to the screen or to a specified file using ggplot functionality.
  ## The user can specify which PCs to plot on the x and y axes with the PCX and PCY options.
  ## Users can specify multiple visual options such as colour, shape and labels of datapoints, 
  ## standard ggplot colours and symbols will be used.
  ## Users can also specify their own colours as a character vector with the palette option.
  ## The function can give the loadings in percentages for the set PCX if returnLoadings is set to TRUE.
  ## The function can return the PCA data to the user by setting returnData to TRUE
  ##
  ## Arguments
  ## x                                      a DESeqDataSet object, preferably vst normalised (DESeqDataSet object)
  ## PCX                                    the PC that will be plotted on the x axis (integer)
  ## PCY                                    the PC that will be plotted on the y axis (integer)
  ## ntop                                   the x highest variance genes to be used in the pca, if higher than the number
  ##                                        of genes in the input object, the maximum number of genes is used (integer)
  ## returnData                             returns the PCA data to the console if TRUE (boolean)
  ## returnLoadings                         returns the loadings of PCX to the console if TRUE (boolean)
  ## palette                                the colours used with the colour option (character vector)
  ## colour                                 the field of column data that is used to colour samples (string)
  ## label                                  the field of column data that is used to label samples (string)
  ## shape                                  the field of column data that is used to set shapes for samples (string)
  ## file                                   the filepath to which the plot should be saved (string)
  ## plotsize                               width and height of the saved plot in cm (vector)
  ##
  ## Examples
  ## PCAplotR(my_DESeq_object)
  ## PCAplotR(my_DESeq_object, file = "/home/user/plots/my_R_plot.png", plotsize = c(10,10))
  ## PCAplotR(my_DESeq_object, PCX = 2, PCY = 3, ntop = 1000)
  ## PCAplotR(my_DESeq_object, PCX = 2, PCY = 3, ntop = 1000, palette = c("#999999", "#ffffff"), colour = "colData_column_1")
  ## PCAplotR(my_DESeq_object, PCX = 2, PCY = 3, ntop = 1000, palette = c(value1="#999999", value2="#00ff00"), colour = "colData_column_1", 
  ## shape = "colData_column_2", label = "colData_column_3")+ ggplot_option1
  
  ## add check for class(DESEqdataset)
  
  if (isFALSE("genefilter" %in% rownames(installed.packages())) || isFALSE("ggplot2" %in% rownames(installed.packages())) || # check if required packages are installed
      isFALSE("ggrepel" %in% rownames(installed.packages())) || isFALSE("DESeq2" %in% rownames(installed.packages()))){
    stop("The packages 'genefilter', 'ggplot2', 'ggrepel' and 'DESeq2' are required to run.")
  }
  if (isFALSE(colour %in% colnames(colData(x)))){                                                         # check if the colour argument is in the column data of the input object
    stop(paste0("the \"colour\" argument should specify columns of coldata(",deparse(substitute(x)),")")) # stop if a wrong argument for "colour" is given
  }
  if (!is.null(label) & isFALSE(label %in% colnames(colData(x)))){                                        # check if the label argument is in the column data of the input object
    stop(paste0("the \"label\" argument should specify columns of coldata(",deparse(substitute(x)),")"))  # stop if a wrong argument for "label" is given
  }
  if (isFALSE(shape %in% colnames(colData(x))) & isFALSE(shQuote(shape) %in% colnames(colData(x)))){                                                          # check if the shape argument is in the column data of the input object
    stop(paste0("the \"shape\" argument should specify columns of coldata(",deparse(substitute(x)),")"))  # stop if a wrong argument for "shape" is given
  }
  if (returnData == T & returnLoadings == T){
    stop("Only one of returnData, returnLoadings may be TRUE at one time")
  }
  if (ntop > nrow(assay(x))){                                                                             # check if ntop is greater than the number of genes in the input
    message("Waring: ntop is greater than the number of rows in the input, using all rows")               #
  }
  if (!is.null(plotsize)){
    if (!is.vector(plotsize) & !length(plotsize) == 2){
      stop("Plotsize must be a vector with 2 elements: width and height, in that order")
    }
  }
  library(genefilter)                                                                                     # set up required libraries
  library(ggplot2)                                                                                        # set up required libraries
  library(ggrepel)                                                                                        # set up required libraries
  library(DESeq2)
  
  rv = rowVars(assay(x))                                                                                  # calculate the variance per row (gene)
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]                                    # select the index numbers of the (ntop) highest variance genes
  pca = prcomp(t(assay(x)[select, ]))                                                                     # perform a pca on the input object's (ntop) highest variance genes
  percentVar = pca$sdev^2/sum(pca$sdev^2)                                                                 # calculate the percentage of variance in each of the principal components
  d = data.frame(assign(paste0("PC", PCX), pca$x[, PCX]), assign(paste0("PC", PCY), pca$x[, PCY]),        # construct a dataframe where column 1 contains the pca values belonging to PCX and column 2 the values belonging to PCY 
                 as.data.frame(colData(x)))                                                               # and further column correspond to the column data of the input object
  
  if (returnLoadings) {                                                                                   # get the loadings for PCX
    ldg = matrix(nrow = nrow(pca$rotation))                                                               # get the factor loadings from prcomp
    ldg[,1] = sapply(as.matrix(pca$rotation[,PCX]), function(val){val^2})                                 # calculate the percentage a gene contributes to the PC
    rownames(ldg)=rownames(pca$rotation)
    ldg = as.matrix(sort(ldg[,1], decreasing = T))                                                        # ouput % contribution of gene to indicated PC
    colnames(ldg) = c(paste0("contribution to PC", PCX, " in %"))
    return(ldg)
  }
  
  if (returnData) {                                                                                       # if the returnData argument is true
    attr(d, "percentVar") = percentVar[PCX:PCY]                                                           # return dataframe d and percentage variance on the PCs
    return(d)                                                                                             # 
  }
  
  plot = ggplot(data = d, aes_string(x = paste0("PC", PCX), y = paste0("PC", PCY),                        # define principal components on axes, points coloured on variable and labeled
                                     color = colour, label = label, shape = shape)) + 
    geom_point(size = 3) +                                                                                # label text size
    theme_light() +                                                                                       #
    xlab(paste0("PC", PCX, ": ", round(percentVar[PCX] * 100), "% variance")) +                           # set x label to have percentage of variance (from pca data input)
    ylab(paste0("PC", PCY, ": ", round(percentVar[PCY] * 100), "% variance")) +                           # set y label to have percentage of variance (from pca data input)
    scale_x_continuous(expand = expand_scale(mult = c(0.1, 0.1))) +                                       # expands the plot area by 10% to the right and left to give space for the labels
    scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.1))) +                                       # expands the plot area by 10% to the top and bottom to give space for the labels
    theme(plot.margin=unit(c(4, 4, 4, 4), "points"))  +                                                   # set white margins around the plot
    ggtitle(
      paste0("PC",PCX," vs. PC",PCY, " for the top ", min(ntop,nrow(assay(x))),                           # set an automatic title based on input parameters
             " highest variance genes of ", deparse(substitute(x))),
      subtitle =  paste0(" coloured on ", colour)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),                                                  # center the title above the plot and increase text size
      plot.subtitle = element_text(hjust = 0.5, size = 12),                                               # center the subtitle above the plot and increase text size
      panel.border = element_rect(fill=NA, size=1),                                                       # draw a border around the graph area
      axis.title = element_text(size = 12),                                                               # increase axis label size
      axis.text = element_text(size = 10),                                                                # increase axis unit size
      legend.background = element_rect(size = 0.7, colour = "grey")                                       # draw a border around the legend(s)
    ) + 
    {if(!is.null(palette))scale_colour_manual(values=palette)}+                                           # set the colour scale for lines and points to the hexadecimal values in a character vector
    {if(!is.null(label))geom_text_repel(size=4, force = 10)} +                                            # set size for label text and repel the labels away from the points
    {if(!is.null(label))ggtitle(paste0("PC",PCX," vs. PC", PCY, " for the top ", min(ntop,nrow(assay(x))),# 
                                       " highest variance genes of ",deparse(substitute(x))),             #
                                subtitle =  paste0(" coloured on ", colour, ", labeled on ", label))}     #                     
  if (is.null(file)){return(plot)}                                                                        # return plot if no output file is defined
  
  if (!is.null(file) & !is.null(plotsize)){
    ggsave(file, width=plotsize[1]/2.54, height=plotsize[2]/2.54)                                         # output to file with width and length in cm  
    message(paste0("Your plot was saved as ", file))                                                      # return the plot location to the console
  }
}
#PCAplotR(vst ,PCX=1, PCY=2, ntop = 500 , palette = cbPalette, label = "name",
#         shape = "dataset", colour= "mutation", file = "/scratch/wpater/Rplot.png", plotsize = c(10,20))+ scale_x_reverse()#+scale_x_reverse()#, file = "/scratch/wpater/Rplot.png")

#cbPalette = c( "20196" = "#999999", NL10 = "#999999", NL2= "#999999", NL3= "#999999", NL4= "#999999",             
#               "-" = "#ffffff", PCLD13= "#ffffff", PCLD14= "#ffffff", PCLD15= "#ffffff", PCLD16= "#ffffff",
#               PCLD2="#00ff00",PCLD8="#ff0000",PCLD9="#ff00ff",PCLD10="#0000ff",PCLD12="#ffa500",PCLD17="#00CC00")   # colour on id
#cbPalette = c(none="#999999", PRKCSH_LOH="#00ff00", PRKCSH="#006400",unknown="#ff0000",SEC63="#ff00ff",PKD1="#0000ff",sol_cyst="#ffa500")   # colour on mutation

#for (i in 1:10) {print(PCAplotR(vst , PCX=i , PCY=i+1 , ntop = 500 , palette = cbPalette, label = "name",shape = "dataset", colour= "mutation"))}


## loadings with gene symbols
#ldg = PCAplotR(vst, PCX = 1, returnLoadings = T)
#ldg = cbind(ldg, sapply(strsplit(row.names(ldg), ".", fixed=T), function(name) name[1]))
#row.names(ldg) = ldg[,2]
#ldg[,2] = mapIds(org.Hs.eg.db, keys=row.names(ldg), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
