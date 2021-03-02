#!/usr/bin/env Rscript

version <- "1.0"

# read in required libraries
library( "ggplot2" )
library( "optparse" )

option_list <- list(

		  make_option( c( "--input"), type="character", default=NA, 
	                  help = "Input table (required)" , metavar="path"),

		  make_option( c( "--output"), type="character", default=NA, 
        	          help = "Output pdf (required)", metavar="path"),

		  make_option( c( "--function_id" ), type="character", default=NA, 
       		          help = "Function id (such as a KEGG ortholog id) (required)" , metavar="string"),

		  make_option( c( "--rank"), type="character", default=NA, 
         	          help = "Taxonomic rank to plot (required). Choices: Phylum, Class, Order, Family, Genus, Species, OTU" , metavar="string"),

		  make_option( c(  "--sort_taxa" ) , type="character" , default="name" ,
			  help = "Plot taxa alphabetically (\"name\") or by function counts (\"count\") in increasing order [default = %default]" ,
			  metavar = "string" ) ,

		  make_option( c(  "--min_OTU_contrib" ) , type="numeric" , default=0 ,
			  help = "Rows in input table with a value of CountContributedByOTU less than this value will be removed [default = %default]" ,
			  metavar = "number" ) ,

		  make_option( c( "--rel_abundance" ) , action = "store_true" , type="logical" , default=FALSE , 
		  	  help = "Flag to indicate that relative abundances should be plotted [default = %default]" ) ,

		  make_option(  c( "--width" )  , type="numeric" , default=7 ,
          	          help = "Manually set pdf width in inches [default = %default]",
			  metavar = "number" ) ,

		  make_option( c(  "--height" ) , type="numeric" , default=7 ,
          	          help = "Manually set pdf height in inches [default = %default]",
 			  metavar = "number" ) , 

		  make_option( c( "--sample_order" ) , type="character" , default=NA , 
		  	  help = "File containing desired plotting order of samples (1 sample per line) [default = %default]" ,
		  	  metavar = "path" ) ,

		  make_option( c(  "-n" , "--no_background_grid" ) , action = "store_true" , type="logical" , default=FALSE , 
		  	  help = "Flag to indicate that the grey background grid should not be plotted [default = %default]" ) ,

		  make_option( c( "--R_default_colours" ) , action = "store_true" , type="logical" , default=FALSE , 
		  	  help = "Flag to indicate that default ggplot2 colours should be used [default = %default]" ) ,

		  make_option( c( "-v" , "--version" ) , action = "store_true" , type="logical" , default=FALSE ,
			  help = "Print out version number and exit", 
			  metavar = "string" ) 

		 )

opt_parser <- OptionParser( option_list=option_list , 
							usage = "%prog [options] --input PATH --output PATH --function_id STRING --rank STRING", )
opt <- parse_args( opt_parser )


# check for version flag
if( opt$version ) {
	cat( "Version" , version , "\n")
	options_tmp <- options(show.error.messages=FALSE) 
    on.exit(options(options_tmp)) 
	stop()
}

# check for required arguments
if ( is.na( opt$input ) ) {
	stop( "Path to input table needs to be specified with --input\nType \"plot_metagenome_contributions.R --help\" for help." )
} else if ( is.na( opt$output )) {
	stop( "Path to output pdf needs to be specified with --output\nType \"plot_metagenome_contributions.R --help\" for help." )
} else if ( is.na( opt$function_id )) {
	stop( "Function of interest needs to be specified with --function_id\nType \"plot_metagenome_contributions.R --help\" for help." )
} else if ( is.na( opt$rank )) {
	stop( "Taxonomic rank needs to be specified with --rank\nType \"plot_metagenome_contributions.R --help\" for help." )
}

# Check that sort_taxa is either "name" or "abundance"
if ( ( opt$sort_taxa != "name" ) & ( opt$sort_taxa != "count" ) ) {
	stop( "--sort_taxa needs to be either name (default) or count, not " , opt$sort_taxa , "\nType \"plot_metagenome_contributions.R --help\" for help.")
}

# Capitalize inputted rank since this would be an easy mistake to make
opt$rank <- paste( toupper(substring(opt$rank, 1,1)) , substring(opt$rank, 2), sep="")

# Set of 12 qualitative colours taken from colorbewer2.org
# These are the colours used unless the "--R_default_colours" flag is used
qual_col <- c( "#a6cee3" ,  "#1f78b4" , "#b2df8a" , "#33a02c" , "#fb9a99" , "#e31a1c" ,
	       "#fdbf6f" , "#ff7f00" , "#cab2d6" , "#6a3d9a" , "#ffff99" , "#b15928" )


return_nearest_lineage <- function( input  , level  ) {

	# Function to return nearest parent level that is classified (along with a flag to
	# distinguish unclassified and unassigned taxa

	taxa_levels <- c( "Kingdom" , "Phylum" , "Class" , "Order" , "Family" , "Genus" , "Species", "OTU" )
	missing_taxa <- c( "p__" , "c__" , "o__" , "f__" , "g__" , "s__" , "", "Unclassified")

	starting_index <- which( taxa_levels == level ) - 1

	vec2return <- vector( mode="character", length=nrow(input) )

	# Loop over each row in df
	for ( i in 1:nrow(input) ) {

		# Most missing taxa are likely "Unclassified", but some could be "Unassigned" (e.g. they are "" rather than "s__")
		missing_type <- "Unclassified"
		if ( input[ i , level ] == "" ) {
			missing_type <- "Unassigned"
		}

		current_index <- starting_index

		while( current_index > 0 ) {

		        # If taxa label is not in missing_taxa vector then it must be a classified label
		        # This label is returned and prefixes the plotted taxa label
			if ( ! ( input[ i , taxa_levels[ current_index ] ] %in% missing_taxa )  ) {
				vec2return[i] <- paste( input[ i , taxa_levels[ current_index ] ] , missing_type , sep=";" )
				break
			}

			current_index <- current_index - 1

		}
	}
	return( vec2return )
}

fix_missing_taxa <- function( table , tax_level ) {

  ### Function to rename missing taxa labels with nearest classified level followed by "Unclassified" or "Unassigned"

  # First check whether specified tax_level is a column in table or not
  if ( ! ( tax_level %in% colnames( table ) ) ) {
    stop( "Error - specified taxa level " , tax_level , " is not a column in input table.")   
  }
  
  # If a taxa label matches any of these strings then that means it's either unclassified or unassigned (in the case of "")
  missing_taxa <- c( "p__" , "c__" , "o__" , "f__" , "g__" , "s__" , "", "Unclassified")
    
  missing_rows <- which( table[ , tax_level ] %in% missing_taxa )

  if ( length( missing_rows ) > 0 ) {
 	
    table[ missing_rows , tax_level ]  <- return_nearest_lineage( 
    	input = table[ missing_rows , , drop = FALSE ] , level = tax_level  )
  }
  
  return( table )
  
}

check_unique_parents <- function( table , level ) {

	# Function to check that each taxa of interest has only 1 unique parent
	# This should catch issues where the same species name is used by 2 different genera

	taxa_levels <- c( "Kingdom" , "Phylum" , "Class" , "Order" , "Family" , "Genus" , "Species")
	missing_taxa <- c( "p__" , "c__" , "o__" , "f__" , "g__" , "s__" , "", "Unclassified")

	higher_levels <- taxa_levels[ 1:(which(taxa_levels %in% level) - 1 ) ]

	parent_taxa_list <- list( c() )

	for ( i in 1:nrow(table)) {

  		parents <- paste( table[i,higher_levels] , collapse=";")
  		taxa_i <- table[i,level]
  
  		# skip row if missing taxa label
  		if ( taxa_i %in% missing_taxa ) {
    		next
  		}
  
  		if ( taxa_i %in% names(parent_taxa_list) ) {
		    	if( ! identical( parent_taxa_list[[taxa_i]] , parents ) ) {
 			     err <- paste( "Error: taxa" , taxa_i , "is a child of both" , parent_taxa_list[[taxa_i]] , "and" , parents , sep=" " )
 			     stop( err )
   			 } 
		} else {
			parent_taxa_list[[taxa_i]] <- parents
  		}
	}
}


# read in file
input <- read.delim( opt$input , stringsAsFactors = FALSE )

# Subset to function
input_subset <- input[ input$Gene == opt$function_id , ]

# Remove rows with value of CountContributedByOTU less than user-specified value
input_subset <- input_subset[ which( input_subset$CountContributedByOTU >= opt$min_OTU_contrib ) , ]

# Check whether subset is empty
if ( nrow(input_subset) == 0 ) {
	stop( "No matches to function " , opt$function_id , " found in Gene column" )
}

# Re-label Unclassified and Unassigned taxa
input_subset_relab <- fix_missing_taxa( input_subset , opt$rank )

# Check that all taxa have only 1 unique parent lineage, IF rank other than OTU
if ( opt$rank != "OTU" ) {
	check_unique_parents( table=input_subset_relab , level=opt$rank )
}

# Convert Sample column to factor and if particular sample order was indicated then set levels in that order
if( is.na(opt$sample_order ) ) {
	# No sample order specified, so just turn Sample column into factor
	input_subset_relab$Sample <- factor( input_subset_relab$Sample )
} else {
	# File with sample order specified so read in and specify factor levels in this order
	sample_order_file <- read.delim( opt$sample_order , header=F , stringsAsFactors=FALSE )
	input_subset_relab$Sample <- factor( input_subset_relab$Sample , levels = sample_order_file$V1 )
}

# If relative abundance is wanted instead then transform data into that format per sample:
if ( opt$rel_abundance ) {

	# get total sums per sample
	sample_sums <- aggregate( input_subset_relab$CountContributedByOTU , 
						by=list(Category=input_subset_relab$Sample) , FUN=sum)

	# loop over each sample and divide each of their counts by their sum and then multiply by 100 
	for ( i in 1:nrow(sample_sums)) {
	  sample <- sample_sums[ i , "Category" ]
 	 input_subset_relab[ which( input_subset_relab$Sample == sample ) , 
 	 			"CountContributedByOTU" ] <- 
 	 			( input_subset_relab[ which( input_subset_relab$Sample == sample ) , 
 	 				"CountContributedByOTU" ] / sample_sums[ i , "x"]) * 100 
	}

}

# make new df with following columns: Sample, CountContributedByOTU and user-specified taxonomic rank
input_subset_relab_focal <- input_subset_relab[ , c( "Sample" , "CountContributedByOTU" , opt$rank ) ]

colnames( input_subset_relab_focal ) <- c("Sample" , "CountContributedByOTU" , "Taxa" )

# sum CountContributedByOTU by taxa and samples
input_subset_relab_focal_sum <- aggregate( CountContributedByOTU ~ Taxa + Sample, data=input_subset_relab_focal , FUN=sum )

# Convert Taxa column to factor.
# Change order of Taxa levels if user specified they wanted them ordered by count
if ( opt$sort_taxa == "name" ) {
	input_subset_relab_focal_sum$Taxa <- factor( input_subset_relab_focal_sum$Taxa )
} else if ( opt$sort_taxa == "count" ) {

	# sum CountContributedByOTU by taxa over all samples (for getting factor level order)
	dataset_wide_aggregate <- aggregate( CountContributedByOTU ~ Taxa, data=input_subset_relab_focal_sum , FUN=sum )
	dataset_wide_aggregate <- dataset_wide_aggregate[ order( dataset_wide_aggregate$CountContributedByOTU , decreasing = FALSE ) , ]
	input_subset_relab_focal_sum$Taxa <- factor( input_subset_relab_focal_sum$Taxa , levels = dataset_wide_aggregate$Taxa )
}

#output plot to pdf 
pdf( opt$output , width = opt$width , height = opt$height )

# if more levels than qualitative colours then just repeat them:
if ( length( levels( input_subset_relab_focal_sum$Taxa ) ) > length( qual_col ) ) {
	fold <- ceiling( length( levels( input_subset_relab_focal_sum$Taxa ) ) / length( qual_col ) )
	qual_col <- rep( qual_col , fold )
} 

ylab_start <- "Abundance of "

# Change y lab to say relative abundance if necessary
if ( opt$rel_abundance) { 
	ylab_start <- "Relative Abundance of "
}

# create stacked bar chart, with x-axis being samples, collapsing to user defined taxonomic rank, and weighted by actual OTU abundances
stacked_plot <- qplot( Sample, data=input_subset_relab_focal_sum, geom="bar",fill=input_subset_relab_focal_sum$Taxa ,
	weight = CountContributedByOTU , xlab = "Samples" , ylab = paste( ylab_start, opt$function_id ) ) +  
# include samples on x-axis that have no counts
	scale_x_discrete( drop = FALSE )
# if default colour flag not set then use fewer, but more discernible set
if ( ! opt$R_default_colours ) {
	stacked_plot <- stacked_plot + scale_fill_manual(values=qual_col) 
}
# remove grid background if specified
if ( opt$no_background_grid ) {
	stacked_plot <- stacked_plot + theme(axis.text.x=element_text(angle = 90, vjust = 0.5) , legend.title=element_blank() , 
		panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), axis.line = element_line(colour = "black"))
} else {
	stacked_plot <- stacked_plot + theme(axis.text.x=element_text(angle = 90, vjust = 0.5) , legend.title=element_blank() )
}

plot(stacked_plot)

graphics.off()
