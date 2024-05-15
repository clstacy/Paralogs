#' library(edgeR)
#' library(dplyr)
#' library(tibble)
#' library(KEGGREST)
#' library(ggkegg)
#' library(AnnotationDbi)
#'
#'
#' #' @title Preprocess Differential Expression Analysis Results
#' #' @description Processes outputs from limma, edgeR, or DESeq2 to create a standardized dataframe.
#' #' @param results Object containing the results of a differential expression analysis.
#' #' @return A dataframe with standardized columns: GeneID, logFC, pValue, adjPValue.
#' #' @export
#' preprocess_DE_results <- function(results) {
#'   # Determine the class of the results object and process accordingly
#'   if (attr(class(results),"package") == "edgeR") {
#'     # edgeR topTags output
#'     data <- topTags(results) %>%
#'       as.data.frame() %>%
#'       dplyr::mutate(GeneID = ORF, logFC = logFC, pValue = PValue, adjPValue = FDR)
#'   } else if (attr(class(results),"package") == "DESeq2") {
#'     # DESeq2 results output
#'     data <- as.data.frame(results) %>%
#'       dplyr::mutate(GeneID = ORF, logFC = log2FoldChange, pValue = pvalue, adjPValue = padj)
#'   } else if (attr(class(results),"package") == "limma") {
#'     # limma results output
#'     data <- as.data.frame(results) %>%
#'       dplyr::mutate(GeneID = ORF, logFC = logFC, pValue = P.Value, adjPValue = adj.P.Val)
#'   } else if (is.data.frame(results)) {
#'     # topTags output from edgeR
#'     data <- as.data.frame(results) %>%
#'       dplyr::mutate(GeneID = rownames(.), logFC = logFC, pValue = pValue, adjPValue = adjPValue)
#'   } else {
#'     stop("Unsupported input type. Please provide an object from edgeR, DESeq2, or limma.")
#'   }
#'
#' # Ensure the output has the correct column names and data types
#'   data <- data %>%
#'     dplyr::select(GeneID, logFC, pValue, adjPValue) %>%
#'     dplyr::arrange(adjPValue)  # Arrange by adjusted P-value for convenience
#'   return(data)
#' }
#'
#'
#'
#'
#' # library(clusterProfiler)
#' # library(org.Sc.sgd.db) # Example for demonstration, ensure this is loaded appropriately
#'
#' #' @title Convert Gene IDs to universal IDs within a DataFrame
#' #' @description Converts gene identifiers in a dataframe to Entrez IDs using the specified organism database.
#' #' @param df Dataframe containing gene identifiers.
#' #' @param gene_id_column Name of the column in df that contains the gene identifiers.
#' #' @param organism_code A string specifying the organism code.
#' #' @param annotation_db_name A string naming the Bioconductor annotation package.
#' #' @param fromType The type of gene ID currently in the dataframe.
#' #' @param toType The target type of gene ID to convert to (default "ENTREZID").
#' #' @return A dataframe with an additional column of converted gene IDs.
#' homogenize_gene_ids <- function(df, gene_id_column, organism_code, annotation_db, fromType, toType = "ENTREZID") {
#'   # Ensure the annotation package is loaded
#'   # if (!requireNamespace(annotation_db_name, quietly = TRUE)) {
#'   #   stop(paste("The package", annotation_db_name, "is not installed. Please install it using BiocManager::install('", annotation_db_name, "').", sep=""))
#'   # }
#'
#'   # Access the organism database object
#'   # OrgDb <- get(annotation_db_name, envir = asNamespace(annotation_db_name))
#'   OrgDb <- org.Sc.sgd.db
#'
#'   # Check if the conversion types are supported
#'   all_keytypes <- columns(OrgDb)
#'   if (!(fromType %in% all_keytypes && toType %in% all_keytypes)) {
#'     stop(print(all_keytypes))
#'     stop("Provided keytypes are not supported by the database.")
#'   }
#'
#'   # Attempt to convert gene IDs
#'   converted_ids <- tryCatch({
#'     bitr(df[[gene_id_column]], fromType = fromType, toType = toType, OrgDb = OrgDb)
#'   }, error = function(e) {
#'     warning("Failed to convert some or all gene IDs: ", e$message)
#'     return(NULL)  # Return NULL if conversion fails
#'   })
#'
#'   # Merge the converted IDs back to the original dataframe
#'   if (!is.null(converted_ids) && "ENTREZID" %in% names(converted_ids)) {
#'     df$ENTREZID <- converted_ids$ENTREZID
#'   } else {
#'     df$ENTREZID <- NA  # Assign NA if no IDs were converted
#'   }
#'
#'   return(df)
#' }
#'
#'
#'
#'
#'
#'
#'
#' #' @title Get Pathway Number
#' #' @description Retrieves the pathway number based on a pathway ID from KEGG data.
#' #' @param kegg_data A dataframe containing KEGG pathway data.
#' #' @param pathway_id A string specifying the pathway ID.
#' #' @return An integer representing the pathway number.
#' get_pathway_number <- function(kegg_data, pathway_id) {
#'   pathway_number <- kegg_data %>%
#'     dplyr::mutate(row_number = dplyr::row_number()) %>%
#'     dplyr::filter(ID == pathway_id) %>%
#'     dplyr::pull(row_number)
#'   return(pathway_number)
#' }
#'
#' #' @title Fetch KEGG Reactions
#' #' @description Retrieves KEGG reaction mappings and formats them for use.
#' #' @return A dataframe containing KEGG reactions with their descriptions.
#' fetch_kegg_reactions <- function() {
#'   if (!exists("KEGG_reactions")) {
#'     KEGG_reactions <- KEGGREST::keggList("reaction")
#'   }
#'     KEGG_reactions_df <- as.data.frame(KEGG_reactions)
#'     colnames(KEGG_reactions_df) <- "long_reaction_description"
#'     KEGG_reactions_df <- KEGG_reactions_df %>%
#'       tibble::rownames_to_column("reaction") %>%
#'       dplyr::mutate(reaction_description = stringr::str_split(
#'         long_reaction_description, ";", simplify = TRUE, n = 2)[,1]) %>%
#'       #dplyr::select(reaction, reaction_description) %>%
#'       dplyr::mutate(reaction = paste0('rn:', reaction))
#'
#'   return(KEGG_reactions_df)
#' }
#'
#'
#' #' @title Prepare KEGG Visualization Data
#' #' @description Prepares data for visualizing a specified KEGG pathway.
#' #' @param kegg_data A dataframe of KEGG pathway data.
#' #' @param pathway_code A string specifying the KEGG pathway code.
#' #' @return A list containing the pathway number and reaction mappings.
#' prepare_kegg_visualization <- function(kegg_data, pathway_code) {
#'   # Get the pathway number
#'   pathway_number <- get_pathway_number(kegg_data, pathway_code)
#'
#'   # Fetch reaction mappings if not already available
#'   KEGG_reactions <- fetch_kegg_reactions()
#'
#'   # Return a list with necessary data for visualization
#'   return(list(pathway_number = pathway_number, KEGG_reactions = KEGG_reactions))
#' }
#'
#'
#' #' @title Create ggkegg Object for KEGG Pathway Visualization
#' #' @description Creates a ggkegg object from KEGG results using specified parameters.
#' #' @param kegg_results A dataframe containing KEGG pathway results.
#' #' @param pathway_number The number representing the pathway in KEGG.
#' #' @param organism_code A string representing the organism code in KEGG (e.g., "sce").
#' #' @return A ggkegg object.
#' create_kegg_ggkegg_object <- function(kegg_results, pathway_code, organism_code = "sce") {
#'   KEGG_data <- kegg_results %>%
#'     ggkegg(
#'       convert_first = FALSE,
#'       convert_collapse = "\n",
#'       pathway_number = pathway_code,
#'       convert_org = c(organism_code),
#'       delete_zero_degree = TRUE,
#'       return_igraph = FALSE
#'     )
#'   return(KEGG_data)
#' }
#'
#'
#'
#' #' @title Process ggkegg Data for Visualization
#' #' @description Processes ggkegg object data to prepare it for visualization, enhancing gene data with KEGG reactions and annotations.
#' #' @param KEGG_data A ggkegg object containing KEGG pathway data.
#' #' @param KEGG_reactions A dataframe containing KEGG reaction descriptions.
#' #' @param annotation_db An AnnotationDbi object for gene annotations.
#' #' @return A dataframe ready for plotting.
#' process_ggkegg_data_for_visualization <- function(KEGG_data, KEGG_reactions, annotation_db) {
#'   graph_data <- KEGG_data$data %>%
#'     filter(type == "gene") %>%
#'     mutate(showname = strsplit(name, " ") %>% str_remove_all("sce:")) %>%
#'     mutate(showname = gsub('c\\(|\\)|"|"|,', '', showname)) %>%
#'     separate_rows(showname, sep = " ") %>%
#'     left_join(rownames_to_column(KEGG_reactions), by = c("showname" = "rowname")) %>%
#'     left_join(KEGG_reactions, by = "reaction") %>%
#'     mutate(gene_name = AnnotationDbi::select(annotation_db, keys = showname, columns = "GENENAME")$GENENAME) %>%
#'     mutate(gene_name = coalesce(gene_name, showname))
#'
#'   return(graph_data)
#' }
#'
#'
#' library(ggplot2)
#' library(ggrepel)
#' library(RColorBrewer)
#' library(patchwork)
#'
#' #' @title Create KEGG Pathway Graphs
#' #' @description Creates comparative graphs for gene expression data on KEGG pathways.
#' #' @param graph_data A dataframe with KEGG pathway data processed for visualization.
#' #' @param pathway_id A string indicating the pathway ID for overlay.
#' #' @param fc_column1 A string naming the first fold change column.
#' #' @param fc_column2 A string naming the second fold change column.
#' #' @param titles A list of two strings for titles of the plots.
#' #' @return A ggplot object displaying the pathway graphs side by side.
#' create_kegg_pathway_graphs <- function(graph_data, pathway_id, fc_column1, fc_column2, titles) {
#'   # Find maximum fold change values for color scale
#'   max_fc <- ceiling(max(c(
#'     abs(graph_data[[fc_column1]]),
#'     abs(graph_data[[fc_column2]])
#'   ), na.rm = TRUE))
#'
#'   # Create the first graph
#'   graph1 <- ggplot(graph_data, aes(x = x, y = y)) +
#'     overlay_raw_map(pathway_id) +
#'     geom_label_repel(
#'       aes(label = gene_name, fill = get(fc_column1)),
#'       box.padding = 0.05,
#'       label.padding = 0.05,
#'       direction = "y",
#'       size = 2,
#'       max.overlaps = 100,
#'       label.r = 0.002,
#'       seed = 123
#'     ) +
#'     theme_void() +
#'     scale_fill_gradientn(
#'       colours = rev(brewer.pal(n = 10, name = "RdYlBu")),
#'       limits = c(-max_fc, max_fc)
#'     ) +
#'     theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
#'     labs(fill = "logFC", title = titles[1])
#'
#'   # Create the second graph
#'   graph2 <- ggplot(graph_data, aes(x = x, y = y)) +
#'     overlay_raw_map(pathway_id) +
#'     geom_label_repel(
#'       aes(label = gene_name, fill = get(fc_column2)),
#'       box.padding = 0.05,
#'       label.padding = 0.05,
#'       direction = "y",
#'       size = 2,
#'       max.overlaps = 100,
#'       label.r = 0.002,
#'       seed = 123
#'     ) +
#'     theme_void() +
#'     guides(fill = FALSE) +
#'     scale_fill_gradientn(
#'       colours = rev(brewer.pal(n = 10, name = "RdYlBu")),
#'       limits = c(-max_fc, max_fc)
#'     ) +
#'     theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
#'     labs(fill = "logFC", title = titles[2])
#'
#'   # Combine graphs side by side
#'   side_by_side_graph <- graph1 + graph2 + plot_annotation(tag_levels = 'A')
#'
#'   return(side_by_side_graph)
#' }
#'
#'
#'
#'
#'
#' # library(AnnotationDbi)
#'
#' #' @title Full Workflow for KEGG Pathway Visualization
#' #' @description Executes a complete workflow from raw DE results to visualizing them on a KEGG pathway.
#' #' @param de_results The differential expression results (from edgeR, limma, DESeq2).
#' #' @param organism_code The organism code for KEGG pathways (e.g., 'sce').
#' #' @param annotation_package The Bioconductor package for annotations relevant to the organism.
#' #' @param keytype The type of key used in the gene IDs (e.g., "SYMBOL", "ENTREZID", "ENSEMBL").
#' #' @param pathway_code The KEGG pathway code to visualize.
#' #' @param fc_column1 The first fold change column name for visualization.
#' #' @param fc_column2 The second fold change column name for visualization.
#' #' @param titles A vector of two strings for the titles of the plots.
#' #' @return A ggplot object displaying the pathway graphs.
#' full_kegg_visualization_workflow <- function(de_results, organism_code, annotation_package, keytype,
#'                                              pathway_code, fc_column1, fc_column2, titles) {
#'   # Step 1: Preprocess DE results to standard format
#'   kegg_results <- preprocess_DE_results(de_results)
#'
#'   # Step 2: Convert Gene IDs to KEGG IDs
#'   kegg_results <- homogenize_gene_ids(kegg_results, "GeneID", organism_code, annotation_package, fromType = keytype)
#'
#'   # Step 3: Fetch KEGG reaction mappings (if needed)
#'   KEGG_reactions <- fetch_kegg_reactions(organism_code = "sce")
#'
#'   return(head(KEGG_reactions))
#'   # # Step 4: Create ggkegg object for specified pathway
#'   # pathway_number <- get_pathway_number(kegg_results, pathway_code)
#'   # KEGG_data <- create_kegg_ggkegg_object(kegg_results, pathway_code, organism_code)
#'   #
#'   # # Step 5: Process ggkegg data for visualization
#'   # graph_data <- process_ggkegg_data_for_visualization(KEGG_data, KEGG_reactions, get(annotation_package))
#'   #
#'   # # Step 6: Create the KEGG pathway graphs
#'   # side_by_side_graph <- create_kegg_pathway_graphs(graph_data, pathway_code, fc_column1, fc_column2, titles)
#'   #
#'   # # Return the final graph
#'   # return(side_by_side_graph)
#' }
#'
#'
#'
#'
#'
#'
#' #'
#' #'
#' #'
#' #' #' @title Process and Map DE Results to KEGG IDs
#' #' #' @description Processes differential expression results and maps the gene IDs to KEGG IDs.
#' #' #' @param results DE analysis results.
#' #' #' @param organism_code KEGG organism code.
#' #' #' @param annotation_package Annotation package for ID conversion.
#' #' #' @param keytype Type of key used in the gene IDs.
#' #' #' @return A dataframe with mapped KEGG IDs.
#' #' process_and_map_to_kegg <- function(results, organism_code, annotation_package, keytype) {
#' #'   # Preprocess the DE results
#' #'   df <- preprocess_DE_results(results)
#' #'
#' #'   # Convert Gene IDs to KEGG IDs
#' #'   df_with_kegg <- convert_dataframe_to_kegg_ids(df, gene_id_column = "GeneID", organism_code, annotation_package, keytype)
#' #'
#' #'   return(df_with_kegg)
#' #' }
#' #'
#' #'
#' #'
#' #'
#' #'
#' #'
#' #'
#' #'
#' #' #' @title Preprocess edgeR topTags Results and Visualize KEGG Pathway
#' #' #' @description Main function to preprocess edgeR topTags results, load KEGG data, and visualize selected pathway for any organism.
#' #' #' @param topTags_result An object returned by edgeR's topTags function or a dataframe with equivalent structure.
#' #' #' @param pathway_to_graph A string indicating the KEGG pathway ID.
#' #' #' @param organism_code A string specifying the organism code in KEGG (e.g., 'sce' for yeast).
#' #' #' @return A plot or graph representation of the KEGG pathway.
#' #' #' @export
#' #' process_and_visualize_kegg <- function(topTags_result, pathway_to_graph = "sce00010", organism_code = "sce") {
#' #'   kegg_results <- preprocess_topTags(topTags_result, organism_code)
#' #'   return(KEGG_visualize_pathway(kegg_results, pathway_to_graph, organism_code))
#' #' }
#' #'
#' #'
#' #'
#' #' # Existing functions here: get_pathway_number, fetch_kegg_reactions, prepare_kegg_data, process_graph_data
#' #'
#' #' # Example of other needed functions for completion
#' #' #' @title Plot KEGG Graph
#' #' #' @description Plots the processed KEGG graph data.
#' #' #' @param graph_data Processed graph data for plotting.
#' #' #' @return A plot object.
#' #' plot_graph <- function(graph_data) {
#' #'   # Placeholder for actual graph plotting logic, possibly using ggraph or similar.
#' #'   print("Graph plotting not implemented.")
#' #' }
#' #'
#' #' # Ensure to include or reference all existing functions adjusted in previous steps.
