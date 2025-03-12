# library(edgeR)
# library(dplyr)
# library(tibble)
# library(KEGGREST)
# library(ggkegg)
# library(AnnotationDbi)
# library(clusterProfiler)


#' @title load Input data (enrichKEGG object)
#' @description Processes outputs from enrichKEGG.
#' @param enrich_results enrichResult Object from enrichKEGG.
#' @return A dataframe with results for subsequent analysis.
#' @export
loadEnrichResults <- function(enrich_results) {
  # Determine the class of the results object and process accordingly
  if (class(enrich_results)[1] == "enrichResult") {
    # edgeR topTags output
    kegg_results <- enrich_results
    return(kegg_results)
  } else {
    stop("Unsupported input type. Please provide an object clusterProfiler::enrichKEGG.")
  }
}

#' @title load Input data (enrichKEGG object)
#' @description Processes outputs from enrichKEGG.
#' @param DE_results gene results object, either edgeR, DESeq2, limma, or data frame.
#' @param FC_col_name A string specifying the column name for log fold change.
#' @return A dataframe with results for subsequent analysis.
#' @export
loadGeneData <- function(DE_results, FC_col_name) {
  # Determine the class of the results object and process accordingly
  if (attr(class(DE_results),"package") == "edgeR") {
    # edgeR topTags output
    DE_results <- edgeR::topTags(DE_results, n = Inf) %>%
      as.data.frame() %>%
        dplyr::mutate(GeneID = ORF, logFC = get(FC_col_name), #logFC2=get(logFC2),
                      pValue = PValue, adjPValue = FDR)
    return(DE_results)
    } else if (attr(class(DE_results),"package") == "DESeq2") {
      # DESeq2 results output
      data <- as.data.frame(DE_results) %>%
        dplyr::mutate(GeneID = ORF, logFC = log2FoldChange, pValue = pvalue, adjPValue = padj)
    } else if (attr(class(DE_results),"package") == "limma") {
      # limma results output
      data <- as.data.frame(DE_results) %>%
        dplyr::mutate(GeneID = ORF, logFC = FC_col_name, pValue = P.Value, adjPValue = adj.P.Val)
    } else if (is.data.frame(DE_results)) {
      # topTags output from edgeR
      data <- as.data.frame(DE_results) %>%
        dplyr::mutate(GeneID = rownames(.), logFC = FC_col_name, pValue = pValue, adjPValue = adjPValue)
  } else {
    stop("Unsupported input type. Please provide an object from edgeR, DESeq2, or limma.")
  }
}

#' @title Fetch KEGG Reactions
#' @description Retrieves KEGG reaction mappings and formats them for use.
#' @return A dataframe containing KEGG reactions with their descriptions.
#' @export
fetchKeggReactions <- function() {
  if (!exists("KEGG_reactions")) {
    KEGG_reactions <- KEGGREST::keggList("reaction")
  }
  KEGG_reactions_df_tmp <- as.data.frame(KEGG_reactions)
  colnames(KEGG_reactions_df_tmp) <- "long_reaction_description"
  KEGG_reactions_df <- KEGG_reactions_df_tmp %>%
    tibble::rownames_to_column("reaction") %>%
    dplyr::mutate(reaction_description = stringr::str_split(
      long_reaction_description, ";", simplify = TRUE, n = 2)[,1]) %>%
    #dplyr::select(reaction, reaction_description) %>%
    dplyr::mutate(reaction = paste0('rn:', reaction))

  return(KEGG_reactions_df)
}


#' @title Get Pathway Number to match row of enrichKEGG to KEGG pathway of interest
#' @description Retrieves the pathway number based on a pathway ID from KEGG data.
#' @param kegg_results An enrichKEGG output from clusterProfiler (provided by user).
#' @param pathway_id A string specifying the pathway ID (provided by user).
#' @return An integer representing the pathway number.
#' @export
getPathwayNumber <- function(kegg_results, pathway_id) {
  pathway_number <- kegg_results %>%
    data.frame() %>%
    dplyr::mutate(row_number = dplyr::row_number()) %>%
    dplyr::filter(ID == pathway_id) %>%
    dplyr::pull(row_number)
  return(pathway_number)
}

#' @title Create ggkegg Object for KEGG Pathway Visualization
#' @description Creates a ggkegg object from KEGG results using specified parameters.
#' @param kegg_results An enrichKEGG output from clusterProfiler (provided by user).
#' @param pathway_id A string specifying the pathway ID (provided by user).
#' @param pathway_number The number representing the pathway in KEGG.
#' @param organism_code A string representing the organism code in KEGG (e.g., "sce").
#' @return A ggkegg object.
#' @export
create_ggkegg <- function(kegg_results, pathway_id, pathway_number, organism_code = "sce") {
  KEGG_data <- # %>%
    ggkegg::ggkegg(kegg_results,
      pid = pathway_id,
      layout = "native",
      return_tbl_graph = TRUE,
      convert_first = FALSE,
      convert_collapse = "\n",
      convert_reaction = FALSE,
      delete_zero_degree = FALSE,
      delete_undefined= FALSE,
      pathway_number = pathway_number,
      convert_org = c(organism_code),
      delete_zero_degree = TRUE,
      numeric_attribute = NULL,
      return_igraph = FALSE,
      node_rect_nudge=0,
      group_rect_nudge=2,
      modulte_type="definition",
      module_definition_type="text"
    )
  return(KEGG_data)
}


#' @title Process ggkegg Data for Visualization
#' @description Processes ggkegg object data to prepare it for visualization, enhancing gene data with KEGG reactions and annotations.
#' @param KEGG_data A ggkegg object containing KEGG pathway data.
#' @param topTags_results A dataframe containing gene results.
#' @param organism_code A string representing the organism code in KEGG (e.g., "sce").
#' @param KEGG_reactions_df A dataframe containing KEGG reaction descriptions.
#' @param annotation_db An AnnotationDbi object for gene annotations.
#' @return A dataframe ready for plotting.
#' @export
process_ggkegg <- function(KEGG_data, topTags_results, KEGG_reactions_df, annotation_db, organism_code = "sce") {
  org_db <- #get(as.character(
    annotation_db#))

  graph_data <- KEGG_data$data %>%
    # filter(type == "gene") %>%
    dplyr::filter(type %in% c("gene", "ortholog")) %>% #pull(name)
    # mutate(showname = strsplit(name, " ") %>% str_remove_all("sce:")) %>%
      # pull(showname)
    dplyr::mutate(showname = #strsplit(name, " ")  %>%
                    stringr::str_remove_all(name, paste0(organism_code,":"))) %>%
    # dplyr::mutate(showname = strsplit(name, " ")) %>%
    # dplyr::mutate(showname = gsub(paste0(organism_code, ":"), "", showname)) %>% pull(showname)
    # dplyr::mutate(showname = gsub('c\\(|\\)|"|"|,', '', showname)) %>%
    tidyr::separate_rows(showname, sep = " ") %>%#pull(showname)
    dplyr::left_join(topTags_results, by = c("showname" = "ORF")) %>%
    dplyr::left_join(KEGG_reactions_df, by = "reaction") #%>%

  gene_name_df <- suppressMessages(
    AnnotationDbi::select(#org_db,
                                        org.Sc.sgd.db,
                                        keys = graph_data$showname,
                                        columns = "GENENAME")
  )

  gene_name <- gene_name_df$GENENAME

  graph_data <- graph_data %>%
    dplyr::mutate(gene_name = gene_name) %>%
    dplyr::mutate(gene_name = dplyr::coalesce(gene_name, showname))

  return(graph_data)
}


#' @title Create KEGG Pathway Graphs
#' @description Creates comparative graphs for gene expression data on KEGG pathways.
#' @param graph_data A dataframe with KEGG pathway data processed for visualization.
#' @param pathway_id A string indicating the pathway ID for overlay.
#' @param fc_column A string naming the first fold change column.
#' @param juke A number for amount of space to separate names. Default is 20.
#' @param title A title for the plot.
#' @return A ggplot object displaying the pathway graphs side by side.
#' @export
createParalogPlots <- function(graph_data, pathway_id,
                                       fc_column, #fc_column2,
                                       juke = 20,
                                       title) {
  # Find maximum fold change values for color scale
  max_fc <- ceiling(max(c(
    abs(graph_data[[fc_column]])
  ), na.rm = TRUE))

  # Create the first graph
  # graph1 <- ggplot(graph_data, aes(x = x, y = y)) +
  #   overlay_raw_map(pathway_id, high_res = FALSE) +
  #   # geom_label(aes(x=x, y=y, label = gene_name, fill=get(fc_column1)),
  #   #            size = 2, color = "black", position=position_fill(vjust = 0.5,
  #   #                                                              reverse = FALSE)) +
  #   geom_label_repel(
  #     aes(label = gene_name, fill = get(fc_column2), group = name),
  #     point.padding = NA,#-0.75,
  #     box.padding = 0.000,
  #     label.padding = 0.15,
  #     #position = position_stack(vjust = -0.5),
  #     direction = "y",
  #     force=0.5,
  #     size = 2,
  #     max.overlaps = 100,
  #     label.r = 0.002,
  #     seed = 123
  #   ) +
  #   theme_void() +
  #   scale_fill_gradientn(
  #     colours = rev(brewer.pal(n = 10, name = "RdYlBu")),
  #     limits = c(-max_fc, max_fc)
  #   ) +
  #   theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  #   labs(fill = "logFC", title = titles[1])

  # Create the second graph
  paralogPlot <-
    graph_data %>%
    group_by(x,y) %>%
    dplyr::mutate(id_in_group = dplyr::row_number(),
           group_size = dplyr::n()) %>%
    dplyr::mutate(tmp.y = dplyr::case_when(
      id_in_group == 1 ~ y,
      id_in_group == 2 ~ y + juke,
      id_in_group == 3 ~ y - juke,
      id_in_group == 4 ~ y,
      id_in_group == 5 ~ y + juke,
      id_in_group == 6 ~ y - juke,
      id_in_group == 7 ~ y,
      id_in_group == 8 ~ y + juke,
      id_in_group == 9 ~ y - juke,
      id_in_group %% 2 == 0 ~ y + (id_in_group-1)*juke,
      id_in_group %% 2 != 0 ~ y - (id_in_group-1)*juke
      # TRUE ~ y
    ),
    tmp.x = dplyr::case_when(
      id_in_group == 1 ~ x,
      id_in_group == 2 ~ x,
      id_in_group == 3 ~ x,
      id_in_group %% 2 == 0 ~ x + 2.5*juke,
      id_in_group %% 2 != 0 ~ x - 2.5*juke
      # TRUE ~ x
    )) %>%
    dplyr::mutate(tmp.y = dplyr::case_when(
      group_size == 2 ~ tmp.y-juke/2,
      TRUE ~ tmp.y)
) %>%
    # adjust for high res image
    dplyr::mutate(tmp.x = tmp.x*2,
           tmp.y = tmp.y*2) %>%
    ggplot(., aes(x = tmp.x, y = tmp.y)) +
    ggkegg::overlay_raw_map(paste0("map",gsub("[^0-9.-]", "", pathway_id)),
                    high_res = TRUE) +
    # geom_label_repel(
    #   aes(label = gene_name, fill = get(fc_column2), group = name),
    #   box.padding = 0.5,
    #   point.padding = -0.1,
    #   label.padding = 0.1,
    #   # position = position_stack(vjust = 0.5),
    #   direction = "y",
    #   nudge_y = .1,
    #   force=0.1,
    #   size = 2,
    #   max.overlaps = 100,
    #   label.r = 0.002,
    #   seed = 123
    # ) +

  geom_tile(aes(x = x*2, y=y*2, height = 2.2*20, width=5.5*20), fill = "white")+
  # geom_label(aes(x=x*2, y=y*2, label = "     ", fill = 0, color = 0),
  #                         label.padding =unit(0.1, "lines"),
  #            label.size = 0,
  #                         max.overlaps=200,
  #                       label.r = unit(0.02, "lines"),
  #                       size = 2) +

  ## this is almost perfect ##
  geom_label(aes(label = gene_name,
                 fill = get(fc_column)),
                 #point.padding = NA,
               label.padding =unit(0.1, "lines"),
               #box.padding = 0.001,
               #max.overlaps=200,
             label.r = unit(0.02, "lines"),
             size = 1.5) +

  ## this is the one I think is better if I use repel...
  # geom_label_repel(
  #   aes(label = gene_name, fill = get(fc_column2)),
  #   point.padding = NA,
  #   box.padding = 0.001,
  #   label.padding = 0.1,
  #   # position = position_stack(vjust = 0.5),
  #   direction = "x",
  #   force=0.5,
  #   force_pull = 5,
  #   size = 2,
  #   max.overlaps = 100,
  #   label.r = 0.002,
  #   seed = 123
  # ) +
    theme_void() +
    # guides(fill = FALSE) +
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")),
      limits = c(-max_fc, max_fc)
    ) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    labs(fill = "logFC", title = title)

  # Combine graphs side by side
  # side_by_side_graph <- #graph1 +
  #   graph2 +
  #   plot_annotation(tag_levels = 'A')

  return(paralogPlot)

  #return(graph2)
}


#' @title Full Workflow for KEGG Pathway Visualization
#' @description Executes a complete workflow from raw DE results to visualizing them on a KEGG pathway.
#' @param enrich_results The enrichKEGG results (from clusterProfiler).
#' @param DE_results The DE results (e.g. from edgeR).
#' @param organism_code The organism code for KEGG pathways (e.g., 'sce').
#' @param annotation_package The Bioconductor package for annotations relevant to the organism.
#' @param keytype The type of key used in the gene IDs (e.g., "SYMBOL", "ENTREZID", "ENSEMBL").
#' @param pathway_code The KEGG pathway code to visualize (e.g. 'sce00020')
#' @param fc_column The first fold change column name for visualization.
#' @param juke A number for amount of space to separate names. Default is 20.
#' @param title A title for the plot.
#' @return A ggplot object displaying the pathway graphs.
#' @export
plotParalogs <- function(enrich_results, DE_results, organism_code="sce",
                         annotation_package = org.Sc.sgd.db,
                         keytype="ORF",
                         pathway_code,
                         fc_column="logFC",
                         juke=20,
                         title=NULL) {

  # Step 1: load KEGG reaction data from keggList
  kegg_results <- loadEnrichResults(enrich_results)
  topTags_results <- loadGeneData(DE_results, fc_column#, fc_column2
                                  )

  # Step 2: Convert Gene IDs to KEGG IDs
  # kegg_results <- homogenize_gene_ids(kegg_results, "GeneID", organism_code, annotation_package, fromType = keytype)
      # don't have to do step 2 anymore if already the result of enrichKEGG (had to be done by user.)

  # Step 3: Fetch KEGG reaction mappings (if needed)
  KEGG_reactions_df <- fetchKeggReactions()

  # return(c(head(KEGG_reactions), head(kegg_results)))
  # Step 4: Create ggkegg object for specified pathway
  pathway_number <- getPathwayNumber(kegg_results, pathway_code)
  KEGG_data <- create_ggkegg(kegg_results, pathway_code, organism_code)
  # return(KEGG_data)

  # # Step 5: Process ggkegg data for visualization
  graph_data <-
    process_ggkegg(KEGG_data,
                   topTags_results,
                   KEGG_reactions_df,
                   # AnnotationDbi::get(
                   annotation_package,
                   organism_code)

  # return(graph_data)

  # Step 6: Create the KEGG pathway graphs
  paralogPlot <- createParalogPlots(graph_data, pathway_code, fc_column, #fc_column2,
                                    juke, title)

  # Return the final graph
  return(paralogPlot)
}

