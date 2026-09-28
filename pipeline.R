# ==============================================================================
# pipeline.R
# Orestias skin microbiota
#
# Analysis pipeline consolidated from:
#   - qiime2R_16S_pipeline.R
#   - Co-ocurrence_network_pipeline.R
#   - integration_picrust2_and_modules.R
#   - ggpicrust2_pipeline.R
#
# QIIME2 preprocessing was performed in Galaxy. This script starts from the
# resulting QIIME2 artifacts and downstream tables.
#
# Input files are expected to remain local and are NOT part of the repository.
# Intermediate objects are saved to local_objects/ for use by figures.R.
# ==============================================================================

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

# ------------------------------------------------------------------------------
# 0. Configuration
# ------------------------------------------------------------------------------

FEATURE_TABLE_QZA <- "clustered_table.qza"
TAXONOMY_QZA      <- "taxonomy_classifications.qza"
METADATA_FILE     <- "samples_table.txt"

PICRUST_ABUNDANCE <- "Orestias_KO_pred_metagenome_unstrat.txt"
PICRUST_CONTRIB   <- "Galaxy107-[KO_pred_metagenome_contrib].tabular"
MODULE_TABLE      <- "Oasc_genus_coocurrence_node_table.csv"
PATHWAY_ABUNDANCE <- "PICRUST2_Orestias_pathway_abundace.txt"

# The original project contains two SPIEC-EASI workflows:
# (1) ASV-level, prevalence-filtered networks
# (2) genus-level networks used for the Cytoscape/module-functional integration.
# Both are preserved, but the ASV-level workflow is optional because it is
# computationally intensive.
RUN_ASV_NETWORK   <- FALSE
RUN_GENUS_NETWORK <- TRUE

dir.create("local_objects", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(qiime2R)
  library(phyloseq)
  library(metacoder)
  library(tidyverse)
  library(ggplot2)
  library(vegan)
  library(agricolae)
  library(ggnewscale)
  library(igraph)
  library(ggraph)
  library(SpiecEasi)
  library(data.table)
  library(KEGGREST)
  library(ade4)
  library(patchwork)
})

# ggpicrust2 is optional: the rest of the pipeline can run without it.
HAS_GGPICRUST2 <- requireNamespace("ggpicrust2", quietly = TRUE)
if (HAS_GGPICRUST2) {
  library(ggpicrust2)
}

# RCy3 is optional because Cytoscape is an external application.
HAS_RCY3 <- requireNamespace("RCy3", quietly = TRUE)
if (HAS_RCY3) {
  library(RCy3)
}

# ------------------------------------------------------------------------------
# 2. Functions
# ------------------------------------------------------------------------------

phyloseq_ntaxa_by_tax <- function(physeq,
                                  TaxRank = "Phylum",
                                  relative = FALSE,
                                  add_meta_data = TRUE) {
  mm <- phyloseq::psmelt(physeq)

  count_otus <- function(z, TaxRank = TaxRank, relative = relative) {
    z <- z[z$Abundance > 0, , drop = FALSE]
    rez <- as.data.frame(table(z[, TaxRank]), stringsAsFactors = FALSE)
    colnames(rez) <- c(TaxRank, "N.OTU")

    if (relative) {
      rez$N.OTU <- with(rez, N.OTU / sum(N.OTU))
    }

    rez
  }

  res <- plyr::ddply(
    .data = mm,
    .variables = "Sample",
    .fun = count_otus,
    TaxRank = TaxRank,
    relative = relative
  )

  if (add_meta_data && !is.null(sample_data(physeq, errorIfNULL = FALSE))) {
    metad <- data.frame(
      Sample = sample_names(physeq),
      sample_data(physeq)
    )
    meta_cols <- setdiff(colnames(metad), "Sample")
    res <- cbind(
      res,
      metad[match(res$Sample, metad$Sample), meta_cols, drop = FALSE]
    )
  }

  rownames(res) <- NULL
  res
}

filter_taxa_prevalence <- function(ps, prev = 0.25) {
  keep <- apply(
    as(otu_table(ps), "matrix"),
    1,
    function(x) sum(x > 0) >= (prev * nsamples(ps))
  )
  prune_taxa(keep, ps)
}

safe_write_csv <- function(x, file) {
  write.csv(x, file = file, row.names = FALSE)
}

# ------------------------------------------------------------------------------
# 3. Import QIIME2 results and metadata
# ------------------------------------------------------------------------------

message("Reading QIIME2 artifacts...")

ps <- qza_to_phyloseq(
  features = FEATURE_TABLE_QZA,
  taxonomy = TAXONOMY_QZA
)

metadata <- read.table(
  METADATA_FILE,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if ("Sample" %in% colnames(metadata)) {
  rownames(metadata) <- metadata$Sample
} else {
  stop("Metadata must contain a column named 'Sample'.")
}

# Keep only metadata samples represented in the feature table.
common_samples <- intersect(sample_names(ps), rownames(metadata))
if (length(common_samples) == 0) {
  stop("No sample IDs are shared between the QIIME2 artifact and metadata.")
}

ps <- prune_samples(common_samples, ps)
metadata <- metadata[common_samples, , drop = FALSE]

sample_data(ps) <- sample_data(metadata)

# ------------------------------------------------------------------------------
# 4. Sequence-depth normalization and bacterial filtering
# ------------------------------------------------------------------------------

median_depth <- median(sample_sums(ps), na.rm = TRUE)

standf <- function(x, t = median_depth) {
  if (sum(x) == 0) return(x)
  round(t * (x / sum(x)))
}

ps_norm <- transform_sample_counts(ps, standf)

if ("Kingdom" %in% rank_names(ps_norm)) {
  bacteria <- subset_taxa(ps_norm, Kingdom == "Bacteria")
} else {
  stop("Taxonomy table does not contain a 'Kingdom' rank.")
}

bacteria <- prune_taxa(taxa_sums(bacteria) > 0, bacteria)

# Preserve sample order used in the original script when available.
preferred_order <- c("A5","A6","A7","A8","C3","C6","C7","C8")
present_order <- preferred_order[preferred_order %in% sample_names(bacteria)]
remaining <- setdiff(sample_names(bacteria), present_order)

sample_data(bacteria)$NewID <- factor(
  sample_names(bacteria),
  levels = c(present_order, remaining)
)

# ------------------------------------------------------------------------------
# 5. Taxonomic composition
# ------------------------------------------------------------------------------

PhylumCounts <- phyloseq_ntaxa_by_tax(bacteria, "Phylum")
ClassCounts  <- phyloseq_ntaxa_by_tax(bacteria, "Class")
OrderCounts  <- phyloseq_ntaxa_by_tax(bacteria, "Order")
FamilyCounts <- phyloseq_ntaxa_by_tax(bacteria, "Family")
GenusCounts  <- phyloseq_ntaxa_by_tax(bacteria, "Genus")

# ------------------------------------------------------------------------------
# 6. Alpha and beta diversity
# ------------------------------------------------------------------------------

alpha_diversity <- estimate_richness(
  bacteria,
  measures = c("Observed", "Shannon", "InvSimpson")
)

alpha_diversity$Sample <- rownames(alpha_diversity)
alpha_diversity <- left_join(
  alpha_diversity,
  metadata %>% rownames_to_column("Sample"),
  by = "Sample"
)

bray_dist <- phyloseq::distance(bacteria, method = "bray")
nmds <- ordinate(bacteria, method = "NMDS", distance = "bray")

metadata_adonis <- as.data.frame(sample_data(bacteria))

permanova <- adonis2(
  bray_dist ~ Host,
  data = metadata_adonis,
  permutations = 999
)

dispersion <- betadisper(
  bray_dist,
  metadata_adonis$Host
)

dispersion_anova <- anova(dispersion)
dispersion_permutation <- permutest(dispersion, permutations = 999)

# Original scripts also evaluated phylum/order richness.
metadata_stats <- metadata_adonis %>%
  rownames_to_column("Sample")

phylum_richness <- PhylumCounts %>%
  count(Sample, name = "Phylum_count")

order_richness <- OrderCounts %>%
  count(Sample, name = "Order_count")

metadata_stats <- metadata_stats %>%
  left_join(phylum_richness, by = "Sample") %>%
  left_join(order_richness, by = "Sample")

if ("Host" %in% colnames(metadata_stats)) {
  phylum_anova <- aov(Phylum_count ~ Host, data = metadata_stats)
  order_anova  <- aov(Order_count ~ Host, data = metadata_stats)

  phylum_tukey <- HSD.test(phylum_anova, "Host", group = TRUE)
  order_tukey  <- HSD.test(order_anova, "Host", group = TRUE)
} else {
  phylum_anova <- order_anova <- NULL
  phylum_tukey <- order_tukey <- NULL
}

# ------------------------------------------------------------------------------
# 7. Metacoder
# ------------------------------------------------------------------------------

metacoder_input <- NULL
meta_obj <- NULL
meta_obj_filt <- NULL
meta_obj_A <- NULL
meta_obj_G <- NULL

# Build a metacoder input directly from the QIIME2-derived phyloseq object.
otu_mat <- as.matrix(otu_table(bacteria))
tax_mat <- as.data.frame(as(tax_table(bacteria), "matrix"))

tax_mat$ASV <- rownames(tax_mat)
otu_df <- as.data.frame(otu_mat)
otu_df$ASV <- rownames(otu_df)

metacoder_input <- left_join(tax_mat, otu_df, by = "ASV")

for (r in c("Kingdom","Phylum","Class","Order","Family","Genus")) {
  if (!r %in% colnames(metacoder_input)) {
    metacoder_input[[r]] <- "Unknown"
  }
  metacoder_input[[r]][is.na(metacoder_input[[r]]) |
                         metacoder_input[[r]] == ""] <- "Unknown"
}

metacoder_input$taxonomy <- paste0(
  "k__", metacoder_input$Kingdom,
  ";p__", metacoder_input$Phylum,
  ";c__", metacoder_input$Class,
  ";o__", metacoder_input$Order,
  ";f__", metacoder_input$Family,
  ";g__", metacoder_input$Genus
)

meta_tax <- metacoder_input %>%
  select(ASV, taxonomy) %>%
  distinct(taxonomy, .keep_all = TRUE) %>%
  column_to_rownames("ASV")

meta_obj <- parse_tax_data(
  meta_tax,
  class_cols = "taxonomy",
  class_sep = ";"
)

meta_obj <- filter_ambiguous_taxa(meta_obj)

sample_cols <- intersect(
  sample_names(bacteria),
  colnames(meta_obj$data$tax_data)
)

meta_obj$data$tax_abund <- calc_taxon_abund(
  meta_obj,
  data = "tax_data",
  cols = sample_cols
)

groups <- metadata[match(sample_cols, rownames(metadata)), "Host"]

if (length(unique(groups)) == 2 && all(!is.na(groups))) {
  meta_obj$data$diff_table <- compare_groups(
    meta_obj,
    data = "tax_abund",
    cols = sample_cols,
    groups = groups
  )
}

# Top-taxon heat tree object.
abund_df <- meta_obj$data$tax_abund %>%
  mutate(total_abund = rowSums(across(-taxon_id)))

top_n_meta <- min(150, nrow(abund_df))
keep_ids <- abund_df %>%
  arrange(desc(total_abund)) %>%
  slice_head(n = top_n_meta) %>%
  pull(taxon_id)

meta_obj_filt <- filter_taxa(
  meta_obj,
  taxon_ids %in% keep_ids,
  subtaxa = FALSE
)

meta_obj_filt$data$tax_abund <- calc_taxon_abund(
  meta_obj_filt,
  data = "tax_data",
  cols = sample_cols
)

if (length(unique(groups)) == 2 && all(!is.na(groups))) {
  meta_obj_filt$data$diff_table <- compare_groups(
    meta_obj_filt,
    data = "tax_abund",
    cols = sample_cols,
    groups = groups
  )

  vals <- meta_obj_filt$data$diff_table$log2_median_ratio
  finite_vals <- vals[is.finite(vals)]

  if (length(finite_vals) > 0) {
    vals[vals == Inf] <- max(finite_vals)
    vals[vals == -Inf] <- min(finite_vals)
    meta_obj_filt$data$diff_table$log2_median_ratio <- vals
  }
}

# Host-specific metacoder objects.
if ("Host" %in% colnames(metadata)) {
  samples_A <- rownames(metadata)[metadata$Host == "O. ascotanensis"]
  samples_G <- rownames(metadata)[metadata$Host == "O. gloriae"]

  if (length(samples_A) > 0) {
    meta_obj_A <- filter_taxa(
      meta_obj,
      taxon_ids %in% keep_ids,
      subtaxa = FALSE
    )
    meta_obj_A$data$tax_abund <- calc_taxon_abund(
      meta_obj_A, data = "tax_data",
      cols = intersect(samples_A, sample_cols)
    )
  }

  if (length(samples_G) > 0) {
    meta_obj_G <- filter_taxa(
      meta_obj,
      taxon_ids %in% keep_ids,
      subtaxa = FALSE
    )
    meta_obj_G$data$tax_abund <- calc_taxon_abund(
      meta_obj_G, data = "tax_data",
      cols = intersect(samples_G, sample_cols)
    )
  }
}

# Alpha diversity from metacoder, retained because it was explicitly present
# in the original workflow.
if (!is.null(meta_obj)) {
  metadata$inv_simp <- diversity(
    meta_obj$data$tax_data[, sample_cols, drop = FALSE],
    index = "invsimpson",
    MARGIN = 2
  )[rownames(metadata)]

  metadata$Shannon <- diversity(
    meta_obj$data$tax_data[, sample_cols, drop = FALSE],
    index = "shannon",
    MARGIN = 2
  )[rownames(metadata)]

  invsimpson_anova <- aov(inv_simp ~ Host, metadata)
  shannon_anova <- aov(Shannon ~ Host, metadata)

  invsimpson_tukey <- HSD.test(invsimpson_anova, "Host", group = TRUE)
  shannon_tukey <- HSD.test(shannon_anova, "Host", group = TRUE)
}

# ------------------------------------------------------------------------------
# 8. SPIEC-EASI co-occurrence networks
# ------------------------------------------------------------------------------

se_asc <- se_glo <- NULL
net_asc <- net_glo <- NULL

# ASV-level workflow from Co-ocurrence_network_pipeline.R.
se_asc_asv <- se_glo_asv <- NULL
net_asc_asv <- net_glo_asv <- NULL

if (RUN_ASV_NETWORK) {
  ps_asc_asv <- subset_samples(bacteria, Host == "O. ascotanensis")
  ps_glo_asv <- subset_samples(bacteria, Host == "O. gloriae")

  ps_asc_asv <- prune_taxa(taxa_sums(ps_asc_asv) > 0, ps_asc_asv)
  ps_glo_asv <- prune_taxa(taxa_sums(ps_glo_asv) > 0, ps_glo_asv)

  ps_asc_asv <- filter_taxa_prevalence(ps_asc_asv, prev = 0.25)
  ps_glo_asv <- filter_taxa_prevalence(ps_glo_asv, prev = 0.25)

  se_asc_asv <- spiec.easi(
    ps_asc_asv, method = "mb",
    lambda.min.ratio = 1e-2,
    nlambda = 20,
    pulsar.params = list(rep.num = 50)
  )

  se_glo_asv <- spiec.easi(
    ps_glo_asv, method = "mb",
    lambda.min.ratio = 1e-2,
    nlambda = 20,
    pulsar.params = list(rep.num = 50)
  )

  net_asc_asv <- adj2igraph(getRefit(se_asc_asv))
  net_glo_asv <- adj2igraph(getRefit(se_glo_asv))
}

# Genus-level workflow from integration_picrust2_and_modules.R.
# This is the network used downstream for module/functional integration.
if (RUN_GENUS_NETWORK) {
  ps_asc <- subset_samples(bacteria, Host == "O. ascotanensis")
  ps_glo <- subset_samples(bacteria, Host == "O. gloriae")

  ps_asc <- prune_taxa(taxa_sums(ps_asc) > 0, ps_asc)
  ps_glo <- prune_taxa(taxa_sums(ps_glo) > 0, ps_glo)

  ps_asc <- filter_taxa_prevalence(ps_asc, prev = 0.25)
  ps_glo <- filter_taxa_prevalence(ps_glo, prev = 0.25)

  ps_asc_filt <- prune_taxa(taxa_sums(ps_asc) > 20, ps_asc)
  ps_glo_filt <- prune_taxa(taxa_sums(ps_glo) > 20, ps_glo)

  ps_asc_genus <- tax_glom(ps_asc_filt, taxrank = "Genus")
  ps_glo_genus <- tax_glom(ps_glo_filt, taxrank = "Genus")

  ps_asc_genus_filt <- ps_asc_genus %>%
    prune_taxa(taxa_sums(.) > 20, .) %>%
    filter_taxa(function(x) sum(x > 0) >= 2, TRUE)

  ps_glo_genus_filt <- ps_glo_genus %>%
    prune_taxa(taxa_sums(.) > 20, .) %>%
    filter_taxa(function(x) sum(x > 0) >= 2, TRUE)

  se_asc <- spiec.easi(
    ps_asc_genus_filt,
    method = "mb",
    lambda.min.ratio = 1e-2,
    nlambda = 20,
    pulsar.params = list(rep.num = 20)
  )

  se_glo <- spiec.easi(
    ps_glo_genus_filt,
    method = "mb",
    lambda.min.ratio = 1e-2,
    nlambda = 15,
    pulsar.params = list(rep.num = 20)
  )

  net_asc <- adj2igraph(getRefit(se_asc))
  net_glo <- adj2igraph(getRefit(se_glo))

  V(net_asc)$name <- taxa_names(ps_asc_genus_filt)
  V(net_glo)$name <- taxa_names(ps_glo_genus_filt)

  tax_asc <- as.data.frame(tax_table(ps_asc_genus_filt))
  tax_glo <- as.data.frame(tax_table(ps_glo_genus_filt))

  V(net_asc)$Phylum <- tax_asc$Phylum
  V(net_asc)$Genus <- tax_asc$Genus
  V(net_glo)$Phylum <- tax_glo$Phylum
  V(net_glo)$Genus <- tax_glo$Genus

  # Export edge/node tables for downstream Cytoscape work.
  edges_asc <- as.data.frame(as_edgelist(net_asc), stringsAsFactors = FALSE)
  colnames(edges_asc) <- c("from", "to")

  edges_glo <- as.data.frame(as_edgelist(net_glo), stringsAsFactors = FALSE)
  colnames(edges_glo) <- c("from", "to")

  nodes_asc <- data.frame(
    id = V(net_asc)$name,
    Phylum = V(net_asc)$Phylum,
    Genus = V(net_asc)$Genus
  )

  nodes_glo <- data.frame(
    id = V(net_glo)$name,
    Phylum = V(net_glo)$Phylum,
    Genus = V(net_glo)$Genus
  )

  safe_write_csv(edges_asc, "local_objects/Oasc_genus_edges.csv")
  safe_write_csv(edges_glo, "local_objects/Ogloriae_genus_edges.csv")
  safe_write_csv(nodes_asc, "local_objects/Oasc_genus_nodes.csv")
  safe_write_csv(nodes_glo, "local_objects/Ogloriae_genus_nodes.csv")
}

# Optional Cytoscape export when Cytoscape is running.
if (HAS_RCY3 && !is.null(net_asc) && !is.null(net_glo)) {
  try({
    if (RCy3::cytoscapePing()) {
      RCy3::createNetworkFromIgraph(
        net_asc,
        title = "O. ascotanensis",
        collection = "Orestias networks"
      )
      RCy3::createNetworkFromIgraph(
        net_glo,
        title = "O. gloriae",
        collection = "Orestias networks"
      )
    }
  }, silent = TRUE)
}

# ------------------------------------------------------------------------------
# 9. PICRUSt2 functional analysis with ggpicrust2
# ------------------------------------------------------------------------------

daa_annotated <- NULL
gsea_kegg <- NULL
gsea_kegg_annotated <- NULL
gsea_go <- NULL
gsea_go_annotated <- NULL
picrust_abundance <- NULL
picrust_metadata <- NULL

if (HAS_GGPICRUST2 && file.exists(PICRUST_ABUNDANCE)) {

  picrust_abundance <- readr::read_delim(
    PICRUST_ABUNDANCE,
    delim = "\t",
    col_names = TRUE,
    trim_ws = TRUE,
    show_col_types = FALSE
  ) %>%
    column_to_rownames("function")

  picrust_metadata <- readr::read_delim(
    METADATA_FILE,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    show_col_types = FALSE
  )

  daa_results <- ggpicrust2::pathway_daa(
    abundance = picrust_abundance,
    metadata = picrust_metadata,
    pathway = "KO",
    group = "Host",
    daa_method = "ALDEx2",
    ko_to_kegg = FALSE,
    p.adjust = "BH",
    select = NULL,
    reference = NULL
  )

  daa_sub_method_results <- daa_results[
    daa_results$method == "ALDEx2_Wilcoxon rank test", ,
    drop = FALSE
  ]

  daa_annotated <- ggpicrust2::pathway_annotation(
    pathway = "KO",
    daa_results_df = daa_sub_method_results
  )

  # KEGG GSEA
  gsea_kegg <- ggpicrust2::pathway_gsea(
    abundance = picrust_abundance,
    metadata = picrust_metadata %>% column_to_rownames("Sample"),
    group = "Host",
    pathway_type = "KEGG",
    method = "fgsea"
  )

  gsea_kegg_annotated <- ggpicrust2::gsea_pathway_annotation(
    gsea_results = gsea_kegg,
    pathway_type = "KEGG"
  )

  # GO biological-process GSEA
  gsea_go <- ggpicrust2::pathway_gsea(
    picrust_abundance,
    metadata = picrust_metadata %>% column_to_rownames("Sample"),
    group = "Host",
    pathway_type = "GO",
    go_category = "BP",
    method = "fgsea"
  )

  gsea_go_annotated <- ggpicrust2::gsea_pathway_annotation(
    gsea_results = gsea_go,
    pathway_type = "GO"
  )
}

# ------------------------------------------------------------------------------
# 10. Integration of network modules with PICRUSt2
# ------------------------------------------------------------------------------

module_func <- NULL
module_func_annot <- NULL
module_L2 <- NULL
module_L3 <- NULL
top_diff <- NULL
module_group <- NULL
ko_df <- NULL

if (
  RUN_GENUS_NETWORK &&
  file.exists(PICRUST_CONTRIB) &&
  file.exists(MODULE_TABLE)
) {

  contrib <- data.table::fread(PICRUST_CONTRIB)
  modules <- read.csv(MODULE_TABLE, check.names = FALSE)

  # Preserve the source-script naming convention.
  if ("X__glayCluster" %in% colnames(modules)) {
    modules$module <- modules$X__glayCluster
  }
  if (!"name" %in% colnames(modules)) {
    stop("Module table must contain a 'name' column.")
  }

  modules$name <- trimws(modules$name)
  contrib$taxon <- trimws(contrib$taxon)

  contrib_mod <- merge(
    contrib,
    modules,
    by.x = "taxon",
    by.y = "name"
  )

  contrib_mod$module <- contrib_mod$X__glayCluster

  if ("function" %in% colnames(contrib_mod)) {
    colnames(contrib_mod)[colnames(contrib_mod) == "function"] <- "KO"
  }

  module_func <- contrib_mod %>%
    group_by(module, KO) %>%
    summarise(
      abundance = sum(taxon_function_abun),
      .groups = "drop"
    ) %>%
    group_by(module) %>%
    mutate(rel_abundance = abundance / sum(abundance)) %>%
    ungroup()

  top_funcs <- module_func %>%
    group_by(module) %>%
    slice_max(rel_abundance, n = 10) %>%
    ungroup()

  # KEGG annotations are queried in chunks as in the original script.
  ko_ids <- unique(module_func$KO)
  chunks <- split(ko_ids, ceiling(seq_along(ko_ids) / 10))

  ko_info_list <- lapply(chunks, function(x) {
    tryCatch(KEGGREST::keggGet(x), error = function(e) NULL)
  })
  ko_info_list <- ko_info_list[!vapply(ko_info_list, is.null, logical(1))]

  if (length(ko_info_list) > 0) {
    ko_info <- do.call(c, ko_info_list)

    ko_df <- data.frame(
      KO = vapply(ko_info, function(x) x$ENTRY %||% NA_character_, character(1)),
      description = vapply(ko_info, function(x) x$NAME %||% NA_character_, character(1))
    )
  }

  ko2path <- KEGGREST::keggLink("pathway", "ko")
  ko2path_df <- data.frame(
    KO = sub("ko:", "", names(ko2path)),
    pathway = sub("path:", "", ko2path)
  )

  module_func_annot <- merge(module_func, ko2path_df, by = "KO")

  path_names <- KEGGREST::keggList("pathway")
  path_df <- data.frame(
    pathway = sub("path:", "", names(path_names)),
    name = unname(path_names)
  )

  module_func_annot <- merge(module_func_annot, path_df, by = "pathway")

  module_func_annot$level2 <- case_when(
    grepl("Glycolysis|Gluconeogenesis|Pentose|Starch|Sugar",
          module_func_annot$name, ignore.case = TRUE) ~ "Carbohydrate metabolism",
    grepl("Amino acid", module_func_annot$name, ignore.case = TRUE) ~ "Amino acid metabolism",
    grepl("Lipid", module_func_annot$name, ignore.case = TRUE) ~ "Lipid metabolism",
    grepl("Energy|Oxidative|Carbon fixation",
          module_func_annot$name, ignore.case = TRUE) ~ "Energy metabolism",
    grepl("Transport|Secretion|ABC",
          module_func_annot$name, ignore.case = TRUE) ~ "Membrane transport",
    grepl("Ribosome|Translation",
          module_func_annot$name, ignore.case = TRUE) ~ "Translation",
    grepl("DNA|Replication|Repair",
          module_func_annot$name, ignore.case = TRUE) ~ "Replication and repair",
    grepl("Transcription", module_func_annot$name, ignore.case = TRUE) ~ "Transcription",
    grepl("Signal", module_func_annot$name, ignore.case = TRUE) ~ "Signal transduction",
    TRUE ~ "Other"
  )

  module_L2 <- module_func_annot %>%
    group_by(module, level2) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    group_by(module) %>%
    mutate(rel_abundance = abundance / sum(abundance)) %>%
    ungroup()

  module_L3 <- module_func_annot %>%
    group_by(module, name) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    group_by(module) %>%
    mutate(rel_abundance = abundance / sum(abundance)) %>%
    ungroup() %>%
    filter(name != "Other")

  module_L3_filt <- module_L3 %>%
    group_by(name) %>%
    filter(sd(rel_abundance, na.rm = TRUE) > 0.001) %>%
    ungroup()

  top_diff <- module_L3_filt %>%
    group_by(module) %>%
    mutate(z = as.numeric(scale(rel_abundance))) %>%
    slice_max(z, n = 5, with_ties = FALSE) %>%
    ungroup() %>%
    distinct(module, name, .keep_all = TRUE)

  top_diff$group <- case_when(
    grepl("Glycolysis|Pentose|Sugar|Carbohydrate", top_diff$name) ~ "Carbohydrate metabolism",
    grepl("Amino acid|Lysine|Histidine|Cysteine", top_diff$name) ~ "Amino acid metabolism",
    grepl("Lipid|Fatty acid|Linoleic", top_diff$name) ~ "Lipid metabolism",
    grepl("Benzoate|Chloro|Degradation", top_diff$name) ~ "Xenobiotic degradation",
    grepl("Transport|Secretion|PTS", top_diff$name) ~ "Membrane transport",
    grepl("Biofilm|Quorum", top_diff$name) ~ "Microbial interaction",
    grepl("Resistance|beta|Vancomycin|CAMP", top_diff$name) ~ "Antibiotic resistance",
    grepl("DNA|Repair|Replication", top_diff$name) ~ "Genetic processes",
    grepl("Energy|Oxidative|Carbon", top_diff$name) ~ "Energy metabolism",
    TRUE ~ "Other"
  )

  module_group <- top_diff %>%
    group_by(module, group) %>%
    summarise(
      abundance = sum(rel_abundance),
      z = mean(z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(module) %>%
    mutate(
      z_scaled = ifelse(
        max(z, na.rm = TRUE) == min(z, na.rm = TRUE),
        1,
        (z - min(z, na.rm = TRUE)) /
          (max(z, na.rm = TRUE) - min(z, na.rm = TRUE))
      )
    ) %>%
    ungroup()

  # Cluster modules by their functional-category profiles.
  mat <- module_group %>%
    select(module, group, z) %>%
    pivot_wider(names_from = group, values_from = z, values_fill = 0)

  if (nrow(mat) > 1) {
    mat_matrix <- as.matrix(mat[, -1, drop = FALSE])
    rownames(mat_matrix) <- mat$module
    hc <- hclust(dist(mat_matrix))
    module_order <- hc$labels[hc$order]
    module_group$module <- factor(module_group$module, levels = module_order)
  }

  top_diff$name_short <- str_replace_all(
    top_diff$name,
    "biosynthesis|metabolism|pathway",
    ""
  )
}

# ------------------------------------------------------------------------------
# 11. Taxonomy / network / function comparisons
# ------------------------------------------------------------------------------

tax_mat_network <- NULL
func_mat <- NULL
mod_mat <- NULL
tax_dist <- func_dist <- mod_dist <- NULL
mantel_tax_func <- mantel_tax_mod <- mantel_func_mod <- NULL
proc_tax_func <- proc_tax_mod <- proc_func_mod <- NULL
coin_tax_func <- coin_tax_mod <- coin_mod_func <- NULL

if (
  RUN_GENUS_NETWORK &&
  !is.null(net_asc) &&
  file.exists(PATHWAY_ABUNDANCE) &&
  !is.null(module_func)
) {

  tax_mat_network <- as.data.frame(t(otu_table(ps_asc_genus_filt)))
  tax_mat_network <- decostand(tax_mat_network, method = "hellinger")

  func_mat <- read.table(
    PATHWAY_ABUNDANCE,
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
  func_mat <- decostand(t(func_mat), method = "hellinger")

  otu <- as.data.frame(otu_table(ps_asc_genus_filt))
  otu$taxon <- rownames(otu)

  modules_df <- modules[, c("id", "module")]
  colnames(modules_df) <- c("taxon", "module")

  otu_mod <- otu %>%
    inner_join(modules_df, by = "taxon") %>%
    group_by(module) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop")

  mod_mat <- as.data.frame(t(otu_mod[, -1, drop = FALSE]))
  colnames(mod_mat) <- otu_mod$module
  mod_mat <- decostand(mod_mat, method = "hellinger")

  common <- Reduce(
    intersect,
    list(
      rownames(tax_mat_network),
      rownames(func_mat),
      rownames(mod_mat)
    )
  )

  tax_sub <- tax_mat_network[common, , drop = FALSE]
  func_sub <- func_mat[common, , drop = FALSE]
  mod_sub <- mod_mat[common, , drop = FALSE]

  tax_dist <- vegdist(tax_sub, method = "bray")
  func_dist <- vegdist(func_sub, method = "bray")
  mod_dist <- vegdist(mod_sub, method = "bray")

  mantel_tax_func <- mantel(tax_dist, func_dist, method = "spearman", permutations = 999)
  mantel_tax_mod  <- mantel(tax_dist, mod_dist, method = "spearman", permutations = 999)
  mantel_func_mod <- mantel(func_dist, mod_dist, method = "spearman", permutations = 999)

  tax_ord  <- cmdscale(tax_dist, k = 2)
  func_ord <- cmdscale(func_dist, k = 2)
  mod_ord  <- cmdscale(mod_dist, k = 2)

  proc_tax_func <- protest(tax_ord, func_ord, permutations = 999)
  proc_tax_mod  <- protest(tax_ord, mod_ord, permutations = 999)
  proc_func_mod <- protest(func_ord, mod_ord, permutations = 999)

  # Co-inertia, retaining the original three pairwise comparisons.
  tax_pca  <- dudi.pca(tax_sub, scannf = FALSE, nf = 2)
  func_pca <- dudi.pca(func_sub, scannf = FALSE, nf = 2)
  mod_pca  <- dudi.pca(mod_sub, scannf = FALSE, nf = 2)

  coin_tax_func <- coinertia(tax_pca, func_pca, scannf = FALSE)
  coin_tax_mod  <- coinertia(tax_pca, mod_pca, scannf = FALSE)
  coin_mod_func <- coinertia(mod_pca, func_pca, scannf = FALSE)
}

# ------------------------------------------------------------------------------
# 12. Save all analysis objects for figures.R
# ------------------------------------------------------------------------------

save(
  ps, ps_norm, bacteria, metadata,
  PhylumCounts, ClassCounts, OrderCounts, FamilyCounts, GenusCounts,
  alpha_diversity, bray_dist, nmds, permanova,
  dispersion, dispersion_anova, dispersion_permutation,
  metadata_stats, phylum_anova, order_anova,
  phylum_tukey, order_tukey,
  metacoder_input, meta_obj, meta_obj_filt, meta_obj_A, meta_obj_G,
  invsimpson_anova, shannon_anova,
  invsimpson_tukey, shannon_tukey,
  se_asc, se_glo, net_asc, net_glo,
  se_asc_asv, se_glo_asv, net_asc_asv, net_glo_asv,
  daa_annotated, picrust_abundance, picrust_metadata,
  gsea_kegg, gsea_kegg_annotated, gsea_go, gsea_go_annotated,
  module_func, module_func_annot, module_L2, module_L3,
  top_diff, module_group, ko_df,
  tax_mat_network, func_mat, mod_mat,
  tax_dist, func_dist, mod_dist,
  mantel_tax_func, mantel_tax_mod, mantel_func_mod,
  proc_tax_func, proc_tax_mod, proc_func_mod,
  coin_tax_func, coin_tax_mod, coin_mod_func,
  file = "local_objects/pipeline_objects.RData"
)

message("Pipeline completed.")
message("Objects saved to local_objects/pipeline_objects.RData")
