# Build a Jumble gene annotation (cancergenes_clinseq, allgenes, allexons) for jumble-reference.R -a <file.RDS>.
#
# biomaRt can no longer connect to Ensembl for hg38: www.ensembl.org redirects its BioMart to a dated archive host,
# and biomaRt's host and archive-list checks fail on those redirects. This script swaps Jumble's two biomaRt helpers for
# plain HTTP queries to the archive BioMart, then runs Jumble's own generate_gene_annotation() unchanged.
#
# Usage (inside the Jumble container, which has the Jumble R package, data.table and httr):
#   singularity exec docker://hydragenetics/jumble:0.5.4 \
#     Rscript make_jumble_annotation.R <genome: hg38|hg19> <output.RDS> [biomart host]
#
# cnv/jumble_annotation_hg38_ensembl116.RDS was made with jumble:0.5.4 and
#   Rscript make_jumble_annotation.R hg38 jumble_annotation_hg38_ensembl116.RDS
# (default host https://jun2026.archive.ensembl.org = Ensembl release 116). Point jumble_reference's
# annotation config at the RDS file; jumble-reference.R accepts "biomart" or the path to such a file.

args <- commandArgs(trailingOnly = TRUE)
genome <- args[1]
out_file <- args[2]
host <- if (length(args) >= 3) args[3] else if (genome == "hg19") "https://grch37.ensembl.org" else "https://jun2026.archive.ensembl.org"

library(data.table)

query_biomart <- function(host, attributes, chrom) {
  xml <- paste0(
    '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>',
    '<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" completionStamp="1">',
    '<Dataset name="hsapiens_gene_ensembl" interface="default">',
    '<Filter name="chromosome_name" value="', chrom, '"/>',
    paste0('<Attribute name="', attributes, '"/>', collapse = ""),
    "</Dataset></Query>"
  )
  for (attempt in 1:5) {
    res <- tryCatch(
      httr::POST(paste0(host, "/biomart/martservice"), body = list(query = xml), encode = "form", httr::timeout(600)),
      error = function(e) NULL
    )
    if (!is.null(res) && httr::status_code(res) == 200) {
      txt <- httr::content(res, as = "text", encoding = "UTF-8")
      # completionStamp="1" makes BioMart end a complete answer with "[success]"; anything else is truncated or an error
      if (grepl("\\[success\\]\\s*$", txt)) {
        txt <- sub("\\[success\\]\\s*$", "", txt)
        if (!nzchar(trimws(txt))) return(data.table())
        dt <- fread(text = txt, header = FALSE, sep = "\t", quote = "", col.names = attributes,
                    colClasses = list(character = intersect(attributes, c("ensembl_gene_id", "external_gene_name",
                                                                           "chromosome_name", "gene_biotype"))))
        return(dt)
      }
    }
    message("  chr ", chrom, ": attempt ", attempt, " failed, retrying")
    Sys.sleep(10 * attempt)
  }
  stop("BioMart query failed for chromosome ", chrom, " at ", host)
}

assignInNamespace("get_ensembl_mart", function(genome, mirror = NULL) host, ns = "Jumble")
assignInNamespace("fetch_ensembl_data", function(mart, genome, attributes, type = "genes") {
  chromosomes <- c(as.character(1:22), "X", "Y")
  message("Fetching ", type, " per chromosome from ", mart)
  rbindlist(lapply(chromosomes, function(chrom) query_biomart(mart, attributes, chrom)), fill = TRUE)
}, ns = "Jumble")

# Jumble's bundled cancer_genes.csv uses some outdated HGNC symbols. GRCh37 Ensembl still knows them, but current
# Ensembl (hg38) does not. For most of them the csv also lists the current name, but FAM46C, HIST1H1C and MRE11A
# would silently drop out of cancergenes_clinseq. Translate old symbols to the current ones before matching and
# drop the resulting duplicates. hugo_symbol then holds the current name.
renamed_symbols <- c(
  FAM46C = "TENT5C", H3F3B = "H3-3B", HIST1H1C = "H1-2", MLL = "KMT2A",
  MLL2 = "KMT2D", MLL3 = "KMT2C", MRE11A = "MRE11", MYCL1 = "MYCL"
)
orig_process_cancer_genes <- Jumble:::process_cancer_genes
assignInNamespace("process_cancer_genes", function(allgenes) {
  if (genome == "hg19") return(orig_process_cancer_genes(allgenes))
  cgc_path <- system.file("extdata", "cancer_genes.csv", package = "Jumble")
  cgenes <- fread(cgc_path)
  old <- cgenes$hugo_symbol %in% names(renamed_symbols)
  cgenes[old, hugo_symbol := renamed_symbols[hugo_symbol]]
  cgenes <- unique(cgenes, by = "hugo_symbol")
  tmp <- tempfile(fileext = ".csv")
  fwrite(cgenes, tmp)
  # Run Jumble's own logic on the translated table by pointing system.file() lookups at the temporary copy
  body_fn <- orig_process_cancer_genes
  environment(body_fn) <- list2env(list(system.file = function(...) tmp), parent = environment(orig_process_cancer_genes))
  body_fn(allgenes)
}, ns = "Jumble")

annot <- Jumble:::generate_gene_annotation(genome = genome)
for (n in names(annot)) message(n, ": ", nrow(annot[[n]]), " rows")
wanted <- fread(system.file("extdata", "cancer_genes.csv", package = "Jumble"))$hugo_symbol
if (genome != "hg19") wanted <- ifelse(wanted %in% names(renamed_symbols), renamed_symbols[wanted], wanted)
missing <- setdiff(wanted, annot$cancergenes_clinseq$hugo_symbol)
message("cancer genes not found in Ensembl: ", if (length(missing)) paste(missing, collapse = " ") else "none")
stopifnot(nrow(annot$allgenes) > 15000, nrow(annot$allexons) > 100000, nrow(annot$cancergenes_clinseq) > 0)
attr(annot, "source") <- paste("Ensembl BioMart", host, "queried", format(Sys.time(), "%Y-%m-%d"))
saveRDS(annot, out_file)
message("Saved ", out_file)
