#!/usr/bin/env Rscript
# site/reference.R: man/*.Rd -> site/reference/<topic>.html fragments and site/reference/index.json.
#
# Run from the package root (R is not run on the author's Mac; the package is rsynced to
# whitebox:Code/emphasis-site/ and the script runs there):
#
#   cd ~/Code/emphasis-site && LC_ALL=C.UTF-8 R_LIBS=~/rlib-emphasis-edfix:~/rlib-emphasis-final \
#     Rscript site/reference.R
#
# One fragment per Rd file, named after the file stem (estimate_rates.Rd -> estimate_rates.html,
# dot-attach_w.Rd -> dot-attach_w.html).  A fragment is the body of tools::Rd2HTML's output with
# the document skeleton removed: no <html>, <head>, <body>, no "topic {emphasis} / R Documentation"
# header table, no <h2> title, no package footer.  It starts at the first <h3> and holds the
# sections Description, Usage, Arguments, Details, Value, Examples, in that order, and every
# custom \section{} of the page (simulate_tree's "Model specification", for example) after
# Details and before Examples.  See Also, Note, References and Author(s) are dropped.  The title
# is in index.json, so the builder writes the heading and the anchor; cross-references between
# topics are written as href="#<topic>" (the file stem), and the builder gives each topic's
# heading id="<topic>".  Inline and displayed equations are \( ... \) and \[ ... \] inside
# <span class="reqn">, for MathJax.
#
# index.json: an array, one object per topic, in the order of the fragments on the page
# (exported functions first, alphabetically; then the S3 print methods; then the package page;
# then the internal topics): topic (file stem), title (the one-line \title), aliases, exported
# (true for a name in NAMESPACE export()), kind ("export", "S3method", "package", "internal"),
# sections (the <h3> headings kept), file (the fragment's path relative to the package root).
# A link to another package's topic (\link[stats]{optim}) points at the R manual on
# stat.ethz.ch.

suppressWarnings(suppressMessages({ library(tools) }))

root <- normalizePath(".")
man_dir <- file.path(root, "man")
out_dir <- file.path(root, "site", "reference")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- NAMESPACE: exports and S3 methods -----------------------------------------------------
ns <- readLines(file.path(root, "NAMESPACE"))
exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
exports <- gsub("[\"`]", "", exports)
s3 <- sub("^S3method\\(([^,]+),([^)]+)\\)$", "\\1.\\2", grep("^S3method\\(", ns, value = TRUE))
s3 <- gsub("[\" ]", "", s3)

# ---- Rd files -----------------------------------------------------------------------------
rd_files <- sort(list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE))
rd_tag <- function(x) { t <- attr(x, "Rd_tag"); if (is.null(t)) "" else t }
rd_text <- function(x) paste(trimws(gsub("\\s+", " ", paste(unlist(x), collapse = ""))), collapse = "")

parsed <- lapply(rd_files, function(f) {
  rd <- tools::parse_Rd(f, encoding = "UTF-8")
  tags <- vapply(rd, rd_tag, "")
  list(file = f,
       topic = sub("\\.Rd$", "", basename(f)),
       rd = rd,
       title = rd_text(rd[[which(tags == "\\title")[1L]]]),
       aliases = vapply(rd[tags == "\\alias"], rd_text, ""))
})
names(parsed) <- vapply(parsed, `[[`, "", "topic")

# Every alias resolves to the anchor of its topic on the one reference page.
# Rd2HTML's resolver expects the form ../../<pkg>/html/<topic>.html; the hrefs are rewritten to
# #<topic> once the fragment is written.
link_of <- function(topic) paste0("../../emphasis/html/", topic, ".html")
links <- unlist(unname(lapply(parsed, function(p) setNames(rep(link_of(p$topic), length(p$aliases)), p$aliases))))
links <- c(links, setNames(link_of(names(parsed)), names(parsed)))
links <- links[!duplicated(names(links))]

# ---- Rd -> fragment -----------------------------------------------------------------------
standard_keep <- c("Description", "Usage", "Arguments", "Details", "Value", "Examples")
standard_drop <- c("See Also", "Note", "References", "Author(s)", "Source", "Format")

to_fragment <- function(p) {
  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp))
  tools::Rd2HTML(p$rd, out = tmp, package = "emphasis", Links = links, no_links = FALSE,
                 standalone = FALSE, texmath = "katex", outputEncoding = "UTF-8")
  h <- readLines(tmp, encoding = "UTF-8", warn = FALSE)
  html <- paste(h, collapse = "\n")
  # Everything up to and including the <h2> title is skeleton (the header table sits above it).
  html <- sub("^.*?</h2>\\s*", "", html, perl = TRUE)
  # The package footer, if Rd2HTML wrote one.
  html <- sub("<hr\\s*/?>\\s*<div[^>]*>\\s*\\[Package.*$", "", html, perl = TRUE)
  html <- sub("<hr\\s*/?>\\s*$", "", html, perl = TRUE)
  # Split at the <h3> section headings.
  pos <- gregexpr("<h3[^>]*>", html, perl = TRUE)[[1L]]
  if (pos[1L] < 0) return(list(html = "", sections = character(0)))
  starts <- as.integer(pos)
  ends <- c(starts[-1L] - 1L, nchar(html))
  chunks <- substring(html, starts, ends)
  heads <- sub("(?s)^<h3[^>]*>(.*?)</h3>.*$", "\\1", chunks, perl = TRUE)
  heads <- gsub("<[^>]+>", "", heads)
  heads <- trimws(heads)
  keep <- !(heads %in% standard_drop)
  chunks <- chunks[keep]; heads <- heads[keep]
  # Order: the standard sections in their fixed order, custom sections between Details and Examples.
  rank <- match(heads, standard_keep)
  custom <- is.na(rank)
  rank[custom] <- match("Details", standard_keep) + 0.5
  ord <- order(rank, seq_along(heads))
  chunks <- chunks[ord]; heads <- heads[ord]
  html <- paste(trimws(chunks), collapse = "\n\n")
  # Rd2HTML writes a within-page link as href="#topic+<alias>"; the page's anchors are the file
  # stems, so each alias is mapped to its topic.
  m <- gregexpr("href=\"#topic\\+[^\"]+\"", html, perl = TRUE)
  hrefs <- regmatches(html, m)[[1L]]
  if (length(hrefs)) {
    alias <- utils::URLdecode(sub("^href=\"#topic\\+([^\"]+)\"$", "\\1", hrefs, perl = TRUE))
    topic <- sub("^\\.\\./\\.\\./emphasis/html/(.*)\\.html$", "\\1", links[alias])
    topic[is.na(topic)] <- alias[is.na(topic)]
    regmatches(html, m)[[1L]] <- paste0("href=\"#", topic, "\"")
  }
  # A link into another package (\link[stats]{optim}) comes out as href="stats.html#topic+optim";
  # it is sent to the R manual.
  html <- gsub("href=\"([A-Za-z0-9.]+)\\.html#topic\\+([^\"]+)\"",
               "href=\"https://stat.ethz.ch/R-manual/R-patched/library/\\1/html/\\2.html\"", html, perl = TRUE)
  # Equations: Rd2HTML's katex mode leaves the TeX in <code class="reqn">, which MathJax skips
  # (code is in its skipHtmlTags).  A displayed equation (\deqn, written inside a centred <p>)
  # becomes \[ ... \] and an inline one (\eqn) \( ... \), both in <span class="reqn">.
  html <- gsub("(?s)(<p style=\"text-align: center;\">\\s*)<code class=\"reqn\">(.*?)</code>",
               "\\1<span class=\"reqn\">\\\\[\\2\\\\]</span>", html, perl = TRUE)
  html <- gsub("(?s)<code class=\"reqn\">(.*?)</code>", "<span class=\"reqn\">\\\\(\\1\\\\)</span>", html, perl = TRUE)
  # No katex script or stylesheet tags may remain (none are expected without standalone).
  html <- gsub("<script[^>]*katex[^>]*>.*?</script>", "", html, perl = TRUE)
  html <- gsub("<link[^>]*katex[^>]*>", "", html, perl = TRUE)
  list(html = paste0(html, "\n"), sections = heads)
}

kind_of <- function(p) {
  if (p$topic == "emphasis-package") "package"
  else if (any(p$aliases %in% exports)) "export"
  else if (any(p$aliases %in% s3)) "S3method" else "internal"
}

index <- list()
for (p in parsed) {
  fr <- to_fragment(p)
  out <- file.path(out_dir, paste0(p$topic, ".html"))
  writeLines(fr$html, out, useBytes = TRUE)
  index[[p$topic]] <- list(topic = p$topic, title = p$title, aliases = I(p$aliases),
                           exported = kind_of(p) == "export", kind = kind_of(p),
                           sections = I(fr$sections),
                           file = file.path("site", "reference", paste0(p$topic, ".html")))
}

# Order: exports (alphabetical), then S3 methods, then internal topics; the package page first
# among the internal ones.
kinds <- vapply(index, `[[`, "", "kind")
topics <- names(index)
ord <- order(match(kinds, c("export", "S3method", "package", "internal")), topics)
index <- unname(index[ord])

# JSON by hand: no dependency beyond base R.
js <- function(x) paste0("\"", gsub("\"", "\\\\\"", gsub("\\\\", "\\\\\\\\", x)), "\"")
js_arr <- function(v) paste0("[", paste(js(v), collapse = ", "), "]")
rows <- vapply(index, function(e) paste0(
  "  {\"topic\": ", js(e$topic),
  ", \"title\": ", js(e$title),
  ", \"aliases\": ", js_arr(e$aliases),
  ", \"exported\": ", if (e$exported) "true" else "false",
  ", \"kind\": ", js(e$kind),
  ", \"sections\": ", js_arr(e$sections),
  ", \"file\": ", js(e$file), "}"), "")
writeLines(c("[", paste(rows, collapse = ",\n"), "]"), file.path(out_dir, "index.json"), useBytes = TRUE)

cat(sprintf("%d fragments written to %s (%d exported, %d S3 methods, 1 package page, %d internal)\n",
            length(index), out_dir, sum(kinds == "export"), sum(kinds == "S3method"), sum(kinds == "internal")))
