# The export list in NAMESPACE and the reference index in `_pkgdown.yml` are two
# lists that must agree, maintained separately. `_pkgdown.yml` is in
# `.Rbuildignore`, so this only runs from the source tree, not under
# `R CMD check` on a built tarball.

pkgdown_root <- function() {
  root <- testthat::test_path("..", "..")
  if (file.exists(file.path(root, "_pkgdown.yml"))) root else NA_character_
}

# Every `\alias{}` in `man/`, mapped to the topic that owns it. The index lists
# topics, not exports, so an export documented under a shared topic --
# `gap_below` under `predicates`, say -- is covered by that topic's entry.
rd_aliases <- function(root) {
  files <- list.files(file.path(root, "man"), pattern = "[.]Rd$", full.names = TRUE)
  out <- lapply(files, function(f) {
    hits <- grep("^\\\\alias\\{", readLines(f, warn = FALSE), value = TRUE)
    aliases <- gsub("^\\\\alias\\{(.*)\\}\\s*$", "\\1", hits)
    stats::setNames(
      rep(sub("[.]Rd$", "", basename(f)), length(aliases)),
      gsub("\\\\", "", aliases)
    )
  })
  unlist(out)
}

index_contents <- function(root) {
  index <- yaml::read_yaml(file.path(root, "_pkgdown.yml"))$reference
  unique(unlist(lapply(index, `[[`, "contents")))
}

test_that("every entry in the pkgdown index is a real documented topic", {
  skip_if_not_installed("yaml")
  root <- pkgdown_root()
  skip_if(is.na(root), "not running from the source tree")

  aliases <- rd_aliases(root)
  expect_setequal(setdiff(index_contents(root), names(aliases)), character())
})

test_that("every export appears somewhere in the pkgdown index", {
  skip_if_not_installed("yaml")
  root <- pkgdown_root()
  skip_if(is.na(root), "not running from the source tree")

  aliases <- rd_aliases(root)
  # A topic covers an export if the index lists the topic that owns its alias.
  covered <- names(aliases)[aliases %in% index_contents(root)]
  expect_setequal(setdiff(getNamespaceExports("ripr"), covered), character())
})

test_that("flat_subnull is internal", {
  expect_false("flat_subnull" %in% getNamespaceExports("ripr"))
})
