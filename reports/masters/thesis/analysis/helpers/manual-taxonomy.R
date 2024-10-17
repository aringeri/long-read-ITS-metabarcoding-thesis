
accepted_synonyms <-
  tibble(alt_family='Pucciniaceae', alt_genus='Puccinia', alt_species='Puccinia psidii', actual_species='Austropuccinia psidii') %>%
  rows_append(
    tibble(alt_family='Debaryomycetaceae', alt_genus='Kurtzmaniella', alt_species=NA, actual_species='Candida boleticola')) %>%
rows_append(
    tibble(alt_family='Debaryomycetaceae', alt_genus='Kurtzmaniella', alt_species=NA, actual_species='Candida zeylanoides'))