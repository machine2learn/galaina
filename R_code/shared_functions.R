update_random_seed <- function(random_seed, parameter = 10) {
  random_seed + parameter
}

version_at_least <- function(pkg, than) {
  as.logical((compareVersion(as.character(packageVersion(pkg)), than) >= 0))
}

print_and_append_to_log <- function(tmp_log_str_v, fileConn) {
  cat(tmp_log_str_v)
  write(paste(tmp_log_str_v, collapse = " "), file = fileConn, append = TRUE)
}

stopGeneric <- function(condition_b, fileConn, msg) {
  if (!(condition_b)) {
    print_and_append_to_log(c(msg, "\n"), fileConn)
    stop(msg)
  }
}

showData <- function(focus_v, fileConn, msg) {
  print_and_append_to_log(c(msg, "\n", sep = ""), fileConn)
  print_and_append_to_log(c(capture.output(focus_v), "\n", sep = "\n"), fileConn)
}

