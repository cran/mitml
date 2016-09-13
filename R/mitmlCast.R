#mitmlCast <- function() {
## cast and uncast variables in (multiply imputed) data sets
#
## purpose 1: transform variables in data frame prior to performing MI (e.g., log transform)
## purpose 2: re-transform variables back to their original scale after MI (e.g., reverse log transform)
## purpose 3: rounding/discretization of continuous data (and reverse operations for induction of noise into discrete classes)
#
## NOTE: possible features
## * save casting expression/formula as attribute specific for each variable
## * allow reverse casting for known transform-inverse-combinations or using a generic algorithm for inverse-finding (if possible)
#
#}

