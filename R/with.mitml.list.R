with.mitml.list <- function(data, expr, ...){
# evaluates an expression for a list of data sets

  expr <- substitute(expr)
  parent <- parent.frame()

  lapply(data, function(x) eval(expr, x, parent))

}
