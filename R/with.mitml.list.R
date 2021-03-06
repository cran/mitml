with.mitml.list <- function(data, expr, include.data = FALSE, ...){
# evaluates an expression for a list of data sets

  expr <- substitute(expr)
  pf <- parent.frame()

  # check include.data argument
  if(is.character(include.data)){
    name.data <- include.data
    include.data <- TRUE
  }else{
    name.data <- "data"
  }

  out <- if(include.data){

    lapply(data, function(d, expr){
      expr[[name.data]] <- substitute(d)
      eval(expr, parent.frame())
    }, expr = expr)

  }else{

    lapply(data, function(d, expr, pf) eval(expr, d, pf), expr = expr, pf = pf)

  }

  class(out) <- c("mitml.result", "list")
  return(out)

}
