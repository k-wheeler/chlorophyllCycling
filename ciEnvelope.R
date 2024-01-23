##' @name ciEnvelope
##' @title ciEnvelope
##' @export 
##' @author Mike Dietze
##' @description plot a confidence/credible interval
##' @param x x-axis data
##' @param ylo y-axis lower confidence bound
##' @param yhi y-axis upper confidence bound
##' @param ... optional graphical parameters
ciEnvelope <- function(x, ylo, yhi, ...) {
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi), ylo[1])), 
          border = NA, ...)
}