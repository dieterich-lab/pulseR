.getIndex <- function(path, example) {
  x <- unlist(example)
  x[] <- seq_along(x)
  index <- .getElement(relist(x, example), path)
  if (length(index) > 1) 
    stop("Not complete path to the element (length of selection > 1)")
  index
}
  
.assignElement <- function(x, i, value) {
  if (length(i) == 1) {
    x[[unlist(i)]] <- value 
    x
  } else {
    x[[unlist(i[1])]] <- Recall(x[[unlist(i[1])]], i[-1], value)
    x 
  }
}

.getElement <- function(x, i) {
  if (length(i) == 1) {
    x[[unlist(i)]] 
  } else {
    Recall(x[[unlist(i[1])]], i[-1])
  }
}