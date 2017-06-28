result

getIndex <- function(path, example) {
  x <- unlist(example)
  x[] <- seq_along(x)
  
}
  
.assignElement <- function(x, i, value) {
  if (length(i) == 1) {
    x[[i]] <- value 
    x
  } else {
    x[[i[1]]] <- Recall(x[[i[1]]], i[-1], value)
    x 
  }
}

.getElement <- function(x, i) {
  if (length(i) == 1) {
    x[[i]] 
  } else {
    Recall(x[[i[1]]], i[-1])
  }
}