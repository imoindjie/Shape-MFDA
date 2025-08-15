
### Takes a list of shape objects, obtained with the last function: shapes= list(shape_1, shape_2, ...)
## The output is an object "shapes", which has also some associated functions 
as.shapes=function(shapes){ 
  ti=attr(shapes[[1]], 'time')
  for(i in 1:length(shapes)){ 
    attr(shapes[[i]], 'time')=NULL
    attr(shapes[[i]], 'class')=NULL
  }
  class(shapes)='shapes'
  attr(shapes, 'time')=ti
  shapes
}


as.mvshapes=function(Xs){
  class(Xs)='mvshapes'
  attr(Xs, 'time')= attr(Xs[[1]], 'time')
  for(k in 1:length(Xs))attr(Xs[[k]], 'time')=NULL
  Xs
}
