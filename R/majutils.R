

#' Length of unique values
#'
#' This function gives you the number of unique values in a vector.
#' @param x Vector of values.
#' @keywords 
#' @export
#' @examples
#' len.unique()
len.unique <- function(x){
  
  return(length(unique(x)))
  
}


#' Residual deviance check for GLM
#'
#' Compares resid dev with chisq having resid.df degrees of freedom.
#' @param myglm A fitted glm model
#' @keywords 
#' @export
#' @examples
#' glm.check.res.dev()
glm.check.res.dev <- function(myglm){
  
  # residual deviance
  dev <- deviance(myglm)
  
  df <- myglm$df.residual
  
  p <- 1 - pchisq(dev, df)
  pr <- residuals(myglm, "pearson")
  
  res1 <- sprintf("ChiSq: %.2f", qchisq(0.95, myglm$df.residual))
  res2 <- sprintf("Dev  : %.2f", deviance(myglm))
  res3 <- sprintf("p    : %.5f", p)
  res4 <- sprintf("Fit  : %s", ifelse(p < 0.05, "poor", "ok"))
  
  return(paste(res1, res2, res3, res4, sep = '\n'))
  
}


#' Codes p-values into significance stars
#'
#' Translates p < 0.05 to * etc.
#' @param x Vector of p-values
#' @keywords 
#' @export
#' @examples
#' corstarsl(x)
corstarsl <- function(x){ 
  require(Hmisc) 
  x <- as.matrix(x) 
  R <- rcorr(x)$r 
  p <- rcorr(x)$P 
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1] 
  
  ## build a new matrix that includes the correlations with their apropriate stars 
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x)) 
  diag(Rnew) <- paste(diag(R), " ", sep="") 
  rownames(Rnew) <- colnames(x) 
  colnames(Rnew) <- paste(colnames(x), "", sep="") 
  
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew) 
  
  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew) 
}


#' Number of levels in a factor vector
#'
#' 
#' @param x Vector of p-values
#' @keywords 
#' @export
#' @examples
#' how_many_levels(x)
how_many_levels <- function(x){
  if (!is.factor(x)){
    "Not applicable"
  } else {
    length(levels(x))
  }
  
}


#' Used by Cbind
#' 
#' @param mydata Data
#' @param rowsneeded Number of padding rows
#' @param first Do first row
#' @keywords 
#' @export
#' @examples
#' padNA(mydata, rowsneeded, first = TRUE)
padNA <- function (mydata, rowsneeded, first = TRUE) 
{
  temp1 = colnames(mydata)
  rowsneeded = rowsneeded - nrow(mydata)
  temp2 = setNames(
    data.frame(matrix(rep(NA, length(temp1) * rowsneeded), 
                      ncol = length(temp1))), temp1)
  if (isTRUE(first)) rbind(mydata, temp2)
  else rbind(temp2, mydata)
}

#' Used by Cbind
#' 
#' @param 
#' @keywords 
#' @export
#' @examples
#' dotnames()
dotnames <- function(...) {
  vnames <- as.list(substitute(list(...)))[-1L]
  vnames <- unlist(lapply(vnames,deparse), FALSE, FALSE)
  vnames
}

#' Version of cbind to use when you have different length data.frames
#'
#' 
#' @param 
#' @keywords 
#' @export
#' @examples
#' Cbind()
Cbind <- function(..., first = TRUE) {
  Names <- dotnames(...)
  datalist <- setNames(list(...), Names)
  
  datalist <- lapply(datalist, as.data.frame)
  
  nrows <- max(sapply(datalist, function(x) 
    ifelse(is.null(dim(x)), length(x), nrow(x))))
  datalist <- lapply(seq_along(datalist), function(x) {
    z <- datalist[[x]]
    if (is.null(dim(z))) {
      z <- setNames(data.frame(z), Names[x])
    } else {
      if (is.null(colnames(z))) {
        colnames(z) <- paste(Names[x], sequence(ncol(z)), sep = "_")
      } else {
        colnames(z) <- paste(Names[x], colnames(z), sep = "_")
      }
    }
    padNA(z, rowsneeded = nrows, first = first)
  })
  as_data_frame(do.call(cbind, datalist))
}


#' Identifies the max sensitivity and specificity cutpoints
#'
#' 
#' @param 
#' @keywords 
#' @export
#' @examples
#' cutpt.max.sens.spec()
cutpt.max.sens.spec <- function(x){
  no <- which.max(x$sensitivities+x$specificities)
  res <- x$thresholds[no]
  res
}

#' Another cut point
#'
#' 
#' @param 
#' @keywords 
#' @export
#' @examples
#' optimal_lr.eta()
optimal_lr.eta=function(x){
  no=which.max(x$res$sens+x$res$spec)[1]
  result=x$res$lr.eta[no]
  result
}

#' Another cut point
#'
#' 
#' @param 
#' @keywords 
#' @export
#' @examples
#' optimal_cutpoint()
optimal_cutpoint=function(x){
  y=optimal_lr.eta(x)
  b0=unname(x$lr$coeff[1])
  b1=unname(x$lr$coeff[2])
  result=(-log(1/y-1)-b0)/b1
  result
}

#' CI for logistic regression model (for odds ratios)
#'
#' 
#' @param 
#' @keywords 
#' @export
#' @examples
#' ci.logistic()
ci.logistic <- function(lm, dp = 3){
  or <- round(exp(coef(lm)), dp)
  ci <- round(exp(confint.default(lm)), dp)
  
  df = data.frame(or = or, lwr = ci[,1], upr = ci[,2])
  
  return(df)
}

#' Replace a value with another value.
#' Useful for datasets where the null value is coded 99 etc.
#'
#' 
#' @param 
#' @keywords 
#' @export
#' @examples
#' replace.me()
replace.me <- function(x, repl = 88){
  y <- ifelse(x == repl, NA, x)
  return(as.character(y))
}




#' How many are missing in my vector x
#'
#' 
#' @param x Vector of p-values
#' @keywords 
#' @export
#' @examples
#' n.miss()
n.miss <- function(x){
  my.stat <- length(x[is.na(x)==T])
  return(my.stat)
}

#' Means and standard deviation reported as string
#'
#' 
#' @param x Vector of p-values
#' @keywords 
#' @export
#' @examples
#' mean.sd()
mean.sd <- function(x, dp = 2){
  
  my.stat <- paste0(round(mean(x, na.rm = T), dp), " (",
                    round(sd(x, na.rm = T), dp),")")
  
  return(my.stat)
}

#' Abbreviatin for length
#'
#' 
#' @param x Vector of values
#' @keywords 
#' @export
#' @examples
#' len()
len <- function(x){
  return(length(x))
}


#' Produces ggplot of missingness pattern
#'
#' 
#' @param x Data frame
#' @keywords 
#' @export
#' @examples
#' ggplot_missing()
ggplot_missing <- function(x){
  
  x %>% 
    is.na %>%
    reshape2::melt(.) %>%
    ggplot(data = .,
           aes(x = Var2,
               y = Var1)) +
    geom_raster(aes(fill = !value)) +
    scale_fill_manual(values = c("black", "white"))+
    #scale_fill_grey(name = "",
    #              labels = c("Present","Missing")) +
    theme(text = element_text(size=10)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(legend.position="none")+
    theme(axis.text.x  = element_text(angle=45, hjust=1, vjust=0.9)) + 
    labs(x = "",
         y = "Participant")
  
}


#' Produces ggplot of missingness pattern and saves to file
#'
#' 
#' @param df Data frame
#' @param filename Filename
#' @keywords 
#' @export
#' @examples
#' ggplot_missing()
plot.missing <- function(df, filename){
  pdf(file=paste0("../../Output/", filename),width=5,height=5)
  print(ggplot_missing(df))
  dev.off()
}







#' returns string w/o leading whitespace
#'
#' 
#' @param x vector of character strings
#' @keywords 
#' @export
#' @examples
#' trim.leading()
trim.leading <- function (x)  sub("^\\s+", "", x)


#' returns string w/o trailing whitespace
#'
#' 
#' @param x vector of character strings
#' @keywords 
#' @export
#' @examples
#' trim.trailing()
trim.trailing <- function (x) sub("\\s+$", "", x)


#' returns string w/o leading or trailing whitespace
#'
#' 
#' @param x vector of character strings
#' @keywords 
#' @export
#' @examples
#' trim()
trim <- function (x) gsub("^\\s+|\\s+$", "", x)




#' returns string up to the first instance of whitespace
#'
#' 
#' @param x vector of character strings
#' @keywords 
#' @export
#' @examples
#' extract.up.to.whitespace()
extract.up.to.whitespace <- function(x) gsub( " .*$", "", x )



