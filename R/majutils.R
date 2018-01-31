


#' Utility function that provides a sorted vector of the
#' values that were not numeric in the vector passed to
#' not_numeric.
#'
#' @param x vector of values
#' @keywords 
#' @export
#' @examples
#' not_numeric()
not_numeric <- function(x, do.vec = F){
  
  df.1 <- data_frame(x = as.character(x))
  df.2 <- suppressWarnings(df.1  %>%
    dplyr::filter(is.na(as.numeric(x))))
  
  if(do.vec){
    return(x = sort(unique(df.2$x)))
  }else{
    return(as_data_frame(x = sort(unique(df.2$x))))
  }
  
  
}


#' Sets up ggplot themes and some other constructor like stuff.
#'
#' @keywords 
#' @export
#' @examples
#' init()
init <- function(){
  
  
  ggplot2::theme_set(theme_bw())
  ggplot2::theme_update(legend.position="bottom")
  ggplot2::theme_update(legend.title=element_blank())
  # See http://ggplot2.tidyverse.org/reference/theme.html
  ggplot2::theme_update(text=element_text(size=12,  family="sans"))
  ggplot2::theme_update(axis.text.x=element_text(size=10,  family="sans"))
  ggplot2::theme_update(axis.text.y=element_text(size=10,  family="sans"))
  f.sep <- .Platform$file.sep
  
}









#' Just reference code.
#' Produces a list of ggplots with different font styles.
#'
#' @keywords emma
#' @export
#' @examples
#' gg_testplot01()
gg_testplot01 <- function(){
  
  
  windowsFonts()
  # $serif
  # [1] "TT Times New Roman"
  # 
  # $sans
  # [1] "TT Arial"
  # 
  # $mono
  # [1] "TT Courier New"
  
  
  library(ggplot2)
  p1 <- ggplot(mtcars, aes(x=wt, y=mpg)) + geom_point() +
    ggtitle("Fuel Efficiency of 32 Cars") +
    xlab("Weight (x1000 lb)") + ylab("Miles per Gallon") +
    theme(text=element_text(size=16,  family="sans"))
  
  p2 <- ggplot(mtcars, aes(x=wt, y=mpg)) + geom_point() +
    ggtitle("Fuel Efficiency of 32 Cars") +
    xlab("Weight (x1000 lb)") + ylab("Miles per Gallon") +
    theme(text=element_text(size=16,  family="serif"))

  
  return(list(p1.sans = p1, p2.serif = p2))
  
}









#' Duplicates a dataframe n times and appends the duplicates
#' to the original one.
#'
#' @param df data frame to duplicate
#' @param times number of copies required
#' @keywords emma
#' @export
#' @examples
#' duplicate_df()
duplicate_df <- function(df, times = 2){
  
  # n <- 10
  # df.new1$subject <- rep(100:109, each = 3)
  
  df.new <- do.call("rbind", replicate(times, df, simplify = FALSE))
  
  return(df.new)
}









#' Produces a basic latex table for glm coefs.
#'
#' @param lm1 glm
#' @param dp dec places to use
#' @param thefilename output filename
#' @keywords emma
#' @export
#' @examples
#' glm_table()
glm_table <- function(lm1, dp = 2,
                      inverselink = exp,
                      outname = 'tmp.tex'){
  
  s1 <- summary(lm1)
  coefs <- s1$coefficients
  
  if (!is.null(inverselink)){
    coefs[,1] <- inverselink(coefs[,1])
  }
  
  ci <- confint(lm1)
  
  if (!is.null(inverselink)){
    ci[,1] <- inverselink(ci[,1])
    ci[,2] <- inverselink(ci[,2])
  }
 # ci <- round(confint(lm1), dp)
  
  tbl <- data_frame(names = rownames(coefs),
                    est = sprintf(paste0("%.", dp, "f"), coefs[,1]), 
                    ci = sprintf(paste0("%.", dp, "f", " to ", "%.", dp, "f"),   ci[,1],   ci[,2]),
                    p = replace.p.values(coefs[,4])) %>%
    sapply(., as.character) %>%
    as_data_frame(.) %>%
    xtable::xtable(.)
  
  sink(outname)
  print(tbl, only.contents=TRUE, include.rownames=F, 
        include.colnames=F, floating=F,
        hline.after=NULL, sanitize.text.function=identity)
  sink()
  
}



#' Converts the p-values to significance stars.
#'
#' @param x vector of p-values
#' @keywords emma
#' @export
#' @examples
#' star.wars()
star.wars <- function(x){
  out <- ifelse(is.na(x), "", ifelse(x < 0.1, ifelse(x < 0.05, ifelse(x < 0.01, "***", "**"), '*'), ""))
  out
}

#' Converts the p-values to the common rounded representations (< 0.05 etc).
#'
#' @param x vector of p-values
#' @param dp dec places to use
#' @keywords emma
#' @export
#' @examples
#' replace.p.values()
replace.p.values <- function(x, dp = 3){
  out <- ifelse(is.na(x), "", ifelse(x < 0.05, ifelse(x < 0.01, "x < 0.01", "x < 0.05"), paste0(round(x, dp))))
  out
}

#' Proportion of entries in a vector equal to 
#' a specified value. Provided as a fraction or
#' percentage.
#'
#' @param x vector of values
#' @param level entry of interest
#' @param dp decimal places
#' @param percent provide in percentage terms
#' @keywords fred
#' @export
#' @examples
#' prop()
prop <- function(x, level, dp = 1, percent = T){
  
  x2 <- as.character(x)
  lvl2 <- as.character(level)
  
  myfreq <- length(x2[x2 == lvl2])
  
  myprop <- myfreq / len(x2)
  
  if (percent){
    myprop <- round(myprop*100, dp)
  } else {
    myprop <- round(myprop, dp)
  }
  
  myprop
}






#' List to data frame (actually tibble). Uses length of 
#' first vector to define number of rows. Uses the names
#' of the vectors as the column names of the data.frame.
#'
#' @param l0 list of vectors all of equal length
#' @keywords cats
#' @export
#' @examples
#' list_to_df()
list_to_df <- function(l0){
  
  df1 <- data.frame(matrix(unlist(l0), nrow=length(l0[[1]]), byrow=F), stringsAsFactors=FALSE)
  names(df1) <- names(l0)
  
  return(as_data_frame(df1))
}



#' Same as the ttesti in stata (for unequal var and with welch approx)
#'
#' looks at x2 - x1
#' @param x1 Vector of 1st sample means
#' @param x2 Vector of 2nd sample means
#' @param s1 Vector of 1st sample sd
#' @param s2 Vector of 2nd sample sd
#' @param n1 Vector of 1st sample size
#' @param n2 Vector of 2nd sample size
#' @keywords cats
#' @export
#' @examples
#' diff_means()
diff_means <- function(x1, x2, s1, s2, n1, n2){
  
  x.diff <- x2 - x1
  
  s.pooled <- pooled_sd(s1, s2, n1, n2)
  
  v <- dof_approx(s1, s2, n1, n2)
  
  
  t.obs <- x.diff / s.pooled
  
  ciMult <- qt(0.95/2 + .5, v) 
  
  ci <- s.pooled * ciMult
  
  lwr <- x.diff - ci
  upr <- x.diff + ci
  
  reslist = list(x.diff = x.diff, 
                 s.pooled = s.pooled,
                 t.obs = t.obs, 
                 t.crit = ciMult, 
                 lwr = lwr,
                 upr = upr)
  
  return(reslist)
  
}



#' Gets pooled standard deviation
#'
#' @param s1 Vector of 1st sample sd
#' @param s2 Vector of 2nd sample sd
#' @param n1 Vector of 1st sample size
#' @param n2 Vector of 2nd sample size
#' @keywords cats
#' @export
#' @examples
#' pooled_sd()
pooled_sd <- function(s1, s2, n1, n2){
  
  s.pooled <- sqrt(((s1^2)/ n1)     +    ((s2^2)/ n2))
  s.pooled
}


#' Welch satterthwaite - degrees of freedom approximation
#'
#' @param s1 Vector of 1st sample sd
#' @param s2 Vector of 2nd sample sd
#' @param n1 Vector of 1st sample size
#' @param n2 Vector of 2nd sample size
#' @keywords cats
#' @export
#' @examples
#' dof_approx()
dof_approx <- function(s1, s2, n1, n2){
  
  var1 <- s1^2
  var2 <- s2^2
  
  n1.2 <- n1^2
  n2.2 <- n2^2
  
  var1.2 <- s1^4
  var2.2 <- s2^4
  
  v1 <- n1 - 1
  v2 <- n2 - 1
  
  top <-    ( (var1/n1) + (var2/n2) )^2
  bottom <- ( (var1.2/(n1.2 * v1))   +   (var2.2/(n2.2 * v2))  )
  
  v <- top / bottom  
  
  return(v)
  
}


#' Descriptive stats as a latex formatted string
#'
#' This function gives you the number of unique values in a vector.
#' @param x Vector of values.
#' @keywords cats
#' @export
#' @examples
#' desc.stat.str()
desc.stat.str <- function(x){
  
  
  n <- sum(!is.na(x))
  mean <- mean(x, na.rm = T)
  sd <- sd(x, na.rm = T)
  min <- min(x, na.rm = T)
  max <- max(x, na.rm = T)
  
  q1 <- as.numeric(quantile(x, probs = 0.25, na.rm = T))
  q2 = as.numeric(median(x, na.rm = T))
  q3 <- as.numeric(quantile(x, probs = 0.75, na.rm = T))
  
  
  res <- sprintf("%s & %.2f & %.2f & %s & %s & %s & %s & %s", n, mean, sd, min, max, q1, q2, q3)
  
  return(res)
  
}






#' Length of unique values
#'
#' This function gives you the number of unique values in a vector.
#' @param df Vector of values.
#' @keywords cats
#' @export
#' @examples
#' col.names()
col.names <- function(df){
  
  return(as.data.frame(names(df)))
  
}


#' Length of unique values
#'
#' This function gives you the number of unique values in a vector.
#' @param x Vector of values.
#' @keywords cats
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
#' @keywords cats
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
#' @keywords cats
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
#' @keywords cats
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
#' @keywords cats
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
#' @keywords cats
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
#' @keywords cats
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
#' @keywords cats
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
#' @keywords cats
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
#' @keywords cats
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
#' @keywords cats
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
#' @keywords cats
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
#' @keywords cats
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
#' @param x Vector of numeric values
#' @keywords cats
#' @export
#' @examples
#' mean_sd()
mean_sd <- function(x, dp = 2){
  
  my.stat <- paste0(round(mean(x, na.rm = T), dp), " (",
                    round(sd(x, na.rm = T), dp),")")
  
  return(my.stat)
}

#' Abbreviatin for length
#'
#' 
#' @param x Vector of values
#' @keywords cats
#' @export
#' @examples
#' len()
len <- function(x){
  return(length(x))
}


#' Produces ggplot of missingness pattern with numbers of 
#' missing at top.
#'
#' 
#' @param x Data frame
#' @keywords cats
#' @export
#' @examples
#' ggplot_missing()
ggplot_missing <- function(x){
  
  require(grid)
  require(ggplot2)
  
  # Rownumber, variable, indicator of missingness
  df.miss <- x %>% 
    is.na %>%
    reshape2::melt(.)
  
  df.sum <- df.miss %>%
    dplyr::group_by(Var2) %>%
    dplyr::summarise(n.miss = sum(value))
  n.recs <- nrow(x)
  
  myylab <- paste0("Participant (n = ", n.recs, ")")

  p <- ggplot(data = df.miss, aes(x = Var2,y = Var1)) +
    geom_raster(aes(fill = !value)) +
    scale_fill_manual(values = c("black", "white"))+
    #scale_fill_grey(name = "",
    #              labels = c("Present","Missing")) +
    theme(text = element_text(size=10)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(legend.position="none")+
    theme(axis.text.x  = element_text(angle=45, hjust=1, vjust=0.9)) + 
    labs(x = "", y = myylab)+
    theme(plot.margin = unit(c(2,1,1,1), "lines"))
  
  # Add number missing at top of plot
  for (i in 1:length(df.sum$Var2))  {
    p <- p + annotation_custom(
      grob = textGrob(label = df.sum$n.miss[i], 
                      hjust = 0, 
                      gp = gpar(cex = 0.7)),
      ymin = n.recs + 10,      # Vertical position of the textGrob
      ymax = n.recs + 10,
      xmin = df.sum$Var2[i],  
      xmax = df.sum$Var2[i])
  }  
  
  
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid.draw(gt)
  
  
  
}


#' Method to call to produce ggplot of missingness 
#' pattern and save to file.
#'
#' 
#' @param df Data frame
#' @param filename file of plot (uses Output directory as prefix)
#' @param showme display plot as well as save
#' @param width width (cm)
#' @param height height (cm)
#' @keywords cats
#' @export
#' @examples
#' plot.missing()
plot.missing <- function(df, filename, showme = T, 
                         width = 5, height = 5){
  pdf(file=paste0("../../Output/", filename),
      width=width,
      height=height)
  print(ggplot_missing(df))
  dev.off()
  
  if(showme){
    print(ggplot_missing(df))
  }
}







#' returns string w/o leading whitespace
#'
#' 
#' @param x vector of character strings
#' @keywords cats
#' @export
#' @examples
#' trim.leading()
trim.leading <- function (x)  sub("^\\s+", "", x)


#' returns string w/o trailing whitespace
#'
#' 
#' @param x vector of character strings
#' @keywords cats
#' @export
#' @examples
#' trim.trailing()
trim.trailing <- function (x) sub("\\s+$", "", x)


#' returns string w/o leading or trailing whitespace
#'
#' 
#' @param x vector of character strings
#' @keywords cats
#' @export
#' @examples
#' trim()
trim <- function (x) gsub("^\\s+|\\s+$", "", x)




#' returns string up to the first instance of whitespace
#'
#' 
#' @param x vector of character strings
#' @keywords cats
#' @export
#' @examples
#' extract.up.to.whitespace()
extract.up.to.whitespace <- function(x) gsub( " .*$", "", x )



