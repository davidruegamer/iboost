% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iboost.R
\name{format_iboost_res}
\alias{format_iboost_res}
\alias{format_iboost}
\title{Function to calculate selective inference based on sampled values}
\usage{
format_iboost_res(res, alpha, method, vT, this_var, computeCI, fac = 2, y)
}
\arguments{
\item{res}{list of sampled values either directly obtained in the \code{\link{iboost}} 
call or when using \code{\link{iboost}} with \code{returnSamples = TRUE}.}

\item{alpha}{numeric; level for confidence interval(s)}

\item{method}{character; method used for sampling}

\item{vT}{list of test vector(s)}

\item{this_var}{numeric; variance to be used for inference}

\item{computeCI}{logical; whether or not to compute confidence intervals}

\item{fac}{numeric; used in inverse search for confidence interval limits}

\item{y}{vector; response vector}
}
\description{
This function calculates p-values and confidence intervals 
based on sampled values and weights obtained by the \code{\link{iboost}} function.
}
\details{
The specific use of this function is to (internally) compute inference based 
on samples in the \code{\link{iboost}} function, to recalculate inference for 
given samples \code{dist} from an \code{iboost} object or to calculate inference in 
the first place when \code{\link{iboost}} has been ran with argument \code{returnSamples = TRUE}. 
In the first and third case, use \code{format_iboost_res}, otherwise \code{format_iboost}.
}
\examples{
if(require("mboost")){

set.seed(0)

n <- 200
x1 <- rnorm(n)
x2 <- rnorm(n) + 0.25 * x1
x3 <- rnorm(n)
eta <- 3 * sin(x1) + x2^2
y <- scale(eta + rnorm(n), scale = FALSE)

spline1 <- bbs(x1, knots = 20, df = 4)
knots.x2 <- quantile(x2, c(0.25, 0.5, 0.75))
spline2 <- bbs(x2, knots = knots.x2, df = 4)
spline3 <- bbs(x3, knots = 20, df = 4)

data <- data.frame(y=y, x1=x1, x2=x2, x3=x3)

mod1 <- mboost(y ~ spline1 + spline2 + spline3,
control=boost_control(mstop = 73), offset = 0, 
data = data)

# calculate p-values and intervals for model with 
# fixed stopping iteration:
# this is done with only B = 100 samples for
# demonstrative purposes and should be increased 
# for actual research questions
res <- iboost(mod1, method = "impsamp", B = 100)

# recalculate inference for different variance or alpha level
format_iboost(res, alpha = 0.1, this_var = var(y)*5)
}
}
