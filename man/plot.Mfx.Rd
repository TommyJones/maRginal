% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mfx_functions.R
\name{plot.Mfx}
\alias{plot.Mfx}
\title{plot method for objects of class \code{Mfx}}
\usage{
\method{plot}{Mfx}(mfx, type = c("response", "derivative"),
  centered = FALSE, ...)
}
\arguments{
\item{mfx}{an object of class \code{Mfx}}

\item{type}{one of either "response" or "derivative" depending on whether you
want to plot ICE curves of the response or derivative function}

\item{centered}{if \code{type = "response"} logical indicating whether or not
to center the ICE curves to be 0 at their inital value}

\item{...}{other parameters to be passed to \code{\link[base]{plot}}}
}
\description{
plot method for objects of class \code{Mfx}
}
\details{
marginal effects are plotted using individual conditional expectation
(ICE) curves in the flavor of \code{\link[ICEbox]{ice}}. The (default) gray
curves are ICE curves. The green dashed curve is the mean curve. The red solid
line is the linear projection given by the marginal effect.

If \code{centered = TRUE}, the curves are centered at their initial values
where \code{x = min(x)}.

If \code{type = "derivative"}, the derivatives of ICE curves are plotted. There
is no \code{centered = TRUE} analogue.

Points are also plotted. If you want a derivative or centered ICE curve, the
actual data points are adjusted and plotted on the lines as well.
}
\examples{
# TO DO
}
