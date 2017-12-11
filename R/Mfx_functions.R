#' Calculate marginal effects and ICE curves
#' @description calculate marginal effects of features on predicted outcomes
#' from an arbitrary prediction model. Also calculates the necessary components
#' to construct individual conditional expectation (ICE) curves.
#'
#' @param object a model object
#' @param X an object of class \code{data.frame} containing data from which
#' marginal effects should be calculated
#' @param pred_fun a function that accepts two arguments corresponding to
#' \code{object} and \code{X}, above. The function should return a vector of
#' predicted responses. If this argument is not passed, the generic
#' \code{predict} function will be tried.
#' @param predictors a character vector of column names of \code{X} for which
#' marginal effects are to be calculated
#' @param max_pts an integer value or \code{NA}. If an integer, then \code{max_pts}
#' number of values will be used to calculate marginal effects and ICE curves.
#' If \code{max_pts} is \code{NA}, then each predictor's unique values will be used.
#' @param min_pts a positive integer. Predictorswith fewer than \code{min_pts}
#' unique values will be evaluated only at their unique values, rather than
#' \code{max_pts} values. If you want to guarantee all points will be evaluated
#' at \code{max_pts} values, set \code{min_pts} to 0.
#' @param dydx_mean logical indicating whether to calculate the marginal effect
#' based on the mean of the derivative of ICE curves. If \code{FALSE}, then
#' the marginal effect is calculated at each actual observation then averaged.
#' @param ... arguments to be passed on to \code{textmineR}'s function
#' code{\link[textmineR]{TmParallelApply}}.
#'
#' @return
#' Returns a \code{list} of class \code{Mfx} or \code{Mfx_list} depending on
#' whether the \code{predictors} argument is of length 1 or greater. If length 1,
#' then the result is of class \code{Mfx}. If greater than 1, the result is of
#' class \code{Mfx_list}, each element of which is of class \code{Mfx}. An object
#' of class \code{Mfx} has the following slots:
#' @slot mfx the calculated marginal effect
#' @slot se the standard error of the calculated marginal effect
#' @slot conf a 95% confidence interval around the effect
#' @slot dy a vector of changes in the predicted value
#' @slot dx a vector of changes in the predictor (set by \code{max_pts})
#' @slot x the predictor vector set by \code{max_pts}
#' @slot yh0 inital prediction corresponding to the first value of \code{x}
#' @slot true_values \code{data.frame} with values for \code{dx}, \code{dy},
#' and predictions for the actual data points.
#' @slot varname the name of the predictor for which a marginal effect has been
#' calculated
#'
#' @details
#' details coming soon
#'
#' @examples
#' # examples coming soon
#'
#'
#' @export
CalcMfx <- function(object, X, pred_fun = predict, predictors = colnames(X),
                    max_pts = 100, min_pts = 0, dydx_mean = FALSE, ...){

  ### Check consistency of inputs ----

  # are min_pts and max_pts within acceptable settings?
  if (min_pts > max_pts) {
    warning("min_pts cannot be greater than max_pts.
              Setting min_pts equal to max_pts.")

    min_pts <- max_pts
  }

  if (is.infinite(max_pts)) {
    stop("max_pts must be a finite value or NA")
  }


  # Is X a data.frame?
  if (class(X) != "data.frame")
    stop("X must be of class data.frame")

  # Do I have any predictors that aren't in X ?
  if (sum(! predictors %in% colnames(X)) > 0)
    stop("predictors must be in colnames(X)")

  # Any non-numeric predictor variables?
  if (length(predictors) == 1) {
    if (class(X[[ predictors ]]) == "factor")
      stop("CalcMfx can only handle numeric predictors. If you have a categorical
           predictor, use Factor2Binary to encode a k class factor into k
           binary variables and re-fit your model.")

    if (! class(X[[ predictors ]]) %in% c("numeric", "integer"))
      stop("CalcMfx can only handle numeric predictors.")

  } else {
    classes <- sapply(X[, predictors], class)

    if (sum(classes == "factor") > 0)
      stop("CalcMfx can only handle numeric predictors. If you have a categorical
           predictor, use Factor2Binary to encode a k class factor into k
           binary variables and re-fit your model.")

    if (sum(classes %in% c("numeric", "integer")) < length(classes))
      stop("CalcMfx can only handle numeric predictors.")

  }


  ### Gets to calculating! ----

  # If you have more than one predictor, do a recursive loop for each one
  if (length(predictors) > 1) {

    result <- textmineR::TmParallelApply(predictors, function(p){
      CalcMfx(object = object, X = X, pred_fun = pred_fun, predictors = p,
              max_pts = max_pts, dydx_mean = dydx_mean)
    }, export = c("object", "X", "pred_fun", "max_pts", "dydx_mean"), ...)

    names(result) <- predictors

    class(result) <- "Mfx_list"

  } else { # if you have only one, do the calculations

    p <- predictors

    # get the sequence to iterate over
    if (is.na(max_pts) | length(unique(X[[ p ]])) <= min_pts) {
      # if you want to use all unique data points or
      # if you have fewer than min_pts unique values

      pts <- sort(unique(X[[ p ]]))

    } else {
      # all points to be evaluated at max_pts values
      pts <- seq(min(X[[ p ]]), max(X[[ p ]]), length.out = max_pts)

    }

    # add a lower bound so we can get non-infinate values for
    # delta x and delta y at the bottom
    pts <- c(pts[ 1 ] - mean(diff(pts)), pts)

    # get predictions for each point in the sequence
    yhat <- lapply(pts, function(point){

      X_new <- X

      X_new[[ p ]] <- point

      yhat <- pred_fun(object, X_new)

      yhat
    })

    yhat <- do.call(rbind, yhat)

    # store a variable used for plot method
    yh0 <- yhat[ 2 , ] # 2 because we drop our added inital value from before

    # get the derivative

    # dy/dx of curves
    dy <- apply(yhat, 2, function(y) c(NA, diff(y)))

    dx <- c(NA, diff(pts))

    dydx <- dy / dx

    # dy/dx of true values
    small_x <- sapply(X[[ p ]], function(x) max(pts[ pts < x ], na.rm = T))

    dx_true <- X[[ p ]] - small_x

    X_new <- X
    X_new[[ p ]] <- small_x

    yhat_true <- pred_fun(object, X)

    dy_true <- yhat_true - pred_fun(object, X_new)

    dydx_true <- dy_true / dx_true

    # get the mfx
    if (dydx_mean) {
      # If you want to use mean values of the curves
      mfx <- mean(dydx[ 2:length(dydx) ], na.rm = TRUE)

      stdd <- sd(dydx[ 2:length(dydx) ], na.rm = TRUE)

      se <- stdd / sqrt(length(dydx) - 1)

      # confidence interval (totally overconfident, does not account for non-linearity)
      conf <- c(mfx - 1.96 * se,
                mfx + 1.96 * se)
    } else {
      # If you only want to consider values at actual data points
      # dy/dx of true values

      mfx <- mean(dydx_true, na.rm = TRUE)

      stdd <- sd(dydx_true, na.rm = TRUE)

      se <- stdd / sqrt(length(dydx_true))

      # confidence interval (totally overconfident, does not account for non-linearity)
      conf <- c(mfx - 1.96 * se,
                mfx + 1.96 * se)
    }


    # return the result
    result <- list(mfx = mfx,
                   se = se,
                   conf = conf,
                   dy = dy[ -1 , ],
                   dx = dx[ -1 ],
                   x = pts[ -1 ],
                   yh0 = yh0,
                   true_values = data.frame(x = X[[ p ]],
                                            dx = dx_true,
                                            y = yhat_true,
                                            dy = dy_true,
                                            stringsAsFactors = FALSE),
                   varname = p)

    class(result) <- "Mfx"

  }

  return(result)

}


#' plot method for objects of class \code{Mfx}
#' @description plot method for objects of class \code{Mfx}
#' @param mfx an object of class \code{Mfx}
#' @param type one of either "response" or "derivative" depending on whether you
#' want to plot ICE curves of the response or derivative function
#' @param centered if \code{type = "response"} logical indicating whether or not
#' to center the ICE curves to be 0 at their inital value
#' @param ... other parameters to be passed to \code{\link[base]{plot}}
#' @details marginal effects are plotted using individual conditional expectation
#' (ICE) curves in the flavor of \code{\link[ICEbox]{ice}}. The (default) gray
#' curves are ICE curves. The green dashed curve is the mean curve. The red solid
#' line is the linear projection given by the marginal effect.
#'
#' If \code{centered = TRUE}, the curves are centered at their initial values
#' where \code{x = min(x)}.
#'
#' If \code{type = "derivative"}, the derivatives of ICE curves are plotted. There
#' is no \code{centered = TRUE} analogue.
#'
#' Points are also plotted. If you want a derivative or centered ICE curve, the
#' actual data points are adjusted and plotted on the lines as well.
#'
#' @examples
#' # TO DO
#' @export
plot.Mfx <- function(mfx, type = c("response", "derivative"), centered = FALSE, ...){

  ### check consistency of inputs ----
  if (! type[ 1 ] %in% c("response", "derivative"))
    stop("type must be one of 'response' or 'derivative'")

  if (class(mfx) != "Mfx")
    stop("mfx must be of class Mfx")

  ### Declare a function to get our plot options in order ----
  DoDotsProcedure <- function(...){

    dots <- list(...)

    if ( ! "ylim" %in% names(dots)) {
      dots$ylim <- range(plotmat, na.rm = TRUE)
    }

    if (! "col" %in% names(dots)) {
      dots$col <- rgb(0,0,0,0.05)
    }

    if (! "ylab" %in% names(dots)) {
      dots$ylab <- ylab
    }

    if (! "xlab" %in% names(dots)) {
      dots$xlab <- mfx$varname
    }

    if ("type" %in% names(dots)) {
      warning("reverting type to type = 'l'")
    }

    dots$type <- "l"

    dots
  }

  ### Do calculations to get the chosen plot ----
  if (type[ 1 ] == "response") {
    plotmat <- rbind(mfx$yh0, mfx$dy[ -1 , ])

    plotmat <- apply(plotmat, 2, cumsum)

    ice_mean <- rowMeans(plotmat)

    plotpoints <- data.frame(x = mfx$true_values$x,
                             y = mfx$true_values$y,
                             stringsAsFactors = FALSE)


    # get intercept for mfx line
    int <- ice_mean[ which.min(abs(mfx$x)) ] - mfx$mfx * mfx$x[ which.min(abs(mfx$x)) ]

    mfx_pred <- mfx$x * mfx$mfx + int

    ylab <- "partial y-hat"

    if (centered) {

      plotpoints$y <- plotpoints$y - plotmat[ 1 , ]

      plotmat <- t(t(plotmat) - plotmat[ 1 , ])

      mfx_pred <- mfx_pred - mfx_pred[ 1 ]

      ylab <- paste(ylab, "(centered)")

    }

  } else { # I assume you mean derivative

    if (centered)
      message("Warning, centered = TRUE has no meaning for derivative. Curve
              will not be centered.")

    plotmat <- mfx$dy / mfx$dx

    mfx_pred <- rep(mfx$mfx, length(mfx$x))

    plotpoints <- data.frame(x = mfx$true_values$x,
                             y = mfx$true_values$dy / mfx$true_values$dx,
                             stringsAsFactors = FALSE)

    ylab <- "derivative of partial y-hat"
  }

  # get our plot options straight
  dots <- DoDotsProcedure(...)

  ### do the plotting ----

  # plot the ICE curves
  dots$x <- mfx$x
  dots$y <- plotmat[ , 1 ]

  do.call(plot, dots)

  for (j in 2:ncol(plotmat)) {
    dots$y <- plotmat[ , j ]
    do.call(lines, dots)
  }

  # plot points of actual values
  dots$x <- plotpoints$x
  dots$y <- plotpoints$y
  dots$pch <- 19
  dots$type <- "p"
  dots$col <- "black"
  do.call(points, dots)

  # plot the curve predicted by the mfx
  dots$type = "l"
  dots$x <- mfx$x
  dots$y <- mfx_pred
  dots$col <- "#fc8d62"
  dots$lty = 1
  dots$lwd = 5

  do.call(lines, dots)


  # plot the average ICE curve
  dots$y <- rowMeans(plotmat)
  dots$col <- "#1b9e77"
  dots$lty = 2
  dots$lwd = 5

  do.call(lines, dots)

}

#' plot method for objects of class \code{Mfx_list}
#' @description plot method for objects of class \code{Mfx_list}
#' @param mfx_list an object of class \code{Mfx_list}
#' @param ask do you want to pause and wait for user input between plots?
#' @param ... other parameters to be passed to \code{plot.Mfx} or
#' \code{\link[base]{plot}}.
#' @details
#' \code{plot.Mfx_list} simply loops over every element of \code{mfx_list} and
#' callsj \code{plot.Mfx} for each element. It pauses between plots depending on
#' the value passed to \code{ask}.
#' @examples
#' # TO DO
#' @export
plot.Mfx_list <- function(mfx_list, ask = TRUE, ...) {
  if (class(mfx_list) != "Mfx_list") {
    check <- unique(sapply(mfx_list, class))

    if (length(check) != 1 | check[ 1 ] != "Mfx")
      stop("mfx_list must be of class Mfx_list")
  }

  # get current par setting
  opar <- par()

  # run plots, pausing between
  for (j in seq_along(mfx_list)) {
    par(ask = ask)
    plot(mfx_list[[ j ]], ...)
  }

  # return par to previous setting
  par(ask = opar$ask)

}

#' summary method for objects of class \code{Mfx_list}
#' @description summary method for objects of class \code{Mfx_list}
#' @param mfx_list an object of class \code{Mfx_list}
#' @param print do you want to print the result to the console?
#' @return a \code{data.frame} with 6 columns:
#' @slot variable the name of the variable or predictor
#' @slot effect the marginal effect of that variable on the outcome
#' @slot se the standard error of the the effect
#' @slot t_stat a t-statistic on the effect with a null hypothesis of
#' \code{effect = 0}
#' @slot 95_conf_low the lower bound of a 95% confidence interval
#' @slot 95_conf_high the upper bound of a 95% confidence interval
#' @examples
#' # TO DO
#' @export
summary.Mfx_list <- function(mfx_list, print = TRUE){
  tab <- do.call(rbind, lapply(mfx_list, function(x){
    data.frame(variable = paste0(rep(x$varname, length(x$mfx)), names(x$mfx)),
               effect = x$mfx,
               stringsAsFactors = FALSE)
  }))

  se <- sapply(mfx_list, function(x) x$se)

  conf <- do.call(rbind, lapply(mfx_list, function(x) x$conf))

  colnames(conf) <- c("95_conf_low", "95_conf_high")

  tab <- cbind(tab, se = se, t_stat = tab$effect / se, conf)

  if(print){
    print(tab)
  } else{
    tab
  }
}

#' print method for objects of class \code{Mfx_list}
#' @description print method for objects of class \code{Mfx_list}
#' @param mfx_list an object of class \code{Mfx_list}
#' @details simply calls \code{summary.Mfx_list} with option \code{print = TRUE}
#' @examples
#' # TO DO
#' @export
print.Mfx_list <- function(mfx_list){
  summary(mfx_list, print = TRUE)
}

#' easily convert a k-level factor variable to k binary variables
#' @description easily convert a k-level factor variable to k binary variables
#' @param y a factor variable with at least 2 levels
#' @return a matrix of k binary columns, one corresponding to each level of
#' the input, \code{y}
#' @examples
#' my_factor <- c("bob", "sharron", "bob", "alice", "alice")
#'
#' my_factor <- as.factor(my_factor)
#'
#' my_binary_matrix <- Factor2Binary(my_factor)
#'
#' my_binary_matrix
#' @export
Factor2Binary <- function(y){
  ynames <- levels(y)

  Y <- data.frame(y = y)

  Y <- model.matrix(~y, Y)

  colnames(Y) <- ynames

  # model.matrix always makes the first column all "1", fix this
  Y[ , 1 ] <- as.numeric(rowSums(Y[ , -1 ]) == 0)


  Y
}
