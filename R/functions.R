
# merge w/ flag -------------------------------------------------------------------

# via mathilde
stata_merge <- function(x, y, name) {
  require(dplyr)
  x <- x %>% mutate(new1 = 1)
  y <- y %>% mutate(new2 = 2)
  df <- merge(x, y, by = name, all = TRUE)
  df <- df %>%
    mutate(new1 = ifelse(is.na(new1), 0, new1)) %>%
    mutate(new2 = ifelse(is.na(new2), 0, new2)) %>%
    mutate(merge_flag = new1 + new2) %>%
    select(-new1,-new2)
  return(df)
}


# ggplot functions for coefplots ------------------------------------------


require(ggplot2)
collidev <- function(data, height = NULL, name, strategy, check.height = TRUE) {
  # Determine height
  if (!is.null(height)) {
    # height set manually
    if (!(all(c("ymin", "ymax") %in% names(data)))) {
      data$ymin <- data$y - height / 2
      data$ymax <- data$y + height / 2
    }
  } else {
    if (!(all(c("ymin", "ymax") %in% names(data)))) {
      data$ymin <- data$y
      data$ymax <- data$y
    }

    # height determined from data, must be floating point constant
    heights <- unique(data$ymax - data$ymin)
    heights <- heights[!is.na(heights)]

    #   # Suppress warning message since it's not reliable
    #     if (!zero_range(range(heights))) {
    #       warning(name, " requires constant height: output may be incorrect",
    #         call. = FALSE)
    #     }
    height <- heights[1]
  }

  # Reorder by x position, relying on stable sort to preserve existing
  # ordering, which may be by group or order.
  data <- data[order(data$ymin), ]

  # Check for overlap
  intervals <- as.numeric(t(unique(data[c("ymin", "ymax")])))
  intervals <- intervals[!is.na(intervals)]

  if (length(unique(intervals)) > 1 & any(diff(scale(intervals)) < -1e-6)) {
    warning(name, " requires non-overlapping y intervals", call. = FALSE)
    # This is where the algorithm from [L. Wilkinson. Dot plots.
    # The American Statistician, 1999.] should be used
  }

  if (!is.null(data$xmax)) {
    plyr::ddply(data, "ymin", strategy, height = height)
  } else if (!is.null(data$x)) {
    data$xmax <- data$x
    data <- plyr::ddply(data, "ymin", strategy, height = height)
    data$x <- data$xmax
    data
  } else {
    stop("Neither x nor xmax defined")
  }
}

# Stack overlapping intervals.
# Assumes that each set has the same horizontal position
pos_stackv <- function(df, height) {
  if (nrow(df) == 1) return(df)

  n <- nrow(df) + 1
  x <- ifelse(is.na(df$x), 0, df$x)
  if (all(is.na(df$y))) {
    heights <- rep(NA, n)
  } else {
    heights <- c(0, cumsum(x))
  }

  df$xmin <- heights[-n]
  df$xmax <- heights[-1]
  df$x <- df$xmax
  df
}

# Stack overlapping intervals and set height to 1.
# Assumes that each set has the same horizontal position.
pos_fillv <- function(df, height) {
  stacked <- pos_stackv(df, height)
  stacked$xmin <- stacked$xmin / max(stacked$xmax)
  stacked$xmax <- stacked$xmax / max(stacked$xmax)
  stacked$x <- stacked$xmax
  stacked
}

# Dodge overlapping interval.
# Assumes that each set has the same horizontal position.
pos_dodgev <- function(df, height) {
  n <- length(unique(df$group))
  if (n == 1) return(df)

  if (!all(c("ymin", "ymax") %in% names(df))) {
    df$ymin <- df$y
    df$ymax <- df$y
  }

  d_height <- max(df$ymax - df$ymin)

  # df <- data.frame(n = c(2:5, 10, 26), div = c(4, 3, 2.666666,  2.5, 2.2, 2.1))
  # ggplot(df, aes(n, div)) + geom_point()

  # Have a new group index from 1 to number of groups.
  # This might be needed if the group numbers in this set don't include all of 1:n
  groupidy <- match(df$group, sort(unique(df$group)))

  # Find the center for each group, then use that to calculate xmin and xmax
  df$y <- df$y + height * ((groupidy - 0.5) / n - .5)
  df$ymin <- df$y - d_height / n / 2
  df$ymax <- df$y + d_height / n / 2

  df
}


#' Adjust position by dodging overlaps to the side. All code written by Jared Lander available from
#' https://github.com/jaredlander/coefplot/blob/master/R/position.r/
#'
#' @inheritParams ggplot2::position_identity
#' @param height Dodging height, when different to the height of the individual
#'   elements. This is useful when you want to align narrow geoms with wider
#'   geoms. See the examples for a use case.
#' @family position adjustments
#' @export
position_dodgev <- function(height = NULL) {
  ggproto(NULL, PositionDodgeV, height = height)
}



PositionDodgeV <- ggproto(`_class` = "PositionDodgeV", `_inherit` = ggplot2::Position,
                          required_aes = "y",
                          height = NULL,
                          setup_params = function(self, data) {
                            if (is.null(data$ymin) && is.null(data$ymax) && is.null(self$height)) {
                              warning("height not defined. Set with `position_dodgev(height = ?)`",
                                      call. = FALSE)
                            }
                            list(height = self$height)
                          },

                          compute_panel = function(data, params, scales) {
                            collidev(data, params$height, "position_dodgev", pos_dodgev, check.height = FALSE)
                          }
)




# extract legend as grob from ggplot --------------------------------------

#' Extract legend from ggplot
#'
#' This function takes a ggplot2 object and extracts the legend as a grob.
#' @param theplot A ggplot2 object.
#' @export
#' @example
#' library(ggplot2)
#' library(grid)
#' data(mtcars)
#' gp = ggplot(mtcars, aes(x = wt, y = mpg, colour = factor(am))) + geom_point()
#' leg = g_legend(gp)
#' grid.draw(leg)
g_legend<-function(theplot){
  tmp <- ggplot_gtable(ggplot_build(theplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}




# get terciles ------------------------------------------------------------

#' Create terciles
#'
#' Create terciles from a numeric variable
#' @param x A numeric variable.
#' @export
#' @example
#' x = runif(10)
#' tercileAssign(x)

tercileAssign = function(x){
  stopifnot(class(x) %in% c("numeric", "integer"))
  x = Hmisc::cut2(x, g = 3, levels.mean = T) %>% as.character %>% as.numeric
  ind.extrema = quantile(x, c(0, 1), na.rm=T)
  x = ifelse(x == ind.extrema[1], 1,
             ifelse(x == ind.extrema[2], 3,
                    ifelse(is.na(x), NA, 2)))
  return(x)
}




# vcovCluster -------------------------------------------------------------

# compute clustered standard errors
# updated to handle missing data automatically
#' Cluster-robust standard errors
#'
#' Computes cluster-robust standard errors, with automatic missing data handling.
#' @param model An estimated regression model from lm()
#' @param cluster A variable that indicates which cluster each observation belongs to.
#' @export
vcovCluster <- function(
  model,
  cluster
)
{
  require(sandwich)
  require(lmtest)

  cluster = as.factor(cluster)

  if (!is.null(model$na.action)){
    omit.rows = model$na.action
    cluster = cluster[-omit.rows]
  }
  if (nrow(model.matrix(model)) != length(cluster)){
    stop("something's not working: cluster variable has different N than model")
  }
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- model$rank

  if(M<50){
    warning("Fewer than 50 clusters, variances may be unreliable (could try block bootstrap instead).")
  }

  dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
  uj  <- na.omit(apply(estfun(model), 2, function(x) tapply(x, cluster, sum)))
  rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
  return(rcse.cov)
}



# check_then_install ------------------------------------------------------


#' Modified install package
#'
#' Checks whether the list of packages is installed. If not, then it installs it.
#' @param pkg A packages to install
#' @param ... arguments passed to install.packages()
#' @export
check_then_install = Vectorize(function(pkg, ...){
  if(pkg %in% rownames(installed.packages())){
    return(paste0(pkg, ' is already installed!'))
  } else {
    install.packages(pkg, ...)
  }
}, "pkg")




# dens_at_grid ------------------------------------------------------------


#' Evaluate density at user-specified grid
#'
#' Estimates the density at a grid of points specified by the user
#' @param point grid of points at which to evaluate the density of data.
#' @param data data with which to estimate the density.
#' @param ... arguments passed to density()

dens_at_grid = Vectorize(function(point, data, ...){
  stopifnot(length(point) == 1)
  density(data, from = point, to = point, n = 1, ...)$y
}, "point")

