#' Cluster stations within a given distance.
#'
#' @param spdf SpatialPointsDataFrame of stations to cluster
#' @param dist clustering distance, i.e., cluster if closer than this
#' @param verbose whether to print the iterations (default = FALSE)
#' 
#' @return SpatialPoints object of (possibly clustered) stations
cluster_points <- function(spdf, dist, verbose = FALSE) {
  require(dplyr)
  require(geosphere)
  require(rgeos)
  require(sp)
  
  # we're going to make horrible assumptions about the structure of spdf. 
  # caveat coder.
  
  # algorithm: for each station, find nearest. if they are within min_dist of
  # each other, find their centroid. repeat with the resulting centroids and 
  # the non-clustered stations. continue until number is stable.
  
  # first get the Cartesian product
  dist.sta.p <- merge(data.frame(sta1_id = spdf$OBJECTID), 
                      data.frame(sta2_id = spdf$OBJECTID),
                      by = NULL)
  dist.sta.p <- dplyr::tbl_df(dist.sta.p)
  
  # remove rows where the station ID is the same
  dist.sta.p <- dplyr::filter(dist.sta.p, sta1_id != sta2_id)
  
  # create the distance matrix
  dist.mtx <- geosphere::distm(spdf)

  # now add distance
  dist.sta.p <- 
    dist.sta.p %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(distance_meters = dist.mtx[sta1_id, sta2_id]) %>%
    dplyr::ungroup()
  
  # now find the minimum distance, i.e., the distance to the nearest station
  dist.sta.nearest <- 
    dist.sta.p %>% 
    dplyr::group_by(sta1_id) %>% 
    dplyr::summarize(dist_to_nearest = min(distance_meters))
  
  dist.sta.nearest <- 
    dist.sta.nearest %>%
    dplyr::inner_join(dist.sta.p, 
                      by = c("sta1_id", "dist_to_nearest" = "distance_meters")) %>%
    dplyr::select(sta1_id, sta2_id, dist_to_nearest)
  
  # now calculate the centroids for those stations closer than min_dist.
  # we want the symmetric matches, where station A is closest to
  # station B and station B is closest to station A.
  sta.cluster <- 
    dist.sta.nearest %>% 
    # the left side of the join consists of those records where distance to 
    # nearest neighbor is less than the minimum, and sta1_id < sta2_id
    dplyr::filter(dist_to_nearest < dist, sta1_id < sta2_id) %>%
    dplyr::select(sta1_id, sta2_id) %>%
    dplyr::inner_join(dist.sta.nearest %>% 
                        # the match is symmetric if we also have a record for 
                        # sta2_id > sta1_id (we already checked for distance)
                        dplyr::filter(sta1_id > sta2_id) %>%
                        dplyr::select(sta2_id, sta1_id),
                      by = c("sta1_id" = "sta2_id",
                             "sta2_id" = "sta1_id"))
  if(verbose) {
    print(paste("clustering", nrow(sta.cluster) * 2, "out of", 
                nrow(spdf), "points"))
  }
  
  i <- 1
  while(i < nrow(sta.cluster) + 1) {
    u <- rgeos::gCentroid(rbind(spdf[spdf$OBJECTID == sta.cluster[i,]$sta1_id,],
                                spdf[spdf$OBJECTID == sta.cluster[i,]$sta2_id,]))
    if(i == 1) {
      x <- u
    } else {
      x <- rbind(x, u)
    }
    i <- i + 1
  }
  
  # now consolidate the new centroids and non-clustered stations
  # the non-clustered station IDs
  sta.nocluster <- 
    unlist(dist.sta.nearest %>% 
             # sta1_id is not in the clustered stations
             dplyr::filter(!(sta1_id %in% c(sta.cluster$sta1_id, 
                                            sta.cluster$sta2_id))) %>%
             dplyr::select(sta1_id),
           use.names = FALSE)
  
  # if we clustered anything, append the centroids of the non-clustered 
  # stations. otherwise, just return the original points.
  new.ctrd <- sp::SpatialPoints(spdf[spdf$OBJECTID %in% sta.nocluster,], 
                                proj4string = spdf@proj4string)
  if(nrow(sta.cluster) > 0) {
    new.ctrd <- rbind(new.ctrd, x)
  }
  
  return(new.ctrd)
}

#' Exact Optimal Allocation of samples to strata.
#' 
#' Perform Exact Optimal Allocation as described in Wright (2014).
#' 
#' @param N.h vector of stratum sizes
#' @param S.h vector of stratum standard deviations
#' @param n total sample size from Neyman Allocation
#' @param eps margin of error
#' @param alpha 1 - alpha is the confidence level
#' @param alg which algorithm to use; currently only II is implemented
#' @return list whose elements are the parameters eps and alpha, 
#'   total sample size n, and stratum sample sizes n.h
#'   
#' @references \url{https://www.census.gov/srd/papers/pdf/rrs2014-07.pdf}
exact_opt_alloc <- function(N.h, S.h, n, eps, alpha, alg = 2) {
  H <- length(N.h)
  if(H != length(S.h)) stop("N_h and S_h must be of same length")
  if(missing(alg)) message("alg not specified, using algorithm II")
  
  # if any of the variances is zero, replace it by 0.5 (we assume
  # that we are estimating a proportion)
  S.h.pos <- ifelse(S.h == 0, 0.5, S.h)
  
  # start by rounding the Neyman sample size down
  n.exactopt <- ceiling(n)
  
  # tolerance for the estimated variance of the estimator
  tol <- (eps / qnorm(alpha / 2, lower.tail = F)) ^ 2
  var.p_hat.str <- Inf
  V_hat.t1 <- (N.h / sum(N.h))^2
  
  while(var.p_hat.str > tol) {
    if(alg == 2) {
      # in algorithm II, we begin by assigning 2 units to each stratum; 
      # then construct an H x (n - 2 * H) priority matrix, where H denotes 
      # the number of strata; then allocate the remaining 
      # n - 2 * H additional samples
      n.addl <- n.exactopt - 2 * H
      
      # create an initial matrix whose (i, j)th entry is N_i * S_i, i.e., 
      # the product of the size and standard deviation of the ith stratum;
      # recall that the elements of data will be recycled to fill the matrix
      # if there are too few of them
      NS <- matrix(data = N.h * S.h.pos, nr = H, nc = n.addl)
      
      # now construct the priority matrix, dividing the (i,j)th term of NS
      # by sqrt((i + 1) * (i + 2))
      x <- seq_len(ncol(NS))
      priority.t <- apply(NS, MARGIN = 1, FUN = function(r) {
        return(r / sqrt((x + 1) * (x + 2)))
      })
      priority <- t(priority.t)
      
      # create a matrix whose entries index the elements of the priority
      # matrix by column
      I <- matrix(seq_along(priority), nr = nrow(priority))
      
      # and create a logical matrix whose (i,j)th element is TRUE if the 
      # corresponding element of the priority matrix should be selected
      sel <- matrix(I %in% order(priority, decreasing = T)[seq_len(n.addl)], 
                    nr = nrow(priority))
      
      # allocate one sample to the jth stratum for each value of j in
      # the row indices of the selection matrix
      y <- table(which(sel, arr.ind = T)[,1])
      
      # recall that every stratum gets at least 2 samples
      alloc <- as.vector(rep(2, H), mode = "integer")
      alloc[as.integer(names(y))] <- y + alloc[as.integer(names(y))]
    }
    
    # compute the estimated variance
    V_hat.h <- (N.h - alloc) * S.h.pos^2 / (N.h * alloc)
    var.p_hat.str <- sum(V_hat.t1 * V_hat.h)
    n.exactopt <- n.exactopt + 1
  }
  
  return(list(eps = eps, alpha = alpha, n = sum(alloc), n.h = alloc))
}

#' Get the centroid ID for each station.
#' 
#' @param sta SpatialPointsDataFrame of stations
#' @param ctrd SpatialPointsDataFrame of centroids
#' 
#' @return data frame whose rows are the OBJECT_ID and centroid ID
#'   of each station (after clustering)
get_ctrd_id <- function(sta, ctrd) {
  require(geosphere)
  
  if(!all(c(class(sta), class(ctrd)) == "SpatialPointsDataFrame")) {
    stop("both sta and ctrd must be of class SpatialPointsDataFrame")
  }
  
  # calculate the distance matrix
  dist.mtx <- geosphere::distm(sta, ctrd)
  
  # for each row, find the minimum distance. then, the column name
  # gives the centroid ID.
  row_min <- apply(dist.mtx, 1, min)
  
  # where are the minimum values?
  idx <- which(dist.mtx == row_min, arr.ind = TRUE)
  
  return(data.frame(OBJECTID = idx[,1], CTRD_ID = idx[,2]))
}

#' Find the maximum radius of a circle whose center is the cluster 
#'   centroid for each of k clusters, such that no two circles overlap.
#' @param spdf SpatialPointsDataFrame of stations to cluster
#' @param max.dist maximum radius length
#' @param min.dist minimum radius length
#' @verbose whether to print the iterations (default = F)
#' 
#' @return matrix whose columns are the numbers of clusters k and 
#'   maximum radius possible for each k
max_radius <- function(spdf, max.dist, min.dist = 0, verbose = FALSE) {
  # take the spdf and determine, for integer values of k in 
  # [min_dist, max_dist], the maximum radius of a circle whose center is
  # the center of the cluster, such that no two circles overlap
  
  # stats:kmeans isn't going to work because we cannot constrain it
  # such that the minimum distance between cluster centers is at least some
  # constant d. so, roll your own!
  
  # the plan is to start with max_dist and see what k pops out, then decrement
  # the distance by 1 until k changes. keep the max value of k and proceed
  # until we hit min_dist, and return the vector of k and the max distance
  # where it held. as the distance decreases, k should increase.
  
  # iterate through the (integer) values of distance
  for(d in max.dist:min.dist) {
    # set difference to the length of the original spdf
    ncl <- length(spdf)
    # do an initial clustering for this d
    spdf.cl <- cluster_points(spdf, d)
    # number of old clusters - number of new clusters
    dif <- ncl - length(spdf.cl)
    # find stable k for this value of d
    while(dif > 0) {
      ncl <- length(spdf.cl)
      spdf.cl <- 
        sp::SpatialPointsDataFrame(spdf.cl, 
                                   data = data.frame(OBJECTID = seq_len(ncl)))
      spdf.cl <- cluster_points(spdf.cl, d)  
      dif <- ncl - length(spdf.cl)
    }
    
    # once we get here, we know that stable k = ncl. we will always return k
    # for max_dist
    if(d == max.dist) {
      x <- c(ncl, d)
      names(x) <- c("k", "dist")
    } else if(ncl > k_old) {
      # the number of clusters has increased, so we know that d is
      # the largest value for k clusters
      x <- rbind(x, c(ncl,d))
    }
    k_old <- ncl
    if(verbose) print(c(k_old, d))
  }
  return(x)
}

#' Neyman allocation of samples to strata.
#' 
#' @param N.h vector of stratum sizes
#' @param S.h vector of stratum standard deviations
#' @param eps margin of error
#' @param alpha 1 - alpha is the confidence level
#' @return list whose elements are the parameters eps and alpha, 
#'   total sample size n, and stratum sample sizes n.h
neyman_alloc <- function(N.h, S.h, eps, alpha) {
  if(length(N.h) != length(S.h)) stop("N.h and S.h must be of same length")
  
  z.alpha.2 <- qnorm(alpha / 2, lower.tail = FALSE)
  
  # if any of the variances is zero, replace it by 0.5 (we assume
  # that we are estimating a proportion)
  S.h.pos <- ifelse(S.h == 0, 0.5, S.h)
  
  # calculate the stratum weights under optimal Neyman allocation
  weights <- N.h * S.h.pos / sum(N.h * S.h.pos)
  
  neyman.sizes <- ((N.h * S.h.pos) ^ 2 / weights) / 
    (sum(N.h * S.h.pos ^ 2) + (sum(N.h) * eps / z.alpha.2) ^ 2)
  
  return(list(eps = eps, alpha = alpha, 
              n = sum(neyman.sizes), n.h = neyman.sizes))
}

#' Calculate the sample mean.
#' 
#' @param data data frame from which to draw samples
#' @param data.str data frame containing stratum data, e.g., stratum size
#' @param seed random seed
#' @stratify whether to use stratify or not, i.e., use SRS (default)
#' 
#' @return sample mean
sample_mean <- function(seed, data, data.str, stratify = FALSE) {
  set.seed(seed)
  
  if(stratify) {
    sampleData <- sapply(data.str$ID, FUN = function(x) {
      numUnits <- data.str$n.h[x]
      sample(data$PROPERTY_CRIMES[data$WARD == x])[1:numUnits]
    })
    
    means <- sapply(sampleData, mean)
    ybar <- sum(means * data.str$Nh) / sum(data.str$Nh)
  } else {
    sampleData <- sample(data$PROPERTY_CRIMES)[1:sum(data.str$n.h)]
    ybar <- mean(sampleData)
  }
  
  return(ybar)
}
