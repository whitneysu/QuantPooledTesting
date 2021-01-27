#' Empirical CDF Calculation for Quantitative Pooled Testing
#'
#' Compute an empirical cumulative distribution function, with
#' results similar to an “ecdf” object.
#'
#' @param v Numeric vector of the observations for calculating ecdf. For quantitative pooling
#' strategies, only CDF on the support of v < cutoff is needed; so \code{v} can be the fraction
#' of observations less than the cutoff and use \code{N} to pass the number of all observations.
#' @param N The number (vector length) of the observations. By default, \code{N = length(v)}.
#' @param cutoff Cutoff value of quantitative assay that defines test positivity. By default,
#' \code{cutoff = max(v)} which calculates the CDF for the entire support of \code{v}.
#' @param ... Arguments to be passed to subsequent methods.
#'
#' @return The function returns a matrix of three columns: The support of (\code{v <= cutoff}),
#' empirical PMF and empirical CDF.
#'
#' @keywords Quantitative Pooled Testing
#'
#' @import stats utils
#' @references
#'
#' Liu T, Hogan JW, Daniels, MJ, Coetzer M, Xu Y, Bove G, et al. Improved HIV-1 Viral Load
#' Monitoring Capacity Using Pooled Testing with Marker-Assisted Deconvolution. Journal of
#' AIDS. 2017;75(5): 580-587.
#'
#' @examples
#'
#' ecdf_pool(round(runif(100, 0, 20)), cutoff = 18)
#'
#' @export

ecdf_pool = function (v, N = length(v), cutoff = max(v), ...){
  v0 = sort(unique(v))
  ecdf.vl0 = cbind(x = v0, pmf = as.vector(tabulate(match(v, v0)))/N)
  ecdf.vl = ecdf.vl0[v0 <= cutoff, ]
  if(length(ecdf.vl)==0)
    ecdf.vl = c(NA, 0)
  dim(ecdf.vl) = c(length(ecdf.vl)/2, 2)
  ecdf.vl = cbind(ecdf.vl, cdf = cumsum(ecdf.vl[, 2]))
  list(out = ecdf.vl, p0 = tail(ecdf.vl[, 3], 1))
}

# 2. NEW FUNCTION BEGINS HERE

#' Convert ecdf_pool output to an \code{ecdf} class
#'
#' @param x Output object from the \code{ecdf_pool()} function
#'
#' @return Output in the \code{ecdf} class
#'
#' @keywords Quantitative Pooled Testing
#'
#' @references
#'
#' Liu T, Hogan JW, Daniels, MJ, Coetzer M, Xu Y, Bove G, et al. Improved HIV-1 Viral Load
#' Monitoring Capacity Using Pooled Testing with Marker-Assisted Deconvolution. Journal of
#' AIDS. 2017;75(5): 580-587.
#'
#' @examples
#'
#' to_ecdf_class(ecdf_pool(c(0,1,2,4,4,4,4,5,0,9,87,12), cutoff = 8))
#' to_ecdf_class(ecdf_pool(round(runif(100, 0, 20)), cutoff = 18))
#'
#' @import stats utils
#' @export

to_ecdf_class = function(x){
  rval <- approxfun(x$out[, 1], x$out[, 3],
                    method = "constant", yleft = 0, yright = x$p0, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  attr(rval, "call") <- sys.call()
  rval
}

# 3. NEW FUNCTION BEGINS HERE

#' Convolution of two empirical distributions
#'
#' The function calculates the convolution of two empirical distributions up to the value of
#' \code{cutoff} (on the sum of two random variables).
#'
#' @param ecdf1 An ecdf that is the matrix output object from the function \code{convol.dens()} itself or from the function
#' \code{ecdf_pool()}. The ecdf can be truncated up by \code{cutoff}.
#' @param ecdf2 An ecdf that is the matrix output object from the function \code{convol.dens()} itself or from the function
#' \code{ecdf_pool()}. The ecdf can be truncated up by \code{cutoff}.
#' @param N1 The size of the empirical support of \code{ecdf1}.
#' @param N2 The size of the empirical support of \code{ecdf2}.
#' @param cutoff Cutoff of the support of the resulting convolution distribution
#'
#' @return The function returns a matrix of three columns: The support of convolution distribution
#' up to \code{cutoff}, empirical PMF and empirical CDF.
#'
#' @keywords Quantitative Pooled Testing, convolution, empirical distribution
#'
#' @import stats utils
#'
#' @references
#'
#' Liu T, Hogan JW, Daniels, MJ, Coetzer M, Xu Y, Bove G, et al. Improved HIV-1 Viral Load
#' Monitoring Capacity Using Pooled Testing with Marker-Assisted Deconvolution. Journal of
#' AIDS. 2017;75(5): 580-587.
#'
#' @examples
#' experiment_1 <- ecdf_pool(c(0,1,2,4,4,4,4,5,0,9,87,12), cutoff = 8)
#' example_ecdf_pool <- ecdf_pool(round(runif(100, 0, 20)), cutoff = 18)
#' convol.dens(example_ecdf_pool$out, example_ecdf_pool$out, cutoff = 18)
#' convol.dens(experiment_1$out, experiment_1$out, cutoff = 8)
#'
#' @export

convol.dens <- function(ecdf1, ecdf2, N1, N2, cutoff) {
  #print((ecdf1))
  if(is.na(ecdf1[1, 1])){
    out = c(NA, 0)
  }
  else{
    NN1 = nrow(ecdf1)
    NN2 = nrow(ecdf2)
    foo1 = kronecker(t(rep(1, NN2)), ecdf1[, 1])
    foo2 = kronecker(rep(1, NN1), t(ecdf2[, 1]))
    foo.sum = foo1 + foo2
    foo1 = kronecker(t(rep(1, NN2)), ecdf1[, 2])
    foo2 = kronecker(rep(1, NN1), t(ecdf2[, 2]))
    foo.freq = foo1 * foo2
    if(all(foo.sum > cutoff)){
      out = c(NA, 0)
    }else{
      out = aggregate( foo.freq[foo.sum <= cutoff], list(foo.sum[foo.sum <= cutoff]), sum)
    }
  }

  dim(out) = c(length(out)/2, 2)
  colnames(out) = c("v", "freq")
  out = cbind(out, cdf = cumsum(out[, 2]))
  list(out = out, p0 = tail(out[, 3], 1))
}


# 4. A New function begins here

#' Intermediate function needed by MP, MiniPooling with Algorithm (MPA), and mMPA
#'
#' This function calculates the ART needed by MPA for a range of pool sizes
#' (from 2 to a pool size specified).
#'
#' @inheritParams ecdf_pool
#' @param max_K Maximum pool size that is considered.
#' @param quietly Logical value; whether print the results to screen or not.
#'
#' @return The function returns a matrix of two columns: The pool size from 1:max_K, the convol_cdf(? question about how to describe this value), and the Average Tests Required by MPA (MPA_ATR) given that pool size
#'
#' @keywords Quantitative Pooled Testing, ATR (Average Tests Required), Intermediate function, Pool Size, Convolution CDF, MPA_ATR (?)
#'
#' @examples
#'
#' foo_atr(c(200, 1500, 1900, 800, 950), max_K = 9, cutoff = 2000, quietly = TRUE)
#' foo_atr(c(200, 1500, 1900, 800, 950), max_K = 6, cutoff = 1000, quietly = TRUE)
#'
#' @import stats utils
#'
#' @export

foo_atr = function (v, N = length(v), max_K = 5, cutoff = 1000, quietly = TRUE)
{
  t1 = ecdf_pool(v, N, cutoff)
  pp = t1$p0
  tt = t1$out

  atr = NA
  for (i in 2:max_K) {
    tt = convol.dens(tt, t1$out, N^(i - 1), N, cutoff)
    pp = c(pp, tt$p0)
    atr0 = 1 + sum(1 - pp[-1])
    if (!quietly)
      cat("With a pool size of k =", i, ", the estimated ATR using MPA is:",
          round(atr0/i * 100, 2), "assays per 100 individuals.\n")
    atr = c(atr, atr0/i*100)
    tt = tt$out
  }
  out = cbind(Pool_Size = 1:max_K, convol_cdf = pp, MPA_ATR = atr)
}


# 5. A New function begins here

#' Average Tests Required (ATR) needed by MiniPooling (MP)
#'
#' This function calculates the ART needed by MP for a range of pool sizes
#' (from 2 to a pool size specified).
#'
#' @inheritParams foo_atr
#'
#' @return The function returns a matrix of two columns: The pool size from 2:max_K and the Average Tests Required by MP given that pool size (MP_ATR)
#'
#' @keywords Quantitative Pooled Testing, ATR (Average Tests Required), Pool Size (2:K), MP_ATR (?)
#'
#' @references
#'
#' Liu T, Hogan JW, Daniels, MJ, Coetzer M, Xu Y, Bove G, et al. Improved HIV-1 Viral Load
#' Monitoring Capacity Using Pooled Testing with Marker-Assisted Deconvolution. Journal of
#' AIDS. 2017;75(5): 580-587.
#'
#' @examples
#'
#' mp_atr(c(200, 1500, 1900, 800, 950), max_K = 9, cutoff = 2000, quietly = TRUE)
#' mp_atr(c(200, 1500, 1900, 800, 950), max_K = 6, cutoff = 1000, quietly = TRUE)
#'
#' @import stats utils
#'
#' @export

mp_atr = function(v, N = length(v), max_K = 5, cutoff = 1000, quietly = TRUE){
  foo = foo_atr(v, N, max_K, cutoff, quietly=T)
  out = cbind(foo, MP_ATR = (1+(1-foo[, 2])*foo[, 1]) / foo[, 1]*100) # columnbind a column called MP_ATR to the output from foo_atr (which has three columns: Pool_Size, convol_cdf, MPA_ATR)
  out[1, 4] = NA #fill coordinate 1, 4 with 'NA'
  out[, 4] = out[, 4]/100 #divide by 100
  out[-1, c(1, 4)] #delete first row and choose the two columns (column 1 (which is pool size) and column 4 (which is MP_ATR))
}



# 6. A New function begins here


#' Average Tests Required (ATR) needed by MiniPooling with Algorithm (MPA)
#'
#' This function calculates the ART needed by MP for a range of pool sizes
#' (from 2 to a pool size specified).
#'
#' @inheritParams foo_atr
#'
#' @return The function returns a matrix of two columns: The pool size from 2:max_K and the Average Tests Required by MPA given that pool size (MPA_ATR)
#'
#' @keywords Quantitative Pooled Testing, ATR (Average Tests Required), Pool Size (2:K), MPA_ATR (?)
#'
#' @references
#'
#' Liu T, Hogan JW, Daniels, MJ, Coetzer M, Xu Y, Bove G, et al. Improved HIV-1 Viral Load
#' Monitoring Capacity Using Pooled Testing with Marker-Assisted Deconvolution. Journal of
#' AIDS. 2017;75(5): 580-587.
#'
#' @examples
#'
#' mpa_atr(c(200, 1500, 1900, 800, 950), max_K = 9, cutoff = 2000, quietly = TRUE)
#' mpa_atr(c(200, 1500, 1900, 800, 950), max_K = 6, cutoff = 1000, quietly = TRUE)
#'
#' @import stats utils
#'
#' @export

mpa_atr = function(v, N = length(v), max_K = 5, cutoff = 1000, quietly = TRUE){
  foo = foo_atr(v, N, max_K, cutoff, quietly)
  out = cbind(foo, MP_ATR = (1+(1-foo[, 2])*foo[, 1]) / foo[, 1]*100)
  out[1, 4] = NA
  out[, 3] = out[, 3]/100
  out[-1, c(1, 3)]
}

# 7. A New function begins here


#' Average Tests Required (ATR) needed by Marker-Assisted Mini-Pooling with Algorithm (mMPA)
#'
#' This function calculates the ATR needed by mMPA for a range of pool sizes
#' (from 2 to a pool size specified).
#'
#' @inheritParams mpa_atr
#' @param s Risk score that is used for ranked testing
#'
#' @return The function returns a matrix of two columns: The pool size from 2:max_K and the Average Tests Required by mMPA given that pool size (mMPA_ATR)
#'
#' @keywords Quantitative Pooled Testing, ATR (Average Tests Required), Pool Size (2:K), mMPA_ATR (?)
#'
#' @references
#'
#' Liu T, Hogan JW, Daniels, MJ, Coetzer M, Xu Y, Bove G, et al. Improved HIV-1 Viral Load
#' Monitoring Capacity Using Pooled Testing with Marker-Assisted Deconvolution. Journal of
#' AIDS. 2017;75(5): 580-587. (?)
#'
#' @examples
#'
#' mmpa_atr(c(200, 1500, 1900, 800, 950), s = c(5, 3, 1, 2, 4))
#'
#' @import stats utils
#'
#' @export

mmpa_atr = function(v, s, N = length(v), max_K = 5, cutoff = 1000, quietly = TRUE){
  nonFailure = (v <= cutoff)
  numNonFailure = sum(nonFailure)

  s = rank(s, ties.method = "random")
  s0 = s[nonFailure]
  v0 = v[nonFailure]

  foo = order(s0)
  s0 = s0[foo]
  v0 = v0[foo]

  PrTGivenS = matrix(NA, max_K-1, length(s0))
  PrTGivenS[, 1] = (c(2:max_K) * v[1]) <= cutoff

  pb <- txtProgressBar(min = 0, max = numNonFailure, style = 3)
  for(i in 2:numNonFailure){
    PrTGivenS[, i] = foo_atr( v0[1:i], N = s0[i], max_K = (max_K-1),
                              cutoff = max(cutoff - v0[i], 50), quietly = quietly )[, 2]
    setTxtProgressBar(pb, i)
  }

  out = NULL
  for(k in 2:max_K){
    cdfs = NULL
    for(j in 2:k){
      cdfs = c(cdfs, PrTGivenS[j-1, ] %*% dbeta(s0/N, j, k+1-j))
    }
    out = c(out, k - sum(cdfs/N))
  }

  cbind(Pool_Size = 2:max_K, mMPA_ATR = out/c(2:max_K))
}



