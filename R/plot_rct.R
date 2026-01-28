#' Get block rectangles as wk::rct
#'
#' @param x A `block` object.
#' @param i Block column indices (0-based). If NULL, all columns.
#' @param j Block row indices (0-based). If NULL, all rows.
#' @param ... Ignored.
#'
#' @return A `wk::rct` vector of block extents.
#' @export
#' @examples
#' \dontrun{
#' b <- block("/path/to/raster.tif")
#' block_rct(b)              # all blocks
#' block_rct(b, 0, 0)        # single block
#' block_rct(b, 0:3, 0:2)    # grid of blocks
#' plot(block_rct(b))        # wk plotting
#' }
block_rct <- new_generic("block_rct", "x")

method(block_rct, block) <- function(x, i = NULL, j = NULL, ...) {
  nb <- x@nblocks

  if (is.null(i)) i <- seq.int(0L, nb[1] - 1L)
  if (is.null(j)) j <- seq.int(0L, nb[2] - 1L)

  pairs <- expand.grid(i = as.integer(i), j = as.integer(j))

  bboxes <- mapply(
    function(ii, jj) block_bbox(x, ii, jj),
    pairs$i, pairs$j,
    SIMPLIFY = FALSE
  )

  # block_bbox returns c(xmin, ymin, xmax, ymax)
  coords <- do.call(rbind, bboxes)

  wk::rct(
    xmin = coords[, 1],
    ymin = coords[, 2],
    xmax = coords[, 3],
    ymax = coords[, 4],
    crs = x@crs
  )
}


#' Get block metadata as data frame
#'
#' @param x A `block` object.
#' @param i Block column indices (0-based). If NULL, all columns.
#' @param j Block row indices (0-based). If NULL, all rows.
#' @param ... Ignored.
#'
#' @return A data.frame with columns: i, j, xoff, yoff, xsize, ysize,
#'   xmin, ymin, xmax, ymax, and geometry (wk::rct).
#' @export
block_data <- new_generic("block_data", "x")

method(block_data, block) <- function(x, i = NULL, j = NULL, ...) {
  nb <- x@nblocks

  if (is.null(i)) i <- seq.int(0L, nb[1] - 1L)
  if (is.null(j)) j <- seq.int(0L, nb[2] - 1L)

  pairs <- expand.grid(i = as.integer(i), j = as.integer(j))
  n <- nrow(pairs)

  # Pre-allocate
  idx <- matrix(0L, nrow = n, ncol = 4)
  bbox <- matrix(0, nrow = n, ncol = 4)

  for (k in seq_len(n)) {
    idx[k, ] <- block_index(x, pairs$i[k], pairs$j[k])
    bbox[k, ] <- block_bbox(x, pairs$i[k], pairs$j[k])
  }

  data.frame(
    i = pairs$i,
    j = pairs$j,
    xoff = idx[, 1],
    yoff = idx[, 2],
    xsize = idx[, 3],
    ysize = idx[, 4],
    xmin = bbox[, 1],
    ymin = bbox[, 2],
    xmax = bbox[, 3],
    ymax = bbox[, 4],
    geometry = wk::rct(
      xmin = bbox[, 1],
      ymin = bbox[, 2],
      xmax = bbox[, 3],
      ymax = bbox[, 4],
      crs = x@crs
    )
  )
}


#' Find blocks intersecting a bounding box
#'
#' @param x A `block` object.
#' @param bbox Numeric vector c(xmin, ymin, xmax, ymax).
#' @param data If TRUE, return full block_data(). If FALSE (default), return wk::rct.
#'
#' @return A `wk::rct` vector or data.frame of intersecting blocks.
#' @export
blocks_in_bbox <- new_generic("blocks_in_bbox", "x")

method(blocks_in_bbox, block) <- function(x, bbox, data = FALSE, ...) {
  if (length(bbox) != 4L) {
    stop("'bbox' must be c(xmin, ymin, xmax, ymax)", call. = FALSE)
  }

  full_bbox <- x@bbox
  bs <- x@blocksize

  res <- x@res
  nb <- x@nblocks

  # Convert bbox to block indices
  i_min <- floor((bbox[1] - full_bbox[1]) / (bs[1] * res[1]))
  i_max <- floor((bbox[3] - full_bbox[1]) / (bs[1] * res[1]))
  j_min <- floor((full_bbox[4] - bbox[4]) / (bs[2] * res[2]))
  j_max <- floor((full_bbox[4] - bbox[2]) / (bs[2] * res[2]))

  # Clamp to valid range
  i_min <- max(0L, as.integer(i_min))
  i_max <- min(nb[1] - 1L, as.integer(i_max))
  j_min <- max(0L, as.integer(j_min))
  j_max <- min(nb[2] - 1L, as.integer(j_max))

  if (i_min > i_max || j_min > j_max) {
    if (data) {
      return(block_data(x, integer(0), integer(0)))
    } else {
      return(wk::rct(crs = x@crs))
    }
  }

  if (data) {
    block_data(x, i_min:i_max, j_min:j_max)
  } else {
    block_rct(x, i_min:i_max, j_min:j_max)
  }
}


#' Find blocks intersecting an extent
#'
#' @param x A `block` object.
#' @param ext Numeric vector c(xmin, xmax, ymin, ymax).
#' @param data If TRUE, return full block_data(). If FALSE (default), return wk::rct.
#'
#' @return A `wk::rct` vector or data.frame of intersecting blocks.
#' @export
blocks_in_ext <- new_generic("blocks_in_ext", "x")

method(blocks_in_ext, block) <- function(x, ext, data = FALSE, ...) {
  if (length(ext) != 4L) {
    stop("'ext' must be c(xmin, xmax, ymin, ymax)", call. = FALSE)
  }
  # Convert ext to bbox
  bbox <- ext[c(1L, 3L, 2L, 4L)]
  blocks_in_bbox(x, bbox, data = data, ...)
}


# Now plot_block becomes a thin wrapper:

#' Plot blocks
#'
#' @param x A `block` object.
#' @param i Block column indices. If NULL, all.
#' @param j Block row indices. If NULL, all.
#' @param add Add to existing plot?
#' @param ... Passed to wk::wk_plot or graphics::rect.
#'
#' @return Invisibly, the wk::rct of plotted blocks.
#' @export
method(plot_block, block) <- function(x, i = NULL, j = NULL, add = FALSE, ...) {
  rcts <- block_rct(x, i, j)

  if (!add) {
    full_bbox <- x@bbox
    plot(NULL, xlim = full_bbox[c(1, 3)], ylim = full_bbox[c(2, 4)],
         asp = 1, xlab = "x", ylab = "y")
  }

  wk::wk_plot(rcts, add = TRUE, ...)
  invisible(rcts)
}
