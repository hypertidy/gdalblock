#' block S7 Class
#'
#' An S7 class for block-based access to GDAL raster datasets.
#'
#' @param dsn A dataset description (file path, URL, or any GDAL-supported URI).
#'
#' @return A `block` object.
#'
#' @details
#' The block class wraps a `gdalraster::GDALRaster` object and provides
#' convenient access to the block (tile) structure of the dataset.
#'
#' ## Properties
#'
#' Access properties using `@`:
#'
#' - `dsn`: The dataset description/URI
#' - `bbox`: Bounding box as c(xmin, ymin, xmax, ymax)
#' - `dimension`: Dimensions as c(ncol, nrow, nbands)
#' - `crs`: Coordinate reference system (WKT string)
#' - `nbands`: Number of bands
#' - `datatype`: Character vector of data types per band
#' - `res`: Resolution as c(xres, yres)
#' - `blocksize`: Standard block size as c(width, height)
#' - `nblocks`: Number of blocks as c(nx, ny)
#' - `blocks`: List with block structure details
#'
#' ## Methods
#'
#' - `read_block(x, i, j, band)`: Read block data as matrix
#' - `block_dim(x, i, j)`: Get block dimensions
#' - `block_bbox(x, i, j)`: Get block bounding box
#' - `block_index(x, i, j)`: Get pixel/line offsets
#' - `block_rct(x, i, j)`: Get blocks as wk::rct
#' - `block_data(x, i, j)`: Get blocks as data.frame with full metadata
#' - `blocks_in_bbox(x, bbox)`: Find blocks intersecting a bbox
#' - `blocks_in_ext(x, ext)`: Find blocks intersecting an extent
#'
#' @export
#' @examples
#' \dontrun{
#' b <- block("/path/to/raster.tif")
#' b@dimension
#' b@blocksize
#' read_block(b, 0, 0)
#' plot(block_rct(b))
#' }

# Generics for methods - must be defined before class and method registrations

#' Read a block of raster data
#'
#' @param x A `block` object.
#' @param i Block column index (0-based).
#' @param j Block row index (0-based).
#' @param band Band number (1-based, default 1).
#'
#' @return A matrix of raster values.
#' @export
read_block <- new_generic("read_block", "x")

#' Get block dimensions
#'
#' @param x A `block` object.
#' @param i Block column index (0-based).
#' @param j Block row index (0-based).
#'
#' @return Integer vector c(ncol, nrow).
#' @export
block_dim <- new_generic("block_dim", "x")

#' Get block bounding box
#'
#' @param x A `block` object.
#' @param i Block column index (0-based).
#' @param j Block row index (0-based).
#'
#' @return Numeric vector c(xmin, ymin, xmax, ymax).
#' @export
block_bbox <- new_generic("block_bbox", "x")

#' Get block pixel/line index
#'
#' @param x A `block` object.
#' @param i Block column index (0-based).
#' @param j Block row index (0-based).
#'
#' @return Integer vector c(xoff, yoff, xsize, ysize).
#' @export
block_index <- new_generic("block_index", "x")

#' Get block rectangles as wk::rct
#'
#' @param x A `block` object.
#' @param i Block column indices (0-based). If NULL, all columns.
#' @param j Block row indices (0-based). If NULL, all rows.
#' @param ... Ignored.
#'
#' @return A `wk::rct` vector of block extents with CRS.
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
#' @examples
#' \dontrun{
#' b <- block("/path/to/raster.tif")
#' block_data(b)              # all blocks with full metadata
#' block_data(b, 0:2, 0:2)    # subset
#' }
block_data <- new_generic("block_data", "x")

#' Find blocks intersecting a bounding box
#'
#' @param x A `block` object.
#' @param bbox Numeric vector c(xmin, ymin, xmax, ymax).
#' @param data If TRUE, return full block_data(). If FALSE (default), return wk::rct.
#' @param ... Ignored.
#'
#' @return A `wk::rct` vector or data.frame of intersecting blocks.
#' @export
#' @examples
#' \dontrun{
#' b <- block("/path/to/raster.tif")
#' blocks_in_bbox(b, c(100, 200, 500, 600))
#' blocks_in_bbox(b, c(100, 200, 500, 600), data = TRUE)
#' }
blocks_in_bbox <- new_generic("blocks_in_bbox", "x")

#' Find blocks intersecting an extent
#'
#' @param x A `block` object.
#' @param ext Numeric vector c(xmin, xmax, ymin, ymax).
#' @param data If TRUE, return full block_data(). If FALSE (default), return wk::rct.
#' @param ... Ignored.
#'
#' @return A `wk::rct` vector or data.frame of intersecting blocks.
#' @export
#' @examples
#' \dontrun{
#' b <- block("/path/to/raster.tif")
#' blocks_in_ext(b, c(100, 500, 200, 600))
#' }
blocks_in_ext <- new_generic("blocks_in_ext", "x")

#' Plot blocks
#'
#' Draw rectangles for specified blocks on a new or existing plot.
#'
#' @param x A `block` object.
#' @param i Block column indices (0-based). If NULL, all columns.
#' @param j Block row indices (0-based). If NULL, all rows.
#' @param add Logical; if `TRUE`, add to existing plot. Default `FALSE`.
#' @param ... Additional arguments passed to [wk::wk_plot()].
#'
#' @return Invisibly returns the `wk::rct` of plotted blocks.
#' @export
#' @examples
#' \dontrun{
#' b <- block("/path/to/raster.tif")
#' plot_block(b)                         # all blocks
#' plot_block(b, 0:2, 0:1, col = "lightblue")
#' }
plot_block <- new_generic("plot_block", "x")

#' Plot blocks within a bounding box
#'
#' Draw rectangles for all blocks that intersect a bounding box.
#'
#' @param x A `block` object.
#' @param bbox Numeric vector c(xmin, ymin, xmax, ymax).
#' @param add Logical; if `TRUE`, add to existing plot. Default `FALSE`.
#' @param plot_bbox Logical; if `TRUE`, also draw the query bbox. Default `TRUE`.
#' @param bbox_args List of arguments passed to [graphics::rect()] for the bbox rectangle.
#' @param ... Additional arguments passed to [wk::wk_plot()] for block rectangles.
#'
#' @return Invisibly returns a list with `rct` (wk::rct) and `data` (data.frame).
#' @export
plot_block_bbox <- new_generic("plot_block_bbox", "x")

#' Plot blocks within an extent
#'
#' Draw rectangles for all blocks that intersect an extent.
#' Uses extent-style coordinate ordering (xmin, xmax, ymin, ymax).
#'
#' @param x A `block` object.
#' @param ext Numeric vector c(xmin, xmax, ymin, ymax).
#' @param add Logical; if `TRUE`, add to existing plot. Default `FALSE`.
#' @param plot_ext Logical; if `TRUE`, also draw the query extent. Default `TRUE`.
#' @param ext_args List of arguments passed to [graphics::rect()] for the extent rectangle.
#' @param ... Additional arguments passed to [wk::wk_plot()] for block rectangles.
#'
#' @return Invisibly returns a list with `rct` (wk::rct) and `data` (data.frame).
#' @export
plot_block_ext <- new_generic("plot_block_ext", "x")


# Class definition
#' @export
block <- new_class(
  "block",
  properties = list(
    ds = new_property(
      class = class_any,
      getter = function(self) self@.ds,
      setter = function(self, value) {
        stop("@ds is read-only", call. = FALSE)
      }
    ),
    dsn = new_property(
      class = class_character,
      getter = function(self) self@.dsn
    ),
    bbox = new_property(
      class = class_numeric,
      getter = function(self) {
        gt <- self@.ds$getGeoTransform()
        dm <- self@dimension
        c(
          xmin = gt[1],
          ymin = gt[4] + dm[2] * gt[6],
          xmax = gt[1] + dm[1] * gt[2],
          ymax = gt[4]
        )
      }
    ),
    dimension = new_property(
      class = class_integer,
      getter = function(self) {
        c(
          ncol = self@.ds$getRasterXSize(),
          nrow = self@.ds$getRasterYSize(),
          nbands = self@.ds$getRasterCount()
        )
      }
    ),
    crs = new_property(
      class = class_character,
      getter = function(self) self@.ds$getProjection()
    ),
    nbands = new_property(
      class = class_integer,
      getter = function(self) self@.ds$getRasterCount()
    ),
    datatype = new_property(
      class = class_character,
      getter = function(self) {
        vapply(
          seq_len(self@nbands),
          function(i) self@.ds$getDataTypeName(i),
          character(1)
        )
      }
    ),
    res = new_property(
      class = class_numeric,
      getter = function(self) {
        gt <- self@.ds$getGeoTransform()
        c(xres = gt[2], yres = abs(gt[6]))
      }
    ),
    blocksize = new_property(
      class = class_integer,
      getter = function(self) {
        bs <- self@.ds$getBlockSize(1L)
        c(width = bs[1], height = bs[2])
      }
    ),
    nblocks = new_property(
      class = class_integer,
      getter = function(self) {
        dm <- self@dimension
        bs <- self@blocksize
        c(
          nx = as.integer(ceiling(dm[1] / bs[1])),
          ny = as.integer(ceiling(dm[2] / bs[2]))
        )
      }
    ),
    blocks = new_property(
      class = class_list,
      getter = function(self) {
        dm <- self@dimension
        bs <- self@blocksize
        nb <- self@nblocks

        # Edge block sizes (may differ from standard)
        edge_width <- dm[1] - (nb[1] - 1L) * bs[1]
        edge_height <- dm[2] - (nb[2] - 1L) * bs[2]

        list(
          size = bs,
          count = nb,
          edge_width = as.integer(edge_width),
          edge_height = as.integer(edge_height)
        )
      }
    ),
    # Private properties (prefixed with .)
    .ds = new_property(class = class_any),
    .dsn = new_property(class = class_character)
  ),
  constructor = function(dsn) {
    ds <- new(gdalraster::GDALRaster, dsn, read_only = TRUE)
    new_object(S7_object(), .ds = ds, .dsn = dsn)
  }
)

method(print, block) <- function(x, ..., width = getOption("width", 80L)) {

  dm <- x@dimension
  bb <- x@bbox
  rs <- x@res
  bs <- x@blocksize
  nb <- x@nblocks
  blk <- x@blocks

  cat("block object\n")

  cat(" DSN      :", x@dsn, "\n")
  cat(" Dim      :", paste(dm, collapse = ", "), "\n")
  cat(" Res      :", sprintf("%.6f, %.6f", rs[1], rs[2]), "\n")
  cat(" Bbox     :", sprintf("%.6f, %.6f, %.6f, %.6f", bb[1], bb[2], bb[3], bb[4]), "\n")
  cat(" Datatype :", paste(x@datatype, collapse = ", "), "\n")
  cat(" Blocks   :", sprintf("%d x %d (%d x %d blocks)", bs[1], bs[2], nb[1], nb[2]), "\n")

  # Show edge block info if different from standard
  if (blk$edge_width != bs[1] || blk$edge_height != bs[2]) {
    cat(" Edges    :", sprintf("right=%d, bottom=%d", blk$edge_width, blk$edge_height), "\n")
  }

  # CRS at bottom with pretty formatting
  crs_wkt <- gdalraster::srs_to_wkt(x@.ds$getSpatialRef(), pretty = TRUE)
  if (!is.null(crs_wkt) && nzchar(crs_wkt)) {
    crs_lines <- strsplit(crs_wkt, "\n")[[1]]
    prefix <- " CRS      : "
    cont_prefix <- "           "  # same width as prefix
    max_content <- width - nchar(prefix)

    for (i in seq_along(crs_lines)) {
      line <- crs_lines[i]
      pfx <- if (i == 1L) prefix else cont_prefix
      # Truncate if too wide
      if (nchar(line) > max_content) {
        line <- paste0(substr(line, 1, max_content - 3), "...")
      }
      cat(pfx, line, "\n", sep = "")
    }
  } else {
    cat(" CRS      : (not set)\n")
  }

  invisible(x)
}


# -----------------------------------------------------------------------------
# Core block accessors
# -----------------------------------------------------------------------------

method(read_block, block) <- function(x, i, j, band = 1L) {
  idx <- block_index(x, i, j)
  x@.ds$read(
    band = as.integer(band),
    xoff = idx["xoff"],
    yoff = idx["yoff"],
    xsize = idx["xsize"],
    ysize = idx["ysize"],
    out_xsize = idx["xsize"],
    out_ysize = idx["ysize"]
  ) |>
    matrix(nrow = idx["xsize"], ncol = idx["ysize"]) |>
    t()
}

method(block_dim, block) <- function(x, i, j) {
  idx <- block_index(x, i, j)
  c(ncol = idx["xsize"], nrow = idx["ysize"])
}

method(block_bbox, block) <- function(x, i, j) {
  idx <- block_index(x, i, j)
  gt <- x@.ds$getGeoTransform()

  xmin <- gt[1] + idx["xoff"] * gt[2]
  xmax <- gt[1] + (idx["xoff"] + idx["xsize"]) * gt[2]
  ymax <- gt[4] + idx["yoff"] * gt[6]
  ymin <- gt[4] + (idx["yoff"] + idx["ysize"]) * gt[6]

  c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)
}

method(block_index, block) <- function(x, i, j) {
  bs <- unname(x@blocksize)
  dm <- unname(x@dimension)
  nb <- unname(x@nblocks)

  i <- as.integer(i)
  j <- as.integer(j)

  if (i < 0L || i >= nb[1]) stop("i out of range [0, ", nb[1] - 1L, "]", call. = FALSE)
  if (j < 0L || j >= nb[2]) stop("j out of range [0, ", nb[2] - 1L, "]", call. = FALSE)

  xoff <- i * bs[1]
  yoff <- j * bs[2]

  # Handle edge blocks
  xsize <- if (i == nb[1] - 1L) dm[1] - xoff else bs[1]
  ysize <- if (j == nb[2] - 1L) dm[2] - yoff else bs[2]

  c(xoff = as.integer(xoff), yoff = as.integer(yoff),
    xsize = as.integer(xsize), ysize = as.integer(ysize))
}


# -----------------------------------------------------------------------------
# Block geometry as data structures
# -----------------------------------------------------------------------------

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


# -----------------------------------------------------------------------------
# Block intersection queries
# -----------------------------------------------------------------------------

#' Convert bbox to block index ranges (internal helper)
#' @noRd
bbox_to_block_range <- function(x, bbox) {
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
  list(
    i_min = max(0L, as.integer(i_min)),
    i_max = min(nb[1] - 1L, as.integer(i_max)),
    j_min = max(0L, as.integer(j_min)),
    j_max = min(nb[2] - 1L, as.integer(j_max))
  )
}

method(blocks_in_bbox, block) <- function(x, bbox, data = FALSE, ...) {
  if (length(bbox) != 4L) {
    stop("'bbox' must be c(xmin, ymin, xmax, ymax)", call. = FALSE)
  }

  rng <- bbox_to_block_range(x, bbox)

  if (rng$i_min > rng$i_max || rng$j_min > rng$j_max) {
    if (data) {
      return(block_data(x, integer(0), integer(0)))
    } else {
      return(wk::rct(crs = x@crs))
    }
  }

  if (data) {
    block_data(x, rng$i_min:rng$i_max, rng$j_min:rng$j_max)
  } else {
    block_rct(x, rng$i_min:rng$i_max, rng$j_min:rng$j_max)
  }
}

method(blocks_in_ext, block) <- function(x, ext, data = FALSE, ...) {
  if (length(ext) != 4L) {
    stop("'ext' must be c(xmin, xmax, ymin, ymax)", call. = FALSE)
  }
  # Convert ext (xmin, xmax, ymin, ymax) to bbox (xmin, ymin, xmax, ymax)
  bbox <- ext[c(1L, 3L, 2L, 4L)]
  blocks_in_bbox(x, bbox, data = data, ...)
}


# -----------------------------------------------------------------------------
# Plotting - thin wrappers around block_rct / blocks_in_*
# -----------------------------------------------------------------------------

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

method(plot_block_bbox, block) <- function(x, bbox, add = FALSE,
                                           plot_bbox = TRUE,
                                           bbox_args = list(border = "red", lwd = 2),
                                           ...) {
  if (length(bbox) != 4L) {
    stop("'bbox' must be c(xmin, ymin, xmax, ymax)", call. = FALSE)
  }
  if (bbox[1] >= bbox[3]) {
    stop("'bbox' xmin must be less than xmax", call. = FALSE)
  }
  if (bbox[2] >= bbox[4]) {
    stop("'bbox' ymin must be less than ymax", call. = FALSE)
  }

  rcts <- blocks_in_bbox(x, bbox, data = FALSE)
  block_df <- blocks_in_bbox(x, bbox, data = TRUE)

  # Setup plot if needed
  if (!add) {
    full_bbox <- x@bbox
    plot(NULL, xlim = full_bbox[c(1, 3)], ylim = full_bbox[c(2, 4)],
         asp = 1, xlab = "x", ylab = "y")
  }

  # Draw block rectangles
  if (length(rcts) > 0L) {
    wk::wk_plot(rcts, add = TRUE, ...)
  }

  # Draw query bbox
  if (plot_bbox) {
    do.call(rect, c(list(xleft = bbox[1], ybottom = bbox[2],
                         xright = bbox[3], ytop = bbox[4]), bbox_args))
  }

  invisible(list(rct = rcts, data = block_df))
}

method(plot_block_ext, block) <- function(x, ext, add = FALSE,
                                          plot_ext = TRUE,
                                          ext_args = list(border = "red", lwd = 2),
                                          ...) {
  if (length(ext) != 4L) {
    stop("'ext' must be c(xmin, xmax, ymin, ymax)", call. = FALSE)
  }
  if (ext[1] >= ext[2]) {
    stop("'ext' xmin must be less than xmax", call. = FALSE)
  }
  if (ext[3] >= ext[4]) {
    stop("'ext' ymin must be less than ymax", call. = FALSE)
  }

  # Convert ext (xmin, xmax, ymin, ymax) to bbox (xmin, ymin, xmax, ymax)
  bbox <- ext[c(1L, 3L, 2L, 4L)]

  plot_block_bbox(x, bbox, add = add, plot_bbox = plot_ext,
                  bbox_args = ext_args, ...)
}
