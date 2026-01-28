#' Block class
#'
#' An S7 class representing GDAL raster block structure.
#'
#' @param dsn Data source name (file path or connection string).
#' @return A `block` object.
#' @export
#' @examples
#' \dontrun{
#' b <- block("/path/to/raster.tif")
#' read_block(b, 0, 0)
#' }

#' @export
read_block <- new_generic("read_block", "x")

#' @export
block_index <- new_generic("block_index", "x")

#' @export
block_bbox <- new_generic("block_bbox", "x")

block <- new_class(
  "gdalblock",
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

method(read_block, block) <- function(x, i, j) {
# ...existing code...
}

method(block_bbox, block) <- function(x, i, j) {
# ...existing code...
}

method(block_index, block) <- function(x, i, j) {
# ...existing code...
}

method(print, block) <- function(x, ..., width = getOption("width", 80L)) {

  dm <- x@dimension
  bb <- x@bbox
  rs <- x@res
  bs <- x@blocksize
  nb <- x@nblocks
  blk <- x@blocks

  cat("gdalblock object\n")
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

#' Plot blocks
#'
#' Draw rectangles for specified blocks on a new or existing plot.
#'
#' @param x A `block` object.
#' @param i Block column index (0-based). Can be a vector, or a 2-column matrix
#'   of (i, j) pairs (in which case `j` is ignored).
#' @param j Block row index (0-based). Can be a vector. If both `i` and `j` are
#'   vectors, they are expanded to a grid via `expand.grid()`.
#' @param add Logical; if `TRUE`, add to existing plot. Default `FALSE`.
#' @param ... Additional arguments passed to [rect()].
#'
#' @return Invisibly returns the bounding boxes of drawn blocks as a matrix.
#' @export
#' @examples
#' \dontrun{
#' b <- block("/path/to/raster.tif")
#' # Plot a single block
#' plot_block(b, 0, 0)
#' # Plot a grid of blocks
#' plot_block(b, 0:2, 0:1, col = "lightblue")
#' # Plot scattered blocks using matrix
#' plot_block(b, cbind(c(0, 2, 4), c(1, 3, 0)), border = "red")
#' }
plot_block <- new_generic("plot_block", "x")

#' Plot blocks within a bounding box
#'
#' Draw rectangles for all blocks that partly or wholly intersect a bounding box.
#'
#' @param x A `gdalblock` object.
#' @param bbox Numeric vector c(xmin, ymin, xmax, ymax).
#' @param add Logical; if `TRUE`, add to existing plot. Default `FALSE`.
#' @param plot_bbox Logical; if `TRUE`, also draw the query bbox. Default `TRUE`.
#' @param bbox_args List of arguments passed to [rect()] for the bbox rectangle.
#' @param ... Additional arguments passed to [rect()] for block rectangles.
#'
#' @return Invisibly returns a list with `blocks` (matrix of i,j pairs) and
#'   `bboxes` (matrix of block bounding boxes).
#' @export
#' @examples
#' \dontrun{
#' b <- gdalblock("/path/to/raster.tif")
#' # Plot blocks intersecting a region
#' plot_block_bbox(b, c(100, 200, 500, 600))
#' # Customize appearance
#' plot_block_bbox(b, c(100, 200, 500, 600),
#'                 col = "lightgray",
#'                 bbox_args = list(border = "red", lwd = 2))
#' }
plot_block_bbox <- new_generic("plot_block_bbox", "x")

#' Plot blocks within an extent
#'
#' Draw rectangles for all blocks that partly or wholly intersect an extent.
#' This is identical to [plot_block_bbox()] but uses extent-style coordinate
#' ordering (xmin, xmax, ymin, ymax) instead of bbox-style (xmin, ymin, xmax, ymax).
#'
#' @param x A `gdalblock` object.
#' @param ext Numeric vector c(xmin, xmax, ymin, ymax).
#' @param add Logical; if `TRUE`, add to existing plot. Default `FALSE`.
#' @param plot_ext Logical; if `TRUE`, also draw the query extent. Default `TRUE`.
#' @param ext_args List of arguments passed to [rect()] for the extent rectangle.
#' @param ... Additional arguments passed to [rect()] for block rectangles.
#'
#' @return Invisibly returns a list with `blocks` (matrix of i,j pairs) and
#'   `bboxes` (matrix of block bounding boxes in xmin, ymin, xmax, ymax order).
#' @export
#' @examples
#' \dontrun{
#' b <- gdalblock("/path/to/raster.tif")
#' # Plot blocks intersecting a region (extent style)
#' plot_block_ext(b, c(100, 500, 200, 600))
#' # Customize appearance
#' plot_block_ext(b, c(100, 500, 200, 600),
#'                col = "lightgray",
#'                ext_args = list(border = "red", lwd = 2))
#' }
plot_block_ext <- new_generic("plot_block_ext", "x")

method(plot_block, block) <- function(x, i, j = NULL, add = FALSE, ...) {
  # Handle matrix input for scattered blocks

if (is.matrix(i)) {
    if (ncol(i) != 2L) {
      stop("Matrix 'i' must have 2 columns (i, j pairs)", call. = FALSE)
    }
    pairs <- i
  } else {
    # Vectorized i and j - expand to grid
    if (is.null(j)) {
      stop("'j' is required when 'i' is not a matrix", call. = FALSE)
    }
    pairs <- as.matrix(expand.grid(i = as.integer(i), j = as.integer(j)))
  }

  # Get bboxes for all blocks
  bboxes <- t(apply(pairs, 1L, function(p) block_bbox(x, p[1L], p[2L])))
  colnames(bboxes) <- c("xmin", "ymin", "xmax", "ymax")

  # Setup plot if needed
  if (!add) {
    full_bbox <- x@bbox
    plot(NULL, xlim = full_bbox[c(1, 3)], ylim = full_bbox[c(2, 4)],
         asp = 1, xlab = "x", ylab = "y")
  }

  # Draw rectangles
  rect(bboxes[, "xmin"], bboxes[, "ymin"], bboxes[, "xmax"], bboxes[, "ymax"], ...)

  invisible(bboxes)
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

  # Get raster properties
  full_bbox <- x@bbox
  bs <- x@blocksize
  res <- x@res
  nb <- x@nblocks

  # Convert bbox to block indices
  # Block i covers pixels from i*bs[1] to (i+1)*bs[1]-1
  # In geo coords: xmin + i*bs[1]*res[1] to xmin + (i+1)*bs[1]*res[1]
  i_min <- floor((bbox[1] - full_bbox[1]) / (bs[1] * res[1]))
  i_max <- floor((bbox[3] - full_bbox[1]) / (bs[1] * res[1]))
  j_min <- floor((full_bbox[4] - bbox[4]) / (bs[2] * res[2]))
  j_max <- floor((full_bbox[4] - bbox[2]) / (bs[2] * res[2]))

  # Clamp to valid range
  i_min <- max(0L, as.integer(i_min))
  i_max <- min(nb[1] - 1L, as.integer(i_max))
  j_min <- max(0L, as.integer(j_min))
  j_max <- min(nb[2] - 1L, as.integer(j_max))

  # Generate block pairs
  if (i_min > i_max || j_min > j_max) {
    pairs <- matrix(integer(0), ncol = 2L)
    bboxes <- matrix(numeric(0), ncol = 4L)
  } else {
    pairs <- as.matrix(expand.grid(i = i_min:i_max, j = j_min:j_max))
    bboxes <- t(apply(pairs, 1L, function(p) block_bbox(x, p[1L], p[2L])))
  }
  colnames(bboxes) <- c("xmin", "ymin", "xmax", "ymax")

  # Setup plot if needed
  if (!add) {
    plot(NULL, xlim = full_bbox[c(1, 3)], ylim = full_bbox[c(2, 4)],
         asp = 1, xlab = "x", ylab = "y")
  }

  # Draw block rectangles
  if (nrow(bboxes) > 0L) {
    rect(bboxes[, "xmin"], bboxes[, "ymin"], bboxes[, "xmax"], bboxes[, "ymax"], ...)
  }

  # Draw query bbox
  if (plot_bbox) {
    do.call(rect, c(list(xleft = bbox[1], ybottom = bbox[2],
                         xright = bbox[3], ytop = bbox[4]), bbox_args))
  }

  invisible(list(blocks = pairs, bboxes = bboxes))
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

  # Delegate to plot_block_bbox

  plot_block_bbox(x, bbox, add = add, plot_bbox = plot_ext,
                  bbox_args = ext_args, ...)
}

