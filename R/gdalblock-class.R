#' gdalblock S7 Class
#'
#' An S7 class for block-based access to GDAL raster datasets.
#'
#' @param dsn A dataset description (file path, URL, or any GDAL-supported URI).
#'
#' @return A `gdalblock` object.
#'
#' @details
#' The gdalblock class wraps a `gdalraster::GDALRaster` object and provides
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
#'
#' @export
#' @examples
#' \dontrun{
#' b <- gdalblock("/path/to/raster.tif")
#' b@dimension
#' b@blocksize
#' read_block(b, 0, 0)
#' }

# Generics for methods - must be defined before class and method registrations
#' Read a block of raster data
#'
#' @param x A `gdalblock` object.
#' @param i Block column index (0-based).
#' @param j Block row index (0-based).
#' @param band Band number (1-based, default 1).
#'
#' @return A matrix of raster values.
#' @export
read_block <- new_generic("read_block", "x")

#' Get block dimensions
#'
#' @param x A `gdalblock` object.
#' @param i Block column index (0-based).
#' @param j Block row index (0-based).
#'
#' @return Integer vector c(ncol, nrow).
#' @export
block_dim <- new_generic("block_dim", "x")

#' Get block bounding box
#'
#' @param x A `gdalblock` object.
#' @param i Block column index (0-based).
#' @param j Block row index (0-based).
#'
#' @return Numeric vector c(xmin, ymin, xmax, ymax).
#' @export
block_bbox <- new_generic("block_bbox", "x")

#' Get block pixel/line index
#'
#' @param x A `gdalblock` object.
#' @param i Block column index (0-based).
#' @param j Block row index (0-based).
#'
#' @return Integer vector c(xoff, yoff, xsize, ysize).
#' @export
block_index <- new_generic("block_index", "x")

# Class definition
#' @export
gdalblock <- new_class(
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

method(print, gdalblock) <- function(x, ..., width = getOption("width", 80L)) {

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

method(read_block, gdalblock) <- function(x, i, j, band = 1L) {
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

method(block_dim, gdalblock) <- function(x, i, j) {
  idx <- block_index(x, i, j)
  c(ncol = idx["xsize"], nrow = idx["ysize"])
}

method(block_bbox, gdalblock) <- function(x, i, j) {
  idx <- block_index(x, i, j)
  gt <- x@.ds$getGeoTransform()

  xmin <- gt[1] + idx["xoff"] * gt[2]
  xmax <- gt[1] + (idx["xoff"] + idx["xsize"]) * gt[2]
  ymax <- gt[4] + idx["yoff"] * gt[6]
  ymin <- gt[4] + (idx["yoff"] + idx["ysize"]) * gt[6]

  c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)
}

method(block_index, gdalblock) <- function(x, i, j) {
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


