test_that("gdalblock constructor works with system file", {
  # Use a GDAL test file if available, skip otherwise
  skip_if_not_installed("gdalraster")

  # Create a temp test file
  tf <- tempfile(fileext = ".tif")
  on.exit(unlink(tf), add = TRUE)

  # Create a simple test raster using gdalraster
  driver <- gdalraster::create(
    format = "GTiff",
    dst_filename = tf,
    xsize = 100,
    ysize = 80,
    nbands = 1,
    dataType = "Int16", return_obj = TRUE
  )

  driver$setGeoTransform(c(0, 1, 0, 80, 0, -1))
  driver$setProjection('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]')
  driver$write(band = 1, xoff = 0, yoff = 0, xsize = 100, ysize = 80, seq_len(8000))
  driver$close()

  b <- gdalblock(tf)

  #expect_s3_class(b, "gdalblock")
  expect_equal(b@dsn, tf)
  expect_equal(b@dimension[["ncol"]], 100L)
  expect_equal(b@dimension[["nrow"]], 80L)
  expect_equal(b@nbands, 1L)
})

test_that("block_index returns correct offsets", {
  skip_if_not_installed("gdalraster")

  tf <- tempfile(fileext = ".tif")
  on.exit(unlink(tf), add = TRUE)

  driver <- gdalraster::create(
    format = "GTiff",
    dst_filename = tf,
    xsize = 100,
    ysize = 80,
    nbands = 1,
    dataType = "Int16",
    options = c("TILED=YES", "BLOCKXSIZE=32", "BLOCKYSIZE=32"), return_obj = TRUE
  )
  driver$setGeoTransform(c(0, 1, 0, 80, 0, -1))
  driver$close()

  b <- gdalblock(tf)

  # First block
  idx <- block_index(b, 0, 0)
  expect_equal(idx[["xoff"]], 0L)
  expect_equal(idx[["yoff"]], 0L)
  expect_equal(idx[["xsize"]], 32L)
  expect_equal(idx[["ysize"]], 32L)

  # Edge block (rightmost)
  idx_edge <- block_index(b, 3, 0)
  expect_equal(idx_edge[["xoff"]], 96L)
  expect_equal(idx_edge[["xsize"]], 4L)  # 100 - 96 = 4
})

test_that("block_bbox returns correct coordinates",
{
  skip_if_not_installed("gdalraster")

  tf <- tempfile(fileext = ".tif")
  on.exit(unlink(tf), add = TRUE)

  driver <- gdalraster::create(
    format = "GTiff",
    dst_filename = tf,
    xsize = 100,
    ysize = 80,
    nbands = 1,
    dataType = "Int16",
    options = c("TILED=YES", "BLOCKXSIZE=32", "BLOCKYSIZE=32"),
    return_obj = TRUE
  )
  driver$setGeoTransform(c(0, 1, 0, 80, 0, -1))
  driver$close()

  b <- gdalblock(tf)

  bb <- block_bbox(b, 0, 0)
  expect_equal(bb[[1]], 0)
  expect_equal(bb[[3]], 32)
  expect_equal(bb[[2]], 48)  # 80 - 32 = 48
  expect_equal(bb[[4]], 80)
})

test_that("read_block returns matrix of correct size", {
  skip_if_not_installed("gdalraster")

  tf <- tempfile(fileext = ".tif")
  on.exit(unlink(tf), add = TRUE)

  driver <- gdalraster::create(
    format = "GTiff",
    dst_filename = tf,
    xsize = 100,
    ysize = 80,
    nbands = 1,
    dataType = "Int16",
    options = c("TILED=YES", "BLOCKXSIZE=32", "BLOCKYSIZE=32"),
    return_obj = TRUE
  )
  driver$setGeoTransform(c(0, 1, 0, 80, 0, -1))
  driver$write(band = 1, xoff = 0, yoff = 0, xsize = 100, ysize = 80, seq_len(8000))
  driver$close()

  b <- gdalblock(tf)

  mat <- read_block(b, 0, 0)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 32)
  expect_equal(ncol(mat), 32)

  # Edge block
  mat_edge <- read_block(b, 3, 2)
  expect_equal(ncol(mat_edge), 4)  # edge width
  expect_equal(nrow(mat_edge), 16)  # edge height (80 - 64 = 16)
})

test_that("out of range blocks throw errors", {
  skip_if_not_installed("gdalraster")

  tf <- tempfile(fileext = ".tif")
  on.exit(unlink(tf), add = TRUE)

  driver <- gdalraster::create(
    format = "GTiff",
    dst_filename = tf,
    xsize = 100,
    ysize = 80,
    nbands = 1,
    dataType = "Int16",
    options = c("TILED=YES", "BLOCKXSIZE=32", "BLOCKYSIZE=32"),
    return_obj = TRUE
  )
  driver$close()

  b <- gdalblock(tf)

  expect_error(block_index(b, -1, 0), "out of range")
  expect_error(block_index(b, 0, -1), "out of range")
  expect_error(block_index(b, 10, 0), "out of range")
  expect_error(block_index(b, 0, 10), "out of range")
})
test_that("gdalblock class exists", {
  expect_true(inherits(gdalblock, "S7_class"))
})
