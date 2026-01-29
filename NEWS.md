# gdalblock (development version)

## New features

* `block_rct()` returns block geometries as `wk::rct` vectors with CRS attached.
  This enables direct plotting via `plot(block_rct(b))` and integration with

  the wk ecosystem.

* `block_data()` returns a data.frame with full block metadata: i, j indices,
  pixel offsets (xoff, yoff, xsize, ysize), geographic extents (xmin, ymin,
  xmax, ymax), and a geometry column as `wk::rct`.

* `blocks_in_bbox()` and `blocks_in_ext()` find blocks intersecting a query

  region, returning either `wk::rct` (default) or full metadata via `data = TRUE`.

* Plotting functions (`plot_block()`, `plot_block_bbox()`, `plot_block_ext()`)
  now use `wk::wk_plot()` internally and return `wk::rct` objects invisibly.

## Dependencies
* Now imports `wk` for geometry representation.
