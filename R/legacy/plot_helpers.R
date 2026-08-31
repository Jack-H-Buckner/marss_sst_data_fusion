#' Bounding box around a set of lon/lat points, with a percentage buffer
#'
#' Buffers by a fraction of the data extent in each direction. Because a degree
#' of longitude shrinks with latitude, a naive percentage buffer produces a
#' visually lopsided map: `equal_km = TRUE` (the default) applies the same
#' buffer distance in kilometres on both axes, converting to degrees separately
#' for longitude and latitude.
#'
#' @param lon,lat Numeric vectors of coordinates (or pass `data` + column names).
#' @param data Optional data frame; if given, `lon` and `lat` are column names.
#' @param buffer Fractional buffer, e.g. 0.1 for 10 percent. Recycled, or give
#'   c(x, y) for different buffers per axis.
#' @param equal_km If TRUE, buffer the same distance in km on both axes.
#' @param min_buffer_km Floor on the buffer, in km. Prevents a degenerate box
#'   when all points are nearly coincident.
#' @param clamp Clip latitude to [-90, 90] and longitude to [-180, 180].
#'
#' @return Named numeric vector: xmin, ymin, xmax, ymax.
bbox_from_points <- function(lon, lat, data = NULL,
                             buffer = 0.10,
                             equal_km = TRUE,
                             min_buffer_km = 2,
                             clamp = TRUE) {
  
  if (!is.null(data)) {
    lon <- data[[deparse(substitute(lon))]]
    lat <- data[[deparse(substitute(lat))]]
  }
  
  ok <- is.finite(lon) & is.finite(lat)
  if (!any(ok)) stop("No finite coordinates supplied.", call. = FALSE)
  if (any(!ok)) message(sum(!ok), " point(s) dropped for non-finite coordinates.")
  lon <- lon[ok]; lat <- lat[ok]
  
  if (max(lon) - min(lon) > 180)
    warning("Longitude span exceeds 180 degrees; points may cross the ",
            "antimeridian and the box will be wrong.", call. = FALSE)
  
  buffer <- rep_len(buffer, 2)          # c(x, y)
  x_rng  <- range(lon); y_rng <- range(lat)
  x_span <- diff(x_rng); y_span <- diff(y_rng)
  
  km_per_deg_lat <- 111.32
  km_per_deg_lon <- 111.32 * cos(mean(y_rng) * pi / 180)
  
  if (equal_km) {
    # one buffer distance in km, from the larger of the two spans
    span_km <- max(x_span * km_per_deg_lon, y_span * km_per_deg_lat)
    buf_km  <- pmax(buffer * span_km, min_buffer_km)
    pad_x   <- buf_km[1] / km_per_deg_lon
    pad_y   <- buf_km[2] / km_per_deg_lat
  } else {
    pad_x <- pmax(buffer[1] * x_span, min_buffer_km / km_per_deg_lon)
    pad_y <- pmax(buffer[2] * y_span, min_buffer_km / km_per_deg_lat)
  }
  
  out <- c(xmin = x_rng[1] - pad_x, ymin = y_rng[1] - pad_y,
           xmax = x_rng[2] + pad_x, ymax = y_rng[2] + pad_y)
  
  if (clamp) {
    out[c("ymin", "ymax")] <- pmin(pmax(out[c("ymin", "ymax")], -90), 90)
    out[c("xmin", "xmax")] <- pmin(pmax(out[c("xmin", "xmax")], -180), 180)
  }
  out
}
