# Version-robust wrappers around scrutiny's GRIM and GRIMMER.
#
# scrutiny changed the interface of grim() and grimmer() after 0.6.1, and the
# two forms are mutually exclusive:
#
#   * scrutiny <= 0.6.1 (currently the CRAN version) takes the reported mean
#     and SD as zero-padded strings, e.g. grim(x = "5.10", n = 28), and reads
#     the reported decimal places off the string. It does not know digits_x.
#   * later versions take them as numbers together with their reported decimal
#     places, e.g. grim(x = 5.1, n = 28, digits_x = 2), and reject a character
#     x outright ("argument is of length zero").
#
# So there is no single call that works against both. Detect which interface
# the installed version exposes and dispatch accordingly, keeping the package
# working against CRAN scrutiny and newer builds alike. The detection is per
# function rather than global, so a release that migrates one before the other
# is still handled.
#
# scrutiny::restore_zeros() builds the zero-padded strings for the older
# interface; it is present, exported and identical in behaviour in both
# generations (restore_zeros(c(5.1, 5.19), width = 2) gives "5.10" "5.19").

# Internal: does this scrutiny entry point take the digits_x-style arguments?
.takes_digits <- function(f) "digits_x" %in% names(formals(f))

# Internal: scrutiny::grim() across both interfaces. `digits` is the number of
# decimal places of the reported mean. Returns scrutiny's verdict unchanged.
.grim_compat <- function(x, n, digits, items = 1, rounding = "up_or_down") {
  if (.takes_digits(scrutiny::grim)) {
    scrutiny::grim(
      x = x,
      n = n,
      digits_x = digits,
      items = items,
      rounding = rounding
    )
  } else {
    scrutiny::grim(
      x = scrutiny::restore_zeros(x, width = digits),
      n = n,
      items = items,
      rounding = rounding
    )
  }
}

# Internal: scrutiny::grimmer() across both interfaces. `digits_x`/`digits_sd`
# are the decimal places of the reported mean and SD. Returns scrutiny's
# verdict unchanged.
.grimmer_compat <- function(
  x,
  sd,
  n,
  digits_x,
  digits_sd,
  items = 1,
  rounding = "up_or_down"
) {
  if (.takes_digits(scrutiny::grimmer)) {
    scrutiny::grimmer(
      x = x,
      sd = sd,
      n = n,
      digits_x = digits_x,
      digits_sd = digits_sd,
      items = items,
      rounding = rounding
    )
  } else {
    scrutiny::grimmer(
      x = scrutiny::restore_zeros(x, width = digits_x),
      sd = scrutiny::restore_zeros(sd, width = digits_sd),
      n = n,
      items = items,
      rounding = rounding
    )
  }
}
