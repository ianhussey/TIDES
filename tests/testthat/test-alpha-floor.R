test_that("a reported positive alpha forces SD > 0 at a whole-number composite mean", {
  # 2-item 0-3 composite (sum window 0..6, affordable), sum mean 3 is a whole number,
  # where the amplified floor vanishes but the Gini envelope does not.
  amped <- sd_bounds(l = 0, u = 3, n = 10, mean = 1.5, granularity = "integer",
                     scoring = "meanscored", n_items = 2)              # alpha-free
  withA <- sd_bounds(l = 0, u = 3, n = 10, mean = 1.5, granularity = "integer",
                     scoring = "meanscored", n_items = 2, alpha = 0.7)
  expect_equal(amped$min_sd, 0)                    # alpha-free floor is 0 here
  expect_gt(withA$min_sd, 0)                       # a positive alpha lifts it off zero
  expect_match(withA$min_rule, "Gini envelope")
})

test_that("the Gini envelope reproduces the article's disproof counterexample", {
  # composite {0, 3,3,3, 4x10}, n=14: SS_S = 15.5, m_max = 2.7468
  env <- strait:::sd_min_alpha_gini(l = 0, u = 4, n = 14, mean = 3.5, m = 2.7468)
  expect_equal(env, sqrt(15.5 / 13), tolerance = 1e-6)
})

test_that("the Gini envelope is a valid lower bound and sharp (brute force)", {
  # exhaustive: n=4 people, k=2 integer items in 1..3
  n <- 4; k <- 2; il <- 1; iu <- 3
  G <- as.matrix(expand.grid(rep(list(il:iu), n * k)))
  person  <- ((seq_len(n * k) - 1) %/% k) + 1
  itemidx <- ((seq_len(n * k) - 1) %%  k) + 1
  comp <- sapply(1:n, function(p) rowSums(G[, person == p, drop = FALSE]))
  compmean <- rowMeans(comp)
  SS_S <- rowSums((comp - compmean)^2)
  V <- rowSums(sapply(1:k, function(it) {
    Xi <- G[, itemidx == it, drop = FALSE]
    apply(Xi, 1, function(r) sum((r - mean(r))^2))
  }))
  m_d  <- ifelse(V > 0, SS_S / V, Inf)
  sd_d <- sqrt(SS_S / (n - 1))
  # validity: envelope at each dataset's own (mean, design factor) never exceeds its SD
  viol <- 0L
  for (i in which(V > 1e-9 & SS_S > 1e-9)) {
    e <- strait:::sd_min_alpha_gini(l = k * il, u = k * iu, n = n,
                                   mean = compmean[i], m = m_d[i])
    if (!is.null(e) && e > sd_d[i] + 1e-9) viol <- viol + 1L
  }
  expect_equal(viol, 0L)
  # sharpness at a whole-number mean
  mu0 <- k * 2
  for (mt in c(1.2, 1.5, 2.0)) {
    cand <- which(abs(compmean - mu0) < 1e-9 & m_d >= mt - 1e-9 & SS_S > 1e-9)
    brute <- min(sd_d[cand])
    e <- strait:::sd_min_alpha_gini(l = k * il, u = k * iu, n = n, mean = mu0, m = mt)
    expect_equal(e, brute, tolerance = 1e-9)
  }
})

test_that("an unaffordable composite window falls back to the proven floor with a note", {
  r <- sd_bounds(l = 1, u = 5, n = 20, mean = 3, granularity = "integer",
                 scoring = "meanscored", n_items = 3, alpha = 0.8)   # W = 12 over budget
  expect_true(r$feasible)
  expect_false(is.na(r$note))
  expect_match(r$note, "not evaluated")
})
