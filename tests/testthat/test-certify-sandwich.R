# The fast certification path (R/certify-sandwich.R). Its whole value is that
# it reaches the same verdicts as a sweeping route for less work, so most of
# what is worth testing is agreement: with brute force where brute force is
# affordable, and with the corridor DP where it is not. The corridor DP is no
# longer shipped -- the fast path settled every cell of a 2.7M-cell sweep
# without exhausting its budget, so the fallback was dead weight in R/ -- and
# now lives in helper-corridor-dp.R purely as this suite's oracle.

brute_states <- function(n, W) {
  grid <- as.matrix(expand.grid(rep(list(0:W), n)))
  data.frame(S = rowSums(grid), Q = rowSums(grid^2))
}

test_that("the sandwich endpoints are the true extremes of the sum of squares", {
  for (n in 2:4) for (W in 1:5) {
    b <- brute_states(n, W)
    for (s in 0:(n * W)) {
      qs <- b$Q[b$S == s]
      expect_identical(min(qs), .q_min_int(s, n))
      expect_identical(max(qs), .q_max_int(s, W))
      # every sample has Q = S (mod 2), which is what makes the parity screen
      # sound rather than merely usually right
      expect_true(all(qs %% 2 == s %% 2))
    }
  }
})

test_that("the search decides every window a small design admits", {
  # exhaustive in both directions: for every sum and every narrow window, the
  # search must agree with enumeration, and any witness it returns must be a
  # real sample hitting the window. Tens of thousands of cases, so the verdicts
  # are accumulated and asserted once rather than one expectation at a time.
  #
  # This is also the only coverage of the lazy ordering's slow path. Children
  # are ordered on a frame's second visit, not its first, and an attainable
  # report is usually reached with no backtracking at all -- so the branch that
  # calls order() is never taken on the happy path. Some of the searches below
  # do take it; validation/validate_sd_bounds_functions.qmd counts how many,
  # rather than a number being asserted here where it would go stale unseen.
  disagreed <- 0L
  bad_witness <- 0L
  for (n in 2:4) for (W in 1:5) {
    b <- brute_states(n, W)
    for (s in 0:(n * W)) {
      ach <- sort(unique(b$Q[b$S == s]))
      for (klo in 0:(n * W * W)) for (khi in klo:min(klo + 2, n * W * W)) {
        got <- .witness_search(W, n, s, klo, khi)
        if (!identical(got$possible, any(ach >= klo & ach <= khi)))
          disagreed <- disagreed + 1L
        if (isTRUE(got$possible)) {
          w <- got$witness
          if (length(w) != n || any(w < 0 | w > W) || sum(w) != s ||
              sum(w^2) < klo || sum(w^2) > khi)
            bad_witness <- bad_witness + 1L
        }
      }
    }
  }
  expect_identical(disagreed, 0L)
  expect_identical(bad_witness, 0L)
})

test_that("deep ladders at tiny n are not smoothed over", {
  # n = 2 on 0-6 with sum 6: the pairs (3,3), (2,4), (1,5), (0,6) give sums of
  # squares 18, 20, 26, 36 -- two wide gaps. Any band logic that assumed the
  # achievable set was dense between floor and ceiling would call the gaps
  # possible, so this is the guard against exactly that.
  ok <- function(lo, hi) .witness_search(6, 2, 6, lo, hi)$possible
  expect_true(ok(18, 18)); expect_true(ok(20, 20))
  expect_true(ok(26, 26)); expect_true(ok(36, 36))
  expect_false(ok(21, 25)); expect_false(ok(27, 35))
  expect_false(ok(22, 22)); expect_false(ok(30, 30))
})

test_that("a binary scale admits only Q = S", {
  # W = 1 makes every value its own square, so the sum of squares is the sum
  for (s in 0:8) {
    expect_true(.witness_search(1, 8, s, s, s)$possible)
    if (s > 0) expect_false(.witness_search(1, 8, s, s - 1, s - 1)$possible)
  }
})

test_that("parity alone settles reports the sandwich admits", {
  # 0-100, n = 100, mean 99.80, sd 0.50: the sum is 9980 and the only integer
  # in the SD's window is 996029, which is odd while the sum is even
  W <- 100; n <- 100
  tg <- .target_states(0, 100, n, 1, 99.80, 0.50, 2, 2, "half_up")
  expect_true(nrow(tg) >= 1)
  expect_null(.k_window(tg$S[1], tg$Q_lo[1], tg$Q_hi[1],
                        tg$lo_incl[1], tg$hi_incl[1], n, W))
  expect_false(.certify_fast(W, n, tg))
})

test_that("the sd = 0 sign trap and the empty candidate range survive", {
  # unround_interval(0, ...) reaches below zero; squaring that endpoint would
  # put a positive floor under the sum of squares and exclude the constant
  # sample, which is the one sample a reported SD of 0 must admit
  W <- 6; n <- 12
  tg <- .target_states(0, 6, n, 1, 3.0, 0.0, 1, 1, "half_up")
  expect_true(.certify_fast(W, n, tg))
  # a mean admitting no integer sum yields no targets at all, so the fast path
  # must report impossible rather than iterating a descending sequence
  tg0 <- .target_states(0, 6, n, 1, 3.51, 1.0, 2, 1, "half_up")
  expect_false(.certify_fast(W, n, tg0))
})

test_that("the fast path agrees with the corridor DP cell for cell", {
  for (cfg in list(c(l = 1, u = 5, n = 9), c(l = 0, u = 6, n = 12),
                   c(l = 1, u = 7, n = 11))) {
    l <- cfg[["l"]]; u <- cfg[["u"]]; n <- cfg[["n"]]; W <- u - l
    gr <- expand.grid(mean = seq(l, u, by = 0.1),
                      sd = seq(0, (u - l) / 2, by = 0.1))
    for (rr in c("half_up", "half_down")) {
      fast <- rep(NA, nrow(gr)); slow <- fast
      for (i in seq_len(nrow(gr))) {
        tg <- .target_states(l, u, n, 1, gr$mean[i], gr$sd[i], 1, 1, rr)
        fast[i] <- .certify_fast(W, n, tg)
        slow[i] <- .attainable_target(W, n, tg)
      }
      expect_false(anyNA(fast))
      expect_identical(fast, slow)
    }
  }
})

test_that("an exhausted search budget errors rather than guessing", {
  # with no fallback route left, a report the search cannot settle must say so
  # rather than return a verdict it has not earned. search_budget = 0 refuses
  # the search everything, which is the only way to provoke this on demand.
  args <- list(l = 0, u = 20, n = 30, digits = 1,
               mean = c(10.0, 9.7, 0.2, 19.9, 3.4),
               sd   = c(5.0, 0.1, 0.5, 0.2, 12.0))
  expect_true(all(do.call(brimmest, args)$possible %in% c(TRUE, FALSE)))
  expect_error(do.call(brimmest, c(args, list(search_budget = 0))),
               "search_budget")
})

test_that("a design the lattice refuses is certified anyway", {
  # 0-63 at n = 50 needs a 3151 x 24801 state table, so only the targeted
  # route can answer at all -- and the fast path answers immediately
  r <- brimmest(l = 0, u = 63, n = 50, mean = 30.42, sd = 12.71, digits = 2)
  expect_true(r$possible)
})

test_that("the search is not bounded by R's evaluation depth", {
  # the tree is one level deep per observation placed, so a recursive search
  # died here long before the node budget bit
  n <- 3000
  S <- n * 5
  k <- .q_min_int(S, n) + 200
  r <- .witness_search(10, n, S, k, k + 40)
  expect_true(r$possible)
  w <- r$witness
  expect_length(w, n)
  expect_equal(sum(w), S)
  expect_true(all(w >= 0 & w <= 10))
  expect_true(sum(w^2) >= k && sum(w^2) <= k + 40)
})
