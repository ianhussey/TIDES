# Fast-path certification for brimmest(): decide one reported tuple by
# arithmetic and a pruned constructive search, rather than by sweeping a state
# space whose size is set by the design.
#
# The full attainable lattice and the corridor DP both answer "which states are
# reachable?" and then read the answer off. That is the right shape for a grid
# of reports and the wrong shape for one. This file answers the membership
# question directly, in three layers of increasing cost.
#
# Everything here works on the shifted integer scale y = multiplier * (x - l),
# so y is a whole number in [0, W] with W = multiplier * (u - l), and
#
#   S = sum(y)      the sample sum, pinned to a few integers by the mean
#   Q = sum(y^2)    the sum of squares, pinned to a window by the SD
#
# Layer 1 -- the sandwich screen, O(1) per candidate sum. For a given S the
# achievable Q values lie between the clustered configuration (q_min_int) and
# the Structure-S configuration (q_max_int). They also share S's parity,
# because y^2 = y (mod 2) makes Q = S (mod 2) for every sample. A window
# meeting no integer of the right parity inside that sandwich is impossible,
# decided in microseconds.
#
# Layer 2 -- a constructive search over non-increasing value sequences, which
# is to say over the partitions of S into at most n parts of size at most W.
# Each node re-applies the layer-1 sandwich to what is left to place, which
# prunes hard: near either wall of the scale the surviving tree is a handful of
# nodes, so those cells (where the closed-form screens leak) are settled
# outright.
#
# Layer 3 -- the same search run to exhaustion. Reaching a leaf is a proof of
# possibility and yields a witness sample; exhausting the tree is a proof of
# impossibility. Only when neither happens within the node budget is a verdict
# withheld.
#
# Two structural choices keep the tree small. Requiring the sequence to be
# non-increasing collapses the n! orderings of a sample to one, and reflecting
# y -> W - y whenever the sum sits above the midpoint means the search always
# builds from the nearer wall, where the partition tree is shallow. A third
# keeps it usable: the search is depth-first on an explicit stack rather than
# recursive, since it is one level deep per observation placed and a sample
# size in the thousands would otherwise exhaust R's evaluation depth long
# before the node budget bit.

# Internal: the smallest sum of squares of n non-negative integers summing to
# S -- the clustered configuration, values split between floor(S/n) and one
# more. The integer form of the sd_min_integer() floor.
q_min_int <- function(sum_total, n) {
  base_value <- sum_total %/% n
  remainder_count <- sum_total - n * base_value

  (n - remainder_count) *
    base_value *
    base_value +
    remainder_count * (base_value + 1) * (base_value + 1)
}

# Internal: the largest sum of squares of integers in [0, cap] summing to S --
# the Structure-S configuration, as many values as possible at cap, one
# remainder, the rest at 0. The integer form of sd_max_structure_s(). A cap of
# 0 forces S = 0, hence a sum of squares of 0.
q_max_int <- function(sum_total, cap) {
  cap <- rep(cap, length.out = length(sum_total))
  out <- numeric(length(sum_total))
  positive <- cap > 0

  if (any(positive)) {
    cap_positive <- cap[positive]
    sum_positive <- sum_total[positive]
    count_at_cap <- sum_positive %/% cap_positive
    remainder <- sum_positive - count_at_cap * cap_positive
    out[positive] <- count_at_cap * cap_positive * cap_positive + remainder^2
  }

  out
}

# Internal: layer 1. Narrow one candidate's Q window to the integers that could
# actually occur -- inside the sandwich and of the right parity. Returns
# c(window_lo, window_hi) or NULL when the candidate is already impossible.
#
# The window arrives as the real interval implied by the reported SD, with
# endpoint inclusion flags, so an endpoint that the rounding rule excludes must
# not be admitted even when it lands exactly on an integer.
k_window <- function(
  sum_total,
  q_lo,
  q_hi,
  lo_incl,
  hi_incl,
  n,
  span,
  tol = TOLERANCE
) {
  window_lo <- ceiling(q_lo - tol)
  if (!lo_incl && dplyr::near(window_lo, q_lo, tol = tol)) {
    window_lo <- window_lo + 1
  }
  window_hi <- floor(q_hi + tol)
  if (!hi_incl && dplyr::near(window_hi, q_hi, tol = tol)) {
    window_hi <- window_hi - 1
  }

  window_lo <- max(window_lo, q_min_int(sum_total, n))
  window_hi <- min(window_hi, q_max_int(sum_total, span))
  if ((window_lo %% 2) != (sum_total %% 2)) {
    window_lo <- window_lo + 1
  }
  if ((window_hi %% 2) != (sum_total %% 2)) {
    window_hi <- window_hi - 1
  }
  if (window_lo > window_hi) {
    return(NULL)
  }

  c(window_lo, window_hi)
}

# Internal: layers 2 and 3. Is there a sample of n integers in [0, span] with
# sum `sum_total` and sum of squares in [window_lo, window_hi]?
#
# Returns a list with `possible` (TRUE, FALSE, or NA when the node budget was
# reached before the tree was exhausted), `witness` (a non-increasing integer
# vector on the shifted scale when possible), and `nodes` (the search cost).
witness_search <- function(
  span,
  n,
  sum_total,
  window_lo,
  window_hi,
  budget = SEARCH_BUDGET
) {
  # Build from the nearer wall. Under y -> span - y the sum becomes
  # n*span - sum_total and the sum of squares shifts by
  # n*span^2 - 2*span*sum_total, so the window travels with it.
  flip <- 2 * sum_total > n * span
  if (flip) {
    shift <- n * span * span - 2 * span * sum_total
    sum_total <- n * span - sum_total
    window_lo <- window_lo + shift
    window_hi <- window_hi + shift
  }

  # The next value to place, given that `remaining` are left, none may exceed
  # `cap`, they must sum to `sum_left`, and their squares must sum into
  # [lo, hi]. Returns the admissible choices and a score for each -- how
  # centrally the surviving window sits inside that child's own sandwich. A
  # target hugging either end of what its subtree can reach is the one
  # likeliest to need backtracking, so the lowest score is the child to try
  # first. NULL when the state is already impossible.
  #
  # When remaining == 1 this is the leaf test: the single remaining value is
  # forced, and it survives the filter only if it lands the sum and the sum of
  # squares exactly. So an admissible choice there is a complete sample.
  #
  # This runs once per node and is the whole cost of the search, so it is
  # written against R's overheads rather than for brevity: no dplyr::if_else()
  # or ifelse(), which evaluate both arms over the full vector; no
  # pmax()/pmin(); and the Structure-S ceiling inlined rather than taken from
  # q_max_int(), whose recycling and mask-scatter are wasted when the cap is
  # already `values`.
  children <- function(remaining, cap, sum_left, lo, hi) {
    # The remaining - 1 values after this one are all <= value, so value must
    # be at least sum_left / remaining; and it cannot exceed the running cap or
    # the sum itself.
    value_hi <- min(cap, sum_left)
    # ceiling(sum_left / remaining), without leaving the integers
    value_lo <- (sum_left + remaining - 1) %/% remaining
    if (value_lo > value_hi) {
      return(NULL)
    }
    values <- value_lo:value_hi

    rest <- remaining - 1
    sum_rest <- sum_left - values
    rest_lo <- lo - values * values
    rest_hi <- hi - values * values
    if (rest == 0) {
      keep <- values[sum_rest == 0 & rest_lo <= 0 & rest_hi >= 0]
      return(
        if (length(keep)) {
          list(value = keep, score = numeric(length(keep)))
        }
      )
    }

    base_value <- sum_rest %/% rest
    remainder_count <- sum_rest - rest * base_value
    q_min <- (rest - remainder_count) *
      base_value *
      base_value +
      remainder_count * (base_value + 1) * (base_value + 1)
    # Structure-S on the remaining values, whose cap is the value just chosen.
    # A zero value can only arise when sum_rest is 0, where the division is
    # undefined and the true ceiling is 0.
    count_at_cap <- sum_rest %/% values
    remainder <- sum_rest - count_at_cap * values
    q_max <- count_at_cap * values * values + remainder * remainder
    q_max[!is.finite(q_max)] <- 0

    child_lo <- rest_lo
    tighter <- q_min > rest_lo
    child_lo[tighter] <- q_min[tighter]
    child_hi <- rest_hi
    tighter <- q_max < rest_hi
    child_hi[tighter] <- q_max[tighter]
    child_lo <- child_lo + ((child_lo %% 2) != (sum_rest %% 2))
    # `sum_rest >= 0` and `sum_rest <= rest * values` need no test: the first
    # is implied by value_hi <= sum_left, the second by
    # value_lo >= sum_left / remaining.
    admissible <- child_lo <= child_hi
    if (!any(admissible)) {
      return(NULL)
    }

    width <- q_max - q_min
    degenerate <- width == 0 # a sandwich of one value: a forced child
    width[degenerate] <- 1
    position <- ((child_lo + child_hi) / 2 - q_min) / width
    position[degenerate] <- 0.5

    list(value = values[admissible], score = abs(position - 0.5)[admissible])
  }

  # Depth-first, on an explicit stack rather than by recursion: the tree is n
  # deep, and a sample size in the thousands would otherwise exhaust R's
  # evaluation depth long before the node budget bit.
  #
  # A frame's children are ordered lazily. Only the best child is wanted on the
  # way down, and the way down is nearly always the whole search -- an
  # attainable report is typically reached in exactly n nodes, with no
  # backtracking at all. So the first visit takes which.min(), and a frame pays
  # for order() only if it is ever returned to. That matters because order()
  # inspects its arguments through match.arg() and two vapply() passes before
  # sorting anything, which costs several times the sort itself at these
  # lengths.
  frame_remaining <- numeric(n)
  frame_sum <- numeric(n)
  frame_lo <- numeric(n)
  frame_hi <- numeric(n)
  frame_at <- integer(n)
  frame_pick <- integer(n)
  frame_sorted <- logical(n)
  frame_children <- vector("list", n)
  frame_scores <- vector("list", n)
  chosen <- numeric(n)

  nodes <- 0
  capped <- FALSE
  hit <- NULL

  root <- children(n, span, sum_total, window_lo, window_hi)
  if (is.null(root)) {
    return(list(possible = FALSE, witness = NULL, nodes = nodes))
  }

  depth <- 1L
  frame_remaining[1] <- n
  frame_sum[1] <- sum_total
  frame_lo[1] <- window_lo
  frame_hi[1] <- window_hi
  frame_at[1] <- 0L
  frame_sorted[1] <- FALSE
  frame_children[[1]] <- root$value
  frame_scores[[1]] <- root$score

  while (depth > 0L) {
    if (nodes >= budget) {
      capped <- TRUE
      break
    }
    choices <- frame_children[[depth]]
    if (!frame_sorted[depth]) {
      # first visit: best child only
      pick <- which.min(frame_scores[[depth]])
      value <- choices[pick]
      frame_pick[depth] <- pick
      frame_sorted[depth] <- TRUE
      frame_at[depth] <- 0L
    } else if (frame_at[depth] == 0L) {
      # returned to: order the rest, once
      scores <- frame_scores[[depth]]
      scores[frame_pick[depth]] <- Inf # the child already tried sorts last
      frame_children[[depth]] <- choices[order(scores)][
        seq_len(length(choices) - 1L)
      ]
      frame_at[depth] <- 1L
      next
    } else {
      at <- frame_at[depth]
      if (at > length(choices)) {
        # this subtree is exhausted
        depth <- depth - 1L
        next
      }
      frame_at[depth] <- at + 1L
      value <- choices[at]
    }

    chosen[depth] <- value
    nodes <- nodes + 1
    if (frame_remaining[depth] == 1) {
      # every value placed, window met
      hit <- chosen
      break
    }

    next_remaining <- frame_remaining[depth] - 1
    next_sum <- frame_sum[depth] - value
    next_lo <- frame_lo[depth] - value * value
    next_hi <- frame_hi[depth] - value * value
    next_children <- children(
      next_remaining,
      value,
      next_sum,
      next_lo,
      next_hi
    )
    if (!is.null(next_children)) {
      depth <- depth + 1L
      frame_remaining[depth] <- next_remaining
      frame_sum[depth] <- next_sum
      frame_lo[depth] <- next_lo
      frame_hi[depth] <- next_hi
      frame_at[depth] <- 0L
      frame_sorted[depth] <- FALSE
      frame_pick[depth] <- 0L
      frame_children[[depth]] <- next_children$value
      frame_scores[[depth]] <- next_children$score
    }
  }

  if (!is.null(hit)) {
    return(list(
      possible = TRUE,
      witness = if (flip) span - hit else hit,
      nodes = nodes
    ))
  }

  list(possible = if (capped) NA else FALSE, witness = NULL, nodes = nodes)
}

# Internal: the fast route for one reported tuple under one rounding rule.
# Takes the (S, Q-window) targets target_states() already builds, so the
# endpoint semantics per rounding rule are shared rather than reimplemented.
#
# Returns TRUE / FALSE / NA, with NA meaning "no verdict within budget".
certify_fast <- function(span, n, targets, budget = SEARCH_BUDGET) {
  if (is.null(targets) || !nrow(targets)) {
    return(FALSE)
  }

  unknown <- FALSE
  for (row_current in seq_len(nrow(targets))) {
    window <- k_window(
      targets$S[row_current],
      targets$Q_lo[row_current],
      targets$Q_hi[row_current],
      targets$lo_incl[row_current],
      targets$hi_incl[row_current],
      n,
      span
    )
    if (is.null(window)) {
      next # layer 1 settles this candidate
    }
    found <- witness_search(
      span,
      n,
      targets$S[row_current],
      window[1],
      window[2],
      budget = budget
    )
    if (isTRUE(found$possible)) {
      return(TRUE)
    }
    if (is.na(found$possible)) {
      unknown <- TRUE
    }
  }

  if (unknown) NA else FALSE
}
