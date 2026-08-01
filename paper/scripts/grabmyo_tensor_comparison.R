#!/usr/bin/env Rscript

# GRABMyo tensor-normal versus every matrix-normal unfolding.
#
# Fixed analysis specification:
#   * GRABMyo v1.1.0, Session 1 only
#   * 43 participants (independent draws)
#   * 16 active gestures
#   * 28 documented active EMG channels (the four U channels are removed)
#   * 7 trials averaged within participant and gesture
#   * 16 equal 312.5 ms movement-phase bins over each five-second trial
#   * log RMS amplitude in physical mV
#
# Each participant produces a 16 gesture x 28 channel x 16 phase tensor.

options(timeout = 3600)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload first: install.packages('pkgload')")
}
if (!requireNamespace("curl", quietly = TRUE)) {
  stop("Install curl first: install.packages('curl')")
}

project_dir <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(project_dir, "DESCRIPTION"))) {
  stop("Run this script from the tensortools package root.")
}
pkgload::load_all(project_dir, quiet = TRUE)

version <- "1.1.0"
session <- 1L
participants <- seq_len(as.integer(Sys.getenv("GRABMYO_N_SUBJECTS", "43")))
gestures <- 1:16
trials <- 1:7
n_phase <- 16L
n_channels_recorded <- 32L
n_samples <- 10240L
sampling_rate <- 2048
process_only <- identical(Sys.getenv("GRABMYO_PROCESS_ONLY", "0"), "1")

data_dir <- file.path(project_dir, "data", "GRABMyo")
raw_dir <- file.path(data_dir, "raw")
subject_tensor_dir <- file.path(data_dir, "subject_tensors_16x28x16")
tensor_file <- file.path(data_dir, "grabmyo_session1_tensors_16x28x16.rds")
result_file <- file.path(data_dir, "grabmyo_tensor_matrix_results.rds")
result_csv <- file.path(data_dir, "grabmyo_tensor_matrix_results.csv")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(subject_tensor_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- paste0(
  "https://physionet-open.s3.amazonaws.com/grabmyo/", version,
  "/Session", session
)

record_stem <- function(participant, gesture, trial) {
  sprintf(
    "session%d_participant%d_gesture%d_trial%d",
    session, participant, gesture, trial
  )
}

record_paths <- function(participant, gesture, trial) {
  stem <- record_stem(participant, gesture, trial)
  participant_dir <- paste0("session", session, "_participant", participant)
  local_dir <- file.path(raw_dir, paste0("Session", session), participant_dir)

  list(
    stem = stem,
    local_dir = local_dir,
    hea = file.path(local_dir, paste0(stem, ".hea")),
    dat = file.path(local_dir, paste0(stem, ".dat")),
    hea_url = paste0(base_url, "/", participant_dir, "/", stem, ".hea"),
    dat_url = paste0(base_url, "/", participant_dir, "/", stem, ".dat")
  )
}

download_one <- function(url, destination, expected_size = NULL,
                         attempts = 4L) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)

  complete <- file.exists(destination) && file.info(destination)$size > 0
  if (!is.null(expected_size)) {
    complete <- complete && file.info(destination)$size == expected_size
  }
  if (complete) return(destination)

  partial <- paste0(destination, ".part")
  for (attempt in seq_len(attempts)) {
    status <- tryCatch(
      suppressWarnings(utils::download.file(
        url,
        destfile = partial,
        mode = "wb",
        method = "libcurl",
        quiet = TRUE
      )),
      error = function(e) 1L
    )

    size_ok <- file.exists(partial) && file.info(partial)$size > 0
    if (!is.null(expected_size)) {
      size_ok <- size_ok && file.info(partial)$size == expected_size
    }

    if (identical(status, 0L) && size_ok) {
      if (file.exists(destination)) unlink(destination)
      if (!file.rename(partial, destination)) {
        stop("Could not finalize ", destination)
      }
      return(destination)
    }
  }

  stop("Download failed: ", url)
}

download_participant <- function(participant) {
  records <- expand.grid(
    gesture = gestures,
    trial = trials,
    KEEP.OUT.ATTRS = FALSE
  )

  jobs <- unlist(lapply(seq_len(nrow(records)), function(i) {
    paths <- record_paths(
      participant,
      records$gesture[i],
      records$trial[i]
    )
    list(
      list(url = paths$hea_url, destination = paths$hea,
           expected_size = NULL),
      list(url = paths$dat_url, destination = paths$dat,
           expected_size = n_samples * n_channels_recorded * 2L)
    )
  }), recursive = FALSE)

  for (attempt in 1:4) {
    complete <- vapply(jobs, function(job) {
      ok <- file.exists(job$destination) &&
        file.info(job$destination)$size > 0
      if (!is.null(job$expected_size)) {
        ok <- ok && file.info(job$destination)$size == job$expected_size
      }
      ok
    }, logical(1))

    if (all(complete)) return(invisible(NULL))

    pending <- jobs[!complete]
    destinations <- vapply(pending, `[[`, character(1), "destination")
    invisible(lapply(unique(dirname(destinations)), dir.create,
                     recursive = TRUE, showWarnings = FALSE))

    result <- curl::multi_download(
      urls = vapply(pending, `[[`, character(1), "url"),
      destfiles = destinations,
      resume = FALSE,
      progress = FALSE,
      multi_timeout = 3600
    )

    bad_http <- is.na(result$status_code) |
      result$status_code < 200L | result$status_code >= 300L
    if (any(bad_http)) {
      unlink(result$destfile[bad_http])
    }
  }

  stop("One or more downloads remain incomplete for participant ", participant)
}

read_grabmyo_record <- function(header_file, data_file) {
  header <- readLines(header_file, warn = FALSE)
  first <- strsplit(trimws(header[1]), "[[:space:]]+")[[1]]

  if (length(first) < 4L || as.integer(first[2]) != n_channels_recorded ||
      as.integer(first[3]) != sampling_rate ||
      as.integer(first[4]) != n_samples) {
    stop("Unexpected WFDB dimensions in ", header_file)
  }

  fields <- strsplit(trimws(header[-1]), "[[:space:]]+")
  gain_fields <- vapply(fields, `[[`, character(1), 3L)
  gains <- as.numeric(sub("\\(.*$", "", gain_fields))
  baselines <- as.numeric(sub(
    "^.*\\(([-+]?[0-9.]+)\\).*$",
    "\\1",
    gain_fields
  ))
  channel_names <- vapply(fields, function(x) x[length(x)], character(1))

  connection <- file(data_file, open = "rb")
  on.exit(close(connection))
  adc <- readBin(
    connection,
    what = integer(),
    n = n_samples * n_channels_recorded,
    size = 2L,
    signed = TRUE,
    endian = "little"
  )
  if (length(adc) != n_samples * n_channels_recorded) {
    stop("Unexpected data length in ", data_file)
  }

  adc <- matrix(adc, nrow = n_samples, ncol = n_channels_recorded,
                byrow = TRUE)
  physical_mv <- sweep(adc, 2L, baselines, "-")
  physical_mv <- sweep(physical_mv, 2L, gains, "/")

  active <- !startsWith(channel_names, "U")
  physical_mv <- physical_mv[, active, drop = FALSE]
  colnames(physical_mv) <- channel_names[active]

  if (ncol(physical_mv) != 28L) {
    stop("Expected 28 active channels in ", header_file)
  }
  physical_mv
}

trial_log_rms <- function(signal, n_phase = 16L) {
  if (nrow(signal) %% n_phase != 0L) {
    stop("The trial length is not divisible by the number of phase bins.")
  }

  samples_per_phase <- nrow(signal) %/% n_phase
  phase <- rep(seq_len(n_phase), each = samples_per_phase)
  rms <- vapply(
    split(seq_len(nrow(signal)), phase),
    function(index) sqrt(colMeans(signal[index, , drop = FALSE]^2)),
    numeric(ncol(signal))
  )
  log(pmax(rms, .Machine$double.xmin))
}

make_participant_tensor <- function(participant) {
  output <- array(
    0,
    dim = c(length(gestures), 28L, n_phase),
    dimnames = list(
      gesture = paste0("gesture", gestures),
      channel = c(paste0("F", 1:16), paste0("W", 1:12)),
      phase = paste0("phase", seq_len(n_phase))
    )
  )

  for (gesture in gestures) {
    trial_features <- vector("list", length(trials))
    for (trial in trials) {
      paths <- record_paths(participant, gesture, trial)
      signal <- read_grabmyo_record(paths$hea, paths$dat)
      trial_features[[trial]] <- trial_log_rms(signal, n_phase)
    }
    output[gesture, , ] <- Reduce(`+`, trial_features) / length(trials)
  }
  output
}

format_time <- function(seconds) {
  seconds <- max(0, round(seconds))
  sprintf(
    "%02d:%02d:%02d",
    seconds %/% 3600,
    (seconds %% 3600) %/% 60,
    seconds %% 60
  )
}

participant_tensors <- vector("list", length(participants))
names(participant_tensors) <- paste0("participant", participants)
started <- Sys.time()

for (i in seq_along(participants)) {
  participant <- participants[i]
  cache_file <- file.path(
    subject_tensor_dir,
    paste0("participant", participant, ".rds")
  )

  if (file.exists(cache_file)) {
    participant_tensors[[i]] <- readRDS(cache_file)
    action <- "cached"
  } else {
    download_participant(participant)
    participant_tensors[[i]] <- make_participant_tensor(participant)
    saveRDS(participant_tensors[[i]], cache_file)
    action <- "processed"
  }

  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  eta <- elapsed / i * (length(participants) - i)
  message(sprintf(
    "GRABMyo %d/%d: participant %d %s | elapsed %s | ETA %s",
    i, length(participants), participant, action,
    format_time(elapsed), format_time(eta)
  ))
}

saveRDS(participant_tensors, tensor_file)
stopifnot(
  length(participant_tensors) == length(participants),
  identical(dim(participant_tensors[[1]]), c(16L, 28L, 16L)),
  all(vapply(participant_tensors, function(x) all(is.finite(x)), logical(1)))
)

if (process_only) {
  message("Processing complete; GRABMYO_PROCESS_ONLY=1, so model fitting is skipped.")
  quit(save = "no", status = 0L)
}

fit_one <- function(label, draws) {
  cache_file <- file.path(data_dir, paste0("fit_logdet_v2_", label, ".rds"))
  if (file.exists(cache_file)) {
    message("Loading cached model: ", label)
    return(readRDS(cache_file))
  }

  message("Fitting model: ", label, " with dimensions ",
          paste(dim(draws[[1]]), collapse = " x "))
  started <- Sys.time()
  fit <- tensor_mle(
    draws,
    model = "normal",
    max_iter = 1000,
    tol = 1e-6,
    quiet = FALSE
  )
  fit$elapsed_seconds <- as.numeric(
    difftime(Sys.time(), started, units = "secs")
  )
  saveRDS(fit, cache_file)
  fit
}

fits <- list()
fits$tensor <- fit_one("tensor", participant_tensors)

for (mode in 1:3) {
  matrix_draws <- lapply(
    participant_tensors,
    matricization,
    k = mode
  )
  fits[[paste0("mat", mode)]] <- fit_one(
    paste0("mat", mode),
    matrix_draws
  )
  rm(matrix_draws)
  gc(FALSE)
}

results <- data.frame(
  model = names(fits),
  dimensions = vapply(
    c(list(participant_tensors), lapply(1:3, function(mode) {
      lapply(participant_tensors[1], matricization, k = mode)
    })),
    function(draws) paste(dim(draws[[1]]), collapse = " x "),
    character(1)
  ),
  loglik = vapply(fits, `[[`, numeric(1), "loglik"),
  k = vapply(fits, `[[`, numeric(1), "k"),
  AIC = vapply(fits, `[[`, numeric(1), "AIC"),
  BIC = vapply(fits, `[[`, numeric(1), "BIC"),
  elapsed_seconds = vapply(fits, `[[`, numeric(1), "elapsed_seconds"),
  row.names = NULL
)
results <- results[order(results$BIC), ]
results$delta_BIC <- results$BIC - min(results$BIC)

saveRDS(list(results = results, fits = fits), result_file)
utils::write.csv(results, result_csv, row.names = FALSE)

print(results, digits = 12, row.names = FALSE)
message("Best model by BIC: ", results$model[1])
message("Results saved to: ", result_csv)
