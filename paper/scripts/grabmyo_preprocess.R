#!/usr/bin/env Rscript

# Generate the processed GRABMyo tensor sample used in tensortools.Rmd.
#
# Fixed processing specification:
#   * GRABMyo v1.1.0, Session 1 only
#   * 43 participants (the independent draws)
#   * 16 active gestures
#   * 28 active EMG channels (the four U channels are removed)
#   * 7 trials averaged within each participant and gesture
#   * 16 equal 312.5 ms movement-phase bins per five-second trial
#   * log RMS amplitude in physical millivolts
#
# The output is a list of 43 tensors. Each participant tensor has dimensions
# 16 gestures x 28 channels x 16 movement-phase bins. Raw files and
# participant-level caches are kept outside the paper directory; the compact
# combined RDS used by the article is written to paper/data/GRABMyo.

options(timeout = 3600)

if (!requireNamespace("curl", quietly = TRUE)) {
  stop("Install curl first: install.packages('curl')")
}

# Locate the package root whether the script is run from the package root,
# paper/, or paper/scripts/.
grabmyo_find_project <- function() {
  candidates <- unique(normalizePath(
    c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "../..")),
    mustWork = FALSE
  ))
  is_project <- vapply(
    candidates,
    function(path) file.exists(file.path(path, "DESCRIPTION")),
    logical(1)
  )
  if (!any(is_project)) {
    stop("Run this script from the tensortools project, paper, or scripts directory.")
  }
  candidates[which(is_project)[1]]
}

grabmyo_project_dir <- grabmyo_find_project()

# Define variables used to read the data.
grabmyo_version <- "1.1.0"
grabmyo_session <- 1L
grabmyo_participants <- 1:43
grabmyo_gestures <- 1:16
grabmyo_trials <- 1:7
grabmyo_n_phase <- 16L
grabmyo_sampling_rate <- 2048L
grabmyo_n_samples <- 10240L # five seconds at 2048 Hz
grabmyo_n_recorded_channels <- 32L

# Large raw files remain outside the package project by default. Set
# `GRABMYO_DATA_DIR` to override this location on another machine.
grabmyo_data_dir <- file.path(grabmyo_project_dir, "data", "GRABMyo")
grabmyo_external_data_dir <- Sys.getenv(
  "GRABMYO_DATA_DIR",
  unset = file.path(dirname(grabmyo_project_dir), "GRABMyo-data")
)
grabmyo_raw_dir <- file.path(grabmyo_external_data_dir, "raw")
grabmyo_subject_dir <- file.path(
  grabmyo_data_dir, "subject_tensors_16x28x16"
)

# The processed sample is small enough to accompany the R Journal article.
grabmyo_article_data_dir <- file.path(
  grabmyo_project_dir, "paper", "data", "GRABMyo"
)
grabmyo_tensor_file <- file.path(
  grabmyo_article_data_dir,
  "grabmyo_session1_tensors_16x28x16.rds"
)

dir.create(grabmyo_raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(grabmyo_subject_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(grabmyo_article_data_dir, recursive = TRUE, showWarnings = FALSE)

# Public URLs for Session 1.
grabmyo_base_url <- paste0(
  "https://physionet-open.s3.amazonaws.com/grabmyo/",
  grabmyo_version, "/Session", grabmyo_session
)

grabmyo_record_paths <- function(participant, gesture, trial) {
  # Pull a specific participant, gesture, and trial.
  stem <- sprintf(
    "session%d_participant%d_gesture%d_trial%d",
    grabmyo_session, participant, gesture, trial
  )
  participant_dir <- paste0(
    "session", grabmyo_session, "_participant", participant
  )
  local_dir <- file.path(
    grabmyo_raw_dir,
    paste0("Session", grabmyo_session),
    participant_dir
  )

  list(
    hea = file.path(local_dir, paste0(stem, ".hea")),
    dat = file.path(local_dir, paste0(stem, ".dat")),
    hea_url = paste0(
      grabmyo_base_url, "/", participant_dir, "/", stem, ".hea"
    ),
    dat_url = paste0(
      grabmyo_base_url, "/", participant_dir, "/", stem, ".dat"
    )
  )
}

grabmyo_download_participant <- function(participant) {
  records <- expand.grid(
    gesture = grabmyo_gestures,
    trial = grabmyo_trials,
    KEEP.OUT.ATTRS = FALSE
  )

  # Each record has one small text header and one binary signal file.
  jobs <- unlist(lapply(seq_len(nrow(records)), function(i) {
    paths <- grabmyo_record_paths(
      participant,
      records$gesture[i],
      records$trial[i]
    )
    list(
      list(
        url = paths$hea_url,
        destination = paths$hea,
        expected_size = NA_real_
      ),
      list(
        url = paths$dat_url,
        destination = paths$dat,
        expected_size = grabmyo_n_samples *
          grabmyo_n_recorded_channels * 2L
      )
    )
  }), recursive = FALSE)

  file_is_complete <- function(job) {
    if (!file.exists(job$destination)) return(FALSE)
    size <- file.info(job$destination)$size
    size > 0 && (is.na(job$expected_size) || size == job$expected_size)
  }

  complete <- vapply(jobs, file_is_complete, logical(1))
  if (all(complete)) return(invisible(NULL))

  # Download every missing file in one concurrent batch. No GRABMyo files
  # failed in the analysis, so no dataset-specific skip list is required.
  pending <- jobs[!complete]
  destinations <- vapply(
    pending, `[[`, character(1), "destination"
  )
  invisible(lapply(
    unique(dirname(destinations)),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
  invisible(curl::multi_download(
    urls = vapply(pending, `[[`, character(1), "url"),
    destfiles = destinations,
    resume = FALSE,
    progress = FALSE,
    multi_timeout = 3600
  ))

  # Keep one integrity check so an interrupted binary file is never analyzed.
  complete <- vapply(jobs, file_is_complete, logical(1))
  if (!all(complete)) {
    stop("Download incomplete for participant ", participant)
  }
  invisible(NULL)
}

grabmyo_read_record <- function(header_file, data_file) {
  header <- readLines(header_file, warn = FALSE) # read header lines
  first_line <- strsplit(trimws(header[1]), "[[:space:]]+")[[1]]

  # Stop if any file has unexpected dimensions.
  if (length(first_line) < 4L ||
      as.integer(first_line[2]) != grabmyo_n_recorded_channels ||
      as.integer(first_line[3]) != grabmyo_sampling_rate ||
      as.integer(first_line[4]) != grabmyo_n_samples) {
    stop("Unexpected WFDB dimensions in ", header_file)
  }

  channel_fields <- strsplit(trimws(header[-1]), "[[:space:]]+")
  gain_fields <- vapply(channel_fields, `[[`, character(1), 3L)
  gains <- as.numeric(sub("\\(.*$", "", gain_fields))
  baselines <- as.numeric(sub(
    "^.*\\(([-+]?[0-9.]+)\\).*$", "\\1", gain_fields
  ))
  channel_names <- vapply(
    channel_fields,
    function(x) x[length(x)],
    character(1)
  )

  connection <- file(data_file, open = "rb")
  on.exit(close(connection))
  adc <- readBin( # read measurement data
    connection,
    what = integer(),
    n = grabmyo_n_samples * grabmyo_n_recorded_channels,
    size = 2L,
    signed = TRUE,
    endian = "little"
  )
  if (length(adc) != grabmyo_n_samples * grabmyo_n_recorded_channels) {
    stop("Unexpected data length in ", data_file)
  }

  # The binary file stores all channels for time point 1, then all channels
  # for time point 2, and so on. Reshape it into time x channel.
  adc <- matrix(
    adc,
    nrow = grabmyo_n_samples,
    ncol = grabmyo_n_recorded_channels,
    byrow = TRUE
  )

  # Convert device ADC counts to millivolts using the baseline and gain stored
  # for each channel in the WFDB header.
  physical_mv <- sweep(adc, 2L, baselines, "-")
  physical_mv <- sweep(physical_mv, 2L, gains, "/")

  active <- !startsWith(channel_names, "U") # remove unused channels
  physical_mv <- physical_mv[, active, drop = FALSE]
  colnames(physical_mv) <- channel_names[active]
  if (ncol(physical_mv) != 28L) {
    stop("Expected 28 active channels in ", header_file)
  }

  physical_mv
}

grabmyo_trial_log_rms <- function(signal, n_phase = grabmyo_n_phase) {
  if (nrow(signal) %% n_phase != 0L) {
    stop("Trial length is not divisible by the number of phase bins.")
  }

  samples_per_phase <- nrow(signal) %/% n_phase # 640 samples = 312.5 ms
  phase <- rep(seq_len(n_phase), each = samples_per_phase)

  # RMS measures the signal strength even though EMG oscillates above and below
  # zero. vapply() returns a channel x phase matrix. The log transformation
  # makes the positive, typically right-skewed RMS amplitudes more symmetric.
  rms <- vapply(
    split(seq_len(nrow(signal)), phase),
    function(index) {
      sqrt(colMeans(signal[index, , drop = FALSE]^2))
    },
    numeric(ncol(signal))
  )
  log(pmax(rms, .Machine$double.xmin))
}

grabmyo_make_participant_tensor <- function(participant) {
  participant_tensor <- array(
    0,
    dim = c(16L, 28L, grabmyo_n_phase),
    dimnames = list(
      gesture = paste0("gesture", grabmyo_gestures),
      channel = c(paste0("F", 1:16), paste0("W", 1:12)),
      phase = paste0("phase", seq_len(grabmyo_n_phase))
    )
  )

  for (gesture in grabmyo_gestures) {
    trial_features <- vector("list", length(grabmyo_trials))

    for (trial in grabmyo_trials) {
      paths <- grabmyo_record_paths(participant, gesture, trial)
      signal <- grabmyo_read_record(paths$hea, paths$dat)
      trial_features[[trial]] <- grabmyo_trial_log_rms(signal)
    }

    # Average seven repeated trials only after applying the same deterministic
    # feature extraction to each trial. Repeats are not IID model inputs.
    participant_tensor[gesture, , ] <-
      Reduce(`+`, trial_features) / length(trial_features)
  }

  participant_tensor
}

grabmyo_format_time <- function(seconds) {
  seconds <- max(0, round(seconds))
  sprintf(
    "%02d:%02d:%02d",
    seconds %/% 3600,
    (seconds %% 3600) %/% 60,
    seconds %% 60
  )
}

# Create or load one cache file per participant. Saving at this granularity
# ensures that a network interruption does not discard processed participants.
grabmyo_tensors <- vector("list", length(grabmyo_participants))
names(grabmyo_tensors) <- paste0("participant", grabmyo_participants)
grabmyo_started <- Sys.time()

for (i in seq_along(grabmyo_participants)) {
  participant <- grabmyo_participants[i]
  participant_cache <- file.path(
    grabmyo_subject_dir,
    paste0("participant", participant, ".rds")
  )

  if (file.exists(participant_cache)) {
    grabmyo_tensors[[i]] <- readRDS(participant_cache)
    action <- "loaded from cache"
  } else {
    grabmyo_download_participant(participant)
    grabmyo_tensors[[i]] <- grabmyo_make_participant_tensor(participant)
    saveRDS(grabmyo_tensors[[i]], participant_cache)
    action <- "downloaded and processed"
  }

  # Report elapsed time and ETA after every participant.
  elapsed <- as.numeric(difftime(
    Sys.time(), grabmyo_started, units = "secs"
  ))
  eta <- elapsed / i * (length(grabmyo_participants) - i)
  message(sprintf(
    "GRABMyo %d/%d: participant %d %s | elapsed %s | ETA %s",
    i, length(grabmyo_participants), participant, action,
    grabmyo_format_time(elapsed), grabmyo_format_time(eta)
  ))
}

# Save the compact processed sample used by the R Journal article. Every list
# element is one independent 16 x 28 x 16 tensor-valued draw.
saveRDS(grabmyo_tensors, grabmyo_tensor_file)
stopifnot(
  length(grabmyo_tensors) == 43L,
  identical(dim(grabmyo_tensors[[1]]), c(16L, 28L, 16L)),
  all(vapply(grabmyo_tensors, function(x) {
    all(is.finite(x))
  }, logical(1)))
)

message("Processed GRABMyo tensors saved to: ", grabmyo_tensor_file)
