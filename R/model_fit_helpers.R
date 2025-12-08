# # ========================================================================= #
# copy_svm2 -------
# # ========================================================================= #

#' Copy the `svm2` methodInfo to the sdm package `methods/sdm` directory
#'
#' This function checks if the file `svm2.R` exists in the `methods/sdm`
#' subdirectory of the installed sdm package. If not, it writes a new file
#' containing the `methodInfo` list for registering the `svm2` method (using
#' `e1071` SVM) in sdm.
#'
#' @details The `methodInfo` list defines the `svm2` modelling method for sdm,
#'   using the `e1071::svm` implementation.
#' @return Returns `invisible(NULL)`. The function is called for its side effect
#'   of writing a file.
#' @author Ahmed El-Gabbas
#' @noRd
#' @keywords internal

copy_svm2 <- function() {

  # Check that sdm R package is installed
  ecokit::check_packages(c("sdm", "e1071"))

  # Target path in the sdm package
  target_dir <- system.file(fs::path("methods", "sdm"), package = "sdm")

  if (target_dir == "") {
    ecokit::stop_ctx(
      paste0(
        "Could not find 'methods/sdm' directory in the 'sdm' package. ", #nolint
        "The package structure may be unexpected or the package may not be ",
        " installed correctly."))
  }

  target_file <- fs::path(target_dir, "svm2.R")

  # Check if file already exists
  if (fs::file_exists(target_file)) {
    return(invisible(NULL))
  }

  # Check if target directory is writeable
  if (file.access(target_dir, 2L) != 0L) {
    ecokit::stop_ctx(
      paste0(
        "Cannot write to the 'sdm' package methods/sdm directory. ", #nolint
        "You do not have write permissions for this directory."))
  }

  # MethodInfo function text to write
  method_text <-
    'methodInfo <- list(
  name = c("svm2", "SVM2", "svm_e1071"),
  packages = "e1071",
  modelTypes = c("pa", "pb", "ab", "n"),
  fitParams = list(
    formula = "standard.formula", data = "sdmDataFrame", v = "sdmVariables"),
  fitSettings = list(kernel = "radial", probability = TRUE),
  fitFunction = function(formula, data, v, ...) {
    x <- sdm:::.getData.sdmMatrix(
      formula, data, normalize = TRUE, frame = v@varInfo$numeric, scale = FALSE)
    y <- sdm:::.getData.sdmY(formula, data)

    # Set class weights for pa/pb (binary response)
    n0 <- sum(y == 0, na.rm = TRUE)
    n1 <- sum(y == 1, na.rm = TRUE)
    max_weight <- 20 # Upper bound for weight
    if (n0 >= n1) {
      # More absences
      class.weights <- c("0" = 1, "1" = min(n0 / n1, max_weight))
    } else {
      # More presences
      class.weights <- c("0" = min(n1 / n0, max_weight), "1" = 1)
    }

    e1071::svm(x = x, y = y, scale = TRUE, class.weights = class.weights, ...)
  },
  settingRules = NULL,
  tuneParams = NULL,
  predictParams = list(
    object = "model", formula = "standard.formula", newx = "sdmDataFrame",
    v = "sdmVariables"),
  predictSettings = list(probability = TRUE),
  predictFunction = function(object, formula, newx, v, ...) {
    newx <- sdm:::.getData.sdmMatrix(
      formula, newx, normalize = TRUE,
      frame = v@varInfo$numeric, scale = FALSE)
    predict(object, newx, ...)
  },

  #------ metadata (optional):
  title = "Support Vector Machines using e1071",
  creator = "Ahmed El-Gabbas"
  )'

  # Write to file, preserving formatting
  writeLines(method_text, target_file)

  invisible(NULL)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# # ========================================================================= #
# reduce_sdm_formulas ------
# # ========================================================================= #

#' Reduce SDM Model Formulas to Base Environment To Reduce File Size
#'
#' This function iterates through the models stored in a `sdm` object, and
#' ensures that all formulas and their associated environments are set to the
#' base environment. This helps to reduce the fitted model's object size, see
#' [here](https://github.com/babaknaimi/sdm/issues/43).
#'
#' @param obj An `sdmModels` object containing fitted SDM models.
#' @return A copy of the input object with all model formulas and terms
#'   environments set to the base environment.
#' @author Ahmed El-Gabbas
#' @noRd
#' @keywords internal

reduce_sdm_formulas <- function(obj) {

  if (is.null(obj) || !inherits(obj, "sdmModels")) {
    ecokit::stop_ctx(
      "Input object must be of class 'sdmModels'.",
      input_object = obj, input_class = class(obj))
  }

  # Make a copy to avoid modifying the original
  obj_copy <- obj

  # Loop through each model using the structure of the SDM object
  for (sp_name in names(obj_copy@models)) {
    for (method_name in names(obj_copy@models[[sp_name]])) {
      for (run_idx in seq_along(obj_copy@models[[sp_name]][[method_name]])) {

        # Check if there's a formula in the call
        model <- obj_copy@models[[sp_name]][[method_name]][[run_idx]]

        if (
          !isS4(model@object) &&
          !is.null(model@object$call$formula) &&
          deparse1(model@object$call$formula) != ".f") {

          # Replace formula in the model call (if it exists and is not ".f")
          formula_text <- deparse1(model@object$call$formula)
          clean_formula <- stats::as.formula(formula_text, env = baseenv())
          model@object$call$formula <- clean_formula

          if (!is.null(model@object$formula)) {
            # Replace formula with clean version
            formula_text <- deparse1(model@object$formula)
            clean_formula <- stats::as.formula(formula_text, env = baseenv())
            model@object$formula <- clean_formula
          }

          # Fix terms environment if it exists
          if (!is.null(model@object$terms)) {
            attr(model@object$terms, ".Environment") <- baseenv()
          }
          # Update the model in the object
          obj_copy@models[[sp_name]][[method_name]][[run_idx]] <- model
        }

      }
    }
  }

  # Force garbage collection
  invisible(gc())

  obj_copy
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# # ========================================================================= #
# extract_sdm_info ------
# # ========================================================================= #

#' Extract SDM Model Information and Evaluation Metrics
#'
#' This function extracts key information, evaluation metrics, variable
#' importance, and response curves from a fitted species distribution model
#' (SDM) object of class `sdmModels`, or from a file path to a saved model. The
#' function works with models fitted using the `sdm` package and provides a tidy
#' summary of model performance for both training and independent test data, as
#' well as variable importance and response curves.
#'
#' @param model An object of class `sdmModels` (fitted using `sdm` package) or a
#'   character path to a saved model file. Required.
#'
#' @return A tibble with the following columns:
#' - `evaluation_training`: A tibble with evaluation metrics for the training
#'   data, including AUC, TSS, Kappa, prevalence, sensitivity, specificity, and
#'   other statistics.
#' - `evaluation_testing`: A tibble with evaluation metrics for the independent
#'   test data, with the same structure as `evaluation_training`.
#' - `variable_importance`: A tibble summarizing the importance of each
#'   predictor variable, including correlation and AUC metrics for the test
#'   data.
#' - `response_curves`: A tibble containing response curves for each predictor
#'   variable, showing the predicted response across the range of values.
#'
#'   If evaluation or variable importance data are missing, the function returns
#'   tibbles with `NA` values for the corresponding metrics. The function also
#'   handles errors gracefully and provides informative messages if extraction
#'   fails.
#'
#' @author Ahmed El-Gabbas
#' @noRd
#' @keywords internal

extract_sdm_info <- function(model = NULL) {

  AUCtest <- corTest <- criteria <- variables <- . <- NULL #nolint

  if (is.null(model)) {
    ecokit::stop_ctx("`model` can not be empty")
  }

  # Validate and load model
  model_fail <- FALSE

  if (!inherits(model, "sdmModels")) {
    if (inherits(model, "character")) {
      if (ecokit::check_data(model, warning = FALSE)) {
        model <- ecokit::load_as(model)
      } else {
        model_fail <- TRUE
      }
    } else {
      model_fail <- TRUE
    }
  }

  if (model_fail) {
    ecokit::stop_ctx(
      "model must be a sdmModels object or a file path",
      model = model, class_model = class(model))
  }
  if (!methods::.hasSlot(model, "setting")) {
    ecokit::stop_ctx(
      "model does not have a 'setting' slot.",
      model = model, class_model = class(model))
  }
  if (!methods::.hasSlot(model@setting, "featureFrame")) {
    ecokit::stop_ctx(
      "model does not have a 'setting@featureFrame' slot.",
      model = model, class_model = class(model))
  }
  if (!methods::.hasSlot(model@setting@featureFrame, "predictors")) {
    ecokit::stop_ctx(
      "model does not have a 'setting@featureFrame@predictors' slot.",
      model = model, class_model = class(model))
  }
  if (!methods::.hasSlot(model@setting, "methods")) {
    ecokit::stop_ctx(
      "model does not have a 'setting@methods' slot.",
      model = model, class_model = class(model))
  }
  if (!methods::.hasSlot(model, "models")) {
    ecokit::stop_ctx(
      "model does not have a 'models' slot.",
      model = model, class_model = class(model))
  }

  # Ensure predictor_names is a character vector
  predictor_names <- as.character(model@setting@featureFrame@predictors)
  if (!is.character(predictor_names)) {
    ecokit::stop_ctx(
      "'model@setting@featureFrame@predictors' must be a character vector.",
      predictors = predictor_names, class_predictors = class(predictor_names))
  }

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Extract model metadata
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  species_model <- tryCatch(
    model@models[[1L]][[1L]][[1L]],
    error = function(e) NULL)

  if (is.null(species_model)) {
    ecokit::stop_ctx("No fitted model found in sdmModels object.")
  }

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Evaluation ----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # columns to keep in evaluation data
  evaluation_sort <- c(
    "TSS", "Kappa", "threshold", "prevalence", "sensitivity", "specificity")

  # Create empty evaluation tibble to handle cases where evaluation is missing
  # or not available
  empty_evaluation_train <- tibble::tibble(
    prevalence_train = NA_real_, auc_train = NA_real_,
    boyce_train = NA_real_, cor_train = NA_real_,
    cor_p_train = NA_real_, deviance_train = NA_real_,
    tss_ess_train = NA_real_, kappa_ess_train = NA_real_,
    threshold_ess_train = NA_real_, prevalence_ess_train = NA_real_,
    sensitivity_ess_train = NA_real_, specificity_ess_train = NA_real_,
    tss_mss_train = NA_real_, kappa_mss_train = NA_real_,
    threshold_mss_train = NA_real_, prevalence_mss_train = NA_real_,
    sensitivity_mss_train = NA_real_, specificity_mss_train = NA_real_)

  empty_evaluation_test <- empty_evaluation_train %>%
    dplyr::rename_with(~ stringr::str_replace(., "train", "test"))

  # Ensure that `evaluation` slot exists
  if (methods::.hasSlot(species_model, "evaluation")) {

    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    ## Evaluation - training ------
    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    if ("training" %in% names(species_model@evaluation)) {

      eval_train <- species_model@evaluation$training

      # `Statistics` slot
      if (methods::.hasSlot(eval_train, "statistics")) {
        eval_train_stats <- eval_train@statistics
        if ("cBoyce" %in% names(eval_train_stats)) {
          cont_boyce_train <- eval_train_stats$cBoyce
        } else {
          cont_boyce_train <- NA_real_
        }
        eval_train_stats <- tibble::tibble(
          prevalence_train = eval_train_stats$Prevalence,
          auc_train = eval_train_stats$AUC,
          boyce_train = cont_boyce_train,
          cor_train = eval_train_stats$COR[1L],
          cor_p_train = eval_train_stats$COR[2L],
          deviance_train = eval_train_stats$Deviance)
      } else {
        # Return empty tibble if no `statistics` slot exists
        eval_train_stats <- tibble::tibble(
          prevalence_train = NA_real_, auc_train = NA_real_,
          boyce_train = NA_real_, cor_train = NA_real_, cor_p_train = NA_real_,
          deviance_train = NA_real_)
      }

      # `threshold_based` slot
      if (methods::.hasSlot(eval_train, "threshold_based")) {
        eval_train_thr <- eval_train@threshold_based
        eval_train_thr_ess <- eval_train_thr %>%
          dplyr::filter(criteria == "sp=se") %>%
          dplyr::select(tidyselect::all_of(evaluation_sort)) %>%
          dplyr::rename_with(~ stringr::str_c(tolower(.), "_ess_train"))
        eval_train_thr_mss <- eval_train_thr %>%
          dplyr::filter(criteria == "max(se+sp)") %>%
          dplyr::select(tidyselect::all_of(evaluation_sort)) %>%
          dplyr::rename_with(~ stringr::str_c(tolower(.), "_mss_train"))
      } else {
        eval_train_thr_ess <- tibble::tibble(
          tss_ess_train = NA_real_, kappa_ess_train = NA_real_,
          threshold_ess_train = NA_real_, prevalence_ess_train = NA_real_,
          sensitivity_ess_train = NA_real_, specificity_ess_train = NA_real_)
        eval_train_thr_mss <- tibble::tibble(
          tss_mss_train = NA_real_, kappa_mss_train = NA_real_,
          threshold_mss_train = NA_real_, prevalence_mss_train = NA_real_,
          sensitivity_mss_train = NA_real_, specificity_mss_train = NA_real_)
      }

      evaluation_training <- dplyr::bind_cols(
        eval_train_stats, eval_train_thr_ess, eval_train_thr_mss)

    } else {
      # no training evaluation data, return empty tibble
      evaluation_training <- empty_evaluation_train
    }

    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Evaluation - testing -----
    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    if ("test.indep" %in% names(species_model@evaluation)) {

      eval_test <- species_model@evaluation$test.indep

      # `Statistics` slot
      if (methods::.hasSlot(eval_test, "statistics")) {
        eval_test_stats <- eval_test@statistics
        if ("cBoyce" %in% names(eval_test_stats)) {
          cont_boyce_test <- eval_test_stats$cBoyce
        } else {
          cont_boyce_test <- NA_real_
        }
        eval_test_stats <- tibble::tibble(
          prevalence_test = eval_test_stats$Prevalence,
          auc_test = eval_test_stats$AUC,
          boyce_test = cont_boyce_test,
          cor_test = eval_test_stats$COR[1L],
          cor_p_test = eval_test_stats$COR[2L],
          deviance_test = eval_test_stats$Deviance)
      } else {
        # Return empty tibble if no `statistics` slot exists
        eval_test_stats <- tibble::tibble(
          prevalence_test = NA_real_, auc_test = NA_real_,
          boyce_test = NA_real_, cor_test = NA_real_, cor_p_test = NA_real_,
          deviance_test = NA_real_)
      }

      # `threshold_based` slot
      if (methods::.hasSlot(eval_test, "threshold_based")) {
        eval_test_thr <- eval_test@threshold_based
        eval_test_thr_ess <- eval_test_thr %>%
          dplyr::filter(criteria == "sp=se") %>%
          dplyr::select(tidyselect::all_of(evaluation_sort)) %>%
          dplyr::rename_with(~ stringr::str_c(tolower(.), "_ess_test"))
        eval_test_thr_mss <- eval_test_thr %>%
          dplyr::filter(criteria == "max(se+sp)") %>%
          dplyr::select(tidyselect::all_of(evaluation_sort)) %>%
          dplyr::rename_with(~ stringr::str_c(tolower(.), "_mss_test"))
      } else {
        eval_test_thr_ess <- tibble::tibble(
          tss_ess_test = NA_real_, kappa_ess_test = NA_real_,
          threshold_ess_test = NA_real_, prevalence_ess_test = NA_real_,
          sensitivity_ess_test = NA_real_, specificity_ess_test = NA_real_)
        eval_test_thr_mss <- tibble::tibble(
          tss_mss_test = NA_real_, kappa_mss_test = NA_real_,
          threshold_mss_test = NA_real_, prevalence_mss_test = NA_real_,
          sensitivity_mss_test = NA_real_, specificity_mss_test = NA_real_)
      }

      evaluation_testing <- dplyr::bind_cols(
        eval_test_stats, eval_test_thr_ess, eval_test_thr_mss)

    } else {
      # no testing evaluation data, return empty tibble
      evaluation_testing <- empty_evaluation_test
    }

  } else {
    evaluation_training <- empty_evaluation_train
    evaluation_testing <- empty_evaluation_test
  }

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Variable importance -----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  empty_var_imp <- dplyr::bind_cols(
    variable = predictor_names, cor_test = NA_real_, auc_test = NA_real_)

  if (methods::.hasSlot(species_model, "varImportance")) {

    variable_importance <- species_model@varImportance

    if ("test.indep" %in% names(variable_importance)) {
      variable_importance <- variable_importance$test.indep

      if (methods::.hasSlot(variable_importance, "varImportance")) {
        variable_importance <- variable_importance@varImportance

        if (is.null(variable_importance) || length(variable_importance) == 0L) {
          variable_importance <- empty_var_imp
        } else {
          variable_importance <- tibble::tibble(variable_importance) %>%
            dplyr::rename(
              variable = variables, cor_test = corTest, auc_test = AUCtest)
        }
      } else {
        variable_importance <- empty_var_imp
      }
    } else {
      variable_importance <- empty_var_imp
    }
  } else {
    variable_importance <- empty_var_imp
  }

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Response curves -----
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Extract response curves
  r_curves <- sdm::getResponseCurve(x = model, id = 1L)

  empty_r_curves <- tibble::tibble(
    variable = predictor_names, x_value = NA_real_, prediction = NA_real_)

  if (methods::.hasSlot(r_curves, "variables") &&
      methods::.hasSlot(r_curves, "response")) {
    if (length(r_curves@variables) == 0L || length(r_curves@response) == 0L) {
      r_curves <- empty_r_curves
    } else {
      r_curves <- purrr::map_dfr(
        .x = r_curves@variables,
        .f = ~{
          r_curves@response[[.x]] %>%
            tibble::tibble() %>%
            stats::setNames(c("x_value", "prediction")) %>%
            dplyr::mutate(variable = .x, .before = 1L)
        })
    }
  } else {
    r_curves <- empty_r_curves
  }

  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # return outputs
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  tibble::tibble(
    evaluation_training = list(evaluation_training),
    evaluation_testing = list(evaluation_testing),
    variable_importance = list(variable_importance),
    response_curves = list(r_curves))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# # ========================================================================= #
# prep_model_data -------
# # ========================================================================= #

prep_model_data <- function(
    cv_fold, dir_model_data,
    model_data, abs_ratio, bias_fix_value, model_n_reps, model_data_r_file,
    bias_as_predictor) {

  model_rep_id <- cv <- species <- NULL

  # Training data - tibble -----
  training_data_file <- fs::path(
    dir_model_data, paste0("training_data_cv_", cv_fold, ".RData"))

  if (ecokit::check_data(training_data_file, warning = FALSE)) {
    training_data <- ecokit::load_as(training_data_file)
  } else {
    training_data <- dplyr::filter(model_data, cv != cv_fold)
    ecokit::save_as(
      object = training_data, object_name = "training_data",
      out_path = training_data_file)
  }

  # Training presences - tibble  ------
  training_pres_file <- fs::path(
    dir_model_data, paste0("training_pres_cv_", cv_fold, ".RData"))

  if (ecokit::check_data(training_pres_file, warning = FALSE)) {
    training_pres <- ecokit::load_as(training_pres_file)
  } else {
    training_pres <- dplyr::filter(training_data, species == 1L)
    ecokit::save_as(
      object = training_pres, object_name = "training_pres",
      out_path = training_pres_file)
  }
  n_training_pres <- nrow(training_pres)
  rm(training_pres, envir = environment())
  invisible(gc())

  # Training absences - tibble ------
  training_abs_file <- fs::path(
    dir_model_data, paste0("training_abs_cv_", cv_fold, ".RData"))
  if (ecokit::check_data(training_abs_file, warning = FALSE)) {
    training_abs <- ecokit::load_as(training_abs_file)
  } else {
    training_abs <- dplyr::filter(training_data, species == 0L)
    ecokit::save_as(
      object = training_abs, object_name = "training_abs",
      out_path = training_abs_file)
  }
  n_training_abs <- nrow(training_abs)

  rm(training_data, training_abs, envir = environment())
  invisible(gc())

  # Training data - raster -----
  training_r_file <- fs::path(
    dir_model_data, paste0("training_r_cv_", cv_fold, ".RData"))
  if (ecokit::check_data(training_r_file, warning = FALSE)) {
    training_r <- ecokit::load_as(training_r_file, unwrap_r = TRUE)
  } else {
    training_r <- ecokit::load_as(model_data_r_file, unwrap_r = TRUE)
    training_r[training_r$cv == cv_fold] <- NA
    training_r <- terra::toMemory(training_r)
    ecokit::save_as(
      object = terra::wrap(training_r),
      object_name = paste0("training_r_cv_", cv_fold),
      out_path = training_r_file)
  }

  # Pseudo-absences -------
  n_pseudo_abs <- pmin(n_training_pres * abs_ratio, n_training_abs)

  training_r_sub <- terra::subset(
    x = training_r, subset = "species", negate = TRUE)
  if (bias_as_predictor) {
    training_r_sub$bias_fixed <- NULL
  }
  training_r_sub <- terra::toMemory(training_r_sub)

  training_r_pres <- training_r
  training_r_pres[training_r_pres$species == 0L] <- NA
  training_r_pres <- terra::subset(
    x = training_r_pres, subset = "species") %>%
    terra::as.points()
  rm(training_r, envir = environment())
  invisible(gc())

  pseudo_abs <- tibble::tibble(model_rep_id = seq_len(model_n_reps)) %>%
    dplyr::mutate(
      pseudo_abs_file = purrr::map_chr(
        .x = model_rep_id,
        .f = function(rep_id) {

          pseudo_abs_file <- fs::path(
            dir_model_data,
            paste0("pseudo_abs_cv_", cv_fold, "_rep_", rep_id, ".RData"))

          if (!ecokit::check_data(pseudo_abs_file, warning = FALSE)) {
            pseudo_abs <- sdm::background(
              x = training_r_sub, n = n_pseudo_abs,
              method = "eDist", sp = training_r_pres) %>%
              tibble::as_tibble()

            ecokit::save_as(
              object = pseudo_abs,
              object_name = paste0(
                "pseudo_abs_cv_", cv_fold, "_rep_", rep_id),
              out_path = pseudo_abs_file)
          }
          pseudo_abs_file
        }
      ))

  rm(training_r_sub, training_r_pres, envir = environment())
  invisible(gc())


  # Testing data -----
  testing_data_file <- fs::path(
    dir_model_data, paste0("testing_data_cv_", cv_fold, ".RData"))
  if (ecokit::check_data(testing_data_file, warning = FALSE)) {
    testing_data <- ecokit::load_as(testing_data_file)
  } else {
    testing_data <- dplyr::filter(model_data, cv == cv_fold)
    # if bias correction is used, add bias_fixed column
    if (bias_as_predictor) {
      testing_data <- dplyr::mutate(
        testing_data, bias_fixed = bias_fix_value)
    }
    ecokit::save_as(
      object = testing_data, object_name = "testing_data",
      out_path = testing_data_file)
  }
  n_test_pres <- sum(testing_data$species == 1L)
  n_test_abs <- sum(testing_data$species == 0L)

  rm(testing_data, envir = environment())
  invisible(gc())


  # Testing raster -----
  testing_r_file <- fs::path(
    dir_model_data, paste0("testing_r_cv_", cv_fold, ".RData"))
  if (!ecokit::check_data(testing_r_file, warning = FALSE)) {
    testing_r <- ecokit::load_as(model_data_r_file, unwrap_r = TRUE)
    testing_r[testing_r$cv == cv_fold] <- NA

    # if bias correction is used, set bias values to fixed value
    if (bias_as_predictor) {
      testing_r$bias_fixed <- terra::classify(
        testing_r$bias, rbind(cbind(-Inf, Inf, bias_fix_value)))
    }
    testing_r <- terra::toMemory(testing_r)

    ecokit::save_as(
      object = terra::wrap(testing_r),
      object_name = paste0("testing_r_cv_", cv_fold),
      out_path = testing_r_file)
    rm(testing_r, envir = environment())
    invisible(gc())
  }

  tibble::tibble(
    cv = cv_fold,
    training_data = training_data_file,
    training_data_r = training_r_file,
    training_pres = training_pres_file,
    n_training_pres = n_training_pres,
    training_abs = training_abs_file,
    n_training_abs = n_training_abs,
    n_model_reps = model_n_reps,
    n_pseudo_abs = n_pseudo_abs,
    pseudo_abs = pseudo_abs,
    testing_data = testing_data_file,
    testing_data_r = testing_r_file,
    n_test_pres = n_test_pres,
    n_test_abs = n_test_abs) %>%
    tidyr::unnest(col = "pseudo_abs")

}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# # ========================================================================= #
# fit_predict -------
# # ========================================================================= #

fit_predict <- function(
    cv_rep_id, bias_fix_value, bias_as_predictor,
    projection_inputs, predictor_names, models_cv, sdm_settings,
    dir_fit, dir_proj_reps, sdm_packages, proj_mask_file) {

  climate_option <- map_paths <- mod_method <- NULL

  cv_rep_data <- dplyr::slice(models_cv, cv_rep_id)
  training_pres <- cv_rep_data$training_pres
  pseudo_abs <- cv_rep_data$pseudo_abs_file
  training_abs <- cv_rep_data$training_abs
  testing_data <- cv_rep_data$testing_data
  sdm_method <- cv_rep_data$sdm_method

  if (sdm_method %in% names(sdm_settings)) {
    sdm_setting <- sdm_settings[sdm_method]
  } else {
    sdm_setting <- NULL
  }

  model_name <- paste0(
    sdm_method, "_cv", cv_rep_data$cv, "_rep", cv_rep_data$model_rep_id)

  model_fit_file <- fs::path(dir_fit, paste0(model_name, "_fit.qs2"))
  model_out_file <- fs::path(dir_fit, paste0(model_name, "_results.qs2"))

  if (ecokit::check_data(
    c(model_out_file, model_fit_file), warning = FALSE)) {
    return(model_out_file)
  }

  # Terra options ------

  terra_temp_dir <- fs::path_temp(
    paste0(model_name, "_", format(Sys.time(), "%H%M%S")))
  fs::dir_create(terra_temp_dir)
  withr::defer(try(fs::dir_delete(terra_temp_dir), silent = TRUE))

  terra::terraOptions(
    # temp dir for terra
    tempdir = terra_temp_dir,
    # fraction of RAM terra may use (0-0.9)
    memfrac = 0.1,
    # (GB) below which mem is assumed available
    memmin = 1L,
    # (GB) cap for terra
    memmax = 10L,
    # silence per-worker progress bars
    progress = 0L,
    todisk = TRUE)

  # Loading packages -------

  # Define packages to load based on sdm_method
  pkg_to_load <- dplyr::filter(sdm_packages, mod_method == sdm_method)
  if (nrow(pkg_to_load) > 0L) {
    pkg_to_load <- c("sdm", pkg_to_load$packages)
  } else {
    pkg_to_load <- "sdm"
  }
  ecokit::load_packages(package_list = pkg_to_load)


  # Model fitting ------

  if (ecokit::check_data(model_fit_file, warning = FALSE)) {

    fitted_m <- ecokit::load_as(model_fit_file)

  } else {

    # Training and testing data -------
    data_pres <- ecokit::load_as(training_pres) %>%
      dplyr::select(tidyselect::all_of(predictor_names)) %>%
      dplyr::mutate(species = 1L, .before = 1L)

    if (sdm_method == "maxent") {
      # Ensure maxent output path is set
      dir_maxent <- fs::path(dir_fit, "maxent_html", model_name)
      if (fs::dir_exists(dir_maxent)) {
        try(fs::dir_delete(dir_maxent), silent = TRUE)
      }
      fs::dir_create(dir_maxent)

      if (is.null(sdm_setting) || !is.list(sdm_setting) ||
          !("maxent" %in% names(sdm_setting))) {
        sdm_setting <- list(maxent = list(path = dir_maxent))
      } else {
        sdm_setting$maxent$path <- dir_maxent
      }

      # Create absence data by combining pseudo-absences and training
      # presences
      data_abs <- ecokit::load_as(training_abs) %>%
        dplyr::bind_rows(data_pres) %>%
        dplyr::select(tidyselect::all_of(predictor_names)) %>%
        dplyr::mutate(species = 0L, .before = 1L)

      # sample a maximum of 1e6 background points (default value) to avoid
      # memory issues
      max_n_bg_maxent <- getOption(
        "onesdm_max_n_bg_maxent", default = 1000000L)
      if (!is.numeric(max_n_bg_maxent)) {
        ecokit::stop_ctx(
          paste0(
            "The `onesdm_max_n_bg_maxent` option must be a single positive ",
            "numeric value."),
          max_n_bg_maxent = max_n_bg_maxent, cat_timestamp = FALSE)
      }
      max_n_bg_maxent <- as.integer(max_n_bg_maxent)
      if (nrow(data_abs) > max_n_bg_maxent) {
        data_abs <- dplyr::slice_sample(data_abs, n = max_n_bg_maxent)
      }

    } else {
      data_abs <- ecokit::load_as(pseudo_abs) %>%
        dplyr::select(tidyselect::all_of(predictor_names)) %>%
        dplyr::mutate(species = 0L, .before = 1L)
    }

    train_data <- dplyr::bind_rows(data_pres, data_abs)
    test_data <- ecokit::load_as(testing_data) %>%
      dplyr::select(tidyselect::all_of(c("species", predictor_names)))
    rm(data_pres, data_abs, envir = environment())
    invisible(gc())

    # Model formula ------
    model_formula <- paste0(
      "species ~ ", paste(predictor_names, collapse = " + ")) %>%
      stats::as.formula(env = baseenv())

    # Model data for sdm package ------
    sdm_data <- sdm::sdmData(
      formula = model_formula,
      train = as.data.frame(train_data), test = as.data.frame(test_data))

    # Model fitting ------
    fitted_m <- ecokit::quietly(
      sdm::sdm(
        formula = model_formula, data = sdm_data,
        methods = sdm_method, modelSettings = sdm_setting))

    # Reduce models objects -------
    fitted_m <- reduce_sdm_formulas(obj = fitted_m)

    rm(
      train_data, test_data, sdm_data, model_formula, envir = environment())
    invisible(gc())

    # delete some not-needed files to save space
    if (sdm_method == "maxent") {
      try({
        fs::path(
          dir_maxent,
          c("presence", "absence", "species_explain.bat", "species.csv")) %>%
          fs::file_delete()
      }, silent = TRUE)
    }

    # Save fitted model ------
    ecokit::save_as(fitted_m, out_path = model_fit_file)

  }

  ## Extract info from fitted model object -----
  extracted_data <- extract_sdm_info(model = fitted_m)

  ## Making projections -----

  projections <- projection_inputs %>%
    dplyr::mutate(
      proj = purrr::map2(
        .x = map_paths, .y = climate_option,
        .f = ~ {

          proj_name <- paste0(model_name, "_", .y)
          proj_file <- fs::path(dir_proj_reps, paste0(proj_name, ".tif"))

          if (ecokit::check_tiff(proj_file, warning = FALSE)) {
            pred2 <- terra::rast(proj_file)
          } else {

            proj_extent_r_0 <- terra::rast(proj_mask_file)
            pred_input_r <- terra::rast(.x) %>%
              terra::crop(proj_extent_r_0) %>%
              terra::mask(proj_extent_r_0)

            if (bias_as_predictor) {
              pred_input_r$bias <- terra::setValues(
                pred_input_r[[1L]], bias_fix_value) %>%
                terra::mask(pred_input_r[[1L]])
            }
            pred2 <- ecokit::quietly(
              predict(object = fitted_m, newdata = pred_input_r)) %>%
              stats::setNames(proj_name)

            rm(pred_input_r, proj_extent_r_0, envir = environment())
            invisible(gc())

            terra::writeRaster(
              x = pred2, filename = proj_file, overwrite = TRUE)
          }

          proj_okay <- terra::global(pred2, range, na.rm = TRUE) %>%
            unlist() %>%
            dplyr::between(-0.0000001, 1.00000001) %>%
            all()

          # Clean up terra temp files as soon as possible to avoid disk space
          # issues
          terra::tmpFiles(current = TRUE, remove = TRUE)

          rm(pred2, envir = environment())
          invisible(gc())

          return(
            tibble::tibble(proj_file = proj_file, proj_okay = proj_okay))
        }
      )) %>%
    tidyr::unnest("proj") %>%
    dplyr::select(-map_paths)

  model_results <- tibble::tibble(
    model_fit_file = model_fit_file,
    extracted_data,
    projections = list(projections))
  ecokit::save_as(object = model_results, out_path = model_out_file)

  rm(fitted_m, envir = environment())
  invisible(gc())

  model_out_file
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# # ========================================================================= #
# summarise_preds_cv -------
# # ========================================================================= #

summarise_preds_cv <- function(cv_id, models_cv_summary, dir_proj_cv) {

  proj_file <- n_model_reps <- climate_model <- climate_scenario <-
    year <- climate_option <- proj_okay <- NULL

  withr::with_envvar(
    new = c(NOAWT = "TRUE"),
    code = {
      suppressPackageStartupMessages(
        ecokit::quietly({

          models_cv_summary_sub <- dplyr::slice(
            models_cv_summary, cv_id) %>%
            dplyr::select(
              tidyselect::all_of(
                c("n_model_reps", "sdm_method", "cv_summary"))) %>%
            tidyr::unnest(cols = "cv_summary") %>%
            dplyr::select(
              -tidyselect::all_of(
                c("evaluation_training", "evaluation_testing",
                  "variable_importance", "response_curves"))) %>%
            tidyr::unnest(cols = "projections") %>%
            dplyr::filter(proj_okay) %>%
            dplyr::select(-tidyselect::all_of("proj_okay")) %>%
            dplyr::group_by(
              n_model_reps, climate_model, climate_scenario,
              year, climate_option) %>%
            dplyr::summarise(proj_file = list(proj_file), .groups = "drop")

          n_reps <- unique(models_cv_summary_sub$n_model_reps)
          gdal_settings <- c("COMPRESS=ZSTD", "ZSTD_LEVEL=22", "TILED=YES")

          if (n_reps > 1L) {

            models_cv_summary_sub <- models_cv_summary_sub %>%
              dplyr::mutate(
                projections = purrr::map(
                  .x = proj_file,
                  .f = function(proj_file) {

                    pred_name <- basename(proj_file[1L]) %>%
                      fs::path_ext_remove() %>%
                      stringr::str_replace("_rep(\\d)+_", "_") #nolint

                    pred_maps <- terra::rast(proj_file)

                    out_file_mean <- fs::path(
                      dir_proj_cv, paste0(pred_name, "_mean.tif"))
                    mean_file_okay <- ecokit::check_tiff(
                      out_file_mean, warning = FALSE)
                    if (!mean_file_okay) {
                      pred_mean <- terra::app(
                        x = pred_maps, fun = mean, na.rm = TRUE) %>%
                        stats::setNames(paste0(pred_name, "_mean"))
                      terra::writeRaster(
                        x = pred_mean, filename = out_file_mean,
                        overwrite = TRUE, gdal = gdal_settings)
                      rm(pred_mean, envir = environment())
                      invisible(gc())
                    }

                    out_file_sd <- fs::path(
                      dir_proj_cv, paste0(pred_name, "_sd.tif"))
                    sd_file_okay <- ecokit::check_tiff(
                      out_file_sd, warning = FALSE)
                    if (!sd_file_okay) {
                      pred_sd <- terra::app(
                        x = pred_maps, fun = stats::sd, na.rm = TRUE) %>%
                        stats::setNames(paste0(pred_name, "_sd"))
                      terra::writeRaster(
                        x = pred_sd, filename = out_file_sd,
                        overwrite = TRUE, gdal = gdal_settings)
                      rm(pred_sd, envir = environment())
                      invisible(gc())
                    }

                    rm(pred_maps, envir = environment())
                    invisible(gc())

                    out_file_cov <- fs::path(
                      dir_proj_cv, paste0(pred_name, "_cov.tif"))
                    cov_file_okay <- ecokit::check_tiff(
                      out_file_cov, warning = FALSE)
                    if (!cov_file_okay) {
                      pred_mean <- terra::rast(out_file_mean)
                      # Replace very small mean values with reasonably small
                      # number to avoid overflow warning
                      pred_mean[pred_mean < 1e-8] <- 1e-8
                      pred_sd <- terra::rast(out_file_sd)
                      pred_cov <- (pred_sd / pred_mean) %>%
                        stats::setNames(paste0(pred_name, "_cov"))
                      terra::writeRaster(
                        x = pred_cov, filename = out_file_cov,
                        overwrite = TRUE, gdal = gdal_settings)
                      rm(
                        pred_mean, pred_sd, pred_cov,
                        envir = environment())
                      invisible(gc())
                    }

                    tibble::tibble(
                      pred_mean = out_file_mean, pred_sd = out_file_sd,
                      pred_cov = out_file_cov)
                  })) %>%
              tidyr::unnest("projections") %>%
              dplyr::select(-tidyselect::all_of("proj_file"))

          } else {

            models_cv_summary_sub <- models_cv_summary_sub %>%
              dplyr::mutate(
                projections = purrr::map(
                  .x = proj_file,
                  .f = function(proj_file) {

                    pred_name <- basename(proj_file[1L]) %>%
                      fs::path_ext_remove() %>%
                      stringr::str_replace("_rep(\\d)+_", "_") #nolint

                    pred_maps <- terra::rast(proj_file)

                    out_file_mean <- fs::path(
                      dir_proj_cv, paste0(pred_name, "_mean.tif"))
                    mean_file_okay <- ecokit::check_tiff(
                      out_file_mean, warning = FALSE)
                    if (!mean_file_okay) {
                      pred_mean <- terra::app(
                        x = pred_maps, fun = mean, na.rm = TRUE) %>%
                        stats::setNames(paste0(pred_name, "_mean"))
                      terra::writeRaster(
                        x = pred_mean, filename = out_file_mean,
                        overwrite = TRUE, gdal = gdal_settings)
                      rm(pred_mean, envir = environment())
                      invisible(gc())
                    }

                    rm(pred_maps, envir = environment())
                    invisible(gc())

                    tibble::tibble(
                      pred_mean = out_file_mean, pred_sd = NA_character_,
                      pred_cov = NA_character_)
                  })) %>%
              tidyr::unnest("projections") %>%
              dplyr::select(-tidyselect::all_of("proj_file"))
          }

          models_cv_summary_sub <- models_cv_summary_sub %>%
            dplyr::select(-tidyselect::all_of("n_model_reps"))

        },
        "The following object is masked from",
        "Attaching package: ", "Loading required package",
        "Loaded ", "This version of "))
    })
  models_cv_summary_sub
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ------

# # ========================================================================= #
# summarise_preds -------
# # ========================================================================= #

summarise_preds <- function(model_id, models_summary, dir_proj) {

  proj_file <- pred_mean <- climate_option <- climate_model <-
    climate_scenario <- year <- NULL

  suppressPackageStartupMessages(
    ecokit::quietly({

      gdal_settings <- c("COMPRESS=ZSTD", "ZSTD_LEVEL=22", "TILED=YES")

      models_summary_sub <- dplyr::slice(models_summary, model_id) %>%
        dplyr::select(tidyselect::all_of(c("sdm_method", "projections"))) %>%
        tidyr::unnest(cols = "projections") %>%
        dplyr::group_by(
          climate_model, climate_scenario, year, climate_option) %>%
        dplyr::summarise(proj_file = list(pred_mean), .groups = "drop") %>%
        dplyr::mutate(
          projections = purrr::map(
            .x = proj_file,
            .f = function(proj_file) {

              pred_name <- basename(proj_file[1L]) %>%
                fs::path_ext_remove() %>%
                stringr::str_replace("_mean", "") %>%
                stringr::str_replace("_cv(\\d)+_", "_") #nolint

              pred_maps <- terra::rast(proj_file)

              out_file_mean <- fs::path(
                dir_proj, paste0(pred_name, "_mean.tif"))
              mean_file_okay <- ecokit::check_tiff(
                out_file_mean, warning = FALSE)
              if (!mean_file_okay) {
                pred_mean <- terra::app(
                  x = pred_maps, fun = mean, na.rm = TRUE) %>%
                  stats::setNames(paste0(pred_name, "_mean"))
                terra::writeRaster(
                  x = pred_mean, filename = out_file_mean,
                  overwrite = TRUE, gdal = gdal_settings)
                rm(pred_mean, envir = environment())
                invisible(gc())
              }

              out_file_sd <- fs::path(dir_proj, paste0(pred_name, "_sd.tif"))
              sd_file_okay <- ecokit::check_tiff(out_file_sd, warning = FALSE)
              if (!sd_file_okay) {
                pred_sd <- terra::app(
                  x = pred_maps, fun = stats::sd, na.rm = TRUE) %>%
                  stats::setNames(paste0(pred_name, "_sd"))
                terra::writeRaster(
                  x = pred_sd, filename = out_file_sd,
                  overwrite = TRUE, gdal = gdal_settings)
                rm(pred_sd, envir = environment())
                invisible(gc())
              }

              rm(pred_maps, envir = environment())
              invisible(gc())

              out_file_cov <- fs::path(dir_proj, paste0(pred_name, "_cov.tif"))
              cov_file_okay <- ecokit::check_tiff(out_file_cov, warning = FALSE)
              if (!cov_file_okay) {
                pred_mean <- terra::rast(out_file_mean)
                # Replace very small mean values with reasonably small
                # number to avoid overflow warning
                pred_mean[pred_mean < 1e-8] <- 1e-8
                pred_sd <- terra::rast(out_file_sd)
                pred_cov <- (pred_sd / pred_mean) %>%
                  stats::setNames(paste0(pred_name, "_cov"))
                terra::writeRaster(
                  x = pred_cov, filename = out_file_cov,
                  overwrite = TRUE, gdal = gdal_settings)
                rm(pred_mean, pred_sd, pred_cov, envir = environment())
                invisible(gc())
              }

              tibble::tibble(
                pred_mean = out_file_mean, pred_sd = out_file_sd,
                pred_cov = out_file_cov)
            })) %>%
        tidyr::unnest("projections") %>%
        dplyr::select(-tidyselect::all_of("proj_file"))

    },
    "The following object is masked from",
    "Attaching package: ", "Loading required package",
    "Loaded ", "This version of "))

  models_summary_sub
}
