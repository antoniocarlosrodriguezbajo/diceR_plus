#' Execute Unsupervised Feature Selection (UFS) Methods in MATLAB
#'
#' @description This function establishes a connection to a running MATLAB server, processes a given dataset
#' using various unsupervised feature selection (UFS) methods, retrieves results, and closes the session.
#'
#' The function dynamically interacts with MATLAB, executing predefined feature selection methods while optionally
#' automating input handling using AutoHotkey.
#'
#' @param data A numeric matrix or data frame representing the dataset to be analyzed.
#' @param UFS_Methods A named list specifying which UFS methods should be executed.
#'        If a method is set to `TRUE`, the AutoHotkey script will send `1` to MATLAB as input for default parameters.
#' @param ahk_executable_path A string specifying the path to the AutoHotkey executable (default: configured via package options).
#' @param ahk_script_path A string specifying the path to the AutoHotkey script (default: configured via package options).
#'
#' @return A named list containing results from each executed UFS method.
#'
#' @references
#' Abedinzadeh Torghabeh, F., Modaresnia, Y., Hosseini, S. A. (2023). *Auto-UFSTool: An Automatic Unsupervised
#' Feature Selection Toolbox for MATLAB*. Journal of AI and Data Mining. doi: [10.22044/jadm.2023.12820.2434](https://doi.org/10.22044/jadm.2023.12820.2434)
#'
#' @examples
#' \dontrun{
#' # Define methods for feature selection
#' methods <- list("InfFS" = FALSE, "CFS" = FALSE, "Laplacian" = TRUE)
#'
#' # Execute UFS methods with default AutoHotkey settings
#' results <- runUFS(Meat$x, methods)
#' print(results)
#'
#' # Custom AutoHotkey paths
#' results <- runUFS(Meat$x, methods,
#'                   ahk_executable_path = "D:/AutoHotkey.exe",
#'                   ahk_script_path = "D:/send_input.ahk")
#' }
#'
#' @export
runUFS <- function(data, UFS_Methods,
                   ahk_executable_path = getOption("ufs.ahk_executable_path"),
                   ahk_script_path = getOption("ufs.ahk_script_path")) {

  # Start MATLAB server
  Matlab$startServer()
  matlab <- Matlab()
  print(matlab)

  # Check if MATLAB server is running
  isOpen <- open(matlab)
  if (!isOpen) {
    print("MATLAB server is not running: waited 30 seconds.")
  }
  print(matlab)

  # Set the dataset variable
  setVariable(matlab, X = data)

  # Evaluate MATLAB command to display variables
  evaluate(matlab, "whos")

  UFS_Results <- list()

  # Loop through methods and execute them in MATLAB
  for (UFS_Method in names(UFS_Methods)) {
    print(UFS_Method)
    send_input <- UFS_Methods[[UFS_Method]]
    cmd <- sprintf('Result = Auto_UFSTool(X, "%s");', UFS_Method)

    # Set random seed in MATLAB
    evaluate(matlab, "rng(42);")

    # Execute AutoHotkey script if required
    if (send_input && !is.null(ahk_executable_path) && !is.null(ahk_script_path)) {
      system(sprintf('cmd /c start /B "%s" "%s"', ahk_executable_path, ahk_script_path))
    }

    # Run the feature selection method in MATLAB
    evaluate(matlab, cmd)
    result <- getVariable(matlab, "Result")
    UFS_Results[[UFS_Method]] <- result
  }

  # Close MATLAB connection
  close(matlab)

  return(UFS_Results)
}

