fix_package_timestamps <- function(pkg_dir = ".") {
  stopifnot(dir.exists(pkg_dir))
  op <- setwd(pkg_dir); on.exit(setwd(op), add = TRUE)

  # 需要排除的目录/文件（不会动这些时间戳）
  exclude_dirs <- c(
    ".git", ".Rproj.user", "packrat", "renv", "revdep",
    "vignettes/.Rhistory", "Meta", ".github", ".vscode"
  )
  exclude_files <- c(
    ".Rhistory", ".RData", ".DS_Store"
  )
  # 任何已存在的检查输出目录也排除
  rcheck_dirs <- list.dirs(".", recursive = FALSE, full.names = TRUE)
  rcheck_dirs <- rcheck_dirs[grepl("\\.Rcheck$", rcheck_dirs)]
  exclude_dirs <- unique(c(exclude_dirs, basename(rcheck_dirs)))

  # 递归列出需要处理的文件
  all_files <- list.files(".", all.files = TRUE, recursive = TRUE, full.names = TRUE, no.. = TRUE)

  keep <- file.info(all_files, extra_cols = FALSE)
  # 排除目录
  drop_dir <- Reduce(`|`, lapply(exclude_dirs, function(p) startsWith(all_files, paste0("./", p))))
  drop_dir[is.na(drop_dir)] <- FALSE
  # 排除单个文件名
  drop_file <- basename(all_files) %in% exclude_files
  # 排除目录本身（只改文件）
  drop_isdir <- keep$isdir %in% TRUE

  target <- all_files[!(drop_dir | drop_file | drop_isdir)]
  if (length(target)) {
    message("Updating timestamps for ", length(target), " files...")
    now <- Sys.time()
    # 分批防止超长向量
    for (idx in split(seq_along(target), ceiling(seq_along(target)/1000))) {
      Sys.setFileTime(target[idx], now)
    }
  } else {
    message("No eligible files to update.")
  }

  # 清理 .Rcheck 目录
  rcheck <- list.dirs(".", recursive = FALSE, full.names = TRUE)
  rcheck <- rcheck[grepl("\\.Rcheck$", rcheck)]
  if (length(rcheck)) {
    message("Removing old check directories: ", paste(basename(rcheck), collapse = ", "))
    unlink(rcheck, recursive = TRUE, force = TRUE)
  }

  invisible(TRUE)
}

# 执行：修正时间戳并重新 build & check
fix_package_timestamps(".")
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Please install.packages('devtools') first.")
}
devtools::build()
devtools::check()
