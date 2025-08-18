# zzz.R
.onAttach <- function(libname, pkgname) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    green <- crayon::green
    blue <- crayon::blue
    magenta <- crayon::magenta
    bold <- crayon::bold

    packageStartupMessage(
      "\n", bold(green("==========================================\n")),
      bold(blue("       🌿 Welcome to phyloR 🌿\n")),
      bold(magenta("  Phylogenetic Analysis Made Simple\n")),
      bold(green("==========================================\n")),
      "Version: ", utils::packageVersion("phyloR"), "\n",
      "Author : Feifei Li & Bo Li\n",
      "Help   : type ", bold(blue("?phyloR")), " to get started.\n"
    )
  } else {
    packageStartupMessage(
      "\n*** Welcome to phyloR ***\n",
      "Version: ", utils::packageVersion("phyloR"), "\n"
    )
  }
}
