check_trfn_nexus <- function (trfn) 
{
  tmptr = try(expr = read.nexus(trfn))
  if (class(tmptr) == "try-error") {
    stoptxt = paste0("STOP ERROR in check_BioGeoBEARS_run() or check_trfn(): There was an error in reading the tree file. You specified the trfn (TRee FileName) as '", 
                     inputs$trfn, "'.  Options to try:\n\n1. Read the error message to see if you can figure out what went wrong.\n\n.2. Try typing, in the R window, 'read.tree(trfn)', and see if this works.\n\n")
    cat("\n\n")
    cat(stoptxt)
    cat("\n\n")
    error_txt = "Error in if (tp[3] != \"\")"
    TF = grepl(pattern = error_txt, x = tmptr) == TRUE
    if (TF == TRUE) {
      txt = paste0("NOTE: Your error message begins with '", 
                   error_txt, "'.\n\nThis error typically is a result of trying to load a NEXUS phylogeny file, when what you need is a Newick phylogeny file. Some advice:\n\n(a) Open the file in a plain-text editor, and see what kind of file it is.\n\n(b) If it is a NEXUS-format phylogeny, read it in with FigTree or read.nexus() in the R package 'ape'.\n\n(c) Then save out to a Newick file (with FigTree) or with write.tree() (with the 'ape' package, in R).\n\n(d) Then use that Newick file as the input for BioGeoBEARS.")
      cat("\n\n")
      cat(txt)
      cat("\n\n")
    }
    return(tmptr)
  }
  return(tmptr)
}
#<bytecode: 0x10a6fbb88>
 # <environment: namespace:BioGeoBEARS>