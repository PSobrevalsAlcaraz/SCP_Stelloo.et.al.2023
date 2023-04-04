
##############
### SCP 2023
##############

### ---- Library
{
  
  
  library(devtools)
  
  load_pack <- function(pack) {
    options("pac_update" = F)
    if(!'pacman' %in% installed.packages()[,'Package']) {
      install.packages('pacman')
    }
    library(pacman)
    p_load(char = pack, character.only = TRUE)
  }
  
  load_pack(c("sva","RegressionLibs","factoextra","randomcoloR","expss","svDialogs","dplyr","scater",'ggbreak','viridis',"seqinr","ggpubr","ggplot2","RColorBrewer","scp","magrittr","tidyr","limma","dplyr","proDA","tcltk","anchors","ggplot2","ggrepel","plotly","DEP","reshape2", "ggthemes"))
  
  
  
  inputs <- function(ini.value){
    fontSub <- tkfont.create(weight="bold", size=10)
    
    xvar <- tclVar(as.character(ini.value))
    
    tt <- tktoplevel()
    tkwm.title(tt,"MedianCV threshold")
    x.entry <- tkentry(tt, textvariable=xvar, width=5)
    
    reset <- function()
    {
      tclvalue(xvar)<-as.character(ini.value)
    }
    
    reset.but <- tkbutton(tt, text="Reset", command=reset)
    
    submit <- function() {
      x <- as.numeric(tclvalue(xvar))
      
      e <- parent.env(environment())
      e$x <- x
      
      tkdestroy(tt)
    }
    submit.but <- tkbutton(tt, text="Submit", command=submit)
    
    
    tkgrid(tklabel(tt,text="Determine MedianCV threshold:",font=fontSub),x.entry,pady = 10, padx =10)
    
    
    tkgrid(submit.but, reset.but,pady = 10, padx =10)
    
    
    
    
    tkwait.window(tt)
    return(c(x))
  }
  
  
  .splitSCE <- function(x, 
                        f) {
    ## Check that f is a factor
    if (is.character(f)) {
      if (length(f) != 1) 
        stop("'f' must be of lenght one")
      if (f %in% colnames(rowData(x))) {
        f <- rowData(x)[, f]
      }
      else if (f %in% colnames(colData(x))) {
        f <- colData(x)[, f]
      }
      else {
        stop("'", f, "' not found in rowData or colData")
      }
      if (!is.factor(f)) 
        f <- factor(f)
    }
    ## Check that the factor matches one of the dimensions
    if (!length(f) %in% dim(x)) 
      stop("length(f) not compatible with dim(x).")
    if (length(f) == nrow(x)) { ## Split along rows
      xl <- lapply(split(rownames(x), f = f), function(i) x[i, ])
    } else { ## Split along columns
      xl <- lapply(split(colnames(x), f = f), function(i) x[, i])
    }
    ## Convert list to an ExperimentList
    do.call(ExperimentList, xl)
  }
  
  
  StatBin2 <- ggproto(
    "StatBin2", 
    StatBin,
    compute_group = function (data, scales, binwidth = NULL, bins = NULL, 
                              center = NULL, boundary = NULL, 
                              closed = c("right", "left"), pad = FALSE, 
                              breaks = NULL, origin = NULL, right = NULL, 
                              drop = NULL, width = NULL) {
      if (!is.null(breaks)) {
        if (!scales$x$is_discrete()) {
          breaks <- scales$x$transform(breaks)
        }
        bins <- ggplot2:::bin_breaks(breaks, closed)
      }
      else if (!is.null(binwidth)) {
        if (is.function(binwidth)) {
          binwidth <- binwidth(data$x)
        }
        bins <- ggplot2:::bin_breaks_width(scales$x$dimension(), binwidth, 
                                           center = center, boundary = boundary, 
                                           closed = closed)
      }
      else {
        bins <- ggplot2:::bin_breaks_bins(scales$x$dimension(), bins, 
                                          center = center, boundary = boundary, 
                                          closed = closed)
      }
      res <- ggplot2:::bin_vector(data$x, bins, weight = data$weight, pad = pad)
      
      # drop 0-count bins completely before returning the dataframe
      res <- res[res$count > 0, ] 
      
      res
    })
  
  
  extract.gene <- function(SCP, fastafile = fasta.file, description = fasta.description){
    fasta <- read.fasta(fastafile,as.string=T, seqtype = "AA",
                        whole.header = TRUE,set.attributes = FALSE)
    fasta <- names(fasta)
    fasta <- gsub(";","",fasta, fixed=T)
    fasta <- gsub(",","",fasta, fixed=T)
    
    fasta <- sub("|",";",fasta,fixed=T)
    fasta <- sub("|",",",fasta,fixed=T)
    fastas <- data.frame(fasta = fasta)
    fastas$Master.Protein.Accessions <- gsub(",.*","",gsub(".*;","",fastas$fasta,perl=T))
    fastas$Gene.name <- trimws(gsub("PE=.*","",gsub(".*GN=","",fastas$fasta,perl=T)))
    fastas$Gene.name[grep(" ",fastas$Gene.name,fixed = T,perl=F)] <- fastas$Master.Protein.Accessions[grep(" ",fastas$Gene.name,fixed = T,perl=F)]
    fastas$Description <- gsub(paste(".*_",description,sep=""),"",fastas$fasta,perl=T)
    fastas$fasta=NULL
    
    scp.name <- data.frame()
    for( r in 1:nrow(SCP)){
      scp.temp <- SCP[r,]
      scp.temp <- data.frame(separate_rows(scp.temp,Master.Protein.Accessions, sep= "; " ))
      scp.temp <- merge(fastas, scp.temp, by = "Master.Protein.Accessions", all.x = F, all.y = T)
      
      genes <- unique(trimws(scp.temp$Gene.name))
      desc <- unique(trimws(scp.temp$Description))
      scp.temp <- collapse_rows(scp.temp, c("Master.Protein.Accessions" ,"Gene.name","Description"),collapse= ";")
      scp.temp$Gene.name <- paste(genes, collapse=";")
      scp.temp$Description <- paste(desc, collapse=";")
      scp.name <- rbind(scp.name,scp.temp)
    }
    scp.name$Gene.name.Unique <- make.unique(scp.name$Gene.name)
    return(scp.name)
  }
  
  random_imp <- function(x) {
    
    m <- mean(x, na.rm = TRUE)
    sdev <- sd(x, na.rm = TRUE)
    n <- sum(is.na(x))
    set.seed(42)
    x[is.na(x)] <- rnorm(n, mean = m, sd = sdev/2)
    x
  }
  
  manual_impute <- function(dat, scale = 0.3, shift = 1.8){
    set.seed(42)
    colnames(dat) <- gsub("-",".",colnames(dat))
    if(is.integer(scale)) scale <- is.numeric(scale)
    if(is.integer(shift)) shift <- is.numeric(shift)
    # Show error if inputs are not the required classes
    
    # Show error if there are no missing values
    if(!any(is.na(dat))) {
      stop("No missing values in '", deparse(substitute(dat)), "'",
           call. = FALSE)
    }
    
    # Get descriptive parameters of the current sample distributions
    stat <- dat %>%
      data.frame() %>%
      tibble::rownames_to_column() %>%
      gather(samples, value, -rowname) %>%
      filter(!is.na(value))  %>%
      group_by(samples) %>%
      summarise(mean = mean(value),
                median = median(value),
                sd = sd(value),
                n = n(),
                infin = nrow(dat) - n)
    # Impute missing values by random draws from a distribution
    # which is left-shifted by parameter 'shift' * sd and scaled by parameter 'scale' * sd.
    for (a in colnames(dat)) {
      dat[is.na(dat[, a]), a] <-
        rnorm(stat$infin[grep(a, stat$samples)[1]],
              mean = stat$median[grep(a, stat$samples)[1]] - shift * stat$sd[grep(a, stat$samples)[1]],
              sd = stat$sd[grep(a, stat$samples)[1]] * scale)
    }
    return(dat)
  }
  
  make.proteins.SCP <- function(PSM){
    psm.file <- PSM
    psm.file = DataFrame(psm.file)
    red <- reduceDataFrame(psm.file,psm.file$Annotated.Sequence, simplify = FALSE)
    diferences <- sapply(1:nrow(red),function(x) length(unique(red$Master.Protein.Accessions[[x]])))
    for( d in which(diferences>1)){
      prots <- unique(psm.file[which(psm.file$Annotated.Sequence==rownames(red[d,])),"Master.Protein.Accessions"])
      for( p in 1:length(prots)){
        psm.file[which((psm.file$Annotated.Sequence==rownames(red[d,])) & (psm.file$Master.Protein.Accessions==prots[[p]])),'Annotated.Sequence'] <-
          paste(rownames(red[d,]),paste(rep("Z",p),collapse=""),sep="")
      }
    }
    psm.file <- data.frame(psm.file)
  }
  
  collapse_rows <- function(data, columns, collapse = ";"){
    initial.order <- colnames(data)
    unique.columns <- colnames(data)[!(colnames(data) %in% columns)]
    new.data <- unique(data[,unique.columns])
    for( rrow in 1:nrow(new.data)){
      temp.data <- data[which(do.call("paste", data[,unique.columns]) %in%  do.call("paste", new.data[rrow,unique.columns])),columns]
      for(x in 1:length(columns)) {
        new.data[rrow,columns[x]] <- paste(temp.data[,x],collapse=";")
      }
    }
    new.data <- new.data[,initial.order]
    return(new.data)
  }
  
  TMTchoose <- function(ini.value){
    
    fontSub <- tkfont.create(weight="bold",size=10)
    
    xvar <- tclVar(as.character(ini.value))
    
    tt <- tktoplevel()
    tkwm.title(tt,"TMT Type")
    tktitle(tt)<-"TMT Type"  
    x.entry <- tkentry(tt, textvariable=xvar, width=5)
    
    reset <- function()
    {
      tclvalue(xvar)<-as.character(ini.value)
    }
    
    reset.but <- tkbutton(tt, text="Reset", command=reset)
    
    submit <- function() {
      x <- as.numeric(tclvalue(xvar))
      
      e <- parent.env(environment())
      e$x <- x
      
      tkdestroy(tt)
    }
    submit.but <- tkbutton(tt, text="Submit", command=submit)
    
    tkgrid(tklabel(tt,text="Determine TMT type:", font=fontSub),x.entry,pady = 10, padx =10)
    
    
    tkgrid(submit.but, reset.but,pady = 10, padx =10)
    
    
    tkwait.window(tt)
    return(c(x))
  }
  
  check.PSM <- function(filename, n.obs){
    
    img <- tclVar()
    #tclimg <- tkimage.create("photo", img, file =  filename)
    
    
    fontTite <- tkfont.create(weight="bold", underline=TRUE)
    fontInfo <- tkfont.create(slant="italic",size=8)
    
    psm.values = list()
    psm.values$n.obs = n.obs
    
    
    n.obs <- tclVar(as.character(n.obs))
    
    
    
    tt <- tktoplevel()
    tkwm.title(tt,"PSMs Threshold")
    n.obs.entry <- tkentry(tt, textvariable = n.obs, width=5)
    tcl("image","create","photo", "imageID", file=filename)
    imgs <- ttklabel(tt, image="imageID", compound="image")
    
    
    reset <- function(){
      tclvalue(n.obs) <- as.character(n.obs)
      
    }
    
    reset.but <- tkbutton(tt, text="Reset", command=reset)
    
    submit <- function(){
      
      psm.values$n.obs <<- as.numeric(tclvalue(n.obs))
      tkdestroy(tt)
    }
    
    
    submit.but <- tkbutton(tt, text = "Sumbit", 
                           command = submit)
    
    
    tkgrid(tklabel(tt,text="Minimum number of PSMs for Sample", font=fontTite),columnspan = 5,pady = 20, padx =10)
    tkgrid(tklabel(tt,text="Minimum number of PSMs:"),n.obs.entry,sticky = "w",pady = 10, padx =10)
    
    tkgrid(imgs,columnspan=5)
    tkgrid(submit.but, reset.but,pady = 10, padx =10)
    
    
    tkwait.window(tt)
    return(psm.values)
    
  }
  
  check.SCR <- function(filename, Top.scr.thr, Low.scr.thr){
    
    img <- tclVar()
    #tclimg <- tkimage.create("photo", img, file =  filename)
    
    
    fontTite <- tkfont.create(weight="bold", underline=TRUE)
    fontInfo <- tkfont.create(slant="italic",size=8)
    
    scr.values = list()
    scr.values$Top.scr.thr = Top.scr.thr
    scr.values$Low.scr.thr = Low.scr.thr
    
    Top.scr.thr <- tclVar(as.character(Top.scr.thr))
    Low.scr.thr <- tclVar(as.character(Low.scr.thr))
    
    
    tt <- tktoplevel()
    tkwm.title(tt,"Median SCR threshold")
    Low.scr.thr.entry <- tkentry(tt, textvariable = Low.scr.thr, width=5)
    Top.scr.thr.entry <- tkentry(tt, textvariable = Top.scr.thr, width=5)
    tcl("image","create","photo", "imageID", file=filename)
    imgs <- ttklabel(tt, image="imageID", compound="image")
    
    
    
    reset <- function(){
      tclvalue(Low.scr.thr) <- as.character(Low.scr.thr)
      tclvalue(Top.scr.thr) <- as.character(Top.scr.thr)
      
    }
    
    reset.but <- tkbutton(tt, text="Reset", command=reset)
    
    submit <- function(){
      
      scr.values$Top.scr.thr <<- as.numeric(tclvalue(Top.scr.thr))
      scr.values$Low.scr.thr  <<- as.numeric(tclvalue(Low.scr.thr))
      tkdestroy(tt)
    }
    
    
    submit.but <- tkbutton(tt, text = "Sumbit", 
                           command = submit)
    
    
    tkgrid(tklabel(tt,text="Median SCR values", font=fontTite),columnspan = 5,pady = 20, padx =10)
    tkgrid(tklabel(tt,text="Please write values for both thresholds:",font=fontInfo,background="white"),
           sticky = "w",pady =0, padx =10)
    tkgrid(tklabel(tt,text="Low Threshold"),Low.scr.thr.entry,tklabel(tt,text="High Threshold"),Top.scr.thr.entry,sticky = "w",pady = 10, padx =10)
    
    tkgrid(imgs,columnspan=5)
    tkgrid(submit.but, reset.but,pady = 10, padx =10)
    
    
    tkwait.window(tt)
    return(scr.values)
    
  }
  
  check.FDR <- function(filename, fdr.thr){
    
    img <- tclVar()
    #tclimg <- tkimage.create("photo", img, file =  filename)
    
    
    fontTite <- tkfont.create(weight="bold", underline=TRUE)
    fontInfo <- tkfont.create(slant="italic",size=8)
    
    fdr.values = list()
    fdr.values$fdr.thr = fdr.thr
    
    
    fdr.thr <- tclVar(as.character(fdr.thr))
    
    tt <- tktoplevel()
    tkwm.title(tt,"FDR threshold")
    fdr.thr.entry <- tkentry(tt, textvariable = fdr.thr, width=5)
    
    tcl("image","create","photo", "imageID", file=filename)
    imgs <- ttklabel(tt, image="imageID", compound="image")
    
    
    
    reset <- function(){
      tclvalue(fdr.thr) <- as.character(fdr.thr)
      
      
    }
    
    reset.but <- tkbutton(tt, text="Reset", command=reset)
    
    submit <- function(){
      
      fdr.values$fdr.thr <<- as.numeric(tclvalue(fdr.thr))
      
      tkdestroy(tt)
    }
    
    
    submit.but <- tkbutton(tt, text = "Sumbit", 
                           command = submit)
    
    
    tkgrid(tklabel(tt,text="FDR values", font=fontTite),columnspan = 5,pady = 20, padx =10)
    tkgrid(tklabel(tt,text="Threshold"),fdr.thr.entry,sticky = "w",pady = 10, padx =10)
    
    tkgrid(imgs,columnspan=5)
    tkgrid(submit.but, reset.but,pady = 10, padx =10)
    
    
    tkwait.window(tt)
    return(fdr.values)
    
  }
  
  
  check.mCVC <- function(filename1, filename2,Top.medianCV, Low.medianCV){
    
    img <- tclVar()
    #tclimg <- tkimage.create("photo", img, file =  filename)
    
    
    fontTite <- tkfont.create(weight="bold", underline=TRUE)
    fontInfo <- tkfont.create(slant="italic",size=8)
    
    mcv.values = list()
    mcv.values$Top.medianCV = Top.medianCV
    mcv.values$Low.medianCV = Low.medianCV
    
    
    Top.medianCV <- tclVar(as.character(Top.medianCV))
    Low.medianCV <- tclVar(as.character(Low.medianCV))
    
    tt <- tktoplevel()
    tkwm.title(tt,"MedianCV threshold")
    Top.medianCV.entry <- tkentry(tt, textvariable = Top.medianCV, width=5)
    Low.medianCV.entry <- tkentry(tt, textvariable = Low.medianCV, width=5)
    
    tcl("image","create","photo", "imageID1", file=filename1)
    imgs1 <- ttklabel(tt, image="imageID1", compound="image")
    tcl("image","create","photo", "imageID2", file=filename2)
    imgs2 <- ttklabel(tt, image="imageID2", compound="image")
    
    
    
    reset <- function(){
      tclvalue(Top.medianCV) <- as.character(Top.medianCV)
      tclvalue(Low.medianCV) <- as.character(Low.medianCV)
      
      
    }
    
    reset.but <- tkbutton(tt, text="Reset", command=reset)
    
    submit <- function(){
      
      mcv.values$Top.medianCV <<- as.numeric(tclvalue(Top.medianCV))
      mcv.values$Low.medianCV <<- as.numeric(tclvalue(Low.medianCV))
      
      tkdestroy(tt)
    }
    
    
    submit.but <- tkbutton(tt, text = "Sumbit", 
                           command = submit)
    
    
    tkgrid(tklabel(tt,text="Median CV", font=fontTite),columnspan = 5,pady = 20, padx =10)
    tkgrid(tklabel(tt,text="Low Threshold"),Low.medianCV.entry,tklabel(tt,text="Top Threshold"),Top.medianCV.entry,sticky = "w",pady = 10, padx =10)
    
    tkgrid(imgs1,imgs2,columnspan=5)
    tkgrid(submit.but, reset.but,pady = 10, padx =10)
    
    
    tkwait.window(tt)
    return(mcv.values)
    
  }
  
  check.mRI <- function(filename1,Top.medianRI, Low.medianRI){
    
    img <- tclVar()
    #tclimg <- tkimage.create("photo", img, file =  filename)
    
    
    fontTite <- tkfont.create(weight="bold", underline=TRUE)
    fontInfo <- tkfont.create(slant="italic",size=8)
    
    mRI.values = list()
    mRI.values$Top.medianRI = Top.medianRI
    mRI.values$Low.medianRI = Low.medianRI
    
    
    Top.medianRI <- tclVar(as.character(Top.medianRI))
    Low.medianRI <- tclVar(as.character(Low.medianRI))
    
    tt <- tktoplevel()
    tkwm.title(tt,"MedianRI threshold")
    Top.medianRI.entry <- tkentry(tt, textvariable = Top.medianRI, width=5)
    Low.medianRI.entry <- tkentry(tt, textvariable = Low.medianRI, width=5)
    
    tcl("image","create","photo", "imageID1", file=filename1)
    imgs1 <- ttklabel(tt, image="imageID1", compound="image")
    
    
    reset <- function(){
      tclvalue(Top.medianRI) <- as.character(Top.medianRI)
      tclvalue(Low.medianRI) <- as.character(Low.medianRI)
      
      
    }
    
    reset.but <- tkbutton(tt, text="Reset", command=reset)
    
    submit <- function(){
      
      mRI.values$Top.medianRI <<- as.numeric(tclvalue(Top.medianRI))
      mRI.values$Low.medianRI <<- as.numeric(tclvalue(Low.medianRI))
      
      tkdestroy(tt)
    }
    
    
    submit.but <- tkbutton(tt, text = "Sumbit", 
                           command = submit)
    
    
    tkgrid(tklabel(tt,text="Median RI", font=fontTite),columnspan = 5,pady = 20, padx =10)
    tkgrid(tklabel(tt,text="Low Threshold"),Low.medianRI.entry,tklabel(tt,text="Top Threshold"),Top.medianRI.entry,sticky = "w",pady = 10, padx =10)
    
    tkgrid(imgs1,columnspan=10)
    tkgrid(submit.but, reset.but,pady = 10, padx =10)
    
    
    tkwait.window(tt)
    return(mRI.values)
    
  }
  
  
  check.qc <- function(filename1, qc.files.remove){
    
    img <- tclVar()
    #tclimg <- tkimage.create("photo", img, file =  filename)
    
    
    fontTite <- tkfont.create(weight="bold", underline=TRUE)
    fontInfo <- tkfont.create(slant="italic",size=8)
    
    QC.values = list()
    QC.values$qc.files.remove = as.character(qc.files.remove)
    
    
    qc.files.remove <- tclVar(as.character(qc.files.remove))
    
    tt <- tktoplevel()
    tkwm.title(tt,"Per cell QC metrics")
    qc.files.remove.entry <- tkentry(tt, textvariable = qc.files.remove, width=10)
    
    tcl("image","create","photo", "imageID1", file=filename1)
    imgs1 <- ttklabel(tt, image="imageID1", compound="image")
    
    
    reset <- function(){
      tclvalue(qc.files.remove) <- as.character(qc.files.remove)
      
    }
    
    reset.but <- tkbutton(tt, text="Reset", command=reset)
    
    submit <- function(){
      
      QC.values$qc.files.remove <<- as.character(tclvalue(qc.files.remove))
      
      tkdestroy(tt)
    }
    
    
    submit.but <- tkbutton(tt, text = "Sumbit", 
                           command = submit)
    
    
    tkgrid(tklabel(tt,text="QC metrics", font=fontTite),columnspan = 5,pady = 20, padx =10)
    tkgrid(tklabel(tt,text="Files to remove"),qc.files.remove.entry,sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="*Leave blank if none, if multiple separate by ','"),sticky = "w",pady = 10, padx =10)
    tkgrid(imgs1,columnspan=10)
    tkgrid(submit.but, reset.but,pady = 10, padx =10)
    
    
    tkwait.window(tt)
    return(QC.values)
    
  }
  
  
  InitialValues <- function(ins.n.obs, ins.Top.scr.thr, ins.Low.scr.thr, ins.booster, 
                            ins.fdr.thr, ins.use.medianCV,ins.Top.medianRI,ins.Low.medianRI,
                            ins.Top.medianCV,ins.Low.medianCV, ins.batch,ins.imputation,
                            ins.extra.run,ins.save.plots,
                            ins.save.tables, ins.experiment, ins.experiment.name, ins.miss.perc,
                            ins.fastafile, ins.description){
    
    Initial.values = list()
    Initial.values$n.obs = ins.n.obs
    Initial.values$Top.scr.thr = ins.Top.scr.thr
    Initial.values$Low.scr.thr = ins.Low.scr.thr
    Initial.values$booster = ins.booster
    Initial.values$fdr.thr = ins.fdr.thr
    Initial.values$use.medianCV = ins.use.medianCV
    Initial.values$Top.medianCV = ins.Top.medianCV
    Initial.values$Low.medianCV = ins.Low.medianCV
    Initial.values$Top.medianRI = ins.Top.medianRI
    Initial.values$Low.medianRI = ins.Low.medianRI
    Initial.values$imputation = ins.imputation
    Initial.values$batch = ins.batch
    Initial.values$extra.run = ins.extra.run
    Initial.values$save.plots = ins.save.plots
    Initial.values$save.tables = ins.save.tables
    Initial.values$experiment = ins.experiment
    Initial.values$experiment.name = ins.experiment.name
    Initial.values$miss.perc = ins.miss.perc
    Initial.values$description = ins.description
    Initial.values$fastafile = ins.fastafile
    
    fontTite <- tkfont.create(weight="bold", underline=TRUE,size=14)
    fontSub <- tkfont.create(weight="bold")
    fontNam <- tkfont.create(weight="bold",size=10)
    fontInfo <- tkfont.create(slant="italic",size=10)
    
    
    experiment.name <- tclVar(as.character(ins.experiment.name))
    n.obs <- tclVar(as.character(ins.n.obs))
    Top.scr.thr <- tclVar(as.character(ins.Top.scr.thr))
    Low.scr.thr <- tclVar(as.character(ins.Low.scr.thr))
    booster <- tclVar(as.character(ins.booster))
    fdr.thr <- tclVar(as.character(ins.fdr.thr))
    use.medianCV <- tclVar(as.character(ins.use.medianCV))
    Top.medianCV <- tclVar(as.character(ins.Top.medianCV))
    Low.medianCV <- tclVar(as.character(ins.Low.medianCV))
    Top.medianRI <- tclVar(as.character(ins.Top.medianRI))
    Low.medianRI <- tclVar(as.character(ins.Low.medianRI))
    batch <- tclVar(as.character(ins.batch))
    description <- tclVar(as.character(ins.description))
    fastafile <- tclVar(as.character(ins.fastafile))
    
    tt <- tktoplevel(background="white")
    imputation1 <- tkradiobutton(tt,background="white")
    imputation0 <- tkradiobutton(tt,background="white")
    imputation <- tclVar(ins.imputation)
    extra.run1 <- tkradiobutton(tt,background="white")
    extra.run0 <- tkradiobutton(tt,background="white")
    extra.run <- tclVar(ins.extra.run)
    save.plots1 <- tkradiobutton(tt,background="white")
    save.plots0 <- tkradiobutton(tt,background="white")
    save.plots <- tclVar(ins.save.plots)
    save.tables1 <- tkradiobutton(tt,background="white")
    save.tables0 <- tkradiobutton(tt,background="white")
    save.tables <- tclVar(ins.save.tables)
    experiment1 <- tkradiobutton(tt,background="white")
    experiment0 <- tkradiobutton(tt,background="white")
    experiment <- tclVar(ins.experiment)
    miss.perc1 <- tkradiobutton(tt,background="white")
    miss.perc0 <- tkradiobutton(tt,background="white")
    miss.perc <- tclVar(ins.miss.perc)
    
    
    
    tkwm.title(tt,"SCP Values")
    experiment.name.entry <- tkentry(tt, textvariable = experiment.name, font=fontNam,width=40)
    n.obs.entry <- tkentry(tt, textvariable = n.obs, width=5)
    Low.scr.thr.entry <- tkentry(tt, textvariable = Low.scr.thr, width=5)
    Top.scr.thr.entry <- tkentry(tt, textvariable = Top.scr.thr, width=5)
    booster.entry <- tkentry(tt, textvariable = booster, width= 25)
    fdr.thr.entry <- tkentry(tt, textvariable = fdr.thr, width=5)
    use.medianCV.entry <- tkentry(tt, textvariable = use.medianCV)
    Top.medianCV.entry <- tkentry(tt, textvariable = Top.medianCV, width=5)
    Low.medianCV.entry <- tkentry(tt, textvariable = Low.medianCV, width=5)
    Top.medianRI.entry <- tkentry(tt, textvariable = Top.medianRI, width=5)
    Low.medianRI.entry <- tkentry(tt, textvariable = Low.medianRI, width=5)
    description.entry <- tkentry(tt, textvariable = description, width= 25)
    fastafile.entry <- tkentry(tt, textvariable = fastafile, width= 40)
    batch.entry <- tkentry(tt, textvariable = batch)
    tkconfigure(imputation1,variable=imputation,value=1)
    tkconfigure(imputation0,variable=imputation,value=0)
    tkconfigure(extra.run1,variable=extra.run,value=1)
    tkconfigure(extra.run0,variable=extra.run,value=0)
    tkconfigure(save.plots1,variable=save.plots,value=1)
    tkconfigure(save.plots0,variable=save.plots,value=0)
    tkconfigure(save.tables1,variable=save.tables,value=1)
    tkconfigure(save.tables0,variable=save.tables,value=0)
    tkconfigure(experiment1,variable=experiment,value=1)
    tkconfigure(experiment0,variable=experiment,value=0)
    tkconfigure(miss.perc1,variable=miss.perc,value=1)
    tkconfigure(miss.perc0,variable=miss.perc,value=0)
    
    
    reset <- function(){
      tclvalue(experiment.name) <- as.character(tclvalue(ins.experiment.name))
      tclvalue(n.obs) <- as.character(tclvalue(ins.n.obs))
      tclvalue(Low.scr.thr) <- as.character(tclvalue(ins.Low.scr.thr))
      tclvalue(Top.scr.thr) <- as.character(tclvalue(ins.Top.scr.thr))
      tclvalue(booster) <- as.character(tclvalue(ins.booster))
      tclvalue(fdr.thr) <- as.character(tclvalue(ins.fdr.thr))
      tclvalue(use.medianCV) <- as.character(tclvalue(ins.use.medianCV))
      tclvalue(Top.medianCV) <- as.character(tclvalue(ins.Top.medianCV))
      tclvalue(Low.medianCV) <- as.character(tclvalue(ins.Low.medianCV))
      tclvalue(Top.medianRI) <- as.character(tclvalue(ins.Top.medianRI))
      tclvalue(Low.medianRI) <- as.character(tclvalue(ins.Low.medianRI))
      tclvalue(batch) <- as.character(tclvalue(ins.batch))
      tclvalue(fastafile) <- as.character(tclvalue(ins.fastafile))
      tclvalue(description) <- as.character(tclvalue(ins.description))
      tclvalue(imputation) <- tclvalue(ins.imputation)
      tclvalue(extra.run) <- tclvalue(ins.extra.run)
      tclvalue(save.plots) <- tclvalue(ins.save.plots)
      tclvalue(save.tables) <- tclvalue(ins.save.tables)
      tclvalue(miss.perc) <- tclvalue(ins.miss.perc)
      tclvalue(experiment) <- tclvalue(ins.experiment)
      
    }
    
    reset.but <- tkbutton(tt, text="Reset to original values", command=reset)
    
    clean <- function(){
      tclvalue(experiment.name) <- ''
      tclvalue(n.obs) <- ''
      tclvalue(Low.scr.thr) <- ''
      tclvalue(Top.scr.thr) <- ''
      tclvalue(booster) <- 'Abundance.126'
      tclvalue(fdr.thr) <- 0.01
      tclvalue(use.medianCV) <- 'Channels'
      tclvalue(Top.medianCV) <- ''
      tclvalue(Low.medianCV) <- ''
      tclvalue(Top.medianRI) <- ''
      tclvalue(Low.medianRI) <- ''
      tclvalue(batch) <- 'Both'
      tclvalue(description) <- ins.description
      tclvalue(fastafile) <- ins.fastafile
      tclvalue(imputation) <- 0
      tclvalue(extra.run) <- 1
      tclvalue(save.plots) <- 1
      tclvalue(save.tables) <- 1
      tclvalue(experiment) <- 0
      tclvalue(miss.perc) <- 0
      
    }
    
    clean.but <- tkbutton(tt, text="Clean unnecessary data", command=clean)
    
    submit <- function(){
      
      Initial.values$experiment.name <<- as.character(tclvalue(experiment.name))
      Initial.values$n.obs <<- as.numeric(tclvalue(n.obs))
      Initial.values$Top.scr.thr <<- as.numeric(tclvalue(Top.scr.thr))
      Initial.values$Low.scr.thr  <<- as.numeric(tclvalue(Low.scr.thr))
      Initial.values$booster <<- as.character(tclvalue(booster))
      Initial.values$fdr.thr <<- as.numeric(tclvalue(fdr.thr))
      Initial.values$use.medianCV <<- as.character(tclvalue(use.medianCV))
      Initial.values$Top.medianCV <<- as.numeric(tclvalue(Top.medianCV))
      Initial.values$Low.medianCV <<- as.numeric(tclvalue(Low.medianCV))
      Initial.values$Top.medianRI <<- as.numeric(tclvalue(Top.medianRI))
      Initial.values$Low.medianRI <<- as.numeric(tclvalue(Low.medianRI))
      Initial.values$imputation <<- as.character(tclvalue(imputation))
      Initial.values$batch <<- as.character(tclvalue(batch))
      Initial.values$description <<- as.character(tclvalue(description))
      Initial.values$fastafile <<- as.character(tclvalue(fastafile))
      Initial.values$extra.run <<- as.numeric(tclvalue(extra.run))
      Initial.values$save.plots <<- as.numeric(tclvalue(save.plots))
      Initial.values$save.tables <<- as.numeric(tclvalue(save.tables))
      Initial.values$miss.perc <<- as.numeric(tclvalue(miss.perc))
      Initial.values$experiment <<- as.numeric(tclvalue(experiment))
      
      
      tkdestroy(tt)
    }
    
    submit.but <- tkbutton(tt, text="Submit", command=submit,font=fontNam)
    
    tkgrid(tklabel(tt,text="SCP Initial Values", font=fontTite,background="white"),columnspan = 5,pady = 12, padx =10) 
    tkgrid(tklabel(tt,text="Experiment name:", font=fontSub,background="white"),experiment.name.entry, sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="Minimum number of PSMs for Sample:", font=fontSub,background="white"),n.obs.entry, sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="Median SCR", font=fontSub,background="white"),tklabel(tt,text="Low Threshold",background="white"),Low.scr.thr.entry,tklabel(tt,text="High Threshold",background="white"),Top.scr.thr.entry,sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="Booster:", font=fontSub,background="white"),booster.entry, sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="FDR Threshold:", font=fontSub,background="white"),fdr.thr.entry, sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="Median RI", font=fontSub,background="white"),tklabel(tt,text="Low Threshold",background="white"),Low.medianRI.entry,tklabel(tt,text="Top Threshold",background="white"),Top.medianRI.entry,sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="Median CV", font=fontSub,background="white"),tklabel(tt,text="Low Threshold",background="white"),Low.medianCV.entry,tklabel(tt,text="Top Threshold",background="white"),Top.medianCV.entry,sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="How to impute?", font=fontSub,background="white"), tklabel(tt,text="Normal Distribution",background="white"),imputation1,tklabel(tt,text="KNN",background="white"),imputation0,sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="Batch effect", font=fontSub,background="white"), sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="Options to choose from: Channels, Samples or Both. Experiment is also an option.",background="white"),batch.entry,sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="Filter missingness threshold selection?", font=fontSub,background="white"), tklabel(tt,text="From booster",background="white"),miss.perc1,tklabel(tt,text="Manual",background="white"),miss.perc0,sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="Make Extra Run?", font=fontSub,background="white"), tklabel(tt,text="Yes",background="white"),extra.run1,tklabel(tt,text="No",background="white"),extra.run0,sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="Save Plots?", font=fontSub,background="white"), tklabel(tt,text="Yes",background="white"),save.plots1,tklabel(tt,text="No",background="white"),save.plots0,sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="Save Tables?", font=fontSub,background="white"), tklabel(tt,text="Yes",background="white"),save.tables1,tklabel(tt,text="No",background="white"),save.tables0,sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="Multiple Experiments?", font=fontSub,background="white"), tklabel(tt,text="Yes",background="white"),experiment1,tklabel(tt,text="No",background="white"),experiment0,sticky = "w",pady = 10, padx =10)
    tkgrid(tklabel(tt,text="Fasta file for gene extraction:", font=fontSub,background="white"),fastafile.entry,tklabel(tt,text="Description for gene name:", font=fontSub,background="white"),description.entry, sticky = "w",pady = 10, padx =10)
    
    tkgrid(submit.but, clean.but, reset.but,pady = 10, padx =10)
    
    tkgrid(tklabel(tt,text="*Minimum # of PSMs, Median SCR, FDR Threshold, Median RI and Median CV can be left empty, to be decided further on.",font=fontInfo,background="white"), sticky = "w",
           columnspan = 5)
    tkgrid(tklabel(tt,text="*For Multiple experiments, a new column should be provided on the InputFiles file named 'Experiment'.",font=fontInfo,background="white"),
           sticky = "w", columnspan = 5)
    
    
    
    tkwait.window(tt)
    return(Initial.values)
  }
  
  
  
  IntensityDistribution <- function(sci,experiment, InputFiles, title ){
    
    df.abun <- data.frame()
    plts <- list()
    for( i in 1:length(sci)){
      a <- assay(sci[[i]]); a<- data.frame(a)
      a$File.ID <- names(sci)[i]
      if(i != 1){
        colnames(df.abun) <- colnames(a) }
      df.abun <- rbind(df.abun,a)
    }
    df.abun <- expss::na_if(df.abun,0)
    
    if(experiment){
      df.abun$Experiment <- sapply(df.abun$File.ID,function(x) InputFiles[which(InputFiles$File.ID == x),"Experiment"][[1]])
      df.abun$Experiment <- factor(df.abun$Experiment )
      df.abun <- df.abun[order(df.abun$Experiment),]
    }
    
    df.abun.m <- melt(df.abun)
    df.abun.m$variable <- sapply(df.abun.m$variable, function(x) strsplit(as.character(x),split = ".",fixed=T,perl=F)[[1]][-1])
    
    p <- list()
    p[[1]] <- ggplot(df.abun.m,aes(x = log10(value), fill = File.ID)) +
      facet_wrap(~variable, scales = "free_y")+
      geom_histogram(position = "stack", bins=100, stat = StatBin2,size=1)+
      theme_classic() +
      scale_fill_manual(values = mako(length(unique(df.abun.m$File.ID))))+
      xlab("log2 intensities") + ylab("Number of PSMs") + labs(title = "Intensties distribution",
                                                               subtitle = title) +
      theme(plot.subtitle=element_text(face="italic", color="black")) +
      guides( fill = guide_legend(ncol=1))
    
    if(experiment){
      p[[2]]  <- ggplot(df.abun.m,aes(x = log10(value), fill = Experiment)) +
        facet_wrap(~variable,scales = "free_y")+
        geom_histogram(position = "stack", bins=100, stat = StatBin2,size=1)+
        theme_classic() +
        scale_fill_manual(values = mako(length(unique(df.abun.m$Experiment))))+
        xlab("log2 intensities") + ylab("Number of PSMs") + labs(title = "Intensties distribution",
                                                                 subtitle = title) +
        theme(plot.subtitle=element_text(face="italic", color="black")) +
        guides(fill  = guide_legend(ncol =1))
    }
    return(p)
    
  }
  
  
}

### ---- Select files

select.psm <- tk_choose.files(caption = "Select PSM file")
if(endsWith(select.psm,".csv")){
  PSM <- read.csv(select.psm,header=T,sep=",")
}else{
  PSM <- read.csv(select.psm,header=T,sep="\t")
}


select.InputFiles <- tk_choose.files(caption = "Select original InputFiles file")
dir <- dirname(select.InputFiles)
InputFiles <- read.table(select.InputFiles,header=T,sep="\t")
org.InputFiles <- InputFiles

### ---- Preprocess

mess <- tk_messageBox(type= "yesno",message = "Are the files already prepared for SCP analysis? \nIf not, we will first proceed to prepare them. ", caption = "Preprocess")

if(mess == "no"){
  
  ### Prepare PSM
   PSM <- separate_rows(PSM, Master.Protein.Accessions, sep=";")
  PSM$Master.Protein.Accessions <- trimws(PSM$Master.Protein.Accessions)
  
  PSM.g <- PSM %>% group_by(File.ID, Master.Protein.Accessions, Annotated.Sequence,Percolator.PEP) %>% 
    summarise_at(grep("Abundance.",colnames(PSM), value=T), median, na.rm=T) 
  
  PSM.m <- make.proteins.SCP(data.frame(PSM.g))
  PSM.m <- PSM.m[,-grep("127C",colnames(PSM.m))]
  colnames(PSM.m) <- c("File.ID", "Master.Protein.Accessions", "Annotated.Sequence", 
                       "Percolator.PEP", "Abundance.126", "Abundance.127N", "Abundance.128N", 
                       "Abundance.128C", "Abundance.129N", "Abundance.129C", "Abundance.130N", 
                       "Abundance.130C", "Abundance.131N", "Abundance.131C", "Abundance.132N", 
                       "Abundance.132C", "Abundance.133N", "Abundance.133C", "Abundance.134N", 
                       "Abundance.134C", "Abundance.135N")
  PSM <- PSM.m
  write.table(PSM.m, sub(".txt",'_forSCP.txt',select.psm), quote=F, row.names=F, sep= "\t")
  
  
  ### Prepare TMT
  
  TMTtype <- TMTchoose(ini.value = 18)
  
  
  f16 <- sort(rep(InputFiles$File.ID,TMTtype))
  TMTs <- c("Abundance.126","Abundance.127N","Abundance.127C","Abundance.128N","Abundance.128C",
            "Abundance.129N","Abundance.129C","Abundance.130N","Abundance.130C","Abundance.131N","Abundance.131C","Abundance.132N",
            "Abundance.132C","Abundance.133N","Abundance.133C","Abundance.134N", "Abundance.134C", "Abundance.135N")
  t16 <- rep(TMTs[1:TMTtype], length(InputFiles$File.ID))
  InputFiles <- data.frame(File.ID = f16, Channels = t16)
  InputFiles <- InputFiles[-grep("127C",InputFiles$Channels),]
  InputFiles <- unique(InputFiles)
  
  if(any(grepl("Experiment", colnames(org.InputFiles)))){
    InputFiles$Experiment <- c(sapply(InputFiles$File.ID, function(x) org.InputFiles[org.InputFiles$File.ID == x, 'Experiment'] ))
  }
  #a <- merge(InputFiles, new.InputFiles, by = "File.ID")
  
  write.table(InputFiles, sub(".txt",'_forSCP.txt',select.InputFiles), quote=F, row.names=F, sep= "\t")
}

# Dinamic values

init <- InitialValues(
  ins.n.obs = 150, ins.Top.scr.thr = 0.5, ins.Low.scr.thr = 0, ins.booster= "Abundance.126", 
  ins.fdr.thr = 0.01, ins.use.medianCV= "Channels",
  ins.Top.medianRI= 0, ins.Low.medianRI = -3,
  ins.Top.medianCV= 0.7, ins.Low.medianCV = 0.1, ins.batch= "Both",ins.imputation= 1,
  ins.extra.run= T, ins.save.plots= T,
  ins.save.tables= T, ins.experiment =T, ins.experiment.name= "SCP", ins.miss.perc= T,
  ins.fastafile = "F:\\fasta\\HUMAN_03012022.fasta", ins.description = "HUMAN")

n.observations <- init$n.obs
Top.scr.thr <- init$Top.scr.thr
Low.scr.thr <- init$Low.scr.thr
booster <-  init$booster
fdr.thr <- init$fdr.thr
use.medianCV <- init$use.medianCV
Top.medianCV <- init$Top.medianCV
Low.medianCV <- init$Low.medianCV
Top.medianRI <- init$Top.medianRI
Low.medianRI <- init$Low.medianRI
qc.files.remove <- ""
batch <- init$batch
imputation <- ifelse(init$imputation==0,"KNN","r-norm")
extra.run <- init$extra.run
save.plots <- init$save.plots
save.tables <- init$save.tables
miss.perc <- init$miss.perc 
experiment <- init$experiment
experiment.name <- init$experiment.name
fasta.file <- init$fastafile
fasta.description <- init$description

# Static values
pattern.channels = "127N|128|129|130|131|132|133|134|135"
file.name.1 = paste("Proteins_",experiment.name,"_Unfiltered_ProteinLevel.txt",sep="")
file.name.2 = paste("Proteins_",experiment.name,"_FilteredImputed.txt",sep="")
file.name.dred = paste(experiment.name,"_DimRedData.txt",sep="")


# Prepare input files

new.PSM <- PSM
new.InputFiles <- InputFiles
new.PSM$Abundance.TotalMedianIntensity <- NA
new.PSM[,grep("Abun",colnames(new.PSM))] <- expss::na_if(new.PSM[,grep("Abun",colnames(new.PSM))],0)

if(!experiment){
  new.InputFiles$Experiment <- experiment.name
}

for(i in unique(new.PSM$File.ID)){
  temp <- new.PSM[new.PSM$File.ID==i,]
  temp <- colMedians(as.matrix(temp[,grep("Abun",colnames(temp))][-grep("Abundance.126|Abundance.127C",colnames(temp[,grep("Abun",colnames(temp))]))]),na.rm=T)
  new.PSM$Abundance.TotalMedianIntensity[new.PSM$File.ID==i] <- sum(temp, na.rm=T)
  
  new.PSM$Experiment[new.PSM$File.ID==i] <-  new.InputFiles$Experiment[new.InputFiles$File.ID==i][1]
}

unique.InputFiles <- unique(new.InputFiles[,c(1,3)])
unique.InputFiles <- cbind(unique.InputFiles$File.ID,"Abundance.TotalMedianIntensity",unique.InputFiles$Experiment)
colnames(unique.InputFiles) <- colnames(new.InputFiles)
new.InputFiles <- rbind(new.InputFiles, unique.InputFiles)
new.InputFiles <- unique(new.InputFiles)
new.InputFiles <- new.InputFiles[order(new.InputFiles$File.ID),]
new.InputFiles$Experiment[new.InputFiles$Channels == "Abundance.126"] <- "Booster"
new.InputFiles$Experiment[new.InputFiles$Channels == "Abundance.127C"] <- "Empty"



######################
### SCP Analysis
######################

{
  ## -----------------------------------------------------------------------------
  scp <- readSCP(featureData = new.PSM,
                   colData = new.InputFiles,
                   channelCol = "Channels",
                   batchCol = "File.ID",
                   removeEmptyCols = TRUE)

  color.files <- mako(length(unique(new.InputFiles$File.ID)))
  color.tmt <- mako(length(unique(new.InputFiles$Channels)))
  color.experiments <- c("#440154FF",distinctColorPalette(length(unique(new.InputFiles$Experiment))))
  
  ## ----overview-----------------------------------------------------------------
  scp
  
  ## ----zeroIsNA-----------------------------------------------------------------
  scp <- zeroIsNA(scp, i = 1:length(scp))
  print(">>>> zeroIsNA")
  
  ## ----rowDataNames-------------------------------------------------------------
  rowDataNames(scp)
  
  ## ----dims---------------------------------------------------------------------
  print(">>>>>>>>>>>>> Inital dimensions: ")
  print(dims(scp))
  
  if( save.plots){
    print(">>>> Plot: PSM x Sample")
    amount <- melt(data.frame(dims(scp))[1,])
    if(experiment){
      amount$Experiment <- sapply(amount$variable,function(x) new.InputFiles[which(new.InputFiles$File.ID == x),"Experiment"][[1]])
      amount <- amount[order(amount$Experiment),]
      amount$variable <- factor(amount$variable, levels = c(amount$variable),
                                labels = c(amount$variable))
    }else{amount$Experiment <- as.factor(0)}
    
    psm.abun <- new.PSM[,grepl("Abund",colnames(new.PSM)),]
    psm.abun <- expss::na_if(psm.abun,0)
    psm.abun.med <- colMedians(as.matrix(log2(psm.abun)),na.rm=T)
    psm.abun.sd <- colSds(as.matrix(log2(psm.abun)),na.rm=T)
    
    psm.abun <- data.frame(variable=colnames(psm.abun), meds = psm.abun.med, sds = psm.abun.sd)
    
    psm.abun$variable <- sub("Abundance.",'',psm.abun$variable)
    
    
    psm.sum <- new.PSM[,grepl("Abund",colnames(new.PSM)),]
    psm.sum <- expss::na_if(psm.sum,0)
    psm.sum <- (colSums(as.matrix((psm.sum)),na.rm=T)/colSums(as.matrix((psm.sum)),na.rm=T)[1])*100
    psm.sum <- melt(psm.sum)
    
    psm.sum$variable <- sub("Abundance.",'',row.names(psm.sum))
    
    
    psm.plt <- ggplot(amount, aes(x= variable, y= value, color=Experiment,fill= Experiment))+
      geom_histogram(alpha=0.4,  stat = "identity",size=1)+
      theme_classic() +
      scale_fill_manual(values=color.experiments)+ scale_color_manual(values=color.experiments)+
      xlab("Sample") + ylab("Number of PSMs") + labs(title = "Amount of PSMs",
                                                     subtitle = "Initial data with selected threshold \n") +
      theme(plot.subtitle=element_text(face="italic", color="black"),
            axis.text.x=element_text(angle=90,hjust=1))
    if(!experiment){
      psm.plt <- psm.plt + theme(legend.position = "none")
    }
    
    if(is.na(n.observations)){
      
      png(file = paste(dir,'\\1_PSMs.png',sep =''),
          width = 600,
          height =450)
      print(psm.plt + theme(legend.position = 'none'))
      dev.off()
      
      
      psm.values <- check.PSM(paste(dir,'\\1_PSMs.png',sep =''), n.obs= n.observations)
      n.observations <- psm.values$n.obs
      
    }
    
    psm.plt <- psm.plt + geom_hline(yintercept = n.observations, linetype="dashed", size=1) 
    
    
    pdf(file = paste(dir,'\\1_Initial_PSM_x_Sample_',experiment.name,'.pdf',sep =''),
        width = 10,
        height = 10, paper="a4r", onefile = T)
    
    if(experiment){print(psm.plt)
    }else{psm.plt <- psm.plt + theme(legend.position = "none")
    print(psm.plt)}
    
    
    print(ggplot(psm.abun[-nrow(psm.abun),], aes(as.factor(variable), meds, fill = as.factor(variable)))+
            geom_bar(stat = "identity")+
            geom_errorbar(aes(ymin=meds-sds, ymax=meds+sds), width=.2,
                          position=position_dodge(.9)) +
            theme_classic() +
            scale_fill_manual(values = color.tmt)+
            theme_minimal() + theme(legend.position = "none")+
            xlab("TMT") + ylab("Intensity (log2)") + labs(title = "Median abundance by TMT"))
    
    dev.off()
    
    
    
  }
  
  ## ----filter_assays------------------------------------------------------------
  # All of our files have more than 150 observations
  
  keepAssay <- dims(scp)[1, ] > n.observations
  scp <- scp[, , keepAssay]
  
  if( save.plots){
    print(">>>> Plot: PSM x Sample")
    amount <- melt(data.frame(dims(scp))[1,])
    if(experiment){
      amount$Experiment <- sapply(amount$variable,function(x) new.InputFiles[which(new.InputFiles$File.ID == x),"Experiment"][[1]])
      amount <- amount[order(amount$Experiment),]
      amount$variable <- factor(amount$variable, levels = c(amount$variable),
                                labels = c(amount$variable))
      
    }else{amount$Experiment <- as.factor(0)}
    
    
    psm.abun <- new.PSM[,grepl("Abund",colnames(new.PSM)),]
    psm.abun <- expss::na_if(psm.abun,0)
    psm.abun.med <- colMedians(as.matrix(log2(psm.abun)),na.rm=T)
    psm.abun.sd <- colSds(as.matrix(log2(psm.abun)),na.rm=T)
    
    psm.abun <- data.frame(variable=colnames(psm.abun), meds = psm.abun.med, sds = psm.abun.sd)
    
    psm.abun$variable <- sub("Abundance.",'',psm.abun$variable)
    
    psm.sum <- new.PSM[,grepl("Abund",colnames(new.PSM)),]
    psm.sum <- expss::na_if(psm.sum,0)
    psm.sum <- (colSums(as.matrix((psm.sum)),na.rm=T)/colSums(as.matrix((psm.sum)),na.rm=T)[1])*100
    psm.sum <- melt(psm.sum)
    
    psm.sum$variable <- sub("Abundance.",'',row.names(psm.sum))
    
    pdf(file = paste(dir,'\\1_Initial_PSM_x_Sample_Cutoff_',experiment.name,'.pdf',sep =''),
        width = 10,
        height = 10, paper="a4r", onefile = T)
    
    p <- ggplot(amount, aes(x= variable, y= value, color=Experiment,fill= Experiment))+
      geom_histogram(alpha=0.4,  stat = "identity",size=1)+
      #geom_hline(yintercept = n.observations, linetype="dashed", size=1) +
      theme_classic() +
      scale_fill_manual(values=color.experiments)+ scale_color_manual(values=color.experiments)+
      xlab("Sample") + ylab("Number of PSMs") + labs(title = "Amount of PSMs",
                                                     subtitle = "Initial data with selected threshold \n",
                                                     color = "Experiment",fill= "Experiment") +
      theme(plot.subtitle=element_text(face="italic", color="black"),
            axis.text.x=element_text(angle=90,hjust=1))
    if(experiment){print(p)
    }else{p <- p + theme(legend.position = "none")
    print(p)}
    
    
    
    
    print(ggplot(psm.abun[-nrow(psm.abun),], aes(as.factor(variable), meds, fill = as.factor(variable)))+
            geom_bar(stat = "identity")+
            geom_errorbar(aes(ymin=meds-sds, ymax=meds+sds), width=.2,
                          position=position_dodge(.9)) +
            theme_classic() +
            scale_fill_manual(values = color.tmt)+
            theme_minimal() + theme(legend.position = "none")+
            xlab("TMT") + ylab("Intensity (log2)") + labs(title = "Median abundance by TMT"))
    
     
    dev.off()
    
    
    
  }
  
  ## ----Distribution---------------------------------------------------------------
  
  if( extra.run){
    p.iDist <- IntensityDistribution(sci = scp, experiment = experiment, InputFiles = new.InputFiles, title = "Initial distribution")
    
    pdf(file = paste(dir,'\\1_Initial_Distribution_',experiment.name,'.pdf',sep =''),
        width = 10,
        height = 12,  onefile = T)
    
    print(p.iDist[[1]])
    if(experiment){print(p.iDist[[2]])}
    dev.off()
    
  }
  
  ## ----computeSCR---------------------------------------------------------------
  
  #The function creates an field .meanSCR and stores it in the rowData of each assay.
  # We use the booster (Carrier/126)
  
  print(">>>> computing SCR")
  
  scp <- computeSCR(scp,
                    i = 1:length(scp),
                    colvar = "Channels",
                    carrierPattern =  booster ,
                    samplePattern = pattern.channels,
                    sampleFUN = "median",
                    rowDataName = ".medSCR")
  
  
  
  ## ----plot_SCR, warning=FALSE, message=FALSE-----------------------------------
  
  if( save.plots){
    
    scr.df <- rbindRowData(scp, i = 1:length(scp)) %>%
      data.frame 
    scr.df$Experiment <- sapply(scr.df$File.ID,function(x) new.InputFiles[which(new.InputFiles$File.ID == x),"Experiment"][[4]])
    
    scr.df$Median <- median(scr.df$.medSCR,na.rm=T)
    scr.df$Mean <- mean(scr.df$.medSCR,na.rm=T)
    
    p.scr <- ggplot(scr.df,aes(x = .medSCR, fill = Experiment)) +
      geom_histogram(bins=500,size=1, poistion = "identity")+
      geom_vline(aes(xintercept = median(scr.df$.medSCR,na.rm=T), color = 'Median'), lty = 2)+
      geom_vline(aes(xintercept = mean(scr.df$.medSCR,na.rm=T), color = 'Mean'),  lty = 2) +
      geom_vline(aes(xintercept = 1/20, color = 'Expected'),  lty = 2) +
      theme_classic() + xlim(0,2)+
      scale_fill_manual(values=color.experiments)+
      scale_color_manual(name = "Statistics", values = c(Median = "red", Mean = "blue", Expected = "darkgreen"))+
      xlab("Median SCR") + ylab("Number of PSMs") + labs(title = "Median SCR",
                                                         subtitle = paste("Median channels / ",booster," \n")) +
      theme(plot.subtitle=element_text(face="italic", color="black"))  +
      guides(fill= guide_legend(ncol=2,byrow=T),color= guide_legend(ncol=3,byrow=F))
    
    
    if(any(is.na(c(Top.scr.thr,Low.scr.thr)))){
      
      png(file = paste(dir,'\\2_MedianSCR.png',sep =''),
          width = 450,
          height =600)
      print(p.scr + theme(legend.position = 'none'))
      dev.off()
      
      
      scr.values <- check.SCR(paste(dir,'\\2_MedianSCR.png',sep =''), Top.scr.thr,Low.scr.thr)
      Top.scr.thr <- scr.values$Top.scr.thr
      Low.scr.thr <-  scr.values$Low.scr.thr
      
      
    }
    
    
    
    pdf(file = paste(dir,"\\2_MedianSCR_",Low.scr.thr,"_",Top.scr.thr,"_",experiment.name,".pdf",sep =''),
        width = 9,
        height = 9,  onefile = T)
    
    print(ggplot(scr.df,aes(x = .medSCR, fill = File.ID)) +
            geom_histogram(bins=500,size=1, poistion = "identity")+
            geom_vline(aes(xintercept = median(scr.df$.medSCR,na.rm=T), color = 'Median'), lty = 2)+
            geom_vline(aes(xintercept = mean(scr.df$.medSCR,na.rm=T), color = 'Mean'),  lty = 2) +
            geom_vline(aes(xintercept = 1/20, color = 'Expected'),  lty = 2) +
            geom_vline(xintercept =Top.scr.thr, linetype="dashed", size=1) +
            geom_vline(xintercept =Low.scr.thr, linetype="dashed", size=1) +
            theme_classic() + xlim(0,2)+
            scale_fill_manual(values=color.files)+
            scale_color_manual(name = "Statistics", values = c(Median = "red", Mean = "blue", Expected = "darkgreen"))+
            xlab("Median SCR") + ylab("Number of PSMs") + labs(title = "Median SCR",
                                                               subtitle = paste("Median channels / ",booster," \n")) +
            theme(plot.subtitle=element_text(face="italic", color="black"))  +
            guides(fill= guide_legend(ncol=2,byrow=T),color= guide_legend(ncol=3,byrow=F)) 
    )
    if(experiment){
      print(ggplot(scr.df,aes(x = .medSCR, fill = Experiment)) +
              geom_histogram(bins=500,size=1, poistion = "identity")+
              geom_vline(aes(xintercept = median(scr.df$.medSCR,na.rm=T), color = 'Median'), lty = 2)+
              geom_vline(aes(xintercept = mean(scr.df$.medSCR,na.rm=T), color = 'Mean'),  lty = 2) +
              geom_vline(aes(xintercept = 1/20, color = 'Expected'),  lty = 2) +
              geom_vline(xintercept =Top.scr.thr, linetype="dashed", size=1) +
              geom_vline(xintercept =Low.scr.thr, linetype="dashed", size=1) +
              theme_classic() + xlim(0,2)+
              scale_fill_manual(values=color.experiments)+
              scale_color_manual(name = "Statistics", values = c(Median = "red", Mean = "blue", Expected = "darkgreen"))+
              xlab("Median SCR") + ylab("Number of PSMs") + labs(title = "Median SCR",
                                                                 subtitle = paste("Median channels / ",booster," \n")) +
              theme(plot.subtitle=element_text(face="italic", color="black"))  +
              guides(fill= guide_legend(ncol=2,byrow=T),color= guide_legend(ncol=3,byrow=F)) 
      )
    }
    
    dev.off()
    
    pdf(file = paste(dir,"\\2_MedianSCR_Threshold_",Low.scr.thr,"_",Top.scr.thr,"_",experiment.name,".pdf",sep =''),
        width = 9,
        height = 9,  onefile = T)
    
    print(ggplot(scr.df,aes(x = .medSCR, fill = File.ID)) +
            geom_histogram(bins=500,size=1, poistion = "identity")+
            geom_vline(aes(xintercept = median(scr.df$.medSCR,na.rm=T), color = 'Median'), lty = 2)+
            geom_vline(aes(xintercept = mean(scr.df$.medSCR,na.rm=T), color = 'Mean'),  lty = 2) +
            geom_vline(aes(xintercept = 1/20, color = 'Expected'),  lty = 2) +
            geom_vline(xintercept =Top.scr.thr, linetype="dashed", size=1) +
            geom_vline(xintercept =Low.scr.thr, linetype="dashed", size=1) +
            theme_classic() + xlim(Low.scr.thr,Top.scr.thr)+
            scale_fill_manual(values=color.files)+
            scale_color_manual(name = "Statistics", values = c(Median = "red", Mean = "blue", Expected = "darkgreen"))+
            xlab("Median SCR") + ylab("Number of PSMs") + labs(title = "Median SCR",
                                                               subtitle = paste("Median channels / ",booster," \n")) +
            theme(plot.subtitle=element_text(face="italic", color="black"))  +
            guides(fill= guide_legend(ncol=2,byrow=T),color= guide_legend(ncol=3,byrow=F)) 
    )
    if(experiment){
      print(ggplot(scr.df,aes(x = .medSCR, fill = Experiment)) +
              geom_histogram(bins=500,size=1, poistion = "identity")+
              geom_vline(aes(xintercept = median(scr.df$.medSCR,na.rm=T), color = 'Median'), lty = 2)+
              geom_vline(aes(xintercept = mean(scr.df$.medSCR,na.rm=T), color = 'Mean'),  lty = 2) +
              geom_vline(aes(xintercept = 1/20, color = 'Expected'),  lty = 2) +
              geom_vline(xintercept =Top.scr.thr, linetype="dashed", size=1) +
              geom_vline(xintercept =Low.scr.thr, linetype="dashed", size=1) +
              theme_classic() + xlim(Low.scr.thr,Top.scr.thr)+
              scale_fill_manual(values=color.experiments)+
              scale_color_manual(name = "Statistics", values = c(Median = "red", Mean = "blue", Expected = "darkgreen"))+
              xlab("Median SCR") + ylab("Number of PSMs") + labs(title = "Median SCR",
                                                                 subtitle = paste("Median channels / ",booster," \n")) +
              theme(plot.subtitle=element_text(face="italic", color="black"))  +
              guides(fill= guide_legend(ncol=2,byrow=T),color= guide_legend(ncol=3,byrow=F)) 
      )
    }
    
    dev.off()
    
    
    
  }
  
  ## ----filter_SCR---------------------------------------------------------------
  
  scp <- filterFeatures(scp,
                        ~ !is.na(.medSCR) &
                          (.medSCR > Low.scr.thr) & (.medSCR < Top.scr.thr))
  
  
  
  ## ----computeFDR---------------------------------------------------------------
  #Filter out PSMs with high false discovery rate
  
  scp <- pep2qvalue(scp,
                    i = 1:length(scp),
                    PEP = "Percolator.PEP",
                    rowDataName = ".FDR")
  
  
  
  print(">>>> FDR")
  
  if( save.plots){ 
    
    fdr.df <- rbindRowData(scp, i = 1:length(scp)) %>%
      data.frame 
    fdr.df$Experiment <- sapply(fdr.df$File.ID,function(x) new.InputFiles[which(new.InputFiles$File.ID == x),"Experiment"][[2]])
    
    
    p.fdr <- ggplot(fdr.df,aes(x = .FDR, fill = File.ID)) +
      geom_histogram(bins=1000,size=1, poistion = "identity")+
      theme_classic() +
      scale_fill_manual(values=color.files)+
      xlab("FDR") + ylab("Number of PSMs") + labs(title = "False Discovery Rate",
                                                  subtitle = "Adjusted pvalue from PEP\n") +
      theme(plot.subtitle=element_text(face="italic", color="black")) +
      guides(fill= guide_legend(ncol=2,byrow=TRUE))
    
    
    if(is.na(fdr.thr)){
      
      png(file = paste(dir,'\\3_FDR.png',sep =''),
          width = 450,
          height = 600)
      print(p.fdr + theme(legend.position = 'none'))
      dev.off()
      
      
      fdr.values <- check.FDR(paste(dir,'\\3_FDR.png',sep =''), fdr.thr)
      fdr.thr <- fdr.values$fdr.thr
      
    }
    
    p.fdr <- p.fdr +   geom_vline(xintercept =fdr.thr, linetype="dashed", size=1)
    
    pdf(file = paste(dir,'\\3_FDR_',fdr.thr,'_',experiment.name,'.pdf',sep =''),
        width = 8,
        height = 8,  onefile = T)
    
    print(p.fdr)
    if(experiment){
      print(ggplot(fdr.df,aes(x = .FDR, fill = Experiment)) +
              geom_histogram(bins=1000,size=1, poistion = "identity")+
              theme_classic() +
              scale_fill_manual(values=color.experiments)+
              xlab("FDR") + ylab("Number of PSMs") + labs(title = "False Discovery Rate",
                                                          subtitle = "Adjusted pvalue from PEP\n") +
              theme(plot.subtitle=element_text(face="italic", color="black")) +
              guides(fill= guide_legend(ncol=2,byrow=TRUE)) +  
              geom_vline(xintercept =fdr.thr, linetype="dashed", size=1))
      
    }
    
    dev.off()
    
  }
  
  ## ----filter_FDR---------------------------------------------------------------
  scp <- filterFeatures(scp,
                        ~ (.FDR) < fdr.thr)
  
  ## ----divideByReference--------------------------------------------------------
  #Relative reporter ion intensity
  #We here divide all columns (using the regular expression wildcard .) by the reference channel (Reference).
  
  scp <- divideByReference(scp,
                           i = 1:length(scp),
                           colvar = "Channels",
                           samplePattern = ".",
                           refPattern = "Abundance.Median126")
  
  
  print(">>>> divideByReference")
  
  ## ----Distribution---------------------------------------------------------------
  
  if( extra.run){
    p.dbrDist <- IntensityDistribution(sci = scp, experiment = experiment, InputFiles = new.InputFiles, title = "Divide By Reference")
    
    pdf(file = paste(dir,'\\1_DivideByReference_Distribution_',experiment.name,'.pdf',sep =''),
        width = 10,
        height = 12, onefile = T)
    
    print(p.dbrDist[[1]])
    if(experiment){print(p.dbrDist[[2]])}
    dev.off()
    
  }
  
  ## ----aggregate_peptides, message = FALSE--------------------------------------
  scp <- aggregateFeaturesOverAssays(scp,
                                     i = 1:length(scp),
                                     fcol = "Annotated.Sequence",
                                     name = paste0("peptides_", names(scp)),
                                     fun = matrixStats::colMedians, na.rm = TRUE)
  
  
  print(">>>> Agregating peptides")
  ## ----show_agg_peptides--------------------------------------------------------
  print(">>>>>>>>>>>>> Added peptides scp: ")
  print(scp)
  
  ## ----joinAssays---------------------------------------------------------------
  scp <- joinAssays(scp,
                    i = ((length(scp)/2)+1):length(scp),
                    name = "peptides")
  print(">>>> Join peptides")
  ## ----show_join----------------------------------------------------------------
  print(">>>>>>>>>>>>> Joint peptides scp: ")
  
  print(scp)
  if( save.plots){
    print(">>>> PLOT: Peptides x Sample")
    amount <- melt(data.frame(dims(scp))[1,])
    amount <- amount[grep("peptides",amount$variable),]
    amount$variable <- sub("peptides_","", amount$variable)
    amount <- amount[-grep("peptides",amount$variable),]
    
    if(experiment){
      amount$Experiment <- sapply(amount$variable,function(x) InputFiles[which(InputFiles$File.ID == x),"Experiment"][[1]])
      amount <- amount[order(amount$Experiment),]
      amount$variable <- factor(amount$variable, levels = c(amount$variable),
                                labels = c(amount$variable))
      
    }else{amount$Experiment <- as.factor(0)}
    
    
    pdf(file = paste(dir,'\\1_Initial_Peptides_x_Sample_',experiment.name,'.pdf',sep =''),
        width = 10,
        height = 10, paper="a4r", onefile = T)
    
    pep.p <- ggplot(amount, aes(x= variable, y= value, color = as.factor(Experiment), fill = as.factor(Experiment)))+
      geom_histogram(alpha=0.4,  stat = "identity",size=1)+
      #geom_hline(yintercept = mean(amount$value), linetype="dashed", size=1) +
      theme_classic() +
      scale_fill_manual(values=color.experiments)+ scale_color_manual(values=color.experiments)+
      xlab("Sample") + ylab("Number of peptides") + labs(title = "Amount of peptides",
                                                         subtitle = "Average line \n",
                                                         color = "Experiment",fill= "Experiment") +
      theme(plot.subtitle=element_text(face="italic", color="black"),
            axis.text.x=element_text(angle=90,hjust=1))
    if(experiment){print(pep.p)
    }else{pep.p <- pep.p + theme(legend.position = "none")
    print(pep.p)}
    
    
    dev.off()
    
    
    
  }
  
  # Filter single-cells
  ## ----transferColDataToAssay---------------------------------------------------
  colData(scp[["peptides"]])
  #scp <- transferColDataToAssay(scp, "peptides")
  colData(scp[["peptides"]])
  
  ## ----subset_single_cells------------------------------------------------------
  sce <- scp[["peptides"]]
  sce <- sce[, grep(pattern.channels, colnames(sce), fixed=F, perl=T)]
  
  if(extra.run){
    scp126 <- scp
    sce126 <- scp126[["peptides"]]
    sce126 <- sce126[, grep("126", colnames(sce126), fixed=F, perl=T)]
    
  }
  
  ## ----addAssay_filter1---------------------------------------------------------
  scp %>%
    addAssay(y = sce,
             name = "peptides_filter1") %>%
    addAssayLinkOneToOne(from = "peptides",
                         to = "peptides_filter1") ->
    scp
  
  if(extra.run){
    scp126 %>%
      addAssay(y = sce126,
               name = "peptides_filter") %>%
      addAssayLinkOneToOne(from = "peptides",
                           to = "peptides_filter") ->
      scp126
  }
  
  
  ## ---- Compute Median Reporter Intensity ------------------------------------
  
  if( save.plots){
    
    medians <- colSums(assay(scp[["peptides"]]), na.rm = TRUE)
    mcri.df <- getWithColData(scp, "peptides") %>%
      colData %>%
      data.frame
    scp$MedianRI <- log10(medians)
    
    mcri.df$File.ID <- as.factor(sub("Abundance.","",mcri.df$File.ID))
    mcri.df$Channels <- as.factor(sub("Abundance.","",mcri.df$Channels))
    mcri.df$MedianRI <- log10(medians)
    mcri.df <- mcri.df[mcri.df$Channels != "Median126",]
    #mcv.df <- mcv.df[mcv.df$Channels != "126",]
    
    
    
    p.ri.c <- ggplot(mcri.df, aes(x = MedianRI, y = Channels)) +
      geom_point(position = position_jitterdodge(jitter.width=1.4,seed = 42),
                 aes(color=Channels ),
                 alpha = 1) +
      geom_boxplot(notch=TRUE,lwd=0.6,
                   notchwidth = 0.8,alpha= 0.8,fill = "white", width=0.4,outlier.alpha = 0, outlier.shape = NA, outlier.color = NA) +
      theme_clean() +
      xlab(expression(MedianRI~~(log[10])))+ ylab("") + labs(title = "MedianRI per Channels")+
      scale_color_manual(values = color.tmt)+
      scale_fill_manual(values = color.tmt) +
      theme(legend.position = "none",legend.title = element_text(size=10),
            legend.text = element_text(size=8), axis.line.x = element_line(size=2),
            axis.text = element_text(size=12), axis.title = element_text(size=13),
            axis.line.y = element_line(size=2), axis.ticks = element_blank()) + coord_flip()+
      guides(fill = "none", color = guide_legend(nrow=3,byrow=TRUE))
    
    
    if(any(is.na(c(Top.medianRI,Low.medianRI)))){
      
      png(file = paste(dir,'\\4_MedianRI_Channel.png',sep =''),
          width = 700,
          height =500)
      print(p.ri.c + theme(legend.position = 'none'))
      dev.off()
      
      
      mri.values <- check.mRI(filename1 = paste(dir,'\\4_MedianRI_Channel.png',sep =''),
                              Low.medianRI=Low.medianRI,Top.medianRI = Top.medianRI)
      Top.medianRI <- mri.values$Top.medianRI
      Low.medianRI <- mri.values$Low.medianRI
      
    }
    
    
    pdf(file = paste(dir,'\\4_MedianRI_Channels_',experiment.name,'.pdf',sep =''),
        width = 20,
        height = 15, paper="a4r")
    print(ggplot(mcri.df, aes(x = MedianRI, y = factor(Channels))) +
            geom_boxplot(alpha=0.07,aes(fill = Channels), outlier.alpha = 0, outlier.shape = NA, outlier.color = NA) +
            geom_point(data =mcri.df[which(mcri.df$MedianRI<=Top.medianRI & mcri.df$MedianRI>=Low.medianRI),],position = position_jitter(seed = 42),aes(color=File.ID))+
            geom_point(data =mcri.df[which(mcri.df$MedianRI>Top.medianRI | mcri.df$MedianRI<Low.medianRI),], pch = 8,aes(color=File.ID))+
            geom_vline(xintercept = Top.medianRI, linetype="dashed") +
            geom_vline(xintercept = Low.medianRI, linetype="dashed") +
            theme_clean() +
            xlab(expression(MedianRI~~(log[10]))) + ylab("TMT Channel") + labs(title = "MedianRI per Channel \n")+
            scale_color_manual(values =color.files)+
            scale_fill_manual(values = color.files) +
            theme(legend.position = "bottom",legend.title = element_text(size=10),
                  legend.text = element_text(size=8)) + coord_flip()+
            guides(fill = "none", color = guide_legend(nrow=4,byrow=TRUE)))
    if(experiment){
      print(ggplot(mcri.df, aes(x = MedianRI, y = Experiment)) +
              geom_boxplot(alpha=0.07,aes(fill = Experiment), outlier.alpha = 0, outlier.shape = NA, outlier.color = NA) +
              geom_point(data =mcri.df[which(mcri.df$MedianRI<=Top.medianRI & mcri.df$MedianRI>=Low.medianRI),],position = position_jitter(seed = 42),aes(color=File.ID))+
              geom_point(data =mcri.df[which(mcri.df$MedianRI>Top.medianRI | mcri.df$MedianRI<Low.medianRI),], pch = 8,aes(color=File.ID))+
              geom_vline(xintercept = Top.medianRI, linetype="dashed") +
              geom_vline(xintercept = Low.medianRI, linetype="dashed") +  theme_clean() +
              xlab(expression(MedianRI~~(log[10])))+ ylab("") + labs(title = "MedianRI per Experiment")+
              scale_color_manual(values = color.files)+
              scale_fill_manual(values = color.files) +
              theme(legend.position = "none",legend.title = element_text(size=10),
                    legend.text = element_text(size=8), axis.line.x = element_line(size=2),
                    axis.text = element_text(size=12), axis.title = element_text(size=13),
                    axis.line.y = element_line(size=2), axis.ticks = element_blank()) + coord_flip()+
              guides(fill = "none", color = guide_legend(nrow=3,byrow=TRUE)))
      
      
      
      
      
    }
    dev.off()
    
    
    
  }
  
  ## ----computeMedianCV----------------------------------------------------------
  print(">>>>>>>>>>>>> MedianCV: ")
  scp <- medianCVperCell(scp,
                         i = "peptides",
                         groupBy = "Master.Protein.Accessions",
                         nobs = 3, 
                         norm = "div.median",
                         na.rm = TRUE,
                         colDataName = "MedianCVs")
  
  
  ## ----plot_medianCV, message = FALSE, warning = FALSE--------------------------
  
  if( save.plots){
    
    mcv.df <- getWithColData(scp, "peptides") %>%
      colData %>%
      data.frame
    
    
    mcv.df$File.ID <- as.factor(sub("Abundance.","",mcv.df$File.ID))
    mcv.df$Channels <- as.factor(sub("Abundance.","",mcv.df$Channels))
    mcv.df <- mcv.df[mcv.df$Channels != "Median126",]
    
    mcv.df <- na.omit(mcv.df)
    
    
    if(any(is.na(c(Top.medianCV,Low.medianCV)))){
      
      p.mcvC <- print(ggplot(mcv.df, aes(x = MedianCVs, y = factor(Channels))) +
                        geom_boxplot(alpha=0.07,aes(fill = Channels), outlier.alpha = 0, outlier.shape = NA, outlier.color = NA) +
                        geom_point(position = position_jitter(seed = 42),aes(color=File.ID))+
                        theme_clean() +
                        xlab("Median CV")+ ylab("TMT Channel") + labs(title = "Median CV per Channel \n")+
                        scale_color_manual(values = color.files)+
                        scale_fill_manual(values = color.files) +
                        theme(legend.position = "bottom",legend.title = element_text(size=10),
                              legend.text = element_text(size=8)) + coord_flip()+
                        guides(fill = "none", color = guide_legend(nrow=3,byrow=TRUE)))
      
      p.mcvS <- print(ggplot(mcv.df, aes(x = MedianCVs, y = factor(File.ID))) +
                        geom_boxplot(alpha=0.07,aes(fill = File.ID), outlier.alpha = 0) +
                        geom_point(position = position_jitter(seed = 42),aes(color=Channels))+
                        theme_clean() +
                        xlab("Median CV")+ ylab("Samples") + labs(title = "Median CV per Samples \n")+
                        scale_color_manual(values = color.tmt)+
                        scale_fill_manual(values = color.files) +
                        theme(legend.position = "bottom",legend.title = element_text(size=10),
                              legend.text = element_text(size=8),
                              axis.text.x=element_text(angle=90,hjust=1)) + coord_flip()+
                        guides(fill = "none", color = guide_legend(nrow=3,byrow=TRUE)))
      
      
      png(file = paste(dir,'\\4_MedianCV_Channel.png',sep =''),
          width = 700,
          height =500)
      print(p.mcvC + theme(legend.position = 'none'))
      dev.off()
      
      
      png(file = paste(dir,'\\4_MedianCV_Samples.png',sep =''),
          width = 700,
          height =500)
      print(p.mcvS + theme(legend.position = 'none'))
      dev.off()
      
      
      mcv.values <- check.mCVC(filename1 = paste(dir,'\\4_MedianCV_Channel.png',sep =''),
                               filename2 = paste(dir,'\\4_MedianCV_Samples.png',sep =''),
                               Low.medianCV=Low.medianCV,Top.medianCV = Top.medianCV)
      Top.medianCV <- mcv.values$Top.medianCV
      Low.medianCV <- mcv.values$Low.medianCV
      
      
    }
    
    
    
    
    
    pdf(file = paste(dir,'\\4_MedianCV_',Low.medianCV,"_",Top.medianCV,'_Channels_',experiment.name,'.pdf',sep =''),
        width = 20,
        height = 15, paper="a4r")
    
    print(ggplot(mcv.df, aes(x = MedianCVs, y = factor(Channels))) +
            geom_boxplot(alpha=0.07,aes(fill = Channels), outlier.alpha = 0, outlier.shape = NA, outlier.color = NA) +
            geom_point(data =mcv.df[which(mcv.df$MedianCVs<=Top.medianCV & mcv.df$MedianCVs>=Low.medianCV),],position = position_jitter(seed = 42),aes(color=File.ID))+
            geom_point(data =mcv.df[which(mcv.df$MedianCVs>Top.medianCV | mcv.df$MedianCVs<Low.medianCV),], pch = 8,aes(color=File.ID))+
            geom_vline(xintercept = Top.medianCV, linetype="dashed") +
            geom_vline(xintercept = Low.medianCV, linetype="dashed") +
            theme_clean() +
            xlab("Median CV")+ ylab("TMT Channel") + labs(title = "Median CV per Channel \n")+
            scale_color_manual(values = color.files)+
            scale_fill_manual(values = color.files) +
            theme(legend.position = "bottom",legend.title = element_text(size=10),
                  legend.text = element_text(size=8)) + coord_flip()+
            guides(fill = "none", color = guide_legend(nrow=4,byrow=TRUE)))
    
    
    dev.off()
    
    
    
    
    
    pdf(file = paste(dir,'\\4_MedianCV_',Low.medianCV,"_",Top.medianCV,'_Samples_',experiment.name,'.pdf',sep =''),
        width = 20,
        height = 15, paper="a4r")
    
    print(ggplot(mcv.df, aes(x = MedianCVs, y = factor(File.ID))) +
            geom_boxplot(alpha=0.07,aes(fill = File.ID), outlier.alpha = 0) +
            geom_point(data =mcv.df[which(mcv.df$MedianCVs<=Top.medianCV & mcv.df$MedianCVs>=Low.medianCV),],position = position_jitter(seed = 42),aes(color=Channels))+
            geom_point(data =mcv.df[which(mcv.df$MedianCVs>Top.medianCV | mcv.df$MedianCVs<Low.medianCV),],pch = 8,aes(color=Channels))+
            geom_vline(xintercept = Top.medianCV, linetype="dashed") +
            geom_vline(xintercept = Low.medianCV, linetype="dashed") +
            theme_clean() +
            xlab("Median CV")+ ylab("Samples") + labs(title = "Median CV per Samples \n")+
            scale_color_manual(values = color.tmt)+
            scale_fill_manual(values = color.files) +
            theme(legend.position = "bottom",legend.title = element_text(size=10),
                  legend.text = element_text(size=8),
                  axis.text.x=element_text(angle=90,hjust=1)) + coord_flip()+
            guides(fill = "none", color = guide_legend(nrow=3,byrow=TRUE)))
    dev.off()
    
    if(experiment){
      pdf(file = paste(dir,'\\4_MedianCV_',Low.medianCV,"_",Top.medianCV,'_Experiment_',experiment.name,'.pdf',sep =''),
          width = 20,
          height = 15, paper="a4r")
      
      print(ggplot(mcv.df, aes(x = MedianCVs, y = factor(Experiment))) +
              geom_boxplot(alpha=0.07,aes(fill = Experiment), outlier.alpha = 0) +
              geom_point(data =mcv.df[which(mcv.df$MedianCVs<=Top.medianCV & mcv.df$MedianCVs>=Low.medianCV),],position = position_jitter(seed = 42),aes(color=Channels))+
              geom_point(data =mcv.df[which(mcv.df$MedianCVs>Top.medianCV | mcv.df$MedianCVs<Low.medianCV),],pch = 8,aes(color=Channels))+
              geom_vline(xintercept = Top.medianCV, linetype="dashed") +
              geom_vline(xintercept = Low.medianCV, linetype="dashed") +
              theme_clean() +
              xlab("Median CV")+ ylab("Samples") + labs(title = "Median CV per Experiment \n")+
              scale_color_manual(values = color.files)+
              scale_fill_manual(values = color.files) +
              theme(legend.position = "bottom",legend.title = element_text(size=10),
                    legend.text = element_text(size=8),
                    axis.text.x=element_text(angle=90,hjust=1)) + coord_flip()+
              guides(fill = "none", color = guide_legend(nrow=3,byrow=TRUE)))
      dev.off()
    }
    
    
  }
  
  
  ## ----create_filter2-----------------------------------------------------------
  
  keepSample <- (scp$MedianCVs < Top.medianCV) &(scp$MedianCVs > Low.medianCV) &
    (scp$MedianRI < Top.medianRI) &(scp$MedianRI > Low.medianRI) &
    grepl(pattern.channels, scp$Channels, fixed=F, perl=T )
  
  
  
  keepSample[is.na(keepSample)] <- FALSE
  
  
  if(extra.run){
    
    keepSample126 <- grepl("126", colnames(scp126[["peptides_filter"]]), fixed=F, perl=T )
    
    keepSample126[is.na(keepSample126)] <- FALSE
    
    ## ----apply_filter126------------------------------------------------------------
    sce126 <- scp126[["peptides_filter"]]
    sce126 <- sce126[, keepSample126]
    
  }
  
  ## ----apply_filter2------------------------------------------------------------
  sce <- scp[["peptides"]]
  sce <- sce[,keepSample ]
  
  ## ----add_filter2--------------------------------------------------------------
  addAssay(scp,
           y = sce,
           name = "peptides_filter2") %>%
    addAssayLinkOneToOne(from = "peptides",
                         to = "peptides_filter2") ->
    scp
  
  if(extra.run){
    addAssay(scp126,
             y = sce126,
             name = "peptides_filter2") %>%
      addAssayLinkOneToOne(from = "peptides_filter",
                           to = "peptides_filter2") ->
      scp126
    
  }
  
  ## ----normalize_scale----------------------------------------------------------
  
  scp <- sweep(scp, 
               i = "peptides_filter2",
               MARGIN = 2,
               FUN = "/",
               STATS = colMedians(assay(scp[["peptides_filter2"]]), na.rm = TRUE),
               name = "peptides_norm_col")
  ## Divide rows by mean
  scp <- sweep(scp,
               i = "peptides_norm_col",
               MARGIN = 1,
               FUN = "/",
               STATS = rowMeans(assay(scp[["peptides_norm_col"]]),  na.rm = TRUE),
               name = "peptides_norm")
  
  
  
  if(extra.run){
    scp126 %>%
      ## Divide columns by median
      sweep(i = "peptides_filter",
            MARGIN = 2,
            FUN = "/",
            STATS = colMedians(assay(scp126[["peptides_filter"]]),
                               na.rm = TRUE),
            name = "peptides_norm_col") %>%
      ## Divide rows by mean
      sweep(i = "peptides_norm_col",
            MARGIN = 1,
            FUN = "/",
            STATS = rowMeans(assay(.[["peptides_norm_col"]]),
                             na.rm = TRUE),
            name = "peptides_norm") ->
      scp126
  }
  ## ----show_sweep---------------------------------------------------------------
  scp
  
  
  ## ----logTransform-------------------------------------------------------------
  scp <- logTransform(scp,
                      base = 2,
                      i = "peptides_norm",
                      name = "peptides_log")
  
  
  ## ----aggregate_proteins-------------------------------------------------------
  scp <- aggregateFeatures(scp,
                           i = "peptides_log",
                           name = "proteins",
                           fcol = "Master.Protein.Accessions",
                           fun = matrixStats::colMedians, na.rm = TRUE)
  
  if(extra.run){
    
    scp126 <- logTransform(scp126,
                           base = 2,
                           i = "peptides_norm",
                           name = "peptides_log")
    
    ## ----aggregate_proteins-------------------------------------------------------
    scp126 <- aggregateFeatures(scp126,
                                i = "peptides_log",
                                name = "proteins",
                                fcol = "Master.Protein.Accessions",
                                fun = matrixStats::colMedians, na.rm = TRUE)
    
  }
  
  ## ----show_agg_proteins--------------------------------------------------------
  scp
  if(save.plots){
    prot <- aggregateFeaturesOverAssays(scp,
                                        i = grep("peptides_F",melt(data.frame(dims(scp))[1,])$variable),
                                        fcol = "Master.Protein.Accessions",
                                        name = paste0("proteins_", names(scp)[grep("peptides_F",melt(data.frame(dims(scp))[1,])$variable)]),
                                        fun = matrixStats::colMedians, na.rm = TRUE)
    
    amount <- melt(data.frame(dims(prot))[1,])
    amount <- amount[grep("proteins",amount$variable),]
    amount$variable <- sub("proteins_peptides_","", amount$variable)
    amount <- amount[-grep("proteins",amount$variable),]
    if(experiment){
      amount$Experiment <- sapply(amount$variable,function(x) new.InputFiles[which(new.InputFiles$File.ID == x),"Experiment"][[2]])
      amount <- amount[order(amount$Experiment),]
      amount$variable <- factor(amount$variable, levels = c(amount$variable),
                                labels = c(amount$variable))
      
    }else{
      amount$Experiment <- as.factor("0")
    }
    
    
    prots.plt <- ggplot(amount, aes(x= variable, y= value,color=as.factor(Experiment),fill=as.factor(Experiment)))+
      geom_histogram(alpha=0.4,  stat = "identity",size=1)+
      #geom_hline(yintercept = mean(amount$value), linetype="dashed", size=1) +
      theme_classic() +
      scale_color_manual(values=color.experiments)+scale_fill_manual(values=color.experiments)+
      xlab("Sample") + ylab("Number of proteins") + labs(title = "Amount of proteins",
                                                         subtitle = "Average line \n",
                                                         color = "Experiment",fill= "Experiment") +
      theme(plot.subtitle=element_text(face="italic", color="black"),
            axis.text.x=element_text(angle=90,hjust=1))
    
    pdf(file = paste(dir,'\\1_Initial_Proteins_x_Sample_',experiment.name,'.pdf',sep =''),
        width = 10,
        height = 10, paper="a4r", onefile = T)
    
    
    if(experiment){print(prots.plt)
    }else{prots.plt <- prots.plt + theme(legend.position = "none")
    print(prots.plt)}
    
    
    dev.off()
    
  }
  
  ## ----normalize_center---------------------------------------------------------
  scp %>%
    ## Center columns with median
    sweep(i = "proteins",
          MARGIN = 2,
          FUN = "-",
          STATS = colMedians(assay(scp[["proteins"]]),
                             na.rm = TRUE),
          name = "proteins_norm_col") %>%
    ## Center rows with mean
    sweep(i = "proteins_norm_col",
          MARGIN = 1,
          FUN = "-",
          STATS = rowMeans(assay(.[["proteins_norm_col"]]),
                           na.rm = TRUE),
          name = "proteins_norm1") ->
    scp
  
  
  ## Center columns with median
  
  
  if(extra.run){
    scp126 %>%
      ## Center columns with median
      sweep(i = "proteins",
            MARGIN = 2,
            FUN = "-",
            STATS = colMedians(assay(scp126[["proteins"]]),
                               na.rm = TRUE),
            name = "proteins_norm_col") %>%
      ## Center rows with mean
      sweep(i = "proteins_norm_col",
            MARGIN = 1,
            FUN = "-",
            STATS = rowMeans(assay(.[["proteins_norm_col"]]),
                             na.rm = TRUE),
            name = "proteins_norm") ->
      scp126
    
  }
  
  ## --- QC2 -------------------------------------------------------------------
  
  if(save.plots){
    
    qc <- perCellQCMetrics(scp[["proteins_norm1"]],assay.type=2)
    qc$sum <- log2(qc$sum)
    qc <- data.frame(qc, colData(scp[["proteins_norm1"]])[rownames(qc), ])
    qc$File.ID <- as.character(sapply(rownames(qc), function(x) (strsplit(as.character(x),"Abundance.")[[1]][1])))
    qc$Channels <- as.character(sapply(rownames(qc), function(x) (strsplit(as.character(x),"Abundance.")[[1]][2])))
    if(experiment){
      qc$Experiment <- sapply(qc$File.ID,function(x) new.InputFiles[which(new.InputFiles$File.ID == x),"Experiment"][[2]])
      qc <- qc[order(qc$Experiment),]
      qc$File.ID <- factor(qc$File.ID, levels = c(qc$File.ID),
                           labels = c(qc$File.ID))
      
    }else{
      qc$Experiment <- as.factor("0")
    }
    
    text.pos.file <- qc %>% group_by(File.ID) %>% summarise_at(c("sum","detected"), median)
    text.pos.tmt <- qc %>% group_by(Channels) %>% summarise_at(c("sum","detected"), median)
    text.pos.exp <- qc %>% group_by(Experiment) %>% summarise_at(c("sum","detected"), median)
    
    p.qc <- ggplot(qc, aes(x = (sum), y = detected, color = File.ID, label = File.ID)) +
      geom_point(alpha = 0.5)+
      theme_clean() +
      geom_text(data = text.pos.file, aes(x=sum,y=detected,color = File.ID)) +
      xlab("Counts per Cell  (log2)")+ ylab("Valid observations") + labs(title = "Per cell QC metrics \n")+
      scale_color_manual(values = color.files)+
      scale_fill_manual(values = color.files) +
      theme(legend.position = "bottom",legend.title = element_text(size=10),
            legend.text = element_text(size=8)) + 
      guides(fill = "none", color = guide_legend(nrow=3,byrow=TRUE))
    
    png(file = paste(dir,'\\4_QCmetrics_Samples.png',sep =''),
        width = 700,
        height = 500)
    print(p.qc + theme(legend.position = 'none'))
    dev.off()
    
    
    qc.values <- check.qc(filename1 = paste(dir,'\\4_QCmetrics_Samples.png',sep =''),
                          qc.files.remove = qc.files.remove)
    
    if(qc.values == "NA" | length(qc.values) == 0){
      qc.values <- toupper(trimws(strsplit(qc.values$qc.files.remove,split=",")[[1]]))
      qc.values = NA
    }
    
    qc$Keep <- 0
    qc$Keep[(qc$File.ID %in% qc.values)] <- 1
    
    pdf(file = paste(dir,'\\4_QCmetrics_',experiment.name,'.pdf',sep =''),
        width = 10,
        height = 10, paper="a4r", onefile = T)
    
    print(ggplot() +
            geom_point(data = qc, aes(x = (sum), y = detected, color = Channels,shape = as.factor(Keep )),
                       alpha = 0.5,size=2)+
            scale_shape_manual(values = c(19,4))+
            theme_clean() +
            geom_text(data = text.pos.tmt, aes(x=sum,y=detected,color = Channels,label = Channels)) +
            xlab("Counts per Cell  (log2)")+ ylab("Valid observations") + 
            labs(title = "Per cell QC metrics \n")+
            guides(shape = "none") +
            scale_color_manual(values = color.tmt)+
            scale_fill_manual(values = color.tmt) +
            theme(legend.position = "bottom",legend.title = element_text(size=10),
                  legend.text = element_text(size=8)) + 
            guides(fill = "none", color = guide_legend(nrow=3,byrow=TRUE)))
    
    print(ggplot() +
            geom_point(data = qc, aes(x = (sum), y = detected, color = File.ID, shape = as.factor(Keep )),alpha = 0.5)+
            theme_clean() +
            scale_shape_manual(values = c(19,4))+
            geom_text(data = text.pos.file, aes(x=sum,y=detected,color = File.ID, label = File.ID)) +
            xlab("Counts per Cell  (log2)")+ ylab("Valid observations") + labs(title = "Per cell QC metrics \n")+
            scale_color_manual(values = color.files)+
            scale_fill_manual(values = color.files) +
            theme(legend.position = "bottom",legend.title = element_text(size=10),
                  legend.text = element_text(size=8)) + 
            guides(fill = "none", color = guide_legend(nrow=3,byrow=TRUE)))
    
    if(experiment){
      print(ggplot(qc, aes(x = (sum), y = detected, color = Experiment, label = Experiment,shape = as.factor(Keep ))) +
              geom_point(alpha = 0.5,size=2)+
              theme_clean() +
              scale_shape_manual(values = c(19,4))+
              geom_text(data = text.pos.exp, aes(x=sum,y=detected,color = Experiment)) +
              xlab("Counts per Cell (log2)")+ ylab("Valid observations") + labs(title = "Per cell QC metrics \n")+
              scale_color_manual(values = color.experiments)+
              scale_fill_manual(values = color.experiments) +
              theme(legend.position = "bottom",legend.title = element_text(size=10),
                    legend.text = element_text(size=8)) + 
              guides(fill = "none", color = guide_legend(nrow=3,byrow=TRUE)))
    }
    
    dev.off()
    
  }
  
  if(!is.na(qc.values)){
    
    ## ----add_filter3-------------------------------------------------------------
    
    keepSample <- !(unlist(sapply(colnames(scp[["proteins_norm1"]]),function(x) strsplit(x,"Abundance.")[[1]][1])) %in% qc.values)
    
    ## ----apply_filter2------------------------------------------------------------
    sce <- scp[["proteins_norm1"]]
    sce <- sce[,keepSample ]
    
    ## ----add_filter2--------------------------------------------------------------
    addAssay(scp,
             y = sce,
             name = "proteins_norm") %>%
      addAssayLinkOneToOne(from = "proteins_norm1",
                           to = "proteins_norm") ->
      scp
    
  }else{
    sce <- scp[["proteins_norm1"]]
    
    ## ----add_filter2--------------------------------------------------------------
    addAssay(scp,
             y = sce,
             name = "proteins_norm") %>%
      addAssayLinkOneToOne(from = "proteins_norm1",
                           to = "proteins_norm") ->
      scp
  }
  
  ## ----missingness--------------------------------------------------------------
  scp[["proteins_norm"]] %>%
    assay %>%
    is.na %>%
    mean
  
  ## ----Save first file unfiltered--------------------------------------------------------------
  print(">>>>>>>>>>>>> File Unfiltered: ")
  print(scp)
  
  if(save.tables){
    unfiltered.proteins <- data.frame(assay(scp[["proteins_norm"]]))
    unfiltered.proteins$Master.Protein.Accessions <- row.names(unfiltered.proteins)
    unfiltered.proteins <- unfiltered.proteins[,c(ncol(unfiltered.proteins),(1:ncol(unfiltered.proteins)-1))]
    
    unfiltered.proteins <- extract.gene(unfiltered.proteins)
    unfiltered.proteins <- unfiltered.proteins[,c(1:3,ncol(unfiltered.proteins),4:(ncol(unfiltered.proteins)-1))]
    
    write.table(unfiltered.proteins, paste(dir,"\\",file.name.1,sep=""),sep="\t",row.names=F)
    
    if(batch == "" | is.na(batch)){
      Selector <- c("Samples","Channels","Experiment","Both")
      batch <- dlg_list(
        choices = Selector,
        preselect = NULL,
        multiple = FALSE,
        title = "Batch correction group"
      )
    }
    
    if(batch=="Channels"){
      batches <- as.character(sapply(colnames(unfiltered.proteins[,-c(1:4)]), function(x) strsplit(x,".",fixed=T)[[1]][2]))
      unfiltered.proteins[,-c(1:4)] <- removeBatchEffect(unfiltered.proteins[,-c(1:4)], batch=as.factor(batches))
      
    }else if (batch =="Samples"){
      batches <- as.character(sapply(colnames(unfiltered.proteins[,-c(1:4)]), function(x) strsplit(x,"Abundance.",fixed=T)[[1]][1]))
      unfiltered.proteins[,-c(1:4)] <- removeBatchEffect(unfiltered.proteins[,-c(1:4)], batch=as.factor(batches))
      
    }else if(batch == "Experiment"){
      
      batches <- as.character(sapply(colnames(unfiltered.proteins[,-c(1:4)]), function(x) strsplit(x,"Abundance.",fixed=T)[[1]][1]))
      batches <- as.character(sapply(batches,function(x) new.InputFiles$Experiment[new.InputFiles$File.ID==x][1]))
      unfiltered.proteins[,-c(1:4)] <- removeBatchEffect(unfiltered.proteins[,-c(1:4)], batch=as.factor(batches))
      
    }else if(batch == "Both"){
      batches <- as.character(sapply(colnames(unfiltered.proteins[,-c(1:4)]), function(x) strsplit(x,"Abundance.",fixed=T)[[1]][1]))
      
      batches2 <- as.character(sapply(colnames(unfiltered.proteins[,-c(1:4)]), function(x) strsplit(x,".",fixed=T)[[1]][2]))
      
      unfiltered.proteins[,-c(1:4)] <- removeBatchEffect(unfiltered.proteins[,-c(1:4)], batch=as.factor(batches),batch2 = as.factor(batches2 ))
      
    }
    
    write.table(unfiltered.proteins, paste(dir,"\\",paste("Proteins_",experiment.name,"_Unfiltered_ProteinLevel_Batchcorrected.txt",sep=""),sep=""),sep="\t",row.names=F)
    
  }
  
  ## ----FilterNA proteinLevel--------------------------------------------------------------
  #Remove proteins with >30% missingness 126 corresponds to 42 in total
  print(">>>>>>>>>>>>> Missingness computation: ")
  
  if(extra.run){
    scp126n <- scp126
    miss126 <- data.frame(miss126 = c(100,99,95,90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,5,1))
    miss126$pep <- NA;
    for( i in 1:nrow(miss126)){
      if(i == 1){
        miss126$pep[1] <- nrow(scp126n[["proteins_norm"]])
      }else{
        j = (miss126$miss126[i])/100
        scp126n <- filterNA(scp126n,
                            i = "proteins_norm",
                            pNA = j)
        
        miss126$pep[i] <- nrow(scp126n[["proteins_norm"]])
      }
    }
    
    amount.prot <- miss126$pep[miss126$miss126 == 30]
    
    pdf(file = paste(dir,'\\5_Missingness_Boosters_',experiment.name,'.pdf',sep =''),
        width = 20,
        height = 15, paper="a4r")
    
    plot(miss126, xlab="Missingness allowed %", ylab="Number of proteins", pch=4, lwd=2)
    axis(1, at=seq(0,100,by=10),labels=seq(0,100,by=10),las=1, font= 1)
    abline(v = 30, col="azure4", lty=2,lwd=1)
    lines(x = 0:30, y = rep(amount.prot,30+1), col="azure4", lty=2,lwd=1)
    points(miss126, pch=4, lwd=2)
    abline(h=miss126$pep[1]/2, col="darkturquoise", lwd=2)
    abline(h=miss126$pep[1]/3, col="darkturquoise", lty=2,lwd=1.5)
    abline(h=(miss126$pep[1]/3)*2, col="darkturquoise", lty=2,lwd=1.5)
    title(expression(italic('Multi Boosters missing estimation')),adj=0)
    dev.off()
    
    
    ## ----add_filter2--------------------------------------------------------------
    
    scpn <- scp
    miss <- data.frame(miss= c(100,99,95,90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,5,1))
    miss$pep <- NA;
    for( i in 1:nrow(miss)){
      if(i == 1){
        miss$pep[1] <- nrow(scpn[["proteins_norm"]])
      }else{
        j = (miss$miss[i])/100
        scpn <- filterNA(scpn,
                         i = "proteins_norm",
                         pNA = j)
        
        miss$pep[i] <- nrow(scpn[["proteins_norm"]])
      }
    }
    
    
    
    
    pdf(file = paste(dir,'\\5_Missingness_SC_',experiment.name,'.pdf',sep =''),
        width = 20,
        height = 15, paper="a4r")
    
    plot(miss, xlab="Missingness allowed %", ylab="Number of proteins", pch=4, lwd=2)
    axis(1, at=seq(0,100,by=10),labels=seq(0,100,by=10),las=1, font= 1)
    abline(v = miss$miss[which.min(abs(miss$pep - amount.prot)) ], col="azure4", lty=2,lwd=1)
    lines(x = 0:miss$miss[which.min(abs(miss$pep - amount.prot)) ], y = rep(amount.prot,miss$miss[which.min(abs(miss$pep - amount.prot)) ]+1), col="azure4", lty=2,lwd=1)
    points(miss, pch=4, lwd=2)
    abline(h=miss$pep[1]/2, col="darkturquoise", lwd=2)
    abline(h=miss$pep[1]/3, col="darkturquoise", lty=2,lwd=1.5)
    abline(h=(miss$pep[1]/3)*2, col="darkturquoise", lty=2,lwd=1.5)
    title(expression(italic("Single Cells missing estimation")),adj=0)
    dev.off()
    
  }
  
  if((miss.perc == 0) | (!extra.run)){
    
    percNA =  as.numeric(dlgInput("Threshold selected for the allowed missingness \n(ex. Allowance of 70% missingness)",default = "0.7" , Sys.info()["user"])$res)
    
    
  }else{
    
    percNA = miss$miss[which.min(abs(miss$pep - amount.prot)) ]/100
    
  }
  
  scp <- QFeatures::filterNA(scp,
                             i = "proteins_norm",
                             pNA = percNA)
  
  if((percNA != miss$miss[which.min(abs(miss$pep - amount.prot)) ]/100) & (save.plots)){
    
    pdf(file = paste(dir,'\\5_Missingness_SC_',experiment.name,'_RealCutoff.pdf',sep =''),
        width = 20,
        height = 15, paper="a4r")
    
    plot(miss, xlab="Missingness allowed %", ylab="Number of proteins", pch=4, lwd=2)
    axis(1, at=seq(0,100,by=10),labels=seq(0,100,by=10),las=1, font= 1)
    abline(v = miss$miss[which.min(abs(miss$pep - amount.prot)) ], col="azure4", lty=2,lwd=1)
    abline(v = percNA*100, col="darkred", lty=2,lwd=1)
    lines(x = 0:miss$miss[which.min(abs(miss$pep - amount.prot)) ], y = rep(amount.prot,miss$miss[which.min(abs(miss$pep - amount.prot)) ]+1), col="azure4", lty=2,lwd=1)
    points(miss, pch=4, lwd=2)
    abline(h=miss$pep[1]/2, col="darkturquoise", lwd=2)
    abline(h=miss$pep[1]/3, col="darkturquoise", lty=2,lwd=1.5)
    abline(h=(miss$pep[1]/3)*2, col="darkturquoise", lty=2,lwd=1.5)
    title(expression(italic("Single Cells missing estimation")),adj=0)
    dev.off()
    
  }
  
  ## ----impute-------------------------------------------------------------------
  print(">>>>>>>>>>>>> Remaining data Imputation: ")
  
  impute <- data.frame(assay(scp[['proteins_norm']]))
  
  impute <- melt(impute)
  impute$tmt <- sapply(impute$variable, function(x) strsplit(as.character(x),".",fixed=T, perl = F)[[1]][2])
  impute$f <- "Observed"
  
  
  if(imputation=="KNN"){
    scp <- QFeatures::impute(scp,
                             i = "proteins_norm",
                             method = "knn",
                             k = 3, rowmax = 1, colmax = 1,
                             maxp = Inf, rng.seed = 42)
  }else{
    scp <- QFeatures::impute(scp,
                             i = "proteins_norm",
                             method = "FUN", FUN=random_imp)
  }
  
  
  
  if(save.plots){
    
    imputed <- data.frame(assay(scp[['proteins_norm']]))
    
    imputed <- melt(imputed)
    imputed$tmt <- sapply(imputed$variable, function(x) strsplit(as.character(x),".",fixed=T)[[1]][2])
    
    
    imp <- merge(imputed, impute, by = c('variable', 'value','tmt'),all=T)
    imp$f[is.na(imp$f)] <- "Imputed"
    imp <- imp[!is.na(imp$value),]
    
    
    ## Barplot
    
    tb <- table(is.na(impute$value))
    tb <- data.frame("PSMs"= c(as.numeric(tb[1]),as.numeric(tb[2])), "is_na" = c("Observed   ","Imputed   "))
    
    p0 <- ggplot(tb, aes(x=is_na,y=PSMs)) +
      geom_histogram(stat="identity",aes(color = is_na, fill= is_na)) +
      theme_classic2() + xlab("") + ylab("Amount of observations") +
      scale_color_manual(values = c("black","black"))+
      scale_fill_manual(values = c("#87CEFF", "#5D7B93")) +
      theme(legend.position = "none", plot.title = element_text(size=7),
            axis.text.x = element_text(angle=30, vjust=0.8)) + labs(title="")
    
    d <- ggplot(imp[imp$f=="Observed",], aes(x=value)) +
      geom_density(alpha=0.8,aes(fill=factor(f),y=..count..),size=0.7) +
      geom_density(data= imp[imp$f=="Imputed",],aes(x=value,fill=factor(f),y=..count..),alpha=0.7, linetype="dashed" ,size=0.7)+
      theme_classic2() +
      scale_fill_manual(values = c("#87CEFF","#5D7B93"),name = NULL) +
      xlab("")+ ylab("") + labs(title = "Imputation distribution \n") +
      theme(legend.position=c(-0.5,1),
            legend.direction="horizontal",
            axis.title.x = element_blank(),
            plot.title = element_text(size=12),
            axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.line.y = element_blank(),
            plot.margin=unit(c(0.5,1,0.5,0.5), "cm"))
    
    b <-  ggplot(impute, aes(x = value, y = factor(tmt))) +
      geom_boxplot(width = 0.1,alpha = 0.2,aes(fill = tmt),  color ='lightgray', outlier.colour = "lightgray") +
      geom_violin(data =imp,alpha=0,aes(color=tmt))+
      theme_classic2() +
      xlab("Intensities")+ ylab("TMT Channel") +
      scale_color_manual(values = color.tmt,guide="none")+
      scale_fill_manual(values = color.tmt,guide="none")  +
      theme(plot.margin=unit(c(0.5,1,0.5,0.5), "cm"),  plot.title = element_blank())
    
    pdf(file = paste(dir,'\\6_Imputation_',imputation,'_',experiment.name,'.pdf',sep =''),
        width = 20,
        height = 15, paper="a4r")
    
    print(ggarrange(NULL,d, p0, b, 
                    ncol = 2, nrow = 2,  align = "v", 
                    widths = c(1, 2), heights = c(2, 5),
                    common.legend = F, hjust=-0.8))
    
    
    dev.off()
  }
  
  
  ## ----missingness_imputed------------------------------------------------------
  scp[["proteins_norm"]] %>%
    assay %>%
    is.na %>%
    mean
  
  ## ----transferColDataToAssay_proteins------------------------------------------
  #scp <- transferColDataToAssay(scp, i = "proteins_norm")
  sce <- scp[["proteins_norm"]]
  
  # ## ----prepare_batch_correction-------------------------------------------------

  print(">>>>>>>>>>>>> Batch correction: ")
  
  if(batch == "" | is.na(batch)){
    Selector <- c("Samples","Channels","Experiment","Both")
    batch <- dlg_list(
      choices = Selector,
      preselect = NULL,
      multiple = FALSE,
      title = "Batch correction group"
    )
  }
  
  if(batch=="Channels"){
    batches <- as.character(sapply(colnames(assay(sce)), function(x) strsplit(x,".",fixed=T)[[1]][2]))
    assay(sce) <- removeBatchEffect(assay(sce), batch=as.factor(batches))
    
  }else if (batch =="Samples"){
    batches <- as.character(sapply(colnames(assay(sce)), function(x) strsplit(x,"Abundance.",fixed=T)[[1]][1]))
    assay(sce) <- removeBatchEffect(assay(sce), batch=as.factor(batches))
    
  }else if(batch == "Experiment"){
    
    batches <- as.character(sapply(colnames(assay(sce)), function(x) strsplit(x,"Abundance.",fixed=T)[[1]][1]))
    batches <- as.character(sapply(batches,function(x) new.InputFiles$Experiment[new.InputFiles$File.ID==x][1]))
    assay(sce) <- removeBatchEffect(assay(sce), batch=as.factor(batches))
    
  }else if(batch == "Both"){
    batches <- as.character(sapply(colnames(assay(sce)), function(x) strsplit(x,"Abundance.",fixed=T)[[1]][1]))
    
    batches2 <- as.character(sapply(colnames(assay(sce)), function(x) strsplit(x,".",fixed=T)[[1]][2]))
    
    assay(sce) <- removeBatchEffect(assay(sce), batch=as.factor(batches),batch2 = as.factor(batches2 ))
    
  }

  ## ----add_batch_correction-----------------------------------------------------
  addAssay(scp,
           y = sce,
           name = "proteins_batch") %>%
    addAssayLinkOneToOne(from = "proteins_norm",
                         to = "proteins_batch") ->
    scp
  
  ## ----load_scater--------------------------------------------------------------
  
  if(save.plots & ((batch =="Samples") | (batch == "Both"))){
    batch.bf <- assay(scp[["proteins_norm"]])
    batch.bf <- melt(batch.bf)
    
    batch.bf$batched <- sapply(batch.bf[,2], function(x) strsplit(as.character(x),"Abundance",fixed=T)[[1]][1])
    
    
    
    batch.after <-assay(scp[["proteins_batch"]]) 
    batch.after <- melt(batch.after)
    
    batch.after$batched <- sapply(batch.after[,2], function(x) strsplit(as.character(x),"Abundance",fixed=T)[[1]][1])
    
    
    bf <- ggplot(batch.bf, aes(x = value, y = batched)) +
      geom_boxplot(alpha=0.4,aes(fill = batched), outlier.shape = 1, outlier.color = "#5a5a5a",outlier.fill = "white",outlier.alpha = 0.4) +
      theme_classic2() +
      xlab("Intensities")+ ylab("Batch correction grouping \n") + labs(title = "Original distribution \n")+
      scale_color_manual(values =  color.files)+
      scale_fill_manual(values =  color.files) +
      theme(legend.position = "none",plot.margin=unit(c(0.5,0,0.5,0.5), "cm")) 
    
    
    af <- ggplot(batch.after, aes(x = value, y = batched)) +
      geom_boxplot(alpha=0.4,aes(fill = batched),outlier.shape = 1,outlier.color = "#5a5a5a",
                   outlier.fill = "white",outlier.alpha = 0.4) +
      theme_classic2() +
      xlab("Intensities")+ ylab("") + labs(title = "Batch corrected distribution \n")+
      scale_color_manual(values =  color.files)+
      scale_fill_manual(values =  color.files) +
      theme(legend.position = "none",plot.margin=unit(c(0.5,0.5,0.5,0), "cm"),
            axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.line.y = element_blank())   
    
    
    pdf(file = paste(dir,'\\7_Batch_correction_Samples_',experiment.name,'.pdf',sep =''),
        width = 10,
        height = 10, onefile=T)
    print(ggplot(batch.bf, aes(x = value)) +
            facet_grid(rows = batched ~ ., switch  = "y")+
            geom_density(alpha=0.5,aes(fill=factor(batched),color=factor(batched),y=..count..),trim=T,size=0.7)+ 
            geom_density(data = batch.after, alpha=0,aes(fill=NA, y=..count..),color = 'black',trim=T,size=0.5)+ 
            geom_density(data = batch.after, alpha=0,aes(fill=NA, y=..count..),color = 'black',trim=T,size=1, lty = 2)+ 
            theme_classic2() +
            xlab("Intensities")+ ylab("Batch correction grouping \n") + labs(title = "Original vs Batch corrected distribution\n")+
            scale_color_manual(values = color.files)+
            scale_fill_manual(values = color.files) +
            theme(legend.position = "none",plot.margin=unit(c(0.5,0,0.5,0.5), "cm"), 
                  axis.text.y.left = element_blank(), axis.ticks.y = element_blank(),
                  strip.text.y.left = element_text(angle = 0)) 
    )
    print(ggplot(batch.after, aes(y = value, x = factor(batched))) +
            geom_violin()+
            geom_boxplot(alpha=0.07,aes(fill = batched),width=0.2, outlier.alpha = 0, outlier.shape = NA, outlier.color = NA) +
            theme_clean() +
            ylab("Intensities")+ xlab("Batch correction grouping \n") + labs(title = "Batch corrected distribution \n")+
            scale_color_manual(values = color.files)+
            scale_fill_manual(values = color.files) +
            theme(legend.position = "bottom",legend.title = element_text(size=10),
                  legend.text = element_text(size=8)) + coord_flip()+
            guides(fill = "none", color = guide_legend(nrow=3,byrow=TRUE)))
    dev.off()
    
  }
  
  if(save.plots & ((batch =="Channels") | (batch == "Both"))){
    batch.bf <- assay(scp[["proteins_norm"]])
    batch.bf <- melt(batch.bf)
    
    batch.bf$batched <- sapply(batch.bf[,2], function(x) strsplit(as.character(x),".",fixed=T)[[1]][2])
    batch.bf$batched <- as.factor(batch.bf$batched)
    
    
    batch.after <-assay(scp[["proteins_batch"]]) 
    batch.after <- melt(batch.after)
    batch.after$batched <- sapply(batch.after[,2], function(x) strsplit(as.character(x),".",fixed=T)[[1]][2])
    
    
    
    bf <- ggplot(batch.bf, aes(x = value)) +
      geom_density(alpha=0.5,aes(fill=factor(batched),color=factor(batched),y=..count..),trim=T,size=0.7)+ 
      geom_density(data = batch.after, alpha=0,aes(fill=NA, y=..count..),color = 'black',trim=T,size=1, lty = 3)+ 
      facet_grid(rows = batched ~ ., switch  = "y")+
      theme_classic2() +
      xlab("Intensities")+ ylab("Batch correction grouping \n") + labs(title = "Original distribution \n")+
      scale_color_manual(values = color.tmt)+
      scale_fill_manual(values = color.tmt) +
      theme(legend.position = "none",plot.margin=unit(c(0.5,0,0.5,0.5), "cm"), 
            axis.text.y.left = element_blank(), axis.ticks.y = element_blank()) 
    
    
    af <- ggplot(batch.after, aes(x = value, y = batched)) +
      geom_boxplot(alpha=0.4,aes(fill = batched),outlier.shape = 1,outlier.color = "#5a5a5a",
                   outlier.fill = "white",outlier.alpha = 0.4) +
      theme_classic2() +
      xlab("Intensities")+ ylab("") + labs(title = "Batch corrected distribution \n")+
      scale_color_manual(values = color.tmt)+
      scale_fill_manual(values = color.tmt) +
      theme(legend.position = "none",plot.margin=unit(c(0.5,0.5,0.5,0), "cm"),
            axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.line.y = element_blank())   
    
    
    pdf(file = paste(dir,'\\7_Batch_correction_Channels_',experiment.name,'.pdf',sep =''),
        width = 10,
        height = 10, onefile=T)
    print(ggplot(batch.bf, aes(x = value)) +
            facet_grid(rows = batched ~ ., switch  = "y")+
            geom_density(alpha=0.5,aes(fill=factor(batched),color=factor(batched),y=..count..),trim=T,size=0.7)+ 
            geom_density(data = batch.after, alpha=0,aes(fill=NA, y=..count..),color = 'black',trim=T,size=0.5)+ 
            geom_density(data = batch.after, alpha=0,aes(fill=NA, y=..count..),color = 'black',trim=T,size=1, lty = 2)+ 
            theme_classic2() +
            xlab("Intensities")+ ylab("Batch correction grouping \n") + labs(title = "Original vs Batch corrected distribution\n")+
            scale_color_manual(values = color.tmt)+
            scale_fill_manual(values = color.tmt) +
            theme(legend.position = "none",plot.margin=unit(c(0.5,0,0.5,0.5), "cm"), 
                  axis.text.y.left = element_blank(), axis.ticks.y = element_blank(),
                  strip.text.y.left = element_text(angle = 0)) 
    )
    print(ggplot(batch.after, aes(y = value, x = factor(batched))) +
            geom_violin()+
            geom_boxplot(alpha=0.07,aes(fill = batched),width=0.2, outlier.alpha = 0, outlier.shape = NA, outlier.color = NA) +
            theme_clean() +
            ylab("Intensities")+ xlab("Batch correction grouping \n") + labs(title = "Batch corrected distribution \n")+
            scale_color_manual(values = color.tmt)+
            scale_fill_manual(values = color.tmt) +
            theme(legend.position = "bottom",legend.title = element_text(size=10),
                  legend.text = element_text(size=8)) + coord_flip()+
            guides(fill = "none", color = guide_legend(nrow=3,byrow=TRUE)))
    dev.off()
    
  }
  
  if(save.plots & experiment ){
    batch.bf <- assay(scp[["proteins_norm"]])
    batch.bf <- melt(batch.bf)
    
    batches <- as.character(sapply(batch.bf[,2], function(x) strsplit(as.character(x),"Abundance.",fixed=T)[[1]][1]))
    batch.bf$batched  <- as.character(sapply(batches,function(x) new.InputFiles$Experiment[new.InputFiles$File.ID==x][1]))
    
    
    
    batch.after <-assay(scp[["proteins_batch"]]) 
    batch.after <- melt(batch.after)
    
    batches <- as.character(sapply(batch.after[,2], function(x) strsplit(as.character(x),"Abundance.",fixed=T)[[1]][1]))
    batch.after$batched  <- as.character(sapply(batches,function(x) new.InputFiles$Experiment[new.InputFiles$File.ID==x][1]))
    
    
    
    
    bf <- ggplot(batch.bf, aes(x = value, y = batched)) +
      geom_boxplot(alpha=0.4,aes(fill = batched), outlier.shape = 1, outlier.color = "#5a5a5a",outlier.fill = "white",outlier.alpha = 0.4) +
      theme_classic2() +
      xlab("Intensities")+ ylab("Batch correction grouping \n") + labs(title = "Original distribution",
                                                                       subtitle = ifelse((batch =="Experiment"),"Batch corrected by Experiment \n", "Not used for batch correction \n"))+
      scale_color_manual(values = color.experiments)+
      scale_fill_manual(values = color.experiments) +
      theme(legend.position = "none",plot.margin=unit(c(0.5,0,0.5,0.5), "cm"), plot.subtitle =element_text(face="italic", color="black") ) 
    
    
    af <- ggplot(batch.after, aes(x = value, y = batched)) +
      geom_boxplot(alpha=0.4,aes(fill = batched),outlier.shape = 1,outlier.color = "#5a5a5a",
                   outlier.fill = "white",outlier.alpha = 0.4) +
      theme_classic2() +
      xlab("Intensities")+ ylab("") + labs(title = "Batch corrected distribution ", 
                                           subtitle = ifelse((batch =="Experiment"),"Batch corrected by Experiment \n", "Not used for batch correction \n"))+
      scale_color_manual(values = color.experiments)+
      scale_fill_manual(values = color.experiments) +
      theme(legend.position = "none",plot.margin=unit(c(0.5,0.5,0.5,0), "cm"),
            axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.line.y = element_blank(), plot.subtitle =element_text(face="italic", color="black"))   
    
    
    pdf(file = paste(dir,'\\7_Batch_correction_Experiment_',experiment.name,'.pdf',sep =''),
        width = 10,
        height = 10, onefile=T)
    print(ggplot(batch.bf, aes(x = value)) +
            facet_grid(rows = batched ~ ., switch  = "y")+
            geom_density(alpha=0.5,aes(fill=factor(batched),color=factor(batched),y=..count..),trim=T,size=0.7)+ 
            geom_density(data = batch.after, alpha=0,aes(fill=NA, y=..count..),color = 'black',trim=T,size=0.5)+ 
            geom_density(data = batch.after, alpha=0,aes(fill=NA, y=..count..),color = 'black',trim=T,size=1, lty = 2)+ 
            theme_classic2() +
            xlab("Intensities")+ ylab("Batch correction grouping \n") + labs(title = "Original vs Batch corrected distribution\n")+
            scale_color_manual(values = color.experiments)+
            scale_fill_manual(values = color.experiments) +
            theme(legend.position = "none",plot.margin=unit(c(0.5,0,0.5,0.5), "cm"), 
                  axis.text.y.left = element_blank(), axis.ticks.y = element_blank(),
                  strip.text.y.left = element_text(angle = 0)) 
    )
    print(ggplot(batch.after, aes(y = value, x = factor(batched))) +
            geom_violin()+
            geom_boxplot(alpha=0.07,aes(fill = batched),width=0.2, outlier.alpha = 0, outlier.shape = NA, outlier.color = NA) +
            theme_clean() +
            ylab("Intensities")+ xlab("Experiment") + labs(title = "Batch corrected distribution \n")+
            scale_color_manual(values = color.experiments)+
            scale_fill_manual(values = color.experiments) +
            theme(legend.position = "bottom",legend.title = element_text(size=10),
                  legend.text = element_text(size=8)) + coord_flip()+
            guides(fill = "none", color = guide_legend(nrow=3,byrow=TRUE)))
    dev.off()
    
  }
  
  
  ## ---Save_Param----------------------------------------------------------------
  print(">>>>>>>>>>>>> Save parameters: ")
  { parameters <- data.frame("Variable" = c("Date and Time", "dir", "experiment.name","file.name.dred", "file.name.2", "file.name.1", 
                                            "pattern.channels", "n.observations", "Top.scr.thr", "Low.scr.thr", 
                                            "booster", "fdr.thr",  "Top.medianRI", "Low.medianRI","use.medianCV", "Top.medianCV", "Low.medianCV", 
                                            "qc.values", "missing.percentage", "imputation", "batch", "extra.run", "save.plots", 
                                            "save.tables", "experiment", "fasta.file") ,
                             
                             "Value"= c( as.character(Sys.time()), dir, experiment.name,file.name.dred, file.name.2,file.name.1,
                                         pattern.channels,n.observations,Top.scr.thr,Low.scr.thr,booster,fdr.thr,
                                         Top.medianRI,Low.medianRI,use.medianCV,Top.medianCV,Low.medianCV,qc.values, percNA,imputation,
                                         batch,extra.run,save.plots,save.tables,experiment,fasta.file),
                             
                             "Description"= c("Date of creation", "Directory to find all files","Experiment name", "Addition for DimensionReduction Files", 
                                              "Name of the filteres and imputed file", "Name of the unfiltered file", 
                                              "Pattern for TMT channels", "Minimum number data x sample", "Top threshold of median SCR", 
                                              "Low threshold fo median SCR", "Name of the booster column", 
                                              "Threshold for FDR", "Top threshold of median Reporter Intensity", 
                                              "Low threshold of median Reporter Intensity", "Which value used for median CV", "Top threshold of median CV", 
                                              "Low threshold of median CV", "Files removed with QCmetrics","Percentage of allowed missingness", 
                                              "Type of imputation", "What to correct for batch effect", "Do extra analysies", 
                                              "Save the plots", "Save the tables", "Are there multiple experiments", "Fasta file used to Gene annotation"
                             ))
    
    write.table(parameters, paste(dir,"\\ParametersData_",experiment.name,".txt",sep=""),sep="\t",row.names=F)
  }
  
  
  ## ---Add Genenames-------------------------------------------------------------
  print(">>>>>>>>>>>>> Adding gene names: ")
  
  if(save.tables ){
    
    SCP <- data.frame(assay(scp[["proteins_batch"]]))
    SCP$Master.Protein.Accessions <- row.names(SCP)
    SCP <- SCP[,c(ncol(SCP),(1:ncol(SCP)-1))]
    
    
    scp.name <- extract.gene(SCP)
    
    write.table(scp.name, paste(dir,"\\",sub(".txt","_GeneNames.txt",file.name.2),sep=""),sep="\t",row.names=F)
    
  }
}

###############################
### Downstream analysis Sorted
###############################

### ---- Sorted experiment
select.prot <- tk_choose.files(caption = "Select Sorted filtered and imputed file")
prot.sorted <- read.csv(select.prot,header=T,sep="\t")

dir <- dirname(select.prot)

### ---- ANOVA
df.anova <- data.frame(gene = character(), signif = numeric())

for(gene in sort(colnames(prot.sorted)[-ncol(prot.sorted)])){
  df.gene <- prot.sorted[,c(gene,"Group")]
  colnames(df.gene)[1] <- "gene.name"
  res.aov <- aov(gene.name ~ as.factor(Group), data = df.gene)
  if(summary(res.aov)[[1]][["Pr(>F)"]][1] < 0.01){
    temp <- data.frame( gene = gene, signif = summary(res.aov)[[1]][["Pr(>F)"]][1])
    df.anova <- rbind(df.anova, temp)
  }
}

# Extract significant anova
df.significant <- prot.sorted[,c(which(colnames(prot.sorted) %in% df.anova$gene))]
prot.sorted <- data.frame(df.significant)

### ---- PCA
res.pca <- prcomp(prot.sorted,  scale = T)

res.ind <- get_pca_ind(res.pca)
pos <- data.frame(res.ind$coord)
pos$TMT <- sapply(row.names(pos), function(x) strsplit(as.character(x),".",fixed=T, perl = F)[[1]][2])
pos$TMT <- sapply(pos$TMT, function(x) strsplit(as.character(x),"_",fixed=T, perl = F)[[1]][1])
pos$Files <- sapply(row.names(pos), function(x) strsplit(as.character(x),"Abundance.",fixed=T, perl = F)[[1]][1])
pos$CellType <- sapply(row.names(pos), function(x) strsplit(as.character(x),"_",fixed=T, perl = F)[[1]][2])

cols <- distinctColorPalette(length(unique(pos$Files)))
cole <- c("#90CAF9", "#EF9A9A" ,"#FFE082","#C5E1A5","lightgray")

pdf(file = paste(dir,'\\PCA_prcomp_',experiment.name,'.pdf',sep =''),
    width = 8,
    height = 7, onefile = T)
print(ggplot(pos, aes(x=  Dim.1, y = Dim.2, color = TMT)) + 
        geom_point(size=2.5)+
        scale_color_manual(values = cols)+
        labs(title = "PCA", subtitle = "Colored by TMT") + theme_bw() + 
        theme(plot.subtitle = element_text(face = "italic")))
print(ggplot(pos, aes(x=  Dim.1, y = Dim.2, color = Files)) + 
        geom_point(size=2.5)+
        scale_color_manual(values = cols)+
        labs(title = "PCA", subtitle = "Colored by Files") + theme_bw() + 
        theme(plot.subtitle = element_text(face = "italic")))
print(ggplot(pos, aes(x=  Dim.1, y = Dim.2, color = CellType)) + 
        geom_point(size=2.5)+
        scale_color_manual(values = cole)+
        labs(title = "PCA", subtitle = "Colored by CellType") + theme_bw() + 
        theme(plot.subtitle = element_text(face = "italic")))
dev.off()

res.var <- get_pca_var(res.pca)
pos <- data.frame(res.var$coord)

write.table(pos, paste(dir,'\\PCA_res.var_coord_',experiment.name,'.txt',sep =''),sep="\t",row.names=T)


####################################
#### Downstream Analysis All cells
####################################

### ---- Sorted experiment (Only significant proteins ANOVA BH 0.01)
select.proteins <- tk_choose.files(caption = "Select Sorted file")
Prot.Sorted <- read.csv(select.proteins,header=T,sep="\t")

### ---- Unsorted experiment all proteins
select.proteins <- tk_choose.files(caption = "Select Unsorted file")
Prot.Unsorted <- read.csv(select.proteins,header=T,sep="\t")

dir <- dirname(select.proteins)

### ---- Extract significant genes
filt <- Prot.Unsorted[which(Prot.Unsorted$Gene.name %in% unique(Prot.Sorted$Gene.name)),]
rownames(filt) <- filt$Gene.name.Unique
filt <- filt[,grep("Abund",colnames(filt))]

### ---- Rnorm Imputation
filt.impute.norm <- manual_impute(dat = t(filt))
filt.impute.norm <- t(scale(filt.impute.norm))
filt.impute <- as.data.frame(filt.impute.norm)

### ---- Merge Sorted & Unsorted
filt.impute$Gene.name <- rownames(filt.impute)
proteins <- merge(filt.impute, Prot.Sorted, by = "Gene.name", all.x = T, all.y = F)
proteins <- proteins[rowSums(is.na(proteins)) == 0,]
rownames(proteins) <- make.unique(proteins$Gene.name)
proteins <- proteins[,-1]

### ---- Normalise mean
proteins <- apply(proteins,2,function(x) x - mean(colMeans(as.matrix(x))))

### ---- Get chips for batch correction
batches <- as.character(sapply(colnames(proteins), function(x) strsplit(x,"Abundance.",fixed=T)[[1]][1]))
Chip <- batches
Chip[Chip %in% as.character(c(paste("F",1:12,sep="")))] <- "C1"
Chip[Chip %in% as.character(c(paste("F",13:24,sep="")))] <- "C2"
Chip[Chip %in% as.character(c(paste("F",25:36,sep="")))] <- "C3"
Chip[Chip %in% as.character(c(paste("F",37:49,sep="")))] <- "C4"
Chip[Chip %in% as.character(c(paste("F",50:62,sep="")))] <- "C5"
Chip[Chip %in% as.character(c(paste("F",63:75,sep="")))] <- "C6"

### ---- Batch correction
proteins.batch <- data.frame(t(proteins))

proteins.batch$Group <- as.character(sapply(rownames(proteins.batch), function(x) paste(strsplit(as.character(x),"_")[[1]][2],collapse="_")))
proteins.batch$Group[proteins.batch$Group=="Red"] <- "Endoderm"
proteins.batch$Group[proteins.batch$Group=="Blue"] <- "Ectoderm"
proteins.batch$Group[proteins.batch$Group=="Green"] <- "Mesoderm"
model <- model.matrix(~ Group, data = proteins.batch)

proteins.batch <- ComBat(dat = proteins,
                         batch = Chip,
                         mod = model)

proteins.batch <- data.frame(t(proteins.batch))
write.table(proteins.batch, paste(dir, "\\FilteredImputed_Proteins_", experiment.name,".txt",sep=""), quote=F, row.names=T, sep= "\t")

#### ---- PCA
res.pca <- prcomp(proteins.batch,  scale = T)

#### ---- UMAP analysis

cole <- c("#90CAF9", "#EF9A9A" ,"#FFE082","#C5E1A5","#D5DBDB")

dim = 4
n = 100
data <- scale(t(res.ind$coord[,seq(1,dim,by=1)]))
umpa <- scater::calculateUMAP(data,
                              ncomponents = 2,
                              ntop = Inf,
                              scale = TRUE,
                              n_neighbors = n)
umpa <- data.frame(umpa)
pos <- umpa
pos$TMT <- sapply(row.names(pos), function(x) strsplit(as.character(x),".",fixed=T, perl = F)[[1]][2])
pos$TMT <- sapply(pos$TMT, function(x) strsplit(as.character(x),"_",fixed=T, perl = F)[[1]][1])
pos$Files <- sapply(row.names(pos), function(x) strsplit(as.character(x),"Abundance.",fixed=T, perl = F)[[1]][1])
pos$CellType <- sapply(row.names(pos), function(x) strsplit(as.character(x),"_",fixed=T, perl = F)[[1]][2])
pos$CellType[pos$CellType == "Blue"] <- "Ectoderm" 
pos$CellType[pos$CellType == "Green"] <- "Mesoderm" 
pos$CellType[pos$CellType == "Red"] <- "Endoderm" 
pos$Experiment <- "Sorted"
pos$Experiment[pos$Files %in% as.character(c(paste("F",1:36,sep="")))] <- "Unsorted"
pos$Chip <- NA
pos$Chip[pos$Files %in% as.character(c(paste("F",1:12,sep="")))] <- "C1"
pos$Chip[pos$Files %in% as.character(c(paste("F",13:24,sep="")))] <- "C2"
pos$Chip[pos$Files %in% as.character(c(paste("F",25:36,sep="")))] <- "C3"
pos$Chip[pos$Files %in% as.character(c(paste("F",37:49,sep="")))] <- "C4"
pos$Chip[pos$Files %in% as.character(c(paste("F",50:62,sep="")))] <- "C5"
pos$Chip[pos$Files %in% as.character(c(paste("F",63:75,sep="")))] <- "C6"

pdf(file = paste(dir,"\\UMAP_Neigh",n,"_Dim",dim,"_",experiment.name,".pdf",sep=""),
    width = 9.5, height = 8, onefile = T)

print(ggplot(pos, aes(x=  X1, y = X2, color = Files)) + 
        geom_point(size=2.5,alpha=0.7)+
        scale_color_manual(values = distinctColorPalette(76))+
        labs(title = paste("UMAP | Neighbours:", n, "Dimensions:", dim)) + theme_bw() + 
        xlab("Dimension 1")+ylab("Dimension 2")+
        guides(color=guide_legend(nrow=4,byrow=T))+
        theme(plot.subtitle = element_text(face = "italic"), legend.position = "bottom"))

print(ggplot(pos, aes(x=  X1, y = X2, color = Chip)) + 
        geom_point(size=2.5,alpha=0.7)+
        scale_color_manual(values = distinctColorPalette(6))+
        labs(title = paste("UMAP | Neighbours:", n, "Dimensions:", dim)) + theme_bw() + 
        xlab("Dimension 1")+ylab("Dimension 2")+
        guides(fill=guide_legend(nrow=3,byrow=TRUE))+
        theme(plot.subtitle = element_text(face = "italic"), legend.position = "bottom"))

print(ggplot(pos, aes(x=  X1, y = X2, color = Experiment)) + 
        geom_point(size=2.5,alpha=0.7)+
        scale_color_manual(values = c('#F39C12','#5DADE2'))+
        labs(title = paste("UMAP | Neighbours:", n, "Dimensions:", dim)) + theme_bw() + 
        xlab("Dimension 1")+ylab("Dimension 2")+
        guides(fill=guide_legend(nrow=3,byrow=TRUE))+
        theme(plot.subtitle = element_text(face = "italic"), legend.position = "bottom"))

print(ggplot(pos, aes(x=  X1, y = X2, color = CellType, shape = Experiment)) + 
        geom_point(data = pos[pos$CellType == "Unsorted",], aes(x=  X1, y = X2), size=2.5,alpha=0.7)+
        geom_point(data = pos[pos$CellType != "Unsorted",], aes(x=  X1, y = X2), size=2.5,alpha=0.7)+
        scale_color_manual(values = cole)+
        scale_shape_manual(values = c(8,19))+
        labs(title = paste("UMAP | Neighbours:", n, "Dimensions:", dim)) + theme_bw() + 
        xlab("Dimension 1")+ylab("Dimension 2")+
        guides(fill=guide_legend(nrow=3,byrow=TRUE))+
        theme(plot.subtitle = element_text(face = "italic"), legend.position = "bottom"))

dev.off()


write.table(pos, paste(dir,"\\UMAP_Neigh",n,"_Dim",dim,"_",experiment.name,".txt",sep=""), quote=F, row.names=F, sep= "\t")


#################
####  Randomised
#################

### ---- Sorted experiment (Only significant proteins ANOVA BH 0.01)

select.proteins <- tk_choose.files(caption = "Select Sorted ANOVA file")
Prot.Sorted <- read.csv(select.proteins,header=T,sep="\t")

### ---- Unsorted experiment all proteins
select.proteins <- tk_choose.files(caption = "Select Unsorted file")
Prot.Unsorted <- read.csv(select.proteins,header=T,sep="\t")

dir <- dirname(select.proteins)
experiment.name <- "SCP_Gastruloids_Randomised"

### ---- Extract significant genes
filt <- Prot.Unsorted[which(Prot.Unsorted$Gene.name %in% unique(Prot.Sorted$Gene.name)),]
rownames(filt) <- filt$Gene.name.Unique
filt <- filt[,grep("Abund",colnames(filt))]

### ---- Rnorm Imputation
filt.impute.norm <- manual_impute(dat = t(filt))
filt.impute.norm <- t(scale(filt.impute.norm))
filt.impute <- as.data.frame(filt.impute.norm)

### ---- Merge Sorted & Unsorted
filt.impute$Gene.name <- rownames(filt.impute)
proteins <- merge(filt.impute, Prot.Sorted, by = "Gene.name", all.x = T, all.y = F)
proteins <- proteins[rowSums(is.na(proteins)) == 0,]
rownames(proteins) <- make.unique(proteins$Gene.name)
proteins <- proteins[,-1]

### ----  RANDOMISE
intensities <- as.numeric(as.matrix(proteins))
intensities <- intensities[!is.na(intensities)]
for( i in colnames(proteins)){
  proteins[,i] <- sample(intensities, nrow(proteins) )
}
rr <- rownames(proteins)
proteins <- apply(proteins, 2, as.numeric)
rownames(proteins) <- rr

### ---- Normalise mean
proteins <- apply(proteins,2,function(x) x - mean(colMeans(as.matrix(x))))

### ---- Get chips for batch correction
batches <- as.character(sapply(colnames(proteins), function(x) strsplit(x,"Abundance.",fixed=T)[[1]][1]))
Chip <- batches
Chip[Chip %in% as.character(c(paste("F",1:12,sep="")))] <- "C1"
Chip[Chip %in% as.character(c(paste("F",13:24,sep="")))] <- "C2"
Chip[Chip %in% as.character(c(paste("F",25:36,sep="")))] <- "C3"
Chip[Chip %in% as.character(c(paste("F",37:49,sep="")))] <- "C4"
Chip[Chip %in% as.character(c(paste("F",50:62,sep="")))] <- "C5"
Chip[Chip %in% as.character(c(paste("F",63:75,sep="")))] <- "C6"

### ---- Batch correction
proteins.batch <- data.frame(t(proteins))

proteins.batch$Group <- as.character(sapply(rownames(proteins.batch), function(x) paste(strsplit(as.character(x),"_")[[1]][2],collapse="_")))
proteins.batch$Group[proteins.batch$Group=="Red"] <- "Endoderm"
proteins.batch$Group[proteins.batch$Group=="Blue"] <- "Ectoderm"
proteins.batch$Group[proteins.batch$Group=="Green"] <- "Mesoderm"
model <- model.matrix(~ Group, data = proteins.batch)

proteins.batch <- ComBat(dat = proteins,
                         batch = Chip,
                         mod = model)

proteins.batch <- data.frame(t(proteins.batch))
write.table(proteins.batch, paste(dir, "\\FilteredImputed_Proteins_", experiment.name,".txt",sep=""), quote=F, row.names=T, sep= "\t")

### ---- PCA
res.pca <- prcomp(proteins.batch[,-ncol(proteins.batch)],  scale = T)


res.ind <- get_pca_ind(res.pca)
pos <- data.frame(res.ind$coord)
pos$TMT <- sapply(row.names(pos), function(x) strsplit(as.character(x),".",fixed=T, perl = F)[[1]][2])
pos$TMT <- sapply(pos$TMT, function(x) strsplit(as.character(x),"_",fixed=T, perl = F)[[1]][1])
pos$Files <- sapply(row.names(pos), function(x) strsplit(as.character(x),"Abundance.",fixed=T, perl = F)[[1]][1])
pos$CellType <- sapply(row.names(pos), function(x) strsplit(as.character(x),"_",fixed=T, perl = F)[[1]][2])
pos$CellType[pos$CellType == "Blue"] <- "Ectoderm" 
pos$CellType[pos$CellType == "Green"] <- "Mesoderm" 
pos$CellType[pos$CellType == "Red"] <- "Endoderm" 
pos$Experiment <- "Exp56"
pos$Experiment[pos$Files %in% as.character(c(paste("F",1:36,sep="")))] <- "Exp55"
pos$Chip <- NA
pos$Chip[pos$Files %in% as.character(c(paste("F",1:12,sep="")))] <- "C1"
pos$Chip[pos$Files %in% as.character(c(paste("F",13:24,sep="")))] <- "C2"
pos$Chip[pos$Files %in% as.character(c(paste("F",25:36,sep="")))] <- "C3"
pos$Chip[pos$Files %in% as.character(c(paste("F",37:49,sep="")))] <- "C4"
pos$Chip[pos$Files %in% as.character(c(paste("F",50:62,sep="")))] <- "C5"
pos$Chip[pos$Files %in% as.character(c(paste("F",63:75,sep="")))] <- "C6"

cols <- distinctColorPalette(length(unique(pos$Files)))
cole <- c("#90CAF9", "#EF9A9A" ,"#FFE082","#C5E1A5","#D5DBDB")


pdf(file = paste(dir,'\\PCA_prcomp_',experiment.name,'.pdf',sep =''),
    width = 8,
    height = 7, onefile = T)
print(ggplot(pos, aes(x=  Dim.1, y = Dim.2, color = TMT)) + 
        geom_point(size=2.5)+
        scale_color_manual(values = viridis(17))+
        labs(title = "PCA", subtitle = "Colored by TMT") + theme_bw() + 
        theme(plot.subtitle = element_text(face = "italic")))
print(ggplot(pos, aes(x=  Dim.1, y = Dim.2, color = Files)) + 
        geom_point(size=2.5)+
        scale_color_manual(values = cols)+
        labs(title = "PCA", subtitle = "Colored by Files") + theme_bw() + 
        theme(plot.subtitle = element_text(face = "italic")))
print(ggplot(pos, aes(x=  Dim.1, y = Dim.2, color = CellType,shape = Experiment)) + 
        geom_point(size=2.5)+
        scale_color_manual(values = cole)+
        scale_shape_manual(values = c(8,19))+
        labs(title = "PCA", subtitle = "Colored by CellType") + theme_bw() + 
        theme(plot.subtitle = element_text(face = "italic")))
print(ggplot(pos, aes(x=  Dim.1, y = Dim.2, color = Chip)) + 
        geom_point(size=2.5)+
        scale_color_manual(values = distinctColorPalette(6))+
        labs(title = "PCA", subtitle = "Colored by Chip") + theme_bw() + 
        theme(plot.subtitle = element_text(face = "italic")))
print(ggplot(pos, aes(x=  Dim.1, y = Dim.2, color = Experiment)) + 
        geom_point(size=2.5)+
        scale_color_manual(values = cols)+
        labs(title = "PCA", subtitle = "Colored by Experiment") + theme_bw() + 
        theme(plot.subtitle = element_text(face = "italic")))
dev.off()

write.table(pos, paste(dir, "\\PCA_", experiment.name,".txt",sep=""), quote=F, row.names=F, sep= "\t")


### ---- UMAP

dim = 4
n = 100
pdf(file = paste(dir,"\\UMAP_Neigh",n,"_Dim",dim,"_",experiment.name,".pdf",sep=""),
    width = 9.5,
    height = 8, onefile = T)

data <- scale(t(res.ind$coord[,seq(1,dim,by=1)]))
umpa <- scater::calculateUMAP(data,
                              ncomponents = 2,
                              ntop = Inf,
                              scale = TRUE,
                              n_neighbors = n)
umpa <- data.frame(umpa)
pos <- umpa
pos$TMT <- sapply(row.names(pos), function(x) strsplit(as.character(x),".",fixed=T, perl = F)[[1]][2])
pos$TMT <- sapply(pos$TMT, function(x) strsplit(as.character(x),"_",fixed=T, perl = F)[[1]][1])
pos$Files <- sapply(row.names(pos), function(x) strsplit(as.character(x),"Abundance.",fixed=T, perl = F)[[1]][1])
pos$CellType <- sapply(row.names(pos), function(x) strsplit(as.character(x),"_",fixed=T, perl = F)[[1]][2])
pos$CellType[pos$CellType == "Blue"] <- "Ectoderm" 
pos$CellType[pos$CellType == "Green"] <- "Mesoderm" 
pos$CellType[pos$CellType == "Red"] <- "Endoderm" 
pos$Experiment <- "Exp56"
pos$Experiment[pos$Files %in% as.character(c(paste("F",1:36,sep="")))] <- "Exp55"
pos$Chip <- NA
pos$Chip[pos$Files %in% as.character(c(paste("F",1:12,sep="")))] <- "C1"
pos$Chip[pos$Files %in% as.character(c(paste("F",13:24,sep="")))] <- "C2"
pos$Chip[pos$Files %in% as.character(c(paste("F",25:36,sep="")))] <- "C3"
pos$Chip[pos$Files %in% as.character(c(paste("F",37:49,sep="")))] <- "C4"
pos$Chip[pos$Files %in% as.character(c(paste("F",50:62,sep="")))] <- "C5"
pos$Chip[pos$Files %in% as.character(c(paste("F",63:75,sep="")))] <- "C6"

print(ggplot(pos, aes(x=  X1, y = X2, color = Files)) + 
        geom_point(size=2.5,alpha=0.7)+
        scale_color_manual(values = distinctColorPalette(76))+
        labs(title = paste("UMAP | Neighbours:", n, "Dimensions:", dim)) + theme_bw() + 
        xlab("Dimension 1")+ylab("Dimension 2")+
        guides(color=guide_legend(nrow=4,byrow=T))+
        theme(plot.subtitle = element_text(face = "italic"), legend.position = "bottom"))


print(ggplot(pos, aes(x=  X1, y = X2, color = Chip)) + 
        geom_point(size=2.5,alpha=0.7)+
        scale_color_manual(values = distinctColorPalette(6))+
        labs(title = paste("UMAP | Neighbours:", n, "Dimensions:", dim)) + theme_bw() + 
        xlab("Dimension 1")+ylab("Dimension 2")+
        guides(fill=guide_legend(nrow=3,byrow=TRUE))+
        theme(plot.subtitle = element_text(face = "italic"), legend.position = "bottom"))

print(ggplot(pos, aes(x=  X1, y = X2, color = Experiment)) + 
        geom_point(size=2.5,alpha=0.7)+
        scale_color_manual(values = c('#F39C12','#5DADE2'))+
        labs(title = paste("UMAP | Neighbours:", n, "Dimensions:", dim)) + theme_bw() + 
        xlab("Dimension 1")+ylab("Dimension 2")+
        guides(fill=guide_legend(nrow=3,byrow=TRUE))+
        theme(plot.subtitle = element_text(face = "italic"), legend.position = "bottom"))

print(ggplot(pos, aes(x=  X1, y = X2, color = CellType, shape = Experiment)) + 
        geom_point(data = pos[pos$CellType == "Unsorted",], aes(x=  X1, y = X2), size=2.5,alpha=0.7)+
        geom_point(data = pos[pos$CellType != "Unsorted",], aes(x=  X1, y = X2), size=2.5,alpha=0.7)+
        scale_color_manual(values = cole)+
        scale_shape_manual(values = c(8,19))+
        labs(title = paste("UMAP | Neighbours:", n, "Dimensions:", dim)) + theme_bw() + 
        xlab("Dimension 1")+ylab("Dimension 2")+
        guides(fill=guide_legend(nrow=3,byrow=TRUE))+
        theme(plot.subtitle = element_text(face = "italic"), legend.position = "bottom"))

dev.off()


write.table(pos, paste(dir,"\\UMAP_Neigh",n,"_Dim",dim,"_",experiment.name,".txt",sep=""), quote=F, row.names=F, sep= "\t")



