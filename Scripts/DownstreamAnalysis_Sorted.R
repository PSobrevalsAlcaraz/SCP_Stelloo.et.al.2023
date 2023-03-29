#################################
####  Downstream Analysis Sorted
#################################

library(factoextra)
library(randomcoloR)
library(viridis)
library("RegressionLibs")
library(tcltk2)
library(sva)

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

### ---- Sorted experiment all proteins
select.prot <- tk_choose.files(caption = "Select Sorted filtered and imputed file")
prot.sorted <- read.csv(select.prot,header=T,sep="\t")

dir <- dirname(select.prot)
experiment.name <- "SCP_Sorted"

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

########################
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
