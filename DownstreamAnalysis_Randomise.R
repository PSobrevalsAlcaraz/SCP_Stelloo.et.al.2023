#################
####  Randomised
#################


library(factoextra)
library(randomcoloR)
library("RegressionLibs")
library(viridis)
library(dplyr)
library(tidyr)
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




########################
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


############# MAKE UMAPs

cole <- c("#90CAF9", "#EF9A9A" ,"#FFE082","#C5E1A5","#D5DBDB")

dim=4
n=100
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

