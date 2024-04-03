/mnt/c/Users/hrzon/Desktop/STUDIA/SEM8/MBI/MBI-Bioinformatics-course/4. Analiza danych sekwencyjnych człowieka


library(data.table)
library(parallel)
library(RCurl)
library(gdata)
library(matrixStats)
library(DNAcopy)
library(GenomicRanges)
library(Rsubread)
library(WES.1KG.WUGSC)
library(CODEX)
# set working directory to workDir
workDir <- "/mnt/c/Users/hrzon/Desktop/STUDIA/SEM8/MBI/MBI-Bioinformatics-course/4. Analiza danych sekwencyjnych człowieka/"
setwd(workDir)
#set number of available cores
cores <- 4


dirPath <- system.file("extdata", package = "WES.1KG.WUGSC")
bamFile <- list.files(dirPath, pattern = '*.bam$')
bamdir <- file.path(dirPath, bamFile)
sampname <- as.matrix(read.table(file.path(dirPath, "sampname")))
bedFile <- file.path(dirPath, "chr22_400_to_500.bed")
chr <- 22
bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile,
sampname = sampname, projectname = "CODEX_demo", chr)
bamdir <- bambedObj$bamdir; sampname <- bambedObj$sampname
ref <- bambedObj$ref; projectname <- bambedObj$projectname; chr <- bambedObj$chr
coverageObj <- getcoverage(bambedObj, mapqthres = 20)
Y <- coverageObj$Y; readlength <- coverageObj$readlength

Y_ac <- apply(Y, 2,function(x){(100*x)/width(ref)})
colnames(Y_ac) <- sampname

summary(apply(Y_ac, 2, median))

4.


projectname <- "TGP99"
outputDir <- paste0(workDir,"codex_output/")
dir.create(outputDir)

cfiles <- dir(paste0(workDir, "data/coverage/"), "*bam.coverage*")
cdf <- rbindlist(mclapply(cfiles, function(f){
print(f);
df <- fread(paste0(workDir, "data/coverage/",f));
df$SampleName <- strsplit(f, "\\.")[[1]][1];df
},mc.cores=cores))
colnames(cdf) <- c("Chr", "Start", "Stop", "ReadCount", "SampleName")
cdf <- cdf[order(cdf$SampleName, cdf$Chr, cdf$Start, cdf$Stop),]
dim(cdf)
#[1] 465498 5
head(cdf)
# Chr Start Stop ReadCount SampleName
#1: 20 68319 68439 103 NA06985
#2: 20 76611 77091 129 NA06985
#3: 20 123208 123358 69 NA06985
#4: 20 125995 126389 105 NA06985
#5: 20 138119 138269 37 NA06985
#6: 20 139359 139719 156 NA06985
bedFile <- paste0(workDir, "data/bed/20130108.exome.targets.bed")
sampname <- unique(cdf$SampleName)
chr <- "20"
targetsChr <- cdf[which(cdf$Chr==chr & cdf$SampleName == cdf$SampleName[1]),
c("Chr", "Start", "Stop")]
selChr <- cdf[which(cdf$Chr==chr),]
Y <- t(do.call(rbind,lapply(sampname,
function(s){selChr$ReadCount [which(selChr$SampleName == s)]})))
colnames(Y) <- sampname
rownames(Y) <- 1:nrow(Y)
dim(Y)
#[1] 4702 99
dim(targetsChr)
#[1] 4702 3
ref <- IRanges(start = targetsChr$Start, end = targetsChr$Stop)
gc <- getgc(chr, ref)
mapp <- getmapp(chr, ref)

mapp_thresh <- 0.9 # remove exons with mapability < 0.9
cov_thresh_from <- 20 # remove exons covered by less than 20 reads
cov_thresh_to <- 4000 # remove exons covered by more than 4000 reads
length_thresh_from <- 20 # remove exons of size < 20
length_thresh_to <- 2000 # remove exons of size > 2000
gc_thresh_from <- 20
gc_thresh_to <- 80 
qcObj <- qc(Y, sampname, chr, ref, mapp, gc,
cov_thresh = c(cov_thresh_from, cov_thresh_to),
length_thresh = c(length_thresh_from, length_thresh_to),
mapp_thresh = mapp_thresh,
gc_thresh = c(gc_thresh_from, gc_thresh_to))
Y_qc <- qcObj$Y_qc; sampname_qc <- qcObj$sampname_qc; gc_qc <- qcObj$gc_qc
mapp_qc <- qcObj$mapp_qc; ref_qc <- qcObj$ref_qc; qcmat <- qcObj$qcmat


normObj <- normalize(Y_qc, gc_qc, K = 1:9)
Yhat <- normObj$Yhat; AIC <- normObj$AIC; BIC <- normObj$BIC
RSS <- normObj$RSS; K <- normObj$K
optK <- choiceofK(AIC, BIC, RSS, K,
filename = paste(projectname, "_", chr, "_choiceofK", ".pdf", sep = ""))
finalcall <- CODEX::segment(Y_qc, Yhat, optK = optK,
K = K, sampname_qc,
ref_qc, chr, lmax = 200,
mode = "integer")
finalcall <- data.frame(finalcall, stringsAsFactors=F)
finalcall$targetCount <- as.numeric(finalcall$ed_exon) - as.numeric(finalcall$st_exon)




plotCall <- function(calls, i, Y_qc, Yhat_opt){
startIdx <- as.numeric(calls$st_exon[i])
stopIdx <- as.numeric(calls$ed_exon[i])
sampleName <- calls$sample_name[i]
wd <- 20
startPos <- max(1,(startIdx-wd))
stopPos <- min((stopIdx+wd), nrow(Y_qc))
selQC <- Y_qc[startPos:stopPos,]
selQC[selQC ==0] <- 0.00001
selYhat <- Yhat_opt[startPos:stopPos,]
png(file="cnv.png")
matplot(matrix(rep(startPos:stopPos, ncol(selQC)),
ncol=ncol(selQC)), log(selQC/selYhat,2),
type="l",lty=1, col="dimgrey", lwd=1,
xlab="exon nr", ylab="logratio(Y/Yhat)")
lines(startPos:stopPos,log( selQC[,sampleName]/ selYhat[,sampleName],2), lwd=3, col="red")
dev.off()
}

cnvId <- 4 # indeks zmiany dla której zostanie sporządzony wykres
plotCall(finalcall, cnvId, Y_qc, Yhat[[optK]])