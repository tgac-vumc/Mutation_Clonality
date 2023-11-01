# Code edited from https://rdrr.io/bioc/Clonality/src/R/get.mutation.frequencies.r
Calculate_TCGA_frequencies <- function(xmut.ids,tcga.cancer.type){
    load(url("https://github.com/IOstrovnaya/MutFreq/blob/master/freqdata.RData?raw=true"))
    if(tcga.cancer.type == 'NSCLC'){
        tcga.cancer.type <- c('LUAD','LUSC')
    }
    f<-freqdata[match(xmut.ids,rownames(freqdata)), tcga.cancer.type]    
    if(!is.null(dim(f))){
        f <- rowSums(f)
        names(f)<-xmut.ids
        f[is.na(f)]<-0
        freq<-(f+1)/(sum(freqdata[1, tcga.cancer.type]+1))
    } else{
        names(f)<-xmut.ids
        f[is.na(f)]<-0
        freq<-(f+1)/(freqdata[1, tcga.cancer.type]+1)
    }
    names(freq)<-rownames(mutation_matrix)
    freq
}
