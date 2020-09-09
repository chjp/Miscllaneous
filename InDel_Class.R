suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg19"))

to_CT <- function(nucleotide){ # convert GA to CT
  if(nucleotide=="G"){
    return("C")
    } else if (nucleotide=="A"){
      return("T")
    }else{
      return(nucleotide)
    }
}

getSeqFrom <- function(chr, start, end) # fetch sequence from UCSC
{
  ret = as.character(getSeq(x=BSgenome.Hsapiens.UCSC.hg19, names=chr, start=start, end=end, strand="+", as.character=TRUE))
  return(invisible(ret))
}

checkRepLen <- function(InDel, next50) # check repeat length in the downstream 50bp
{
  ret = 0
  Len = nchar(InDel)
  for (i in 0:floor(50/Len-1)) {
    j = i*Len + 1
    if (InDel == substr(next50, j, j+Len-1)) {
      ret = ret + 1
    }else{break}
  }
  return(ret)
}

checkHomLen <- function(deleted, next50) # check deletion micro-homology length
{
  ret = 0
  for (i in 1:nchar(deleted)) {
    if (substr(deleted, 1, i) == substr(next50, 1, i)) {
      ret = i
    }
  }
  return(ret)
}

indel = read.table("./exam/question4/inputIndel.txt",
                   sep = "\t", header = T)
indel$ref = as.character(indel$ref)
indel$alt = as.character(indel$alt)
indel$chrom = as.character(indel$chrom)
indel$specimen = as.character(indel$specimen)
indel$anno = 0

for (i in 1:nrow(indel)) {
  chr = indel[i,"chrom"]
  pos = indel[i, "pos"]
  ref = indel[i,"ref"]
  alt = indel[i,"alt"]
  sample = indel[i,"specimen"]
  start = as.numeric(pos)+1
  if( abs(nchar(ref)-nchar(alt)) == 1){ # 1bp InDel
    if (nchar(ref)-nchar(alt)==1){ # 1bp Del
      type = "DEL"
      deleted = substr(ref, 2, 3)
      nextn = getSeqFrom(chr = chr, start = start + 1, end = start + 50) # start+1 cause 1bp deletion
      repLen = checkRepLen(InDel = deleted, next50 = nextn)
      Nt = to_CT(deleted)
      if(repLen>=5){repLen="5+"}
      indel$anno[i] = paste0(type,"_", Nt,"_1", "_", repLen)  
    } else {                      # 1bp Ins
      type = "INS"
      inserted = substr(alt, 2, nchar(alt))
      nextn = getSeqFrom(chr = chr, start = start, end = start + 50)
      repLen = checkRepLen(InDel = inserted, next50 = nextn)
      Nt = to_CT(inserted)
      if(repLen>=5){repLen="5+"}
      indel$anno[i] = paste0(type,"_", Nt,"_1", "_", repLen)  
    }
  }else if (nchar(ref) < nchar(alt) ) { # >=2bp Ins
    type = "INS"
    inserted = substr(alt, 2, nchar(alt))
    InLen=nchar(alt)-1
    if(InLen>=5){InLen="5+"}
    nextn = getSeqFrom(chr = chr, start = start, end = start + 50)
    repLen = checkRepLen(InDel = inserted, next50 = nextn)
    if(repLen>=5){repLen="5+"}
    indel$anno[i] = paste0(type,"_repeats","_",InLen, "_", repLen)  #INS_repeats_5+_3
  }else if ( nchar(ref) > nchar(alt) ) { # >=2bp Del
    type = "DEL"
    deleted = substr(ref, 2, nchar(ref))
    DelLen = nchar(ref)-1
    nextn = getSeqFrom(chr = chr, start = start + DelLen, end = start + DelLen + 50) # start+1 cause 1bp deletion
    repLen = checkRepLen(InDel = deleted, next50 = nextn)
    if(DelLen>=5){ReportLen="5+"}else{ReportLen=DelLen}
    if(repLen>0){  #ignor MH if it is a del-repeat, although upstream may be MH
      if(repLen>=5){repLen="5+"}
      indel$anno[i] = paste0(type,"_repeats_", ReportLen,"_", repLen)
      next
    }
    
    prevn = getSeqFrom(chr = chr, start = start - DelLen, end = start - 1)
    nextn = getSeqFrom(chr = chr, start = start + DelLen, end = start + 2*DelLen - 1)
    MH1=checkHomLen(deleted = deleted, next50 = prevn)
    MH2=checkHomLen(deleted = deleted, next50 = nextn)
    MH=MH1+MH2
    if(MH>=5){MH="5+"}
    if(MH==0){
      indel$anno[i] = paste0(type,"_repeats_", ReportLen,"_0")
    }else{
      indel$anno[i] = paste0(type,"_MH_", ReportLen,"_",MH)
    }
  }
}

InDelSummary = matrix(0, nrow = length(unique(indel$anno)), ncol = length(unique(indel$specimen)))
rownames(InDelSummary) = as.character(sort(unique(indel$anno)))
colnames(InDelSummary) = as.character(unique(indel$specimen))

for (i in 1:nrow(indel)){
  samp = indel$specimen[i]
  class = indel$anno[i]
  InDelSummary[class,samp] = InDelSummary[class,samp] + 1
}

write.table(InDelSummary, "./InDel.tsv",
            sep = "\t", quote = F)

