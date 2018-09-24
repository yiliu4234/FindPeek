require(data.table)

#' Find peak areas in location series
#'
#' Find peak or signal areas in location series without noise
#'
#' The first method is to judge the region as a peak according to the continuous rise and fall.
#' The second method is to give a background value (CutMethod = T). If it is higher than the background
#'  value continuously, the region is considered as a peak. The background value can be given manually or
#'  automatically generated according to the maximum value.
#'
#' @param value Numeric vector. Value that you want to find peek
#' @param ConDistCut Integer. Cutoff of continuely point. If CutMethod is TRUE, this argument means how
#'  many continuely point above cutoff.If CutMethod is FALSE, this argument means how many continuely rise point,
#' @param AdmitCut Integer. Cutoff of admitted number of miscontinuely point in continuely point.
#' @param Shape Numeric. Should be between 0 to 1, select how high the peek should be. If Cutmethod is TRUE, then
#'  the CutMethodValue equal to max(value) * Shape. Else
#' @param CutMethod Logical.Use cut method instead of increase method
#' @param CutMethodValue Numeric. Cut method using value, default is NULL, using shape
#' @return List. Including index of peek, upper boundary , lower boundary.
#' @examples
#' test <- data.table(Chr=paste0("chr", rep(1, length=250)),Start=as.integer(seq(1, 100, length.out = 250)), End=as.integer(seq(5, 250, length.out = 250)), Value=c(rep(0.01, 50),dnorm(sort(rnorm(100,0,3)),0,3), dnorm(sort(rnorm(100)))))
#' Peek_dt <- FindPeek(test$Value)
#' Start_End <- data.table(Start = Peek_dt$Peek - Peek_dt$ContDistLeft + 1, End = Peek_dt$Peek + Peek_dt$ContDistRight -1, Peek = Peek_dt$Peek)
#'
#' # generate wig
#' wig_dt <- data.table(Chr= test[Start_End$Start, Chr], Start = test[Start_End$Start,  Start], End = test[Start_End$End, End], Value=test[Start_End$Peek, Value])
#'
#' # remove cross chromosome and ramdom
#' wig_dt <- wig_dt[End > Start,][!Chr %like% "rand" ]
#' @export FindPeek
#' @import data.table
FindPeek <- function(value, ConDistCut = 6, AdmitCut = 2,  Shape = 0.1, CutMethod = T, CutMethodValue = NULL){
  # print(paste0("Start at ",date()))
  if(CutMethod){
    #Cutoff method according to value
    if(is.null(CutMethodValue)){
     ValueCut <-  max(value)*Shape
    }else{
      ValueCut <- CutMethodValue
    }
     Diff <- value >= ValueCut
  }else{
  # increase filter
  Diff <- c(1, diff(value))
  }
  DiffPosi <- Diff > 0
  if(AdmitCut){
    # change Admitcut continue pos to continue
    BreakPoint <- which(diff(DiffPosi)!=0)
    ContDist <- c(Inf,diff(BreakPoint))
    ToAdmitIdx <- ContDist <= AdmitCut
    ToAdmit <- BreakPoint[ToAdmitIdx]
    ToAdmitConDist <- ContDist[ToAdmitIdx]
    if(length(ToAdmit) > 0){
      idx_dt <- data.table(idx=1:length(ToAdmit), Start=ToAdmit - ToAdmitConDist +1, End = ToAdmit)
      idx_pos <- idx_dt[,.(pos=Start:End),by=.(idx, Start)][,.(pos, replacevalue=DiffPosi[Start -1], replacevalue2=Diff[Start - 1]), by=.(idx, Start)]
      DiffPosi[idx_pos$pos] <- idx_pos$replacevalue
      Diff[idx_pos$pos] <- idx_pos$replacevalue2
    }
  }

  DiffPosi[length(DiffPosi)] <- Inf # add Inf  to define end position as break point
  BreakPoint <- which(diff(DiffPosi)!=0)
  if(CutMethod == F){
    ContDistLeft <- c(-Inf, diff(BreakPoint)) # -inf as first point ContDist
    ContDistRight <- c(diff(BreakPoint) + 1, -Inf) # -inf as last point ContDist
    BreakPointFilt <- ContDistLeft > ConDistCut & ContDistRight > ConDistCut &  Diff[BreakPoint] >= 0
    BreakPoint <- BreakPoint[BreakPointFilt]
    ContDistLeft <- ContDistLeft[BreakPointFilt]
    ContDistRight <- ContDistRight[BreakPointFilt]
  }else{
    if(length(BreakPoint)%%2!=0){BreakPoint <- c(BreakPoint, length(value))}
    idx_dt <- data.table(Start = BreakPoint[seq(1,length(BreakPoint), by = 2)], End = BreakPoint[seq(2,length(BreakPoint), by = 2)], idx = 1:(length(BreakPoint)%/%2))
    idx_dt[,BreakPoint:=c(Start:End)[which.max(value[Start:End])], by=idx]
    BreakPoint <- idx_dt$BreakPoint
    ContDistLeft <- idx_dt$BreakPoint - idx_dt$Start
    ContDistRight <- idx_dt$End -idx_dt$BreakPoint + 1
  }
  # print(paste0("Finish first step at ", date()))

  if(CutMethod == F & Shape > 0){
    # shape according to Peek value
    ValueCut <-  max(value)*Shape
    BreakPointFilt <- value[BreakPoint] >= ValueCut
    BreakPoint <- BreakPoint[BreakPointFilt]
    ContDistLeft <- ContDistLeft[BreakPointFilt]
    ContDistRight <- ContDistRight[BreakPointFilt]
  }
  # print(paste0("Finish all step at ", date()))
  return(list(Peek =  BreakPoint, ContDistLeft = ContDistLeft, ContDistRight = ContDistRight))
}


