###### Functions used in the BFAST01 turning points exercise
### March 31st, 2020
### Paulo Negri Bernardino

bfast01classify <- function(object, alpha=0.05, pct_stable=NULL, typology=c("standard", "drylands")) { 
  ## output array 
  out <- rep(NA,8)
  names(out) <- c("flag_type","flag_significance","p_segment1","p_segment2",
                  "pct_segment1","pct_segment2","flag_pct_stable","flag_subtype")
  
  ## Segment and break point parameters
  object.zoo <- as.zoo(object) # data series
  ## Determine regression object (take first class from model )
  reg = class(object$model[[2]])[1]
  ## if break, list segment and break point parameters (p$..)
  if(object$breaks != 0) {
    ToB <- as.numeric(object$breakpoints[[1]])  # time of break
    s1 <- object$model[[2]]$coefficients[3] # slope segment 1
    s2 <- object$model[[2]]$coefficients[4] # slope segment 2
    m <- as.numeric(object.zoo$trend[ToB+1]) - as.numeric(object.zoo$trend[ToB]) # magnitude of abrupt change
  } 
  
  ## ANOVA and PCTCHANGE      
  for (segment in 1:(object$breaks+1)) {
    # subset zoo object for segment
    date.start <- if(segment==1)  object$data$time[1] else object$data$time[ToB+1]
    date.end <- if(segment==2 || object$breaks==0) object$data$time[nrow(object$data)] else object$data$time[ToB]
    object.zoo.subset <- window(object.zoo, start=date.start, end=date.end)
    # Anova
    segment.anova <- anova(lm((object.zoo.subset$response-object.zoo.subset$season)~time(object.zoo.subset))) 
    # linear model of deseasonalized trend versus time
    out[segment+2] <- segment.anova$Pr[1]
    # PctChange
    obs.start <- if(segment==1)  1 else ToB+1
    obs.end <- if(segment==2 || object$breaks==0) nrow(object$data) else ToB
    if(object.zoo$trend[[obs.end]] / object.zoo$trend[[obs.start]] > 0){
      segment.pctchange <- 
        ( (object.zoo$trend[[obs.end]] / object.zoo$trend[[obs.start]])^(1/(date.end-date.start)) -1) * 100      
    } else {
      if(object.zoo$trend[[obs.start]] < object.zoo$trend[[obs.end]]){
        value.start <- object.zoo$trend[[obs.start]] + 2 * abs(object.zoo$trend[[obs.start]])
        value.end <- object.zoo$trend[[obs.end]] + 2 * abs(object.zoo$trend[[obs.start]])
        segment.pctchange <- ( (value.end / value.start)^(1/(date.end-date.start)) -1) * 100               
      } else {
        value.start <- object.zoo$trend[[obs.start]] + 2 * abs(object.zoo$trend[[obs.end]])
        value.end <- object.zoo$trend[[obs.end]] + 2 * abs(object.zoo$trend[[obs.end]])
        segment.pctchange <- ( (value.end / value.start)^(1/(date.end-date.start)) -1) * 100                           
      }
    }
    out[segment+4] <- segment.pctchange
  }
  
  ## CLASSIFICATION 
  ## Standard Typo
  if(typology=="standard"){
    ## monotonic if no break
    if(object$breaks == 0) {     
      slope <- object$model[[1]]$coefficients[2] # slope
      if(slope > 0) out[1] <- 1
      if(slope < 0) out[1] <- 2
    }else{
      ## classes with break
      # with break, but still monotonic
      if(s1 > 0 && s2 > 0 && m > 0) out[1] <- 3
      if(s1 < 0 && s2 < 0 && m < 0) out[1] <- 4
      # interrupted gradual change (setback or boost)
      if(s1 > 0 && s2 > 0 && m < 0) out[1] <- 5
      if(s1 < 0 && s2 < 0 && m > 0) out[1] <- 6
      # trend reversal (greening to browning v.v.)
      if(s1 > 0 && s2 < 0) out[1] <- 7
      if(s1 < 0 && s2 > 0) out[1] <- 8
    }
  }
  
  # Drylands Typo
  if(typology=="drylands"){
    ## monotonic if no break
    if(object$breaks == 0) {
      if(out[3] > alpha){
        out[1] <- 0 #fluctuating/no change
      }else{
        slope <- object$model[[1]]$coefficients[2] # slope
        if(slope > 0) out[1] <- 1 # stable increase
        if(slope < 0) out[1] <- 2 # stable decrease  
      }
    } else {
      ## classes with break
      if(out[3] > alpha && out[4] > alpha){
        out[1] <- 0 #fluctuating/no change
      }else{
        # with break, but still same direction
        if(s1 > 0 && s2 > 0) out[1] <- 3 # interrupted increase
        if(s1 < 0 && s2 < 0) out[1] <- 4 # interrupted decrease
        # trend reversal (greening to browning v.v.)
        if(s1 > 0 && s2 < 0) out[1] <- 5 # negative reversal
        if(s1 < 0 && s2 > 0) out[1] <- 6 # positive reversal 
      }
    }
    ## Sub-classification
    if(object$breaks != 0 && out[1] != 0) { #if there's a break and at least one sign. trend
      if( (s1 > 0 && s2 > 0) || (s1 < 0 && s2 < 0)){
        # non mono trends (sub-types 1 and 2, "slowing down" and "accelerating")
        if( (s1 > 0 && s2 > 0) || (s1 < 0 && s2 < 0) ){
          if(abs(s1) > abs(s2)) out[8] <- 1 # slowing down
          if(abs(s1) < abs(s2)) out[8] <- 2 # accelerating
        }
      }else if( (s1 < 0 && s2 > 0) || (s1 > 0 && s2 < 0)){
        # reversal trend (sub-types 3 and 4, "transition" and "complete")
        if( (out[3] <= alpha && out[4] > alpha) || (out[3] > alpha && out[4] <= alpha) ) out[8] <- 3 # transition
        if(out[3] <= alpha && out[4] <= alpha) out[8] <- 4 # complete
      }
    }
  }
  
  ## Segment significance flag
  # code: 0 = both segments significant (or no break and significant), 
  # 1 = only first segment significant, 
  # 2 = only 2nd segment significant, 
  # 3 = both segments insignificant (or no break and not significant)
  
  # no break
  if(object$breaks == 0) {     
    if(out[3] <= alpha) out[2] <- 0
    if(out[3] > alpha) out[2] <- 3
    # with break
  } else {
    if(out[3] <= alpha && out[4] <= alpha) out[2] <- 0
    if(out[3] <= alpha && out[4] > alpha) out[2] <- 1
    if(out[3] > alpha && out[4] <= alpha) out[2] <- 2
    if(out[3] > alpha && out[4] > alpha) out[2] <- 3
  }   
  
  ## Segment stability flag
  # code: 0 = both segments beyond stable (or no break and not stable), 
  # 1 = only first segment beyond stable, 
  # 2 = only 2nd segment beyond stable, 
  # 3 = both segments stable (or no break and stable)
  
  if(!is.null(pct_stable)) {
    # no break
    if(object$breaks == 0) {     
      if(abs(out[5]) > pct_stable) out[7] <- 0
      if(abs(out[5]) <= pct_stable) out[7] <- 3
      # with break
    } else {
      if(abs(out[5]) > pct_stable && abs(out[6]) > pct_stable) out[7] <- 0
      if(abs(out[5]) > pct_stable && abs(out[6]) <= pct_stable) out[7] <- 1
      if(abs(out[5]) <= pct_stable && abs(out[6]) > pct_stable) out[7] <- 2
      if(abs(out[5]) <= pct_stable && abs(out[6]) <= pct_stable) out[7] <- 3
    }        
  } else {
    out[7] <- NA
  }
  
  ## Organizing how the output should be displayed
  if(typology=="standard"){
    return(as.data.frame(t(out))[1:7])
  }else if(typology=="drylands"){
    return(as.data.frame(t(out[c(1,8,2:7)])))
  }
}