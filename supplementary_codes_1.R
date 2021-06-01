#'##############################################################################
#' 3-segmented leaf curve
#'############################################################################## 


leaf.curve = function(.endKc, # vector containing end kc values for init, midle and end
                      .lenghPhase, # vector contanign lenght of each phenological stage
                      .namePhase=NULL,
                      .dae,
                      .digits=3,
                      .stage=NULL
){
  
  
  stage.by.dae <- function(.dae=NULL, .breaks=NULL, .names=NULL){
    if(is.null(.dae)) stop(".dae is missing")
    if(is.null(.breaks)) .breaks <-seq(from=1-min(.dae),to=max(.dae)+10,by=10)
    if(is.null(.names))  .names  <-paste0("Interval_",.breaks)
    #  .breaks <- c(.breaks,Inf)
    pstage = cut(x = .dae,breaks=.breaks,right = T)
    levels(pstage) = .names
    return(pstage)
  }
  
  
  #'======================================================================================#
  #' model assumptions:
  #' based on FAO-56 recommendations to create kc curves using interpolation
  #'======================================================================================#
  #' 1. Input parameters check
  #'======================================================================================#
  #' if lenghPhase is missing, the function will crash
  if(is.null(.lenghPhase)){stop("missing lenghPhase vector")}
  .lenghPhase = c(0,.lenghPhase)
  #'--------------------------------------------------------------------------------------#
  #' if dae is missing, the function will crash
  if(is.null(.dae)){}
  if(is.null(.namePhase)){.namePhase=c("initial","crop development","midle season","late season")}
  #'--------------------------------------------------------------------------------------#
  #' if .stage is null, .stage is computed using stage.by.dae function
  if(is.null(.stage)){
    .stage = stage.by.dae(.dae = .dae,.breaks = .lenghPhase,.names = .namePhase)}
  
  #'======================================================================================#
  #' 2. Computing Kc curve
  #'======================================================================================#
  Kc =c(.endKc[1],.endKc[1],.endKc[2],.endKc[2],.endKc[3],.endKc[3])
  Lp = .lenghPhase
  df = data.frame(kc=1,stage = .stage)
  
  suppressWarnings(for(i in 1:length(.namePhase)){
    
    if(i == 1){df$kc[df$stage %in% .namePhase[1]] = Kc[1]}
    if(i == 3){df$kc[df$stage %in% .namePhase[3]] = Kc[3]}
    if(i == 2 | i == 4){
      .rate =(Kc[i+1]-Kc[i])/(Lp[i+1]-Lp[i])
      df$kc[df$stage %in% .namePhase[i]] = round(seq(from = Kc[i], 
                                                     to = Kc[i+1],by=.rate),.digits)}
  })
  df$stage = as.character(df$stage)
  df$stage[is.na(df$stage)] = "No Crop"
  return(data.frame(dae=.dae, df))}


#'###############################################################################
# apparent photosynthetic radiation intercepted by the canopy (aPAR)
#'###############################################################################

aPAR = function(Srad, # solar radiation
                fraction=.5, # fraction of solar radition assumed as PAR
                K=.5,       # canopy light extinction coefficient
                LAI         # leaf area index
                ){
  PAR = Srad*fraction
  return(PAR*(1-exp(-K*LAI)))}



panel_EnvAssembly=function(ECs,bottom.size=0.3,left.size=2.5,order.row=T,order.col=T){
  if (!requireNamespace('superheat')) utils::install.packages("superheat")
  if (!requireNamespace('viridis')) utils::install.packages("viridis")
  superheat(ECs,
            pretty.order.rows = order.row,
            pretty.order.cols = order.col,
            grid.vline.col = "transparent",
            grid.hline.col = "transparent",
            #row.dendrogram = T,
            legend.width = 4,
            left.label.size = 0.3,
            left.label.text.size=left.size,
            bottom.label.text.size = bottom.size,
            bottom.label.size = 0.3,
            bottom.label.text.angle = 90,
            #heat.pal = viridis::inferno(100),
            heat.pal = viridis::magma(100),
            legend.text.size = 10,
            #   X.text = round(as.matrix(a),1),X.text.col="white",
            legend.height=0.08
  )
  
}
