
# df prep
dfprep<-function(df1, focalcol,phylo_covar_mat,spatial_covar_mat, Fcols01){
  df1<-df1[!is.na(df1[,focalcol]),]
  colnames(df1)[colnames(df1)%in%focalcol]<-"focalcol"
  df1<-df1[,colnames(df1)%in%c("v0.2.socName",Fcols01, "focalcol")]
  df1<-df1[complete.cases(df1),]
  
  df1$focalcol<-as.factor(df1$focalcol)
  
  # scale
  if("AnimalHusbandry"%in%Fcols01) df1$AnimalHusbandry<-scale(df1$AnimalHusbandry)[,1]
  if("FarmingPropensity"%in%Fcols01) df1$FarmingPropensity<-scale(df1$FarmingPropensity)[,1]
  if("MeanAnnualTemperature"%in%Fcols01) df1$MeanAnnualTemperature<-scale(df1$MeanAnnualTemperature)[,1]
  if("NumberParagraphs"%in%Fcols01) df1$NumberParagraphs<-scale(df1$NumberParagraphs)[,1]

  df1->df_extended
  
  # include matrices
  tokeep<-intersect(rownames(phylo_covar_mat), df1$v0.2.socName) 
  df1<-df1[df1$v0.2.socName%in%tokeep,]
  df1<-df1[match(tokeep, df1$v0.2.socName),]
  spatial_covar_mat1<-spatial_covar_mat[tokeep,tokeep]
  phylo_covar_mat1<-phylo_covar_mat[tokeep,tokeep]
  
  toret<-list(df1, df_extended, phylo_covar_mat1, spatial_covar_mat1)
  names(toret)<-c("df1","df_extended" ,"phylo_covar_mat1", "spatial_covar_mat1")
  return(toret)  
}