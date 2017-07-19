############################################################################
# LM  model function to run models and extract stats for logk or logAlpha  #
# 								           #
# ms:       statements of models in the form that lm() accepts		   #
# dataset:  which dataset to use					   #
# varNames: row names of Regions					   #
############################################################################

regDataLM <- function(ms,dataset,varNames) {
        nn<-length(varNames)
        model.results<-lapply(ms, function(x) {
        foo<-summary(lm(as.formula(x), data=dataset))
        # takes stats for first covariate, specifically
        return(foo$coefficients[2,])
       	}) 
	# convert results to data frame
        gamOut1<-data.frame(matrix(unlist(model.results), nrow=nn, byrow=T))
        names(gamOut1) <-c("coef","se","tval","p")
        row.names(gamOut1)<-varNames
        # multiple comparison testing
        gamOut1[,5]<-p.adjust(gamOut1[,4], method = "fdr", n = nn)
        gamOut1[,6]<-ms
        names(gamOut1)[5:6]<-c("fdr.p","model")
        return(gamOut1)
}


