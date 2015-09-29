###########################################################################
# CoxRegression.GP.ez.R (ver.1)                            May 20, 2008
#
#      Cox regression using clinical data and .cls file
#        Yujin Hoshida (Broad Institute)
###########################################################################
parseCmdLine <- function(...)
{
	input.clin.data.filename= ""
    input.cls.filename="NA"
    output.file=""
    response.variable="NA"
    variable.continuous="NA"
    variable.category="NA"
    variable.interaction.terms="NA"
    input.subgroup="NA"
    variable.selection="none"

	args <- list(...)
	for(i in 1:length(args)) {
		flag <- substring(args[[i]], 0, 4)
		value <- substring(args[[i]], 5, nchar(args[[i]]))

		if(value == "")
        {
            next;
        }
		if(flag=='-inp') {
			input.surv.data.filename = value
		} else if(flag=='-cls') {
			input.cls.filename = value
		} else if(flag=='-out') {
			output.file = value
		} else if(flag=='-tim') {
          			time = value
        }else if(flag=='-sta') {
			status = value
		} else if(flag=='-con') {
			variable.continuous <- value
		} else if(flag=='-cat') {
			variable.category <- value
		} else if(flag=='-int') {
			variable.interaction.terms <- value
		} else if(flag=='-str') {
			strata <- value
		} else if(flag=='-sub') {
          			input.subgroup <- value
        }else if(flag=='-sel') {
			variable.selection <- value
        }else{
            stop(paste("unknown flag ", flag, " value ", value, sep=""), .call=FALSE)
        }
	}

	CoxRegression.GP.ez(input.surv.data.filename=input.surv.data.filename,
                            input.cls.filename=input.cls.filename,
                            output.file=output.file,
                            time=time,
                            status=status,
                            variable.continuous=variable.continuous,
                            variable.category=variable.category,
                            variable.interaction.terms=variable.interaction.terms,
                            strata=strata,
                            input.subgroup=input.subgroup,
                            variable.selection=variable.selection
    )
}

CoxRegression.GP.ez <- function(
  input.surv.data.filename,
  input.cls.filename="NA",

  output.file="CoxRegression_result",

  time="time",
  status="status",

  variable.continuous="NA",
  variable.category="NA",
  variable.interaction.terms="NA",

  strata="NA",
  input.subgroup="NA", #    e.g. "variable:2" only 1 subgroup can be selected

  variable.selection="none"    # "both","forward","backward"

  )
{

  # Advanced parameters
  
  tie.data.handling="efron"  # "breslow", "exact"
  robust.variance="F" # jackknife estimate of the variance. Cannot be selected when variable selection is performed

  test.proportional.hazardness="no"  # no=no;yeswoplot=yes (wo/ plot);yeswithplot=yes (w/ plot), Proportinal hazards tests using Schoenfeld residuals
  
  plot.residuals="no"  # no=no;martingale=Martingale;deviance=Deviance;score=Score;schoenfeld=Schoenfeld;dfbeta=dfbeta;dfbetas=dfbetas;scaledsch=Scaledsch;partial=Partial
  plot.to="NA"     # variable assigned to x-axis. variable not included in the model
  residuals.file="residuals.png"

  # read input files

  surv.data<<-read.surv.data(input.surv.data.filename)
  if (input.cls.filename!="NA"){
    cls<-read.cls(input.cls.filename)
    cls.label<-read.cls.label(input.cls.filename)

    num.cls<-length(table(cls))

    if(length(surv.data[,1])!=length(cls)){
      stop("### Survival data and class file don't match! ###")
    }
    surv.data<<-cbind(surv.data,cls)
  }

  # load packages

  library(splines)
  library(survival)
  library(MASS)
  
  # variable conversion

  if (variable.continuous!="NA"){
    variable.continuous<-gsub(",","+",variable.continuous)
  }else{
    variable.continuous<-""
  }
  if (variable.category!="NA"){
    variable.category<-convert.category(variable.category)
  }else{
    variable.category<-""
  }
  if (variable.interaction.terms!="NA"){
    variable.interaction.terms<-gsub(",","+",variable.interaction.terms)
  }else{
    variable.interaction.terms<-""
  }
  if (strata!="NA"){
    strata<-paste("strata(",strata,")",sep="")
  }else{
    strata<-""
  }

  combined.variables<-paste(variable.continuous,variable.category,variable.interaction.terms,strata,sep="+")
  combined.variables<-gsub("\\++","+",combined.variables)
  combined.variables<-gsub("\\++","+",combined.variables)
  combined.variables<-sub("\\+$","",combined.variables)
  combined.variables<-sub("^\\+","",combined.variables)

  # Cox fit

  if (input.subgroup=="NA"){
    coxph.command<-paste("coxph(Surv(",time,",",status,")~",combined.variables,",data=surv.data,method=\"",tie.data.handling,"\",robust=",robust.variance,")",sep="")
  }else{
    subgroup<-unlist(strsplit(input.subgroup,":"))
    subgroup.variable<-subgroup[1]
    subgroup.member<-subgroup[2]

    coxph.command<-paste("coxph(Surv(",time,",",status,")~",combined.variables,",data=surv.data,method=\"",tie.data.handling,"\",robust=",robust.variance,",subset=",subgroup.variable,"==",subgroup.member,")",sep="")
  }

  coxph.fit<-eval(parse(text=coxph.command))
  fit.vector<-unlist(coxph.fit)

  # variable selection

  if (variable.selection!="none"){
    if (robust.variance=="T"){
      stop("### Robust variance cannot be used in variable selection! ###")
    }
    coxph.fit<-stepAIC(coxph.fit,direction=variable.selection,trace=0)
    fit.vector<-unlist(coxph.fit)
    # null model?
    if (names(fit.vector[1])=="loglik"){
      stop("### Null model! No variable is selected! ###")
    }  # end null model?
  }

  # number of plots

  num.plot<-0
  for (n in 1:length(fit.vector)){
    if (regexpr("^coefficients.",names(fit.vector[n]))!=-1){
      num.plot<-num.plot+1
    }
  }

  if (num.plot==1){
    num.row.col.plot<-1
  }else{
    num.row.col.plot<-as.integer(sqrt(num.plot)+1)
  }

  # Text output

  output.file <- paste(output.file,".txt",sep="")
  
  if (tie.data.handling=="efron"){
    method<-"tie data handling:  Efron approximation"
  }
  if (tie.data.handling=="breslow"){
    method<-"tie data handling:  Breslow approximation"
  }
  if (tie.data.handling=="exact"){
    method<-"tie data handling:  Exact method"
  }

  if (variable.selection=="none"){
    selection<-"variable selection:  none"
  }
  if (variable.selection=="both"){
    selection<-"variable selection:  Stepwise selection"
  }
  if (variable.selection=="forward"){
    selection<-"variable selection:  Forward selection"
  }
  if (variable.selection=="backward"){
    selection<-"variable selection:  Backward selection"
  }

  write.table("# Cox regression results",output.file,quote=F,sep="\t",row.names=F,col.names=F)
  write.table("",output.file,quote=F,sep="\t",row.names=F,col.names=F,append=T)
  write.table(method,output.file,quote=F,sep="\t",row.names=F,col.names=F,append=T)
  write.table(selection,output.file,quote=F,sep="\t",row.names=F,col.names=F,append=T)
  write.table("",output.file,quote=F,sep="\t",row.names=F,col.names=F,append=T)

  write.table(capture.output(summary(coxph.fit)),output.file,quote=F,sep="\t",row.names=F,col.names=F,append=T)

  # Test proportinal hazard assumption

  if (test.proportional.hazardness=="yeswoplot"){
    zph.result<-cox.zph(coxph.fit,global=T)
    write.table("",output.file,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    write.table("# Test for proportional hazard assumption",output.file,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    write.table("",output.file,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    write.table(capture.output(zph.result),output.file,quote=F,sep="\t",row.names=F,col.names=F,append=T)
  }
  if (test.proportional.hazardness=="yeswithplot"){
    zph.result<-cox.zph(coxph.fit,global=T)
    write.table("",output.file,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    write.table("# Test for proportional hazard assumption",output.file,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    write.table("",output.file,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    write.table(capture.output(zph.result),output.file,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    
    if (capabilities("png")==T){
      png("TestProportionalHazardness.png")
      par(mfrow=c(num.row.col.plot,num.row.col.plot))
      plot(zph.result)
      dev.off()
    }
  }

  # Plot residuals

  if (plot.residuals!="no" & capabilities("png")==T){
    residuals<-residuals(coxph.fit,type=plot.residuals)
    if (is.vector(residuals)){
      num.resid<-1
      variable.name<-"residuals"
    }else{
      num.resid<-length(residuals[1,])
      variable.name<-colnames(residuals)
    }
    if (input.subgroup!="NA"){
      x.variable<-eval(parse(text=paste("surv.data$",plot.to,"[surv.data$",subgroup.variable,"==\"",subgroup.member,"\"]",sep="")))
    }
    if (input.subgroup=="NA"){
      x.variable<-eval(parse(text=paste("surv.data$",plot.to,sep="")))
    }

    if (plot.residuals=="martingale"){residuals.title<-"Martingale"}
    if (plot.residuals=="deviance"){residuals.title<-"Deviance (derived from Martingale)"}
    if (plot.residuals=="score"){residuals.title<-"Score"}
    if (plot.residuals=="dfbeta"){residuals.title<-"dfBeta"}
    if (plot.residuals=="dfbetas"){residuals.title<-"dfBeta stded by SE"}
    if (plot.residuals=="partial"){residuals.title<-"Partial"}

    png(residuals.file)

    if (is.numeric(x.variable[1])==FALSE){     # x-axis is character (e.g., case ID)
      if (num.resid==1){
        plot(residuals,type="h",xlab=plot.to,ylab=variable.name,main=residuals.title)
      }else{
        num.row.col.residplot<-as.integer(sqrt(num.resid)+1)
        par(mfrow=c(num.row.col.residplot,num.row.col.residplot))
    
        for (res in 1:num.resid){
          mfg.row<-ceiling(res/num.row.col.residplot)
          mfg.col<-res%%num.row.col.residplot
          if (mfg.col==0){mfg.col<-num.row.col.residplot}
          par(mfg=c(mfg.row,mfg.col))

          if (res==1){
            plot.command<-paste("plot(residuals[,",res,"],type=\"h\",xlab=plot.to,ylab=variable.name[",res,"],main=residuals.title)",sep="")
          }else{
            plot.command<-paste("plot(residuals[,",res,"],type=\"h\",xlab=plot.to,ylab=variable.name[",res,"])",sep="")
          }
          eval(parse(text=plot.command))
        }
      }

    }else{        # x-axis is numeric

      if (num.resid==1){
        plot(x.variable,residuals)
        scatter.smooth(x.variable,residuals,xlab=plot.to,ylab=variable.name,main=residuals.title)
      }else{
        num.row.col.residplot<-as.integer(sqrt(num.resid)+1)
        par(mfrow=c(num.row.col.residplot,num.row.col.residplot))
    
        for (res in 1:num.resid){
          mfg.row<-ceiling(res/num.row.col.residplot)
          mfg.col<-res%%num.row.col.residplot
          if (mfg.col==0){mfg.col<-num.row.col.residplot}
          par(mfg=c(mfg.row,mfg.col))

          if (res==1){
            plot.command<-paste("plot(x.variable,residuals[,",res,"],scatter.smooth(x.variable,residuals[,",res,"],xlab=plot.to,ylab=variable.name[",res,"],main=residuals.title),xlab=plot.to,ylab=variable.name[",res,"],main=residuals.title)",sep="")
          }else{
            plot.command<-paste("plot(x.variable,residuals[,",res,"],scatter.smooth(x.variable,residuals[,",res,"],xlab=plot.to,ylab=variable.name[",res,"]),xlab=plot.to,ylab=variable.name[",res,"])",sep="")
          }
          eval(parse(text=plot.command))
        }
      }
    }   # end x-axis is character or not

    dev.off()
  }

} # end of main



# read clinical data file

read.surv.data<-function(filename)
{
  if (regexpr(".txt$",filename)==-1){
    stop("### surv data should be .txt file! ###")
  }
  surv.data<-read.delim(filename,header=T)
  return(surv.data)
}

# read .cls file

read.cls<-function(filename="NA")
{
  if (regexpr(".cls$",filename)==-1){
    stop("### class data should be .cls file! ###")
  }
  cls<-as.vector(t(read.delim(filename,header=F,sep=" ",skip=2)))
  return(cls)
}

read.cls.label<-function(filename="NA")
{
  cls.label<-scan(filename,what="charactor",sep=" ",nlines=2,quiet=T)
  num.cls<-as.numeric(cls.label[2])
  cls.label<-cls.label[5:(5+num.cls-1)]
  return(cls.label)
}

convert.category<-function(input.variable.category)
{
  splitted.variable<-unlist(strsplit(input.variable.category,","))
  num.cat.variable<-length(splitted.variable)
  variable.category<-""
  for (i in 1:num.cat.variable){
    variable.category<-paste(variable.category,"factor(",splitted.variable[i],")+",sep="")
  }
  return(variable.category)
}
