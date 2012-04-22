# To do:  column chart/mosaic plot.

fncGetConfIntText <- function(estimate, CIinf, CIsup, .faraCI=FALSE) {
	if(.faraCI == TRUE) {
		paste(round(estimate, decimals), " (", round(CIinf, decimals), " - ", round(CIsup, decimals), ")", sep="")
	} else {
		paste(round(estimate, decimals), " (95% CI ", round(CIinf, decimals), " - ", round(CIsup, decimals), ")", sep="")
	}
}

fncNbPercents <- function(.Table, .decimals) {
	.RowTable <- rowPercents(.Table, digits=.decimals)
	cbind(.Table[,1],rowPercents(.Table)[,1])
	for(i in 1:length(.Table[,1])) {
		for(j in 1:length(.Table[1,])) {
			.RowTable[i,j] <- paste(.Table[i,j], " (", .RowTable[i,j], ")", sep="")
			if (i == 1) {colnames(.RowTable)[j]  <- paste(colnames(.RowTable)[j], ": nb. (%)", sep="")}
			if (j == 1) {.RowTable[i,length(.RowTable[1,])] <- paste(.RowTable[i,length(.RowTable[1,])], " (100)", sep="")}
		}
	}
	colnames(.RowTable)[length(.RowTable[1,])] <- "Total (%)"
	.RowTable <- .RowTable[,c(1:length(.Table[1,]),length(.RowTable[1,]))]
	.RowTable
}

fncEBMCrossTab <- function(.table, .x, .y, .xlab, .ylab, .percents, .chisq, .expected, .chisqComp, .fisher, .indicators, .decimals=2) {
library(abind, pos=4)

if (class(.x) != 'factor') {
	.Table = .table
} else {
	.Table <- xtabs(~.x+.y, exclude = c(NA, NaN))
	.Table <<- .Table
}

.TableOriginal <<- .Table

.TableExample <- matrix(c("a", "b", "c", "d"), 2, 2, byrow=TRUE)

.decimals <- as.numeric(.decimals)
decimals <<- .decimals

if (.indicators != "dg") {# to avoid showing chi square test for diagnostic tests
        if (.chisq == '1') {#code from Rcmdr - John Fox
            command <- "chisq.test(.Table, correct=FALSE)"
            #logger(paste(".Test <- ", command, sep=""))
            assign(".Test", justDoIt(command), envir=.GlobalEnv)
            doItAndPrint(".Test")
            if (.expected == '1') doItAndPrint(".Test$expected # Expected Counts")
            warnText <- NULL
            if (0 < (nlt1 <- sum(.Test$expected < 1))) warnText <- paste(nlt1,
                gettextRcmdr("expected frequencies are less than 1"))
            if (0 < (nlt5 <- sum(.Test$expected < 5))) warnText <- paste(warnText, "\n", nlt5,
                gettextRcmdr(" expected frequencies are less than 5"), sep="")
            if (!is.null(warnText)) Message(message=warnText,
                type="warning")
            if (.chisqComp == '1') {
                command <- paste("round(.Test$residuals^2, ", .decimals, ") # Chi-square Components", sep="")
                doItAndPrint(command)
                }
            #logger("remove(.Test)")
            remove(.Test, envir=.GlobalEnv)
            }

	if (.fisher == '1') {#code from Rcmdr - John Fox
		.TestFisher <- fisher.test(.Table)
		print(.TestFisher)
	}
}


	if (.percents == "row") {
		cat("\n# Row Percentages\n")
		command <- "fncNbPercents(.Table, decimals)" # to allow the table to be exported to html through Export plugin
		logger(paste(".TheTable <- ", command, sep=""))
		assign(".TheTable", justDoIt(command), envir=.GlobalEnv)
		doItAndPrint(".TheTable")
	}
        if (.percents == "column") {
		cat("\n# Column Percentages\n")
		#print(t(fncNbPercents(t(.Table), .decimals)))
		command <- "t(fncNbPercents(t(.Table), decimals))"
		logger(paste(".TheTable <- ", command, sep=""))
		assign(".TheTable", justDoIt(command), envir=.GlobalEnv)
		doItAndPrint(".TheTable")
	}
        if (.percents == "total") {
		cat("\n# Percentages of Total\n")
		#print(totPercents(.Table))
		command <- "totPercents(.Table, decimals)"
		logger(paste(".TheTable <- ", command, sep=""))
		assign(".TheTable", justDoIt(command), envir=.GlobalEnv)
		doItAndPrint(".TheTable")
	}


if ((.x == "") || (length(levels(.x))==2 && length(levels(.y))==2) ) {
	require(epiR)

    if (.indicators == "pg") {
	#epi.2by2(dat = .TableOriginal , method = "cohort.count", conf.level = 0.95, units = 1, verbose = FALSE)

	.epi <- epi.2by2(dat = .TableOriginal , method = "cohort.count", conf.level = 0.95, units = 1, verbose = TRUE)

	cat("\n# Notations for calculations\n")
	colnames(.TableExample) <- c("Disease +", "Disease -")
	rownames(.TableExample) <- c("Exposure +", "Exposure -")
	print(.TableExample)
	
	texttest <- paste("\n# Risk difference = ", fncGetConfIntText(round(100*.epi$AR$est,.decimals), round(100*.epi$AR$lower,.decimals), round(100*.epi$AR$upper,.decimals)), " %. Computed using formula: [a / (a + b)] - [c / (c + d)]", sep="")
	cat(texttest)
	texttest <- paste("\n# Relative risk = ", fncGetConfIntText(round(.epi$RR$est,.decimals), round(.epi$RR$lower,.decimals), round(.epi$RR$upper,.decimals)), " %. Computed using formula: [a / (a + b)] / [c / (c + d)]", sep="")
	cat(texttest)
	texttest <- paste("\n# Odds ratio = ", fncGetConfIntText(round(.epi$OR$est,.decimals), round(.epi$OR$lower,.decimals), round(.epi$OR$upper,.decimals)), " . Computed using formula: (a / b) / (c / d)", sep="")
	cat(texttest)
	
	cat("\n# To find more about the results, and about how confidence intervals were computed, type ?epi.2by2 .\n")
    }
        if (.indicators == "th") {
	#epi.2by2(dat = .TableOriginal , method = "cohort.count", conf.level = 0.95, units = 1, verbose = FALSE)
	.epi <- epi.2by2(dat = .TableOriginal , method = "cohort.count", conf.level = 0.95, units = 1, verbose = TRUE)
	
	cat("\n# Notations for calculations\n")
	colnames(.TableExample) <- c("Event +", "Event -")
	rownames(.TableExample) <- c("Treatment", "Control")
	print(.TableExample)

	.ARR.est <-  - .epi$AR$est
	.ARR1 <-  - .epi$AR$lower
	.ARR2 <-  - .epi$AR$upper
	.ARR.lower <- min(.ARR1, .ARR2)
	.ARR.upper <- max(.ARR1, .ARR2)
	
	.NNT.est <- 1/.ARR.est
	.NNT1 <- 1/.ARR.lower
	.NNT2 <- 1/.ARR.upper
	.NNT.lower <- min(.NNT1, .NNT2)
	.NNT.upper <- max(.NNT1, .NNT2)
	if (.ARR.lower < 0) {
		.NNT.lower <- .NNT.upper
		.NNT.upper <- 1/0
	}

	.RR.est <-  .epi$RR$est
	.RR1 <-  .epi$RR$lower
	.RR2 <-  .epi$RR$upper
	.RR.lower <- min(.RR1, .RR2)
	.RR.upper <- max(.RR1, .RR2)

	.RRR.est <- 1 - .RR.est
	.RRR1 <- 1 - .RR.lower
	.RRR2 <- 1 - .RR.upper
	.RRR.lower <- min(.RRR1, .RRR2)
	.RRR.upper <- max(.RRR1, .RRR2)

	texttest <- paste("\n# Absolute risk reduction (ARR) = ", fncGetConfIntText(round(100*.ARR.est,.decimals), round(100*.ARR.lower,.decimals), round(100*.ARR.upper,.decimals)), " %. Computed using formula: [c / (c + d)] - [a / (a + b)] ", sep="")
	cat(texttest)
	texttest <- paste("\n# Relative risk = ", fncGetConfIntText(round(.RR.est,.decimals), round(.RR.lower,.decimals), round(.RR.upper,.decimals)), " %. Computed using formula: [c / (c + d)] / [a / (a + b)]", sep="")
	cat(texttest)
	texttest <- paste("\n# Odds ratio = ", fncGetConfIntText(round(.epi$OR$est,.decimals), round(.epi$OR$lower,.decimals), round(.epi$OR$upper,.decimals)), ". Computed using formula: (a / b) / (c / d)", sep="")
	cat(texttest)
	texttest <- paste("\n# Number needed to treat = ", fncGetConfIntText(round(.NNT.est,.decimals), round(.NNT.lower,.decimals), round(.NNT.upper,.decimals)), ". Computed using formula: 1 / ARR", sep="")
	cat(texttest)	
	texttest <- paste("\n# Relative risk reduction = ", fncGetConfIntText(round(100*.RRR.est,.decimals), round(100*.RRR.lower,.decimals), round(100*.RRR.upper,.decimals)), " %. Computed using formula: { [c / (c + d)] - [a / (a + b)] } / [c / (c + d)] ", sep="")
	cat(texttest)
	
  cat("\n# To find more about the results, and about how confidence intervals were computed, type ?epi.2by2 . The confidence limits for NNT were computed as 1/ARR confidence limits. The confidence limits for RRR were computed as 1 - RR confidence limits.\n")
}
    if (.indicators == "dg") {

    	cat("\n# Notations for calculations\n")
	colnames(.TableExample) <- c("Disease +", "Disease -")
	rownames(.TableExample) <- c("Test +", "Test -")
	print(.TableExample)
 
	#.dd <- epi.tests(a = .TableOriginal[1,1], b = .TableOriginal[1,2], c = .TableOriginal[2,1], d = .TableOriginal[2,2], conf.level = 0.95)

	.dd <- epi.tests(.TableOriginal, conf.level = 0.95, verbose = TRUE)
	texttest <- paste("\n# Sensitivity (Se) = ", fncGetConfIntText(round(100*.dd$se$est,.decimals), round(100*.dd$se$lower,.decimals), round(100*.dd$se$upper,.decimals)), " %. Computed using formula: a / (a + c)", sep="")
	cat(texttest)

	texttest <- paste("\n# Specificity (Sp) = ", fncGetConfIntText(round(100*.dd$sp$est,.decimals), round(100*.dd$sp$lower,.decimals), round(100*.dd$sp$upper,.decimals)), " %. Computed using formula: d / (b + d)", sep="")
	cat(texttest)	
    
	texttest <- paste("\n# Diagnostic acuracy (% of all correct results) = ", fncGetConfIntText(round(100*.dd$diag.acc$est,.decimals), round(100*.dd$diag.acc$lower,.decimals), round(100*.dd$diag.acc$upper,.decimals)), " %. Computed using formula: (a + d) / (a + b + c + d)", sep="")
	cat(texttest)
    
	texttest <- paste("\n# Youden's index = ", fncGetConfIntText(.dd$youden$est, .dd$youden$lower, .dd$youden$upper), ". Computed using formula: Se + Sp - 1", sep="")
	cat(texttest)
    
    	texttest <- paste("\n# Likelihood ratio of a positive test = ", fncGetConfIntText(.dd$plr$est, .dd$plr$lower, .dd$plr$upper), ". Computed using formula: Se / (Sp - 1)", sep="")
    	cat(texttest)
    
    texttest <- paste("\n# Likelihood ratio of a negative test = ", fncGetConfIntText(.dd$nlr$est, .dd$nlr$lower, .dd$nlr$upper), ". Computed using formula: (1 - Se) / Sp", sep="")
    cat(texttest)
    
	texttest <- paste("\n# Positive predictive value = ", fncGetConfIntText(round(100*.dd$ppv$est,.decimals), round(100*.dd$ppv$lower,.decimals), round(100*.dd$ppv$upper,.decimals)), " %. Computed using formula: a / (a + b)", sep="")
	cat(texttest)
  
	texttest <- paste("\n# Negative predictive value = ", fncGetConfIntText(round(100*.dd$npv$est,.decimals), round(100*.dd$npv$lower,.decimals), round(100*.dd$npv$upper,.decimals)), " %. Computed using formula: d / (c + d)", sep="")
	cat(texttest)
    
	texttest <- paste("\n# Number needed to diagnose = ", fncGetConfIntText(.dd$nnd$est, .dd$nnd$lower, .dd$nnd$upper), ". Computed using formula: 1 / [Se - (1 - Sp)]", sep="")
	cat(texttest)	
    
	cat("\n# To find more about the results, and about how confidence intervals were computed, type ?epi.tests .\n")
 
    }
}
#remove(.Table)

#if (.percents != "none") {
#	.rowTable <- rowPercents(.Table, digits=.decimals)
	#.Table <- cbind(.rowTable[,1], .rowTable[,2])
#}

	mosaicplot(t(.Table),color=heat.colors(2),  main="", ylab=.ylab, xlab=.xlab)

#par(xpd=T, mar=par()$mar+c(0,0,0,8))##
#barplot(.Table, xlab=.xlab, ylab=.ylab, beside=F, horiz=F, col=rainbow(2), cex.axis=1, cex=0.8)
#legend(2.75, .Table[1,1]+.Table[2,1], c(colnames(.Table)[1], colnames(.Table)[2]), cex=0.8, fill=rainbow(2));
#legent 7 pt beside = T
#par(mar=c(5, 4, 4, 2) + 0.1)

#remove(.Table)

remove(.x)
remove(.y)
}

#=========================================================================================================================================

fncEBMCrossTabWin <- function(){ #code from Rcmdr - John Fox, with modifications
    initializeDialog(title=gettextRcmdr("Multiple Cross Tab Tests"))
    variablesFrame <- tkframe(top) #
    .numeric <- Numeric()
    yBox <- variableListBox(variablesFrame, Factors(), selectmode="multiple",
        title=gettextRcmdr("Response variable (pick one or more)"))
    xBox <- variableListBox(variablesFrame, Factors(), selectmode="multiple", title=gettextRcmdr("Group variable (pick one or more)"))
    
    onOK <- function(){
        x <- getSelection(xBox)
        y <- getSelection(yBox)
	percents <- as.character(tclvalue(percentsVariable))
        chisq <- tclvalue(chisqTestVariable)
        chisqComp <- tclvalue(chisqComponentsVariable)
        expected <- tclvalue(expFreqVariable)
        fisher <- tclvalue(fisherTestVariable)
	indicators <- as.character(tclvalue(indicatorsVariable))
	digits <- tclvalue(Digitsname)
	
        if (0 == length(y)) {
            errorCondition(recall=fncFastTestT, message=gettextRcmdr("You must select a response variable."))
            return()
            }
        if (0 == length(x)) {
            errorCondition(recall=fncFastTestT, message=gettextRcmdr("No explanatory variables selected."))
            return()
            }

	subset <- tclvalue(subsetVariable) #
	subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) "" else paste(", subset=", subset, sep="") #
        closeDialog()
        .activeDataSet <- ActiveDataSet()

	for (j in 1:length(x)) {
	    for (i in 1:length(y)) {
		doItAndPrint(paste("fncEBMCrossTab(.table='', .y=", .activeDataSet, "$", y[i], ", .x=", .activeDataSet, "$", x[j], ", .ylab='", y[i], "', .xlab='", x[j], "', .percents='", percents, "', .chisq='", chisq, "', .expected='", expected, "', .chisqComp='", chisqComp, "', .fisher='",fisher, "', .indicators='",indicators,"', .decimals='",digits,"')", sep=""))
	    }
        }
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="epi.tests")
    
radioButtons(name="percents",
        buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"),
        values=c("row", "column", "total", "none"), initialValue="row",
        labels=gettextRcmdr(c("Row percentages", "Column percentages", "Percentages of total", "No percentages")), title=gettextRcmdr("Compute Percentages"))  

checkBoxes(frame="testsFrame", boxes=c("chisqTest", "chisqComponents", "expFreq", "fisherTest"), initialValues=c("1", "0", "0", "0"),
        labels=gettextRcmdr(c("Chi-square test of independence", "Components of chi-square statistic",
            "Print expected frequencies", "Fisher's exact test")))
	    
radioButtons(name="indicators",
        buttons=c("pg", "dg", "th"),
        values=c("pg", "dg", "th"), initialValue="pg",
        labels=gettextRcmdr(c("Prognosis", "Diagnosis", "Therapy")), title=gettextRcmdr("Medical indicators"))  

    optionsFrame <- tkframe(top) #
    
    Digitsname <- tclVar("2")
    DigitsEntry <- ttkentry(optionsFrame, width="5", textvariable=Digitsname)

    subsetBox() #

    tkgrid(getFrame(xBox), labelRcmdr(variablesFrame, text="    "), getFrame(yBox), sticky="nw") # 
    tkgrid(variablesFrame, sticky="nw") #
    tkgrid(percentsFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Hypothesis Tests"), fg="blue"), sticky="w")
    tkgrid(testsFrame, sticky="w")   
    tkgrid(indicatorsFrame, sticky="w")  
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Options"), fg="blue"), sticky="w")
    tkgrid(optionsFrame, sticky="nw") #
    tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Digits")), sticky="w")
    tkgrid(DigitsEntry, sticky="w")
    tkgrid(subsetFrame, sticky="w") #
    tkgrid(buttonsFrame, sticky="w") #
    dialogSuffix(rows=6, columns=1) #
    }    

#=========================================================================================================================================

# based on rcmdr code enterTable

enterTableEBMCrossTab <- function(){#code from Rcmdr - John Fox, with modifications
    Library("abind")
    env <- environment()
    initializeDialog(title=gettextRcmdr("Enter Two-Way Table for Evidence Based Medicine medical indicators"))
    outerTableFrame <- tkframe(top)
    nrows <- 2
    ncols <- 2
    assign(".tableFrame", tkframe(outerTableFrame), envir=env)
    #setUpTable <- function(...){
        tkdestroy(get(".tableFrame", envir=env))
        assign(".tableFrame", tkframe(outerTableFrame), envir=env)
        nrows <- 2
        ncols <- nrows
        make.col.names <- "labelRcmdr(.tableFrame, text='')"
        for (j in 1:ncols) {
            col.varname <- paste(".colname.", j, sep="")
            assign(col.varname, tclVar(j), envir=env)
            make.col.names <- paste(make.col.names, ", ", "ttkentry(.tableFrame, width='5', textvariable=",
                    col.varname, ")", sep="")
            }
        eval(parse(text=paste("tkgrid(", make.col.names, ")", sep="")), envir=env)
        for (i in 1:nrows){
            varname <- paste(".tab.", i, ".1", sep="")
            assign(varname, tclVar("") , envir=env)
            row.varname <- paste(".rowname.", i, sep="")
            assign(row.varname, tclVar(i), envir=env)
            make.row <- paste("ttkentry(.tableFrame, width='5', textvariable=",
                row.varname, ")", sep="")
            make.row <- paste(make.row, ", ", "ttkentry(.tableFrame, width='5', textvariable=",
                varname, ")", sep="")
            for (j in 2:ncols){
                varname <- paste(".tab.", i, ".", j, sep="")
                assign(varname, tclVar(""), envir=env)
                make.row <- paste(make.row, ", ", "ttkentry(.tableFrame, width='5', textvariable=",
                    varname, ")", sep="")
                }
            eval(parse(text=paste("tkgrid(", make.row, ")", sep="")), envir=env)
            }
        tkgrid(get(".tableFrame", envir=env), sticky="w")
        #}
    rowColFrame <- tkframe(top)

    onOK <- function(){
        nrows <- 2
        ncols <- nrows
        cell <- 0
        counts <- rep(NA, nrows*ncols)
        row.names <- rep("", nrows)
        col.names <- rep("", ncols)
        chisq <- tclvalue(chisqTestVariable)
        chisqComp <- tclvalue(chisqComponentsVariable)
        expected <- tclvalue(expFreqVariable)
        fisher <- tclvalue(fisherTestVariable)
	indicators <- as.character(tclvalue(indicatorsVariable))
	digits <- tclvalue(Digitsname)
	
        for (i in 1:nrows) row.names[i] <-
            eval(parse(text=paste("tclvalue(", paste(".rowname.", i, sep=""),")", sep="")))
        for (j in 1:ncols) col.names[j] <-
            eval(parse(text=paste("tclvalue(", paste(".colname.", j, sep=""),")", sep="")))
        for (i in 1:nrows){
            for (j in 1:ncols){
                cell <- cell+1
                varname <- paste(".tab.", i, ".", j, sep="")
                counts[cell] <- as.numeric(eval(parse(text=paste("tclvalue(", varname,")", sep=""))))
                }
            }
        counts <- na.omit(counts)
        if (length(counts) != nrows*ncols){
            errorCondition(recall=enterTableEBMCrossTab, message=sprintf(gettextRcmdr("Number of valid entries (%d)\nnot equal to number of rows (%d) * number of columns (%d)."), length(counts), nrows, ncols))
            return()
            }
        if (length(unique(row.names)) != nrows){
            errorCondition(recall=enterTableEBMCrossTab, message=gettextRcmdr("Row names are not unique."))
            return()
            }
        if (length(unique(col.names)) != ncols){
            errorCondition(recall=enterTableEBMCrossTab, message=gettextRcmdr("Column names are not unique."))
            return()
            }
        percents <- as.character(tclvalue(percentsVariable))
	#test <- as.character(tclvalue(testVariable)) #
        closeDialog()
        command <- paste("matrix(c(", paste(counts, collapse=","), "), ", nrows, ", ", ncols,
            ", byrow=TRUE, dimnames = list(c('", row.names[1], "', '", row.names[2], "'), c('", col.names[1], "', '", col.names[2], "')))", sep="")
        assign(".Table", justDoIt(command), envir=.GlobalEnv)
        logger(paste(".Table <- ", command, sep=""))

		#doItAndPrint(paste("fncCrossTab(.Table[,1], .x=.Table[,2], .ylab='', .xlab='', .chisim='", chisim, "', .percents='", percents, "', .fisher='",fisher, "', .lbltest='", linearbylinear, "', .indicators='",indicators,"')", sep=""))
		doItAndPrint(paste("fncEBMCrossTab(.table=.Table, .x='', .y='', .ylab='', .xlab='', .percents='", percents, "', .chisq='", chisq, "', .expected='", expected, "', .chisqComp='", chisqComp, "', .fisher='",fisher, "', .indicators='",indicators,"', .decimals=2)", sep=""))


	tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="epi.2by2")
    radioButtons(name="percents", buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"), values=c("row", "column", "total", "none"),
        initialValue="none", labels=gettextRcmdr(c("Row percentages", "Column percentages",  "Percentages of total", "No percentages")), title=gettextRcmdr("Compute Percentages"))

    optionsFrame <- tkframe(top) #
	
    checkBoxes(frame="testsFrame", boxes=c("chisqTest", "chisqComponents", "expFreq", "fisherTest"), initialValues=c("1", "0", "0", "0"),
        labels=gettextRcmdr(c("Chi-square test of independence", "Components of chi-square statistic",
            "Print expected frequencies", "Fisher's exact test")))

radioButtons(name="indicators",
        buttons=c("pg", "dg", "th"),
        values=c("pg", "dg", "th"), initialValue="pg",
        labels=gettextRcmdr(c("Prognosis", "Diagnosis", "Therapy")), title=gettextRcmdr("Medical indicators"))  

    optionsFrame <- tkframe(top) #
    
    Digitsname <- tclVar("2")
    DigitsEntry <- ttkentry(optionsFrame, width="5", textvariable=Digitsname)

    tkgrid(rowColFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Enter counts:"), fg="blue"), sticky="w")
    tkgrid(outerTableFrame, sticky="w")
    tkgrid(percentsFrame, sticky="w")
    
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Hypothesis Tests"), fg="blue"), sticky="w")
    tkgrid(testsFrame, sticky="w")
    #tkgrid(optionsFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Options"), fg="blue"), sticky="w")
    tkgrid(optionsFrame, sticky="nw") #
    tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Digits")), sticky="w")
    tkgrid(DigitsEntry, sticky="w")
    tkgrid(indicatorsFrame, sticky="w")  
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=8, columns=2)
    }

#=========================================================================================================================================

fncEBMTherapy <- function(){#code from Rcmdr - John Fox with modifications
    initializeDialog(title=gettextRcmdr("Therapy"))
    variablesFrame <- tkframe(top) #
    .numeric <- Numeric()
    yBox <- variableListBox(variablesFrame, Factors(), selectmode="multiple",
        title=gettextRcmdr("Response variable (pick one or more)"))
    xBox <- variableListBox(variablesFrame, Factors(), selectmode="multiple", title=gettextRcmdr("Group variable (pick one or more)"))
    
    onOK <- function(){
        x <- getSelection(xBox)
        y <- getSelection(yBox)
	percents <- as.character(tclvalue(percentsVariable))
        chisq <- tclvalue(chisqTestVariable)
        chisqComp <- tclvalue(chisqComponentsVariable)
        expected <- tclvalue(expFreqVariable)
        fisher <- tclvalue(fisherTestVariable)
	digits <- tclvalue(Digitsname)
	
        if (0 == length(y)) {
            errorCondition(recall=fncEBMTherapy, message=gettextRcmdr("You must select a response variable."))
            return()
            }
        if (0 == length(x)) {
            errorCondition(recall=fncEBMTherapy, message=gettextRcmdr("No explanatory variables selected."))
            return()
            }
	for (j in 1:length(x)) {
		if (length(which(y == x[j])) != 0) {
			errorCondition(recall=fncEBMTherapy, message=gettextRcmdr("Explanatory and response variables must be different."))
			return()
		}
	}

	subset <- tclvalue(subsetVariable) #
	subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) "" else paste(", subset=", subset, sep="") #
        closeDialog()
        .activeDataSet <- ActiveDataSet()

	for (j in 1:length(x)) {
	    for (i in 1:length(y)) {
		doItAndPrint(paste("fncEBMCrossTab(.table='', .y=", .activeDataSet, "$", y[i], ", .x=", .activeDataSet, "$", x[j], ", .ylab='", y[i], "', .xlab='", x[j], "', .percents='", percents, "', .chisq='", chisq, "', .expected='", expected, "', .chisqComp='", chisqComp, "', .fisher='",fisher, "', .indicators='th', .decimals='",digits,"')", sep=""))
	    }
        }
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="epi.2by2")
    
radioButtons(name="percents",
        buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"),
        values=c("row", "column", "total", "none"), initialValue="row",
        labels=gettextRcmdr(c("Row percentages", "Column percentages", "Percentages of total", "No percentages")), title=gettextRcmdr("Compute Percentages"))  

checkBoxes(frame="testsFrame", boxes=c("chisqTest", "chisqComponents", "expFreq", "fisherTest"), initialValues=c("1", "0", "0", "0"),
        labels=gettextRcmdr(c("Chi-square test of independence", "Components of chi-square statistic",
            "Print expected frequencies", "Fisher's exact test")))
	    
    optionsFrame <- tkframe(top) #
    
    Digitsname <- tclVar("2")
    DigitsEntry <- ttkentry(optionsFrame, width="5", textvariable=Digitsname)

    subsetBox() #

    tkgrid(getFrame(xBox), labelRcmdr(variablesFrame, text="    "), getFrame(yBox), sticky="nw") # 
    tkgrid(variablesFrame, sticky="nw") #
    tkgrid(percentsFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Hypothesis Tests"), fg="blue"), sticky="w")
    tkgrid(testsFrame, sticky="w")   
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Options"), fg="blue"), sticky="w")
    tkgrid(optionsFrame, sticky="nw") #
    tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Digits")), sticky="w")
    tkgrid(DigitsEntry, sticky="w")
    tkgrid(subsetFrame, sticky="w") #
    tkgrid(buttonsFrame, sticky="w") #
    dialogSuffix(rows=6, columns=1) #
    }    

#=========================================================================================================================================

fncEBMPrognosis <- function(){#code from Rcmdr - John Fox, with modifications
    initializeDialog(title=gettextRcmdr("Prognosis"))
    variablesFrame <- tkframe(top) #
    .numeric <- Numeric()
    yBox <- variableListBox(variablesFrame, Factors(), selectmode="multiple",
        title=gettextRcmdr("Response variable (pick one or more)"))
    xBox <- variableListBox(variablesFrame, Factors(), selectmode="multiple", title=gettextRcmdr("Group variable (pick one or more)"))
    
    onOK <- function(){
        x <- getSelection(xBox)
        y <- getSelection(yBox)
	percents <- as.character(tclvalue(percentsVariable))
        chisq <- tclvalue(chisqTestVariable)
        chisqComp <- tclvalue(chisqComponentsVariable)
        expected <- tclvalue(expFreqVariable)
        fisher <- tclvalue(fisherTestVariable)
	digits <- tclvalue(Digitsname)
	
        if (0 == length(y)) {
            errorCondition(recall=fncEBMPrognosis, message=gettextRcmdr("You must select a response variable."))
            return()
            }
        if (0 == length(x)) {
            errorCondition(recall=fncEBMPrognosis, message=gettextRcmdr("No explanatory variables selected."))
            return()
            }
	for (j in 1:length(x)) {
		if (length(which(y == x[j])) != 0) {
			errorCondition(recall=fncEBMPrognosis, message=gettextRcmdr("Explanatory and response variables must be different."))
			return()
		}
	}

	subset <- tclvalue(subsetVariable) #
	subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) "" else paste(", subset=", subset, sep="") #
        closeDialog()
        .activeDataSet <- ActiveDataSet()

	for (j in 1:length(x)) {
	    for (i in 1:length(y)) {
		doItAndPrint(paste("fncEBMCrossTab(.table='', .y=", .activeDataSet, "$", y[i], ", .x=", .activeDataSet, "$", x[j], ", .ylab='", y[i], "', .xlab='", x[j], "', .percents='", percents, "', .chisq='", chisq, "', .expected='", expected, "', .chisqComp='", chisqComp, "', .fisher='",fisher, "', .indicators='pg', .decimals='",digits,"')", sep=""))
	    }
        }
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="epi.2by2")
    
radioButtons(name="percents",
        buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"),
        values=c("row", "column", "total", "none"), initialValue="row",
        labels=gettextRcmdr(c("Row percentages", "Column percentages", "Percentages of total", "No percentages")), title=gettextRcmdr("Compute Percentages"))  

checkBoxes(frame="testsFrame", boxes=c("chisqTest", "chisqComponents", "expFreq", "fisherTest"), initialValues=c("1", "0", "0", "0"),
        labels=gettextRcmdr(c("Chi-square test of independence", "Components of chi-square statistic",
            "Print expected frequencies", "Fisher's exact test")))
	    
    optionsFrame <- tkframe(top) #
    
    Digitsname <- tclVar("2")
    DigitsEntry <- ttkentry(optionsFrame, width="5", textvariable=Digitsname)

    subsetBox() #

    tkgrid(getFrame(xBox), labelRcmdr(variablesFrame, text="    "), getFrame(yBox), sticky="nw") # 
    tkgrid(variablesFrame, sticky="nw") #
    tkgrid(percentsFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Hypothesis Tests"), fg="blue"), sticky="w")
    tkgrid(testsFrame, sticky="w")   
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Options"), fg="blue"), sticky="w")
    tkgrid(optionsFrame, sticky="nw") #
    tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Digits")), sticky="w")
    tkgrid(DigitsEntry, sticky="w")
    tkgrid(subsetFrame, sticky="w") #
    tkgrid(buttonsFrame, sticky="w") #
    dialogSuffix(rows=6, columns=1) #
    }    

#=========================================================================================================================================

fncEBMDiagnosis <- function(){#code from Rcmdr - John Fox, with modifications
    initializeDialog(title=gettextRcmdr("Diagnosis and screening"))
    variablesFrame <- tkframe(top) #
    .numeric <- Numeric()
    yBox <- variableListBox(variablesFrame, Factors(), selectmode="multiple",
        title=gettextRcmdr("Response variable (pick one or more)"))
    xBox <- variableListBox(variablesFrame, Factors(), selectmode="multiple", title=gettextRcmdr("Group variable (pick one or more)"))
    
    onOK <- function(){
        x <- getSelection(xBox)
        y <- getSelection(yBox)
	percents <- as.character(tclvalue(percentsVariable))
	digits <- tclvalue(Digitsname)
	
        if (0 == length(y)) {
            errorCondition(recall=fncEBMDiagnosis, message=gettextRcmdr("You must select a response variable."))
            return()
            }
        if (0 == length(x)) {
            errorCondition(recall=fncEBMDiagnosis, message=gettextRcmdr("No explanatory variables selected."))
            return()
            }
	for (j in 1:length(x)) {
		if (length(which(y == x[j])) != 0) {
			errorCondition(recall=fncEBMDiagnosis, message=gettextRcmdr("Explanatory and response variables must be different."))
			return()
		}
	}

	subset <- tclvalue(subsetVariable) #
	subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) "" else paste(", subset=", subset, sep="") #
        closeDialog()
        .activeDataSet <- ActiveDataSet()

	for (j in 1:length(x)) {
	    for (i in 1:length(y)) {
		doItAndPrint(paste("fncEBMCrossTab(.table='', .y=", .activeDataSet, "$", y[i], ", .x=", .activeDataSet, "$", x[j], ", .ylab='", y[i], "', .xlab='", x[j], "', .percents='", percents, "', .chisq='0', .expected='0', .chisqComp='0', .fisher='0', .indicators='dg', .decimals='",digits,"')", sep=""))
	    }
        }
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="epi.tests")
    
radioButtons(name="percents",
        buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"),
        values=c("row", "column", "total", "none"), initialValue="none",
        labels=gettextRcmdr(c("Row percentages", "Column percentages", "Percentages of total", "No percentages")), title=gettextRcmdr("Compute Percentages"))  

    optionsFrame <- tkframe(top) #
    
    Digitsname <- tclVar("2")
    DigitsEntry <- ttkentry(optionsFrame, width="5", textvariable=Digitsname)

    subsetBox() #

    tkgrid(getFrame(xBox), labelRcmdr(variablesFrame, text="    "), getFrame(yBox), sticky="nw") # 
    tkgrid(variablesFrame, sticky="nw") #
    tkgrid(percentsFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Options"), fg="blue"), sticky="w")
    tkgrid(optionsFrame, sticky="nw") #
    tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Digits")), sticky="w")
    tkgrid(DigitsEntry, sticky="w")
    tkgrid(subsetFrame, sticky="w") #
    tkgrid(buttonsFrame, sticky="w") #
    dialogSuffix(rows=6, columns=1) #
    }    

#=========================================================================================================================================

fncEBMPostTest <- function(.pretest, .LR){
	.pretest <- as.numeric(.pretest)
	.LR <- as.numeric(.LR)
	.pretestodds <- .pretest / (1 - .pretest)
	.posttestodds <- .pretestodds * .LR
	.posttest <- .posttestodds / (.posttestodds + 1)
	cat("\n# Post-test probability:")
	.posttest
}

fncEBMPostTestWin <- function(){
    initializeDialog(title=gettextRcmdr("Compute Post-test probability"))
    
    onOK <- function(){
	pretest <- as.numeric(tclvalue(Pretestname))
	LR <- as.numeric(tclvalue(LRname))
	
        if ((is.na(pretest)) || (is.na(LR))) {
            errorCondition(recall=fncEBMPostTestWin, message=gettextRcmdr("Both fields should be filled."))
            return()
            }
	if ((pretest < 0) || (LR < 0)) {
            errorCondition(recall=fncEBMPostTestWin, message=gettextRcmdr("Values must be positive."))
            return()
            }
        if (pretest > 1) {
            errorCondition(recall=fncEBMPostTestWin, message=gettextRcmdr("Pre-test values are from 0 to 1."))
            return()
            }

        closeDialog()

		doItAndPrint(paste("fncEBMPostTest(.pretest=", pretest, ", .LR=", LR, ")", sep=""))

        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="fncEBMPostTest")
    
    optionsFrame <- tkframe(top) #
    
    Pretestname <- tclVar("")
    PretestEntry <- ttkentry(optionsFrame, width="5", textvariable=Pretestname)

    LRname <- tclVar("")
    LREntry <- ttkentry(optionsFrame, width="5", textvariable=LRname)

    subsetBox() #

    tkgrid(labelRcmdr(top, text=gettextRcmdr("Options"), fg="blue"), sticky="w")
    tkgrid(optionsFrame, sticky="nw") #
    tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Pre-test probability")), sticky="w")
    tkgrid(PretestEntry, sticky="w")
    tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Likelihood value")), sticky="w")
    tkgrid(LREntry, sticky="w")
    tkgrid(buttonsFrame, sticky="w") #
    dialogSuffix(rows=6, columns=1) #
    }    

#=========================================================================================================================================