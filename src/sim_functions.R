library(parallel)

global.default <- list(
	mc.cores          = min(detectCores() - 1, 60),
	mc.cores.internal = 1,
	nb.loci           = 10000,
	neutral.loci      = numeric(0), 
	piRNA.loci        = 1:30,
	piRNA.prob        = 30/1000,
	fitness.FUN       = function(n.sel)    exp(-n.sel*0.01),
	regulation.FUN    = function(n.piRNA, n)  if (n.piRNA == 0) 1 else 0,
	u                 = 0.1,
	G                 = 100,
	summary.every     = 10,
	N                 = 1000,
	rep               = 10,
	init.TE.ind       = 1,
	simulator         = NA
)


Simulicron <- function(global) {
	simulicron.path <- "../../Simulicron/src/simulicronalpha/simulicron.py"
	
	simID <- tempfile()
	# translation from parameter names in global to Simulicron parameters
	simulicron.param <- list(
		generations      = global$G, 
		individuals      = global$N, 
		loci             = global$nb.loci,
		selectionPenalty = log(global$fitness.FUN(1)), #Assuming fitness.FUN(0) = 1
		tau              = 1,
		ExcisionRate     = global$u, 
		FrequencyOfInsertion = 1/10,
		Chromosomes      = 30,
		RecombinationRate= 0.499,
		NumberOfInsertions= 10*global$init.TE.ind, # Double check this
		piPercentage      = 100*global$piRNA.prob,
		numberOfPiRNA     = length(global$piRNA.loci),
		piRNASelection   =  if(length(global$neutral.loci) > 0) "False" else "True",
		regulationStr    = 1/global$regulation.FUN(0,1) - 1,
		FileName         = simID
	)
	
	command <- paste0("python3 ", simulicron.path, " ", paste(names(simulicron.param), simulicron.param, sep="=", collapse=" "))
#~ 	print(command)
	system(command, ignore.stdout=TRUE)
			
	dd <- read.table(simID, header=FALSE)
	unlink(simID)
	
	ans <- data.frame(n.tot.mean=dd[,1], n.piRNA.mean=dd[,2])
	ans
}


run.simul <- function(user.global=list(), simulator=c("Simulicron")[1]) {
	global <- global.default
	valid.user.global <- names(user.global)[names(user.global) %in% names(global)]
	if (length(valid.user.global) != length(user.global)) {
		warning(names(user.global)[! names(user.global) %in% global], ": unknown parameter names.")
	}
	global[valid.user.global] <- user.global[valid.user.global]
	if (is.na(global$mc.cores.internal))
		global$mc.cores.internal <- max(1, global$mc.cores %/% global$rep)
	global$notpiRNA.loci <- (1:global$nb.loci)[! 1:global$nb.loci %in% global$piRNA.loci]
	simFUN <- if (simulator == "Simulicron") {
			Simulicron
		} else {
			stop("Unknown simulator: ", simulator)
		}
	mclapply(1:global$rep, function(i) simFUN(global=global), mc.cores=global$mc.cores)
}
