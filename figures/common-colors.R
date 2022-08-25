
library("viridis")

col <- c(
	default = "black", 
	`u+`    = magma(4)[2],
	`u-`    = magma(4)[3], 
	`s+`    = viridis(4)[2],
	`s++`   = viridis(4)[3],
	`k2`    = plasma(5)[3],
	`k5`    = plasma(5)[4],
	'Neutral'="gray", 
	'Del TEs - neutral clusters'="violet", 
	'Del TEs - del clusters'="green", 
	'Copy number regul'="darkblue")

makeTransparent<-function(someColor, alpha=70)
	{ # from https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
		newColor<-col2rgb(someColor)
		apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
		blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
	}
