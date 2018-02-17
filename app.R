### INITIALIZE
options(stringsAsFactors = F)
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

setwd("/srv/shiny-server/homegenome/")

## load expect distr for GRSs in UKBB
ukbbSum = read.table("ukbb_summarizedgrs.tsv", h=T, sep="\t", quote='')
h2 = read.table("ukbb_all_h2univar_results.txt", h=T, sep="\t", quote="")
h2 = h2[,-c(2, 7:14, 19)]
notes = read.table("phenosummary_final_11898_18597.tsv", h=T, sep="\t", quote="")
grs = read.table("grs_everyone.tsv", h=T)

vcf = reactiveValues(raw = NULL)

## load phenotype categories
catDict = read.table("UKBB_category_dict.tsv", h=T, sep="\t", quote="")
catDict = separate_rows(ukbbSum, Categories, sep=",") %>%
	group_by(Categories) %>%
	summarize(numDiagnoses=n()) %>%
	full_join(catDict, ., by=c("codes"="Categories"))

# fix the descriptions for one category
ukbbSum$Description[grepl("cd", ukbbSum$Categories)] = paste("Cause of death:",
		substring(ukbbSum$Description[grepl("cd", ukbbSum$Categories)], first=8))

## set some general options
DTopts = list(pageLength=10, dom="tp", searching=F, ordering=F)
h2opts = list("Observed scale h2"="o", "Liability scale h2"="l", "-log10 P(observed h2>0)"="p")

ui <- fluidPage(
	tags$head(
		tags$style(HTML("hr {border-top: 1px solid #26386C;}
						.multicol {-webkit-column-count: 3; -moz-column-count: 3;}"))
	),
	
	titlePanel("Home Genome"),
	
	h3("Identify yourself"),
	fluidRow(column(4, selectInput("id", "Sample ID:", 1:3)),
			 column(2, br(), actionButton("advancedButton", "Advanced settings >>")),
			 column(2, br(), actionButton("loadButton", "Load sample"))),
	
	conditionalPanel("input.advancedButton%2==1", fluidRow(
			 column(3, radioButtons("h2selection", "reliability measure:", h2opts)),
			 column(3, sliderInput("h2thr", "use only GRSs w/ -log10 P below:", 0, 4, 0)),
			 column(2, radioButtons("refselection", "Reference population:", list("UKBB", "HRC")))
	)),
	
	h3("Here's your data:"),
	
	fluidRow(column(6, plotOutput("grsHist")),
			 column(6, plotOutput("grsScatter", brush="plot_brush"))),
	
	h3("Explore phenotypes"),

	navlistPanel(id = "navpanel", widths = c(3,9),
			tabPanel(value = "z", "Most extreme Z-scores",
					 h4("Phenos w/ extreme z-scores:"),
					 DT::dataTableOutput("sortedbyz")),
			tabPanel(value = "h", "Most reliable GRSs",
					 h4("Phenos w/ reliable GRSs:"),
					 DT::dataTableOutput("sortedbyh")),
			tabPanel(value = "cat", "Browse by category",
					 h4("Select categories to display below:"),
					 tags$div(class="multicol",
					 		 radioButtons("catselector", "Diseases:",
					 				   choiceNames = catDict$categories,
					 				   choiceValues = catDict$codes)),
					 htmlOutput("cattext"),
					 DT::dataTableOutput("filteredcat")),
			tabPanel(value = "br", "Select in plot",
					 h4("Selected points:"),
					 DT::dataTableOutput("brushed")),
			tabPanel(value = "cv", "Clinical variants",
					 h4("Variants associated with clinical consequences:"),
					 htmlOutput("cvtext"),
					 DT::dataTableOutput("cvtable")),
			tabPanel(value = "raw", "Browse raw data",
					 h4("Full table:"),
					 DT::dataTableOutput("rawdata"))
	),
	
	h3("Investigate phenotype"),
	
	div(style="border-style:solid; border-color:#5275BC; border-width:2px; margin:20px",
		fluidRow(
			column(6, plotOutput("phenoplot", height="300px"), style="padding:40px"),
			column(6, align="center", style="padding:50px; padding-top:80px",
			 	   htmlOutput("phenotext"),
				   htmlOutput("phenotext2"),
			 	   htmlOutput("phenotext3")
			)
		),
		fluidRow(style="padding:30px", hr(),
				 h4("Phenotype notes:"),
				 textOutput("phenotext4"))
	)
)

server <- function(input, output) {
	vcf = reactiveValues()
	selphenos = NULL
	selected = reactiveValues(current=data.frame(), selcv=data.frame())
	
	observeEvent(input$loadButton, {
		# read clinvar data and filter GRSs
		cvfile = paste0("clinvar_", input$id, ".tsv")
		vcf$cv = read.table(cvfile, comment.char="", h=T, sep="\t", quote="")
		vcf$grs = filter(grs, id==input$id)

		# merge with ukbb summaries
		if(input$refselection=="UKBB"){
			ukbbSum$GRSmean = ukbbSum$GRSmeanU
			ukbbSum$GRSsd = ukbbSum$GRSsdU
			vcf$cv$maf = vcf$cv$af
			vcf$cv$L = vcf$cv$pU
		} else {
			ukbbSum$GRSmean = ukbbSum$GRSmeanH
			ukbbSum$GRSsd = ukbbSum$GRSsdH
			vcf$cv$maf = vcf$cv$AF_GLOBAL
			vcf$cv$L = vcf$cv$pH
		}
		# categorize likelihood for simple presentation
		vcf$cv = arrange(vcf$cv, L) %>%
				 mutate(LCat = cut(L, c(0, 0.01, 0.05, 0.3, 1),
								include.lowest = T,
								labels = c("very rare", "rare", "unusual", "usual")))

		vcf$grs = left_join(vcf$grs, ukbbSum, by="pheno")
		vcf$grs = mutate(vcf$grs, myQ = pnorm(GRS, GRSmean, GRSsd),
						 myZ = (GRS-GRSmean)/GRSsd)
		
		# attach the currently selected h2 measure
		vcf$h2 = left_join(vcf$grs, h2, by=c("pheno"="phenotype"))
		
		# implement h2 p-value cutoff, if was stricter than 0
		# (NAs are kept only if cutoff left at 0)
		if(input$h2thr>0) vcf$h2 = filter(vcf$h2, -log(h2_p, 10) > input$h2thr)
		
		if(input$h2selection=="o"){
			vcf$h2 = mutate(vcf$h2, h2 = ifelse(is.na(h2_observed), 0, h2_observed))
		} else if(input$h2selection=="l"){
			vcf$h2 = mutate(vcf$h2, h2 = ifelse(is.na(h2_liability), 0, h2_liability))
		} else {
			vcf$h2 = mutate(vcf$h2, h2 = ifelse(is.na(h2_p), 0, -log(h2_p, 10)))
		} 
		
		# categorize myQ for simple presentation
		vcf$h2 = mutate(vcf$h2, myCat = cut(myQ, c(0, 0.15, 0.33, 0.67, 0.85, 1),
							include.lowest = T,
							labels = c("very low", "low", "average", "high", "very high")))
		
		# prepare tables used by navlistPanels
		vcf$sortedbyz = top_n(vcf$h2, 200, abs(myZ)) %>%
			arrange(desc(myZ))
		vcf$sortedbyh = top_n(vcf$h2, 200, h2) %>%
			arrange(desc(h2))
		vcf$filteredcat = vcf$sortedbyz
	})
	
	## top plots
	
	output$grsHist <- renderPlot({
		if(is.null(vcf$grs)) return(ggplot())
		
		ggplot(vcf$grs, aes(x=myZ)) +
			geom_histogram(color="skyblue4", fill="skyblue") + 
			stat_bin(aes(y=..count.., label=ifelse(..count..>0, ..count.., "")),
					 geom="text", vjust=-0.5) +
			scale_y_continuous(expand = c(0,0)) + expand_limits(y=250) +
			xlab("my Z score") + theme_gray() 
	})
	
	output$grsScatter <- renderPlot({
		if(is.null(vcf$h2)) return(ggplot())
		
		# this is for row selection in tables
		if(input$navpanel=="z"){
			selphenos = vcf$sortedbyz$pheno[input$sortedbyz_rows_selected]
		} else if(input$navpanel=="h"){
			selphenos = vcf$sortedbyh$pheno[input$sortedbyh_rows_selected]
		} else if(input$navpanel=="cat"){
			selphenos = vcf$filteredcat$pheno[input$filteredcat_rows_selected]
		} else if(input$navpanel=="br"){
			selphenos = vcf$brushed$pheno[input$brushed_rows_selected]
		} else if(input$navpanel=="raw"){
			selphenos = vcf$h2$pheno[input$rawdata_rows_selected]
		} else if(input$navpanel=="cv"){
			selhrc = vcf$cv$rsHRC[input$cvtable_rows_selected]
			selected$selcv = filter(vcf$cv, rsHRC %in% selhrc)
		} else {
			selphenos = c()
		}
		selected$current = filter(vcf$h2, pheno %in% selphenos)

		ggplot(vcf$h2, aes(x=myZ, y=h2)) +
			geom_point(aes(col=is.na(h2_observed)), size=1) +
			geom_point(data=selected$current, col="orange", size=3) +
			scale_color_manual(values=c("darkblue", "grey50"), guide="none") +
			xlab("my Z score") + ylab("phenotype h2") + theme_bw()
	})
	
	## middle tabsets
	
	output$brushed <- DT::renderDataTable({
		vcf$brushed = brushedPoints(vcf$h2, input$plot_brush)
		
		if(!nrow(vcf$brushed)) return(data.frame())
		
		vcf$brushed[,c("pheno", "Description", "h2", "myQ", "myZ")] %>%
			datatable(rownames=F, selection="single", options=DTopts) %>%
			formatRound(c("h2", "myQ", "myZ"), digits=3)
	})
	
	output$sortedbyz <- DT::renderDataTable({
		if(is.null(vcf$sortedbyz)) return(data.frame())
		
		vcf$sortedbyz[,c("pheno", "Description", "myQ", "myZ")] %>%
			datatable(rownames=F, selection="single", options=DTopts) %>%
			formatRound(c("myQ", "myZ"), digits=3)
	})
	
	output$sortedbyh <- DT::renderDataTable({
		if(is.null(vcf$sortedbyh)) return(data.frame())
		
		vcf$sortedbyh[,c("pheno", "Description", "h2", "myQ", "myZ")] %>%
			datatable(rownames=F, selection="single", options=DTopts) %>%
			formatRound(c("h2", "myQ", "myZ"), digits=3)
	})
	
	output$filteredcat <- DT::renderDataTable({
		if(is.null(vcf$sortedbyz)) return(data.frame())
		
		vcf$filteredcat = separate_rows(vcf$h2, Categories, sep=",") %>%
			filter(Categories == input$catselector) %>%
			arrange(desc(myZ))
		
		vcf$filteredcat[,c("pheno", "Description", "h2", "myQ", "myZ")] %>%
			datatable(rownames=F, selection="single", 
				options = list(pageLength=10, dom="ftp", searching=T, ordering=T)) %>%
			formatRound(c("h2", "myQ", "myZ"), digits=3)
	})
	
	output$cattext <- renderText({
		selcat = catDict[catDict$codes==input$catselector,]
		paste("<br />Selected<font color=\"6CA5CC\" size=\"3\"><b>", selcat$numDiagnoses,
			  "</b></font> diagnoses in category <font color=\"6CA5CC\" size=\"3\">",
			  selcat$categories, "</font><br/>Description:", selcat$definition, "<br />")
	})
	
	output$cvtable <- DT::renderDataTable({
		if(is.null(vcf$cv)) return(data.frame())
		
		vcf$cv[,c("Gene", "Name", "L", "gt", "maf",
				  "Clinical", "PhenotypeList")] %>%
			datatable(rownames=F, selection="single",
				options = list(pageLength=10, dom="ftp", searching=T, ordering=T)) %>%
			formatSignif(c("maf", "L"), 2)
	})
	
	output$cvtext <- renderText({
		sc = selected$selcv
		
		if(nrow(sc)){
			t1 = paste("You have<font color=\"6CA5CC\"><b>", round(sc$ds),
					   "</b></font>copies of allele", sc$ALT, "in variant", sc$rsHRC)
			t2 = paste("This is<font color=\"6CA5CC\"><b>", sc$LCat,
					   "</b></font> - <font color=\"6CA5CC\"><b>", signif(sc$L*100, 2),
					   "%</b></font> of people have this genotype.")
			t3 = paste("It affects the following phenotypes:<br /><font color=\"6CA5CC\">",
					   sc$PhenotypeList)
			paste(t1, t2, t3, sep="<br />")
		}
	})
	
	output$rawdata <- DT::renderDataTable({
		select(vcf$h2, -one_of(c("h2", "myCat", "id", "GRSmean", "GRSsd"))) %>%
			datatable(rownames=F, selection="single",
				options = list(pageLength=15, dom="ftp", searching=T, ordering=T)) %>%
			formatSignif(c("GRS", "GRSmeanU", "GRSsdU", "GRSmeanH", "GRSsdH", "myQ", "myZ",
						   "prevelence", "h2_observed", "h2_observed_se",
						   "h2_liability", "h2_liability_se", "h2_p"), digits=3)
	})
	
	## bottom panel
	
	output$phenoplot <- renderPlot({
		sc = selected$current
		if(nrow(sc)){
			xlims = c(min(sc$GRS-sc$GRSmean, -2*sc$GRSsd),
					  max(sc$GRS-sc$GRSmean, 2*sc$GRSsd)) * 1.1 + sc$GRSmean
			pdfnorm = data.frame(x=seq(xlims[1], xlims[2], length.out=500))
			pdfnorm$y = dnorm(pdfnorm$x, sc$GRSmean, sc$GRSsd)
			myY = dnorm(sc$GRS, sc$GRSmean, sc$GRSsd)
			
			ggplot(pdfnorm, aes(x=x, y=y)) + geom_area(aes(fill=x>sc$GRS)) +
				annotate("text", x=sc$GRS, y=myY, size=10,
						 col="skyblue4", label="'' %down% ''", parse=T, vjust=0) +
				annotate("text", x=sc$GRS, y=myY, size=5,
						 col="skyblue4", label="you", parse=T, vjust=-2.3) +
				xlab("GRS") + ylab("frequency") +
				scale_x_continuous(limits=xlims) +
				scale_fill_manual(values=c("#8CC5EC", "grey60"), guide="none") +
				theme_minimal() + theme(panel.grid.major=element_blank(),
										axis.line.x=element_line(color="black"))
		}
	}, height=300)
	
	output$phenotext <- renderText({
		paste("for phenotype<br /><font color=\"6CA5CC\" size=\"3\">",
			  selected$current$Description, "</font><br /><br />")
	})
	
	output$phenotext2 <- renderText({
		paste("your GRS is<br /><font color=\"6CA5CC\" size=\"5\"><b>",
			  selected$current$myCat, "</b></font><br /><br />")
	})
	
	output$phenotext3 <- renderText({
		paste("<font color=\"6CA5CC\" size=\"4\"><b>",
			  format(selected$current$myQ*100, digits=3), "% </b></font>",
			  "<font size=\"3\">of people have a smaller GRS than you</font>")
	})
	
	output$phenotext4 <- renderText({
		notes$Notes[notes$Field.code==selected$current$pheno]
	})
}

# Run the application 
shinyApp(ui = ui, server = server)


