makeNotes: 
	$(info make clean = delete all output to start over)
	$(info make demux = demultiplex raw MiSeq data with Pheniqs)
	$(info make trim = remove primers and sequencing adapters) 
	$(info make denoise = denoise with DADA2)
	$(info make compile = merge libraries, add metadata and taxonomy)
	$(info make analysis = run a simple example analysis)
	
clean:
	rm -r output
	$(info you have a clean slate)
	
##############################
### demultiplex MiSeq data ###
##############################

muxDir = output/processing/demux
demux: ${muxDir}

${muxDir} ${muxDir}/demux.report.txt: code/demux.config.json 
	mkdir -p ${muxDir} && touch ${muxDir}
	pheniqs mux -R ${muxDir}/demux.report.txt -c $< --base-input data/MiSeq_raw --base-output ${muxDir}/

###############################
### trim primers / adapters ###
###############################

trimDir = output/processing/trim
trim: ${trimDir}

${trimDir} ${trimDir}/summary.tab: code/trim.sh ${muxDir}
	mkdir -p ${trimDir} && touch ${trimDir}
	$< ${muxDir} ${trimDir}

##############
## denoise ###
##############

dadaDir = output/processing/dada
denoise: ${dadaDir}

${dadaDir} ${dadaDir}/dadaSmmary.csv: code/denoise.R ${trimDir}
	mkdir -p ${dadaDir} && touch ${dadaDir}
	Rscript $< ${trimDir} ${dadaDir}

##############################
### compile processed data ###
##############################

compDir = output/processing/compiled/
compile: ${compDir}

${compDir} ${compDir}/phy.rds: code/compile.R ${dadaDir}
	mkdir -p ${compDir} && touch ${compDir}
	Rscript $< 

################
### Analysis ###
################
	
analysis: output/html/analysis.html
	
output/html/analysis.html: code/analysis.Rmd ${compDir}/phy.rds
	mkdir -p output/figs
	mkdir -p output/html
	Rscript -e 'rmarkdown::render("$<", output_dir = "output/html", knit_root_dir = "$(CURDIR)")'
	rm Rplots.pdf


	