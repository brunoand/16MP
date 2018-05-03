#!/usr/bin/env nextflow
	
/**
	Nsilico Another Metagenomic Pipeline (NMP)
	Copyright (C) 2018 	Dr Bruno Andrade 	      
	
*/

version='1.0'
timestamp='20181903'

/**
	Prints version when asked for
*/
if (params.version) {
	System.out.println("")
	System.out.println("Nsilico metagenomic pipeline (NMP) - Version: $version ($timestamp)")
	exit 1
}

/**
	Prints help when asked for
*/

if (params.help) {
	System.out.println("")
	System.out.println("Nsilico metagenomic pipeline (NMP) - Version: $version ($timestamp)")
	System.out.println("")
	System.out.println("Usage")
	System.out.println("   nextflow run -c nextflow.config NMP.nf --reads_R1 R1 --reads_R2 R2 --date date --prefix mysample --outdir path --mode MODE  ")
	System.out.println("                [options] [-with-docker]")
	System.out.println("")
	System.out.println("Mandatory arguments:")
	System.out.println("    --reads_R1   R1      Path for forward reads (library = paired-end) or for all reads (library = single-end")
	System.out.println("    [--reads_R2] R2      Path for reverse reads (library = paired-end)")
	System.out.println("    --prefix   prefix  Prefix used to name the result files")
	System.out.println("    --outdir   path    Output directory (will be outdir/prefix/date)")
	System.out.println("    --mode     <QC|Taxonomic|Complete>")
	System.out.println("Options:")
	System.out.println("    --library <single-end|paired-end>")
	System.out.println("    --dedup         <true|false>   whether to perform de-duplication")
	System.out.println("")
	System.out.println("Container:")
	System.out.println("    Docker image to use with docker is")
	System.out.println("    'docker://bgabriel/nmp'")
    exit 1
}

	
/**
	STEP 0. 
	
	Checks input parameters and (if it does not exists) creates the directory 
	where the results will be stored. 
	Initialises the log file.
	
	The working directory is named after the prefix and the correspondent sample date
        located in the outdir folder. The log file, that will save summary statistics,
        execution time,	and warnings generated during the pipeline execution, will be 
        saved in the working directory as "prefix.log".
*/


//Checking user-defined parameters	
if (params.mode != "QC" && params.mode != "Taxonomic" && params.mode != "Complete") {
	exit 1, "Mode not available. Choose any of <QC, Taxonomic, Complete>"
}	

if (params.library != "paired-end" && params.library != "single-end") { 
	exit 1, "Library layout not available. Choose any of <single-end, paired-end>" 
}   
if (params.Pcoding != 33 && params.Pcoding != 64) { 
	exit 1, "Input quality offset (Pcoding) not available. Choose either 33 (ASCII+33) or 64 (ASCII+64)" 
}
//--reads_R2 can be omitted when the library layout is "single-end"

if (params.library != "single-end" && (params.reads_R2 == "null") ) {
	exit 1, "If dealing with paired-end reads, please set the reads_R2 arguments and if dealing with single-end reads, please set the library argument to 'single-end'"
}


System.out.println(params.outdir)


//Creates working dir
workingpath = params.outdir
workingdir = file(workingpath)
if( !workingdir.exists() ) {
    if( !workingdir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $workingpath"
    } 
}

//Creating other folders
matrixpath = params.outdir + "/Matrix"
QCpath = params.outdir + "/QC"
Fastq_path = params.outdir + "/Fastq"
Binpath = params.outdir + "/Bin"
matrixdir = file(matrixpath)
QC_dir = file(QCpath)
Fastq_dir = file(Fastq_path)
Bindir = file(Binpath)


if( !matrixdir.exists() ) {
    if( !matrixdir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $matrixpath"
    } 
}
if( !QC_dir.exists() ) {
    if( !QC_dir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $QCpath"
    } 
}
if( !Fastq_dir.exists() ) {
    if( !Fastq_dir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $Fastq_path"
 } 
}


//Creates main log file
mylog = file(params.outdir +  "log.txt")

//Logs headers
mylog <<  """---------------------------------------------
|Nsilico Metagenomic Pipeline (NMP) - Version: $version ($timestamp)|
----------------------------------------------

"""
	   
//Fetches information on OS and java versions, including user name
osname = System.getProperty("os.name") //Operating system name
osarch = System.getProperty("os.arch") //Operating system architecture
osversion = System.getProperty("os.version") //Operating system version
osuser = System.getProperty("user.name") //User's account name

javaversion = System.getProperty("java.version") //Java Runtime Environment version
javaVMname = System.getProperty("java.vm.name") //Java Virtual Machine implementation name
javaVMVersion = System.getProperty("java.vm.version") //Java Virtual Machine implementation version

//Gets starting time		
sysdate = new java.util.Date() 
		
//Logs starting time and other information about the run		
mylog << """ 
Analysis starting at $sysdate by user: $osuser
Analysed sample(s): $params.reads_R1 and $params.reads_R2
Results will be saved at $workingdir
New files will be saved using the '$params.prefix' prefix

Analysis mode? $params.mode
Library layout? $params.library	

---------------------------------------------------------------------------------------
									    
	Analysis introspection:            				    
									    
	Operating System:						    
		name:         $osname					    
		architecture: $osarch					    
		version:      $osversion				    
									    
	Java:								    
		version: $javaversion					    
		Java Virtual Machine: $javaVMname ; version: $javaVMVersion 
									    
	nextflow:							    
		version:   $nextflow.version				    
		build:     $nextflow.build				    	
		timestamp: $nextflow.timestamp				    	
									    
	Container:							    
		Docker image: $workflow.container			    
									    
	Analysis environment:						    
									    
		projectDir: $workflow.projectDir			    
		launchDir:  $workflow.launchDir				    
		workingDir: $workflow.workDir				    
									    
		command line: $workflow.commandLine			    
									    
		Run name:   $workflow.runName				    
		Session ID: $workflow.sessionId				    
		profile:    $workflow.profile				    
									    
---------------------------------------------------------------------------------------   
	   
""" 

/**
	Quality Assessment - STEP 1. Read quality assessment using the FastQC software. 
	Multiple  plots are generated to show average phred quality scores and other metrics.
*/

	if (params.library == "paired-end") {
		rawreads = Channel.fromFilePairs(params.reads_R1 + '*R{1,2}*fastq')
	}
	else {
		rawreads = Channel.fromPath(params.reads_R1 + '*.fastq')
	}



// ------------------------------------------------------------------------------
//                             Microbiome barcode profiling                            
// ------------------------------------------------------------------------------

/**
        Microbiome profiling.

        Metaphlan2 will be used to anotate shotgun reads in taxonomic levels, this step will generate 3 outputs, a matrix containing the taxonomic levels, relative abundaces, read counts for each clade, sample size and average genome size in the working file and two others containing only the relative abundance and the read_counts in Matrix/Abundance and Matrix/Read_counts, respectively.
*/
mocktrim = Channel.from("null")

        if (params.library == "paired-end") {
                to_trim = Channel.fromFilePairs(params.reads_R1 + '*R{1,2}*fastq')
        }
        else {
                to_trim = Channel.fromPath(params.reads_R1 + '*.fastq')
        }

process trim {
	publishDir Fastq_dir, mode: 'copy', pattern: "*{_trimmed_R*.fq,html}"
	input:
	set val(id), file(reads1), file(reads2) from to_trim
	file(adapters) from file(params.adapters)
	file(artifacts) from file(params.artifacts)
	file(phix174ill) from file(params.phix174ill)

	output:
	file("*_trimmed*.fq") into todecontaminate
	file("*_trimmed*.fq") into trimmedreads
	file "*_fastqc.html"
	val Fastq_dir into to_bin
	when:
	params.mode == "QC" || params.mode == "Complete"
      script:
      """
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing Quality Control. STEP 2 [Trimming] at \$sysdate\" >>  $mylog
	echo \" \" >>  $mylog
	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
	#Defines command for trimming of adapters and low quality bases

	if [ \"$params.library\" = \"paired-end\" ]; then
		CMD=\"bbduk.sh -Xmx\"\$maxmem\" in=${reads1[0]} in2=${reads1[1]} out=${id}_trimmed_R1_tmp.fq out2=${id}_trimmed_R2_tmp.fq outs=${id}_trimmed_singletons_tmp.fq ktrim=r k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.Quality  minlength=$params.minlength ref=$adapters qin=$params.Pcoding threads=${task.cpus} tbo tpe ow\"
	else
		CMD=\"bbduk.sh -Xmx\"\$maxmem\" in=$reads1 out=`basename $reads1 | sed -r 's/_R.{1,}//'`_trimmed_tmp.fq ktrim=r k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.Quality  minlength=$params.minlength ref=$adapters qin=$params.Pcoding threads=${task.cpus} tbo tpe ow\"
	fi
	#Trims adapters and low quality bases	
	exec \$CMD 2>&1 | tee tmp.log
	
	#Logs some figures about sequences passing trimming
	echo  \"BBduk's trimming stats (trimming adapters and low quality reads): \" >>  $mylog
	sed -n '/Input:/,/Result:/p' tmp.log >>  $mylog
	echo \" \" >>  $mylog			
	if [ \"$params.library\" = \"paired\" ]; then
		unpairedR=\$(wc -l ${id}_trimmed_singletons_tmp.fq | cut -d\" \" -f 1)
		unpairedR=\$((\$unpairedR/4))
		echo  \"\$unpairedR singleton reads whose mate was trimmed shorter preserved\" >>  $mylog
		echo \" \" >>  $mylog
	fi
	
	#Logs version of the software and executed command (BBMap prints on stderr)
	version=\$(bbduk.sh --version 2>&1 >/dev/null | grep \"BBMap version\") 
	echo \"Using bbduk.sh in \$version \" >>  $mylog
	echo \"Using adapters in $adapters \" >>  $mylog
	echo \"Using synthetic contaminants in $params.phix174ill and in $params.artifacts \" >>  $mylog
	
			
	if [ \"$params.library\" = \"paired-end\" ]; then
		unpairedR=\$(wc -l ${id}_trimmed_singletons_tmp.fq | cut -d\" \" -f 1)
		unpairedR=\$((\$unpairedR/4))
		echo  \"\$unpairedR singleton reads whose mate was trimmed shorter preserved\" >>  $mylog
		echo \" \" >>  $mylog
	fi
	#Defines command for removing synthetic contaminants
	if [ \"$params.library\" = \"paired-end\" ]; then
		bbduk.sh -Xmx\"\$maxmem\" in=${id}_trimmed_R1_tmp.fq in2=${id}_trimmed_R2_tmp.fq out=${id}_trimmed_R1.fq out2=${id}_trimmed_R2.fq k=31 ref=$phix174ill,$artifacts qin=$params.Pcoding threads=${task.cpus} ow
		fastqc --quiet --noextract --format fastq --outdir=. --threads ${task.cpus} ${id}_trimmed_R1.fq ${id}_trimmed_R2.fq
	else
		bbduk.sh -Xmx\"\$maxmem\" in=`basename $reads1 | sed -r 's/_R.{1,}//'`_trimmed_tmp.fq out=`basename $reads1 | sed -r 's/_R.{1,}//'`_trimmed.fq k=31 ref=$phix174ill,$artifacts qin=$params.Pcoding threads=${task.cpus} ow
		fastqc --quiet --noextract --format fastq --outdir=. --threads ${task.cpus} ${id}_trimmed_R1.fq
	fi
	#Removes synthetic contaminants and logs some figures (singleton read file, 
	#that exists iif the library layout was 'paired')
	if [ \"$params.library\" = \"paired-end\" ]; then
		bbduk.sh -Xmx\"\$maxmem\" in=${id}_trimmed_singletons_tmp.fq out=${id}_trimmed_singletons.fq k=31 ref=$phix174ill,$artifacts qin=$params.Pcoding threads=${task.cpus} ow
		
		
fi
	
	#Removes tmp files. This avoids adding them to the output channels
	rm -rf ${id}_trimmed*_tmp.fq



	"""
}



	to_make_bin =  Channel.value(Fastq_dir)



// ------------------------------------------------------------------------------
//                    		    Binning
// ------------------------------------------------------------------------------

process binning {
        publishDir workingdir, mode: 'copy', pattern: "*.txt"
        input:
        file(path) from to_make_bin
        
        output:
        file "OTUs.txt" into to_classify
	file "OTUs.txt"
	when:
	params.mode == "Taxonomic"
	script:
	"""
        #Measures execution time
        sysdate=\$(date)
        starttime=\$(date +%s.%N)
        echo \" Performing sequencing binning with DADA2 at \$sysdate\" >>  $mylog
        echo \" \" >>  $mylog
	CMD=\" Rscript --vanilla /storage/raid/home/m991833/16sMP/Pipeline/Scripts/dada3.R $workingdir/$path OTUs.txt \"
	exec \$CMD 2>&1 | tee tmp.log

	#Logs some data about sequence dereplication
	echo  \"Dada2 dereplication stats: \" >> $mylog
	echo \" \" >>  $mylog

	cat tmp.log >> $mylog



	"""

}

process complete_binning {
        publishDir workingdir, mode: 'copy', pattern: "*.txt"
        input:
        file(path) from to_bin.take(1)

        output:
        file "OTUs.txt"

        when:
        params.mode == "Complete"
        script:
        """
        #Measures execution time
        sysdate=\$(date)
        starttime=\$(date +%s.%N)
	echo \" Performing sequencing binning with DADA2  at \$sysdate\"
        echo \" \" >>  $mylog




        CMD=\" Rscript --vanilla /storage/raid/home/m991833/16sMP/Pipeline/Scripts/dada3.R $workingdir/$path OTUs.txt \"
        exec \$CMD 2>&1 | tee tmp.log

        #Logs some data about sequence dereplication
        echo  \"Dada2 dereplication stats: \" >> $mylog
        echo \" \" >>  $mylog



        cat tmp.log >> $mylog



        """

}

process classification {
        publishDir workingdir, mode: 'copy', pattern: "*{tsv,fasta}"
        input:
        file(OTU) from to_classify

        output:
        file "*.tsv"
	//file "*.fasta"

        when:
        params.mode == "Complete" || params.mode == "Taxonomic"
        script:
        """
        #Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	#echo \"Performing Taxonomic classification at \$sysdate\" >>  $mylog
	#echo \" \" >>  $mylog


	echo \"Retrieving OTU table and fasta file of representative sequences \"	
	DadatoOtu.py $OTU OTU.tsv Representatives.fasta
	
	echo \"Converting OTU.txt into OTU.biom \"
	#CMD = \" biom convert -i OTU.txt table.from_txt_json.biom --table-type='OTU table' --to-json -o OTU.biom \"
	#exec \$CMD 2>&1 | tee tmp.log


        """

}



// ------------------------------------------------------------------------------
//                     Quality assessment and visualization                    
// ------------------------------------------------------------------------------


//Creates the channel which performs the QC
toQC = rawreads 

//Process performing all the Quality Assessment
process qualityAssessment {
	
	publishDir  QC_dir, mode: 'copy', pattern: "*.{html,txt}"
	  	
	input:
   	set val(step), file(reads), val(label), val(stem) from toQC
	output:
	file "*_fastqc.html" 
	when:
	params.mode == "QC" | params.mode == "Complete"
   	script:
	"""	
	
	#Logs version of the software and executed command
	version=\$(fastqc --version) 
	fastqc --quiet --noextract --format fastq --outdir=. --threads ${task.cpus} $reads
	
	"""	
}

