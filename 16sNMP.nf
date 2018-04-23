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
	System.out.println("    --mode     <QC|complete>")
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
if (params.mode != "QC" && params.mode != "complete") {
	exit 1, "Mode not available. Choose any of <QC, characterisation, complete>"
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
workingpath = params.outdir + "/" + params.prefix + "/"
workingdir = file(workingpath)
if( !workingdir.exists() ) {
    if( !workingdir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $workingpath"
    } 
}

//Creating other folders
matrixpath = params.outdir + "/" + params.prefix + "/Matrix"
abundancepath = params.outdir + "/" + params.prefix + "/Matrix/Abundance"
read_count_path = params.outdir + "/" + params.prefix + "/Matrix/Read_count"
plotpath = params.outdir + "/" + params.prefix + "/Plots"
matrixdir = file(matrixpath)
abundancedir = file(abundancepath)
read_count_dir = file(read_count_path)
plotdir = file(plotpath)


if( !matrixdir.exists() ) {
    if( !matrixdir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $matrixpath"
    } 
}
if( !abundancedir.exists() ) {
    if( !abundancedir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $abundancepath"
    } 
}
if( !read_count_dir.exists() ) {
    if( !read_count_dir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $read_count_path"
    } 
}
if( !plotdir.exists() ) {
    if( !plotdir.mkdirs() ) {
	exit 1, "Cannot creat working directory: $plotpath"
    }
}


//Creates main log file
mylog = file(params.outdir + "/" + params.prefix + "/" + params.prefix + ".log")

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


	rawreads = Channel.fromPath(params.reads_R1) /**.buffer(size: 2, remainder: true)*/




// ------------------------------------------------------------------------------
//                             Microbiome barcode profiling                            
// ------------------------------------------------------------------------------

/**
        Microbiome profiling.

        Metaphlan2 will be used to anotate shotgun reads in taxonomic levels, this step will generate 3 outputs, a matrix containing the taxonomic levels, relative abundaces, read counts for each clade, sample size and average genome size in the working file and two others containing only the relative abundance and the read_counts in Matrix/Abundance and Matrix/Read_counts, respectively.
*/


process profileTaxa {





	
	input:
	file infile from rawreads
	
	


       

	script:
	"""

#	Rscript --vanilla Scripts/dada3.R $infile $params.reads_R1`basename $infile | sed 's/R1/R2/'` `basename $infile | sed 's/.fastq/.fq.gz/'` `basename $infile | sed 's/.fastq/.fq.gz/'| sed 's/R1/R2/'`
	#echo $infile $params.reads_R1`basename $infile | sed 's/R1/R2/'` `basename $infile | sed 's/.fastq/.fq.gz/'` `basename $infile | sed 's/.fastq/.fq.gz/'| sed 's/R1/R2/'` >> ../../../teste2.txt
	Rscript --vanilla /home/cafofo/Documentos/Pipeline/Scripts/dada3.R $params.reads_R1 $params.outdir

        """
}
