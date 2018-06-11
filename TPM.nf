/**
	Taxonomic Profiling Pipeline for Metabarcoding (TPM)
	Dr Bruno Andrade 	      
	
*/

version='1.0'
timestamp='20181903'

/**
	Prints version when asked for
*/
if (params.version) {
	System.out.println("")
	System.out.println("Taxonomic Profiling Pipeline for Metabarcoding (MT2P) - Version: $version ($timestamp)")
	exit 1
}

/**
	Prints help when asked for
*/

if (params.help) {
	System.out.println("")
	System.out.println("Taxonomic Profiling Pipeline for Metabarcoding (TPM) - Version: $version ($timestamp)")
	System.out.println("")
	System.out.println("Usage")
	System.out.println("   nextflow run -c nextflow.config NMP.nf --in Path --out path --m metadata --mode MODE")
	System.out.println("                [options] [-with-docker]")
	System.out.println("")
	System.out.println("Mandatory arguments:")
	System.out.println("    --input	Input Path	Path for all reads (library = paired-end or single-end")
	System.out.println("    --metadata	metadata	Sample metadata")
	System.out.println("    --outdir	path	Output directory (will be outdir/prefix/date)")
	System.out.println("    --mode	<QC|Classification|Complete>")
	System.out.println("Options:")
	System.out.println("    --l	Library	<single-end|paired-end>")
	System.out.println("")
	System.out.println("Container:")
	System.out.println("    Docker image to use with docker is")
	System.out.println("    'docker://bgabriel/MT2P'")
    exit 1
}

	
/**

	
	Parameters check and creation of folders for the analysis


*/


//Checking user-defined parameters	
if (params.mode != "QC" && params.mode != "Classification" && params.mode != "Complete") {
	exit 1, "Mode not available. Choose any of <QC, Classification, Complete>"
}	

if (params.library != "paired-end" && params.library != "single-end") { 
	exit 1, "Library layout not available. Choose any of <single-end, paired-end>" 
}   

if (params.metadata == null){
	exit 1, "A metadata file is required"
}

//Creates working dir

Working_dir = file(params.outdir)
Matrix_dir = file(params.outdir + "/Matrix")
QC_dir = file(params.outdir + "/QC")
Fastq_dir = file(params.outdir + "/Fastq")
Qiime_dir = file(params.outdir + "/Qiime")
Plot_dir = file(params.outdir + "/Plots")
Log_dir = file(params.outdir + "/Logs")


if( !Working_dir.exists() ) {
    if( !Working_dir.mkdirs() ) 	{
        exit 1, "Cannot create working directory"
    } 
}
if( !Matrix_dir.exists() ) {
    if( !Matrix_dir.mkdirs() ) 	{
        exit 1, "Cannot create Matrix  folder"
    } 
}
if( !QC_dir.exists() ) {
    if( !QC_dir.mkdirs() ) 	{
        exit 1, "Cannot create QC folder"
    } 
}
if( !Fastq_dir.exists() ) {
    if( !Fastq_dir.mkdirs() ) 	{
        exit 1, "Cannot create Fastq folder"
 } 
}
if( !Qiime_dir.exists() ) {
    if( !Qiime_dir.mkdirs() )   {
        exit 1, "Cannot create Qiime folder"
 }
}
if( !Plot_dir.exists() ) {
    if( !Plot_dir.mkdirs() )   {
        exit 1, "Cannot create Plot folder"
 }
}
if( !Log_dir.exists() ) {
    if( !Log_dir.mkdirs() )   {
        exit 1, "Cannot create Log folder"
 }
}


	   


/**
	Quality assessment of raw reads using the FastQC software
*/

	if (params.library == "paired-end") {
		rawreads = Channel.fromFilePairs(params.input + '*R{1,2}*fastq')
	}
	else {
		rawreads = Channel.fromPath(params.input + '*.fastq')
	}



// ------------------------------------------------------------------------------
//                             Trimming and Filtering
// ------------------------------------------------------------------------------

/**
	The trimming process will be performed using the bbduk software

*/


        if (params.library == "paired-end") {
                to_trim = Channel.fromFilePairs(params.input + '*R{1,2}*{fq,fastq}')
        }
        else {
                to_trim = Channel.fromPath(params.input + '*.fastq')
        }

process trim {
	publishDir Fastq_dir, mode: 'copy', pattern: "*_trimmed*{fq,fastq}"
	publishDir QC_dir, mode: 'copy', pattern: "*.html"
	publishDir Log_dir, mode: 'copy', pattern: "Trim_log.txt"
	input:
	set val(id), file(reads1), file(reads2) from to_trim
	file(adapters) from file(params.adapters)
	file(artifacts) from file(params.artifacts)
	file(phix174ill) from file(params.phix174ill)

	output:
	file("*_trimmed*.fq") into todecontaminate
	file("*_trimmed*.fq") into trimmedreads
	file "*_fastqc.html"
	file "Trim_log.txt"
	val Fastq_dir into to_bin
	when:
	params.mode == "QC" || params.mode == "Complete"
      script:
      """
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing Quality Control. STEP 2 [Trimming] at \$sysdate\" >> Trim_log.txt
	echo \" \"
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
	cat Trim_log.txt tmp.log
	"""
}


	//to_make_bin =  Channel.value(Fastq_dir)


// ------------------------------------------------------------------------------
//                    		    Binning
// ------------------------------------------------------------------------------

/**
        Instead of using the binning algorithms found in Qiime1, we decided to use the R package Dada2.

*/


        if (params.mode == "Complete") {
                to_make_bin = to_bin.take(1)
        }
	else { 
		to_make_bin = Channel.value(Fastq_dir)
	}




process binning {
        publishDir Matrix_dir, mode: 'copy', pattern: "*.txt"
	publishDir Log_dir, mode: 'copy', pattern: "Log_binning.txt"
        input:
        file(path) from to_make_bin

        output:
        file "OTUs.txt" into to_classify
	file "Log_binning.txt"
	when:
	params.mode == "Classification" || params.mode == "Complete"
	script:
	"""
        #Measures execution time
        sysdate=\$(date)
        starttime=\$(date +%s.%N)
        echo \" Performing sequencing binning with DADA2 at \$sysdate\" >> Log_binning.txt
        echo \" \" >> Log_binning.txt
	CMD=\" Rscript /opt/Scripts/dada2.R $Working_dir/$path OTUs.txt \"
	exec \$CMD 2>&1 | tee -a Log_binning.txt 
	sed -ri 's/Not all.{1,}//g' Log_binning.txt
	sed -ri 's/L00(1|2)_trimmed_R1.(fq|fastq)//g' OTUs.txt
	

	#cat Log_binning.txt tmp.log
	"""

}


to_Reference = Channel.from(params.table)


// ------------------------------------------------------------------------------
//                             Microbiome barcode profiling
// ------------------------------------------------------------------------------

/**
        Microbiome profiling using Qiime.
        Qiime will be used to profile 16s/18s reads in taxonomic levels, this step will generate several outputs, an OTU table with taxonomic information, a phylogenetic tree and barplots of relative abundance and tables of summarized taxonomic distribution in every taxonomic level.
*/


process classification {
        publishDir Matrix_dir, mode: 'copy', pattern: "*.{tsv,fasta,biom,tre}"
	publishDir Plot_dir, mode: 'copy', pattern: "*.{png,pdf}"
	publishDir Log_dir, mode: 'copy', pattern: "log_Class.txt"

	input:
        file(Input_class) from to_classify
        file(Ref_fasta) from file(params.reference)
        file(Tax_table) from file(params.table)
        file(Meta) from file(params.metadata)

        output:
        file "OTU.biom" into to_diversity
        file "Representatives.tre" into to_tree
        file "OTU.tsv"
        file "Representatives.fasta"
        file "OTU.biom"
        file "Representatives.tre"
        file "*.{png,pdf}"
        file "log_Class.txt"

        when:
        params.mode == "Complete" || params.mode == "Classification"
        script:
        """
        #Measures execution time
        sysdate=\$(date)
        starttime=\$(date +%s.%N)
        echo \"Performing Taxonomic classification at \$sysdate\" >> log_Class.txt
        echo \" \" >>  log_Class.txt
        echo \"Retrieving OTU table and representative sequences \" >> log_Class.txt
        echo \" \" >>  log_Class.txt
        DadatoOtu.py $Input_class OTU.tsv Representatives.fasta
        echo \"Done\" >> log_Class.txt

        echo \"Converting OTU.txt into OTU.biom \" >> log_Class.txt
        biom convert -i OTU.tsv --table-type='OTU table' --to-json -o OTU.biom
        echo \"Done\" >> log_Class.txt

        echo \" \" >> log_Class.txt
        echo \"Performing Taxonomic classification with Qiime \" >> Log_Class.txt
	CMD=\" assign_taxonomy.py -i Representatives.fasta -r $Ref_fasta -t $Tax_table -o $Qiime_dir \"
	exec \$CMD 2>&1 | tee -a Log_Class.txt
        echo \" \" >> Log_Class.txt
        echo \"Aligning representative sequences with Mafft using ${task.cpus} cores\" >> Log_Class.txt
	mafft --thread ${task.cpus} Representatives.fasta > Representatives_aligned.fasta
	echo \"Done \" >> Log_Class.txt
        echo \"Building tree using $params.tree method \"
	CMD=\" make_phylogeny.py -i Representatives_aligned.fasta -t $params.tree -o Representatives.tre \"
        exec \$CMD 2>&1 | tee -a Log_Class.txt
	echo \"Done\" >>Log_Class.txt
	echo \" \" >> Log_Class.txt
        echo \"Generating the full OTU table, summarized tables and barplots \" >> Log_Class.txt
        merge_taxonomy.py OTU.tsv $Qiime_dir/Representatives_tax_assignments.txt OTU_tax.tsv ./ $params.RelP $params.Fig

	"""
}



process diversity {

	publishDir Matrix_dir, mode: 'copy', pattern: "{weighted,unweighted}*.txt"
	publishDir Qiime_dir, mode: 'copy', pattern: "OTU_table_even100.biom"
        publishDir Plot_dir, mode: 'copy', pattern: "*{png,pdf}"
	publishDir Log_dir, mode: 'copy', pattern: "Log_Diversity.txt"

        input:
        file(OTU) from to_diversity
        file(Tree) from to_tree
        file(Metadata) from file(params.metadata)

        output:
        file "*.biom"
        file "*{.pdf,.png}"
	file "Log_Diversity.txt"
	file "{weighted,unweighted}*.txt"

        when:
        params.mode == "Complete" || params.mode == "Classification"
        script:
        """
        #Measures execution time
        sysdate=\$(date)
        starttime=\$(date +%s.%N)
        echo \"Performing Normalization and Diversity analyses at \$sysdate\" >> Log_Diversity.txt
        echo \" \" >>  log_Diversity.txt
        #echo \"Starting single rarefaction with \$params.Depth\" >> Log_Diversity.txt
        single_rarefaction.py -i $OTU -o OTU_table_${params.Depth}.biom -d $params.Depth
	echo \"Done\" >> Log_Diversity.txt
	echo \" \" >>  Log_Diversity.txt
        #echo \"Starting Alpha diversity analysis with \$params.alpha_metrics\" >> Log_Diversity.txt
        alpha_diversity.py -i $OTU -m $params.alpha_metrics -t $Tree -o Alpha_diversity${params.alpha_metrics}.txt
	
        echo \" \" >>  Log_Diversity.txt
        echo \"Starting Beta diversity analysis using Weichted unifrac metrics\" >> Log_Diversity.txt
        beta_diversity.py -i OTU_table_${params.Depth}.biom -m weighted_unifrac -t $Tree -o .
        echo \" \" >>  Log_Diversity.txt
        
	echo \"Starting Beta diversity analysis using Unweighted unifrac metrics\" >> Log_Diversity.txt
        beta_diversity.py -i OTU_table_${params.Depth}.biom -m unweighted_unifrac -t $Tree -o .

        echo \" \" >>  Log_Diversity.txt
        echo \"Estimating PCoA\" >> Log_Diversity.txt
        python3 /opt/Scripts/PCoA.py weighted_unifrac_OTU_table_${params.Depth}.txt $Metadata ./ ${params.Fig}
	python3 /opt/Scripts/PCoA.py unweighted_unifrac_OTU_table_${params.Depth}.txt $Metadata ./ ${params.Fig}

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
