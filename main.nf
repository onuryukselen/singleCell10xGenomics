$HOSTNAME = ""
params.outdir = 'results'  


if (!params.cellBarcodeFile){params.cellBarcodeFile = ""} 
if (!params.mate_split){params.mate_split = ""} 
if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.cellBarcodePattern){params.cellBarcodePattern = ""} 
if (!params.cutoff_for_reads_per_cell){params.cutoff_for_reads_per_cell = ""} 
if (!params.gtf_url){params.gtf_url = ""} 
if (!params.genome_url){params.genome_url = ""} 
if (!params.commondb_url){params.commondb_url = ""} 

Channel.value(params.cellBarcodeFile).set{g_7_cellBarcodeFile_g_72}
Channel.value(params.mate_split).into{g_13_mate_g116_11;g_13_mate_g117_3;g_13_mate_g117_11;g_13_mate_g115_16}
Channel
	.fromFilePairs( params.reads , size: (params.mate != "pair") ? 1 : 2 )
	.ifEmpty { error "Cannot find any read_pairs matching: ${params.reads}" }
	.set{g_27_read_pairs_g_72}

Channel.value(params.mate).into{g_38_mate_g_72;g_38_mate_g_106}
Channel.value(params.cellBarcodePattern).set{g_73_cellBarcodePattern_g_72}
Channel.value(params.cutoff_for_reads_per_cell).into{g_111_cutoff_g122_90;g_111_cutoff_g122_91;g_111_cutoff_g123_90;g_111_cutoff_g123_91;g_111_cutoff_g124_90;g_111_cutoff_g124_91}
g_118_gtf_url_g114_15 = file(params.gtf_url, type: 'any') 
g_119_genome_url_g114_15 = file(params.genome_url, type: 'any') 
Channel.value(params.commondb_url).set{g_120_commondb_url_g114_15}

UMIqualityFilterThreshold = params.UmiExtract_10xGenomics.UMIqualityFilterThreshold
insertSide = params.UmiExtract_10xGenomics.insertSide


if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2000
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "long"
}

process UmiExtract_10xGenomics {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_valid.fastq$/) "valid_fastq/$filename"
}

input:
 val cellBarcode from g_7_cellBarcodeFile_g_72
 val mate from g_38_mate_g_72
 set val(name), file(reads) from g_27_read_pairs_g_72
 val cellBarcodePattern from g_73_cellBarcodePattern_g_72

output:
 set val(name), file("*_valid.fastq")  into g_72_valid_fastq_g_106

script:
readArray = reads.toString().split(' ')
"""
if [ "${mate}" == "pair" ]; then
umi_tools extract --bc-pattern='${cellBarcodePattern}' \
                  --extract-method=regex \
                  --stdin ${readArray[0]} \
                  --stdout out${name}_R1.fastq \
                  --read2-in ${readArray[1]} \
                  --read2-out=out${name}_R2.fastq \
                  --filter-cell-barcode \
                  --whitelist=${cellBarcode} \
				  --quality-filter-threshold=${UMIqualityFilterThreshold} \
				  --quality-encoding=phred33
mv out${name}_${insertSide}.fastq ${name}_valid.fastq
else
umi_tools extract --bc-pattern='${cellBarcodePattern}' \
                  --extract-method=regex \
                  --stdin ${reads} \
                  --stdout ${name}_valid.fastq \
                  --filter-cell-barcode \
                  --whitelist=${cellBarcode} \
				  --quality-filter-threshold= ${UMIqualityFilterThreshold} \
				  --quality-encoding=phred33
fi
sed -i 's/_/:/g' *valid.fastq
"""
}

params.run_Split_Fastq =  "no"  //* @dropdown @options:"yes","no" @show_settings:"SplitFastq"
readsPerFile = params.SplitFastq.readsPerFile
//Since splitFastq operator requires flat file structure, first convert grouped structure to flat, execute splitFastq, and then return back to original grouped structure
//.map(flatPairsClosure).splitFastq(splitFastqParams).map(groupPairsClosure)

//Mapping grouped read structure to flat structure
flatPairsClosure = {row -> if(row[1] instanceof Collection) {
        if (row[1][1]){
            tuple(row[0], file(row[1][0]), file(row[1][1]))
        } else {
            tuple(row[0], file(row[1][0]))
        }
    } else {
        tuple(row[0], file(row[1]))
    }
}

//Mapping flat read structure to grouped read structure
groupPairsClosure = {row -> tuple(row[0], (row[2]) ? [file(row[1]), file(row[2])] : [file(row[1])])}

// if mate of split process different than rest of the pipeline, use "mate_split" as input parameter. Otherwise use default "mate" as input parameter
mateParamName = (params.mate_split) ? "mate_split" : "mate"
splitFastqParams = ""
if (params[mateParamName] != "pair"){
    splitFastqParams = [by: readsPerFile, file:true]
}else {
    splitFastqParams = [by: readsPerFile, pe:true, file:true]
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill
if (!(params.run_Split_Fastq == "yes")){
g_72_valid_fastq_g_106.into{g_106_reads_g115_16; g_106_reads_g116_11; g_106_reads_g117_11}
} else {


process SplitFastq {

input:
 val mate from g_38_mate_g_106
 set val(name), file(reads) from g_72_valid_fastq_g_106.map(flatPairsClosure).splitFastq(splitFastqParams).map(groupPairsClosure)

output:
 set val(name), file("split/*")  into g_106_reads_g115_16, g_106_reads_g116_11, g_106_reads_g117_11

when:
params.run_Split_Fastq == "yes"

script:
"""    
mkdir -p split
mv ${reads} split/.
"""
}
}


params.gtf =  ""  //* @input
params.genome =  ""  //* @input
params.commondb =  ""  //* @input

process Check_and_Build_Module_Check_Genome_GTF {

input:
 file fasta from g_119_genome_url_g114_15
 file downGtf from g_118_gtf_url_g114_15
 val commondb_url from g_120_commondb_url_g114_15

output:
 val "${params.genome}"  into g114_15_genomePath_g114_0, g114_15_genomePath_g114_6, g114_15_genomePath_g114_8, g114_15_genomePath_g114_10, g114_15_genomePath_g114_5, g114_15_genomePath_g114_13
 val "${params.gtf}"  into g114_15_gtfPath_g114_0, g114_15_gtfPath_g114_6, g114_15_gtfPath_g114_8, g114_15_gtfPath_g114_10, g114_15_gtfPath_g114_4, g114_15_gtfPath_g114_13
 val "${params.commondb}"  into g114_15_commondb_path_g114_18

when:
params.run_checkAndBuild == "yes"

script:
gtf_dir  = params.gtf.substring(0, params.gtf.lastIndexOf('/')) 
genome_dir  = params.genome.substring(0, params.genome.lastIndexOf('/')) 
slashCount = commondb_url.count("/")
cutDir = slashCount - 3;

"""
downGenomePath=\$(realpath $fasta)
downGtfPath=\$(realpath $downGtf)
if [ ! -e "${params.genome}" ] ; then
    echo "${params.genome} not found"
    mkdir -p ${genome_dir}
    cp -n \$downGenomePath ${params.genome}
fi
if [ ! -e "${params.gtf}" ] ; then
    echo "${params.gtf} not found"
    mkdir -p ${gtf_dir}
    cp -n \$downGtfPath ${params.gtf}
fi
if [ ! -e "${params.commondb}" ] ; then
    echo "${params.commondb} not found"
    mkdir -p ${params.commondb}
    wget -l inf -nc -nH --cut-dirs=$cutDir -R 'index.html*' -r --no-parent --directory-prefix=${params.commondb} $commondb_url
fi

"""




}


process Check_and_Build_Module_Bowtie_Index {

input:
 val genome from g114_15_genomePath_g114_13
 val gtf from g114_15_gtfPath_g114_13

output:
 val resultDir  into g114_13_genomeIndexPath_g114_18

when:
(params.use_Bowtie_Index == "yes") && (params.run_checkAndBuild == "yes")

script:
bowtie_build_parameters = params.Check_and_Build_Module_Bowtie_Index.bowtie_build_parameters
basedir  = genome.substring(0, genome.lastIndexOf('/'))
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
basename = genome.substring(genome.lastIndexOf('/')+1,genome.lastIndexOf('.'))
filename = genome.substring(genome.lastIndexOf('/')+1,genome.length())
newDirName = "BowtieIndex"
resultDir = basedir.substring(0, basedir.lastIndexOf('/')) +"/"+ newDirName 
tmpResultDir = basedir.substring(0, basedir.lastIndexOf('/')) +"/_tmp_"+ newDirName 

"""
if [ ! -e "${resultDir}/${basename}.rev.2.ebwt" ] ; then
    echo "${resultDir}/${basename}.rev.2.ebwt Bowtie index not found"
    rm -rf $tmpResultDir $resultDir && mkdir -p $tmpResultDir && cd $tmpResultDir
    ln -s ../main/${filename} ${filename}
    bowtie-build ${bowtie_build_parameters} ${filename} ${basename}
    cd .. && mv $tmpResultDir $resultDir
fi
"""




}

params.gtf2bed_path =  ""  //* @input

process Check_and_Build_Module_Check_GTF2BED12 {

input:
 val gtf from g114_15_gtfPath_g114_4


when:
params.run_checkAndBuild == "yes"

script:
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
"""
if [ ! -e "${params.bed}" ] ; then
    echo "${params.bed} not found"
    perl ${params.gtf2bed_path} $gtf > ${params.bed}
fi
"""




}

params.gtf2bed_path =  ""  //* @input

process Check_and_Build_Module_Check_chrom_sizes_and_index {

input:
 val genome from g114_15_genomePath_g114_5


when:
params.run_checkAndBuild == "yes"

script:
genome_dir  = genome.substring(0, genome.lastIndexOf('/')) 
basename_and_path  = genome.substring(0, genome.lastIndexOf('.'))

"""
if [ ! -e "${params.genome_sizes}" ] ; then
    echo "${params.genome_sizes} not found"
    cat ${genome} | awk '\$0 ~ ">" {print c; c=0;printf substr(\$0,2,100) "\\t"; } \$0 !~ ">" {c+=length(\$0);} END { print c; }' > ${basename_and_path}.chrom.sizes
    ##clean first empty line
    sed -i '1{/^\$/d}' ${basename_and_path}.chrom.sizes
fi
"""




}


process Check_and_Build_Module_Check_Build_Rsem_Index {

input:
 val genome from g114_15_genomePath_g114_10
 val gtf from g114_15_gtfPath_g114_10

output:
 val cmdAr  into g114_10_command_g114_11
 val resultDirAr  into g114_10_path_g114_11

when:
(params.use_RSEM_Index == "yes") && (params.run_checkAndBuild == "yes")

script:
create_bowtie_rsem_index = params.Check_and_Build_Module_Check_Build_Rsem_Index.create_bowtie_rsem_index
create_bowtie2_rsem_index = params.Check_and_Build_Module_Check_Build_Rsem_Index.create_bowtie2_rsem_index
create_star_rsem_index = params.Check_and_Build_Module_Check_Build_Rsem_Index.create_star_rsem_index
transcript_to_gene_map = params.Check_and_Build_Module_Check_Build_Rsem_Index.transcript_to_gene_map
RSEM_build_parameters = params.Check_and_Build_Module_Check_Build_Rsem_Index.RSEM_build_parameters

genome_dir  = genome.substring(0, genome.lastIndexOf('/'))
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
basenameGenome = genome.substring(genome.lastIndexOf('/')+1,genome.lastIndexOf('.'))
newDirNameAr = []
cmdAr = []
resultDirAr = []
if (create_bowtie_rsem_index == "true"){ newDirNameAr.push('RSEM_ref_Bowtie') }
if (create_bowtie2_rsem_index == "true"){ newDirNameAr.push('RSEM_ref_Bowtie2') }
if (create_star_rsem_index == "true"){ newDirNameAr.push('RSEM_ref_STAR') }

transcript_to_gene_mapText = ""
if (transcript_to_gene_map?.trim()){
    transcript_to_gene_mapText = "--transcript-to-gene-map " + transcript_to_gene_map
}

for (i = 0; i < newDirNameAr.size(); i++) {
    resultDir = gtf_dir.substring(0, gtf_dir.lastIndexOf('/')) +"/"+ newDirNameAr[i]
    tmpResultDir = gtf_dir.substring(0, gtf_dir.lastIndexOf('/')) +"/_tmp_"+ newDirNameAr[i]
    resultDirAr.push(resultDir)
    cmd = ""
    indexType = ""
    if (newDirNameAr[i] == 'RSEM_ref_Bowtie'){
        indexType = "--bowtie "
        checkFile = "${basenameGenome}.rev.2.ebwt" 
    } else if (newDirNameAr[i] == 'RSEM_ref_Bowtie2'){
        indexType = "--bowtie2 "
        checkFile = "${basenameGenome}.rev.2.bt2" 
    } else if (newDirNameAr[i] == 'RSEM_ref_STAR'){
        indexType = "--star "
        checkFile = "genomeParameters.txt" 
    }
    cmd = "if [ ! -e \"${resultDir}/${checkFile}\" ] ; then rm -rf $tmpResultDir $resultDir && mkdir -p $tmpResultDir && cd $tmpResultDir && rsem-prepare-reference ${RSEM_build_parameters} --gtf ${gtf} ${transcript_to_gene_mapText} ${indexType} ${genome} ${basenameGenome} && cd .. && mv $tmpResultDir $resultDir; fi"
    cmdAr.push(cmd)
}


"""

"""

}


process Check_and_Build_Module_Build_Index_RSEM_run {

input:
 val resultDir from g114_10_path_g114_11.flatten()
 val command from g114_10_command_g114_11.flatten()

output:
 val resultDir  into g114_11_genomeIndexPath

script:
"""    
$command
"""
}


process Check_and_Build_Module_Hisat2_Index {

input:
 val genome from g114_15_genomePath_g114_8
 val gtf from g114_15_gtfPath_g114_8

output:
 val resultDir  into g114_8_genomeIndexPath_g116_11

when:
(params.use_Hisat2_Index == "yes") && (params.run_checkAndBuild == "yes")

script:
hisat2_build_parameters = params.Check_and_Build_Module_Hisat2_Index.hisat2_build_parameters
genome_dir  = genome.substring(0, genome.lastIndexOf('/'))
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
basenameGenome = genome.substring(genome.lastIndexOf('/')+1,genome.lastIndexOf('.'))
basenameGTF = gtf.substring(gtf.lastIndexOf('/')+1,gtf.lastIndexOf('.'))
filename = genome.substring(genome.lastIndexOf('/')+1,genome.length())
newDirName = "Hisat2Index"
resultDir = gtf_dir.substring(0, gtf_dir.lastIndexOf('/')) +"/"+ newDirName 
tmpResultDir = gtf_dir.substring(0, gtf_dir.lastIndexOf('/')) +"/_tmp_"+ newDirName 

extract_splice_sites = "hisat2_extract_splice_sites.py ${gtf} > ${tmpResultDir}/${basenameGTF}.hisat2_splice_sites.txt"
extract_exons = "hisat2_extract_exons.py ${gtf}> ${tmpResultDir}/${basenameGTF}.hisat2_exons.txt"
ss = "--ss ${basenameGTF}.hisat2_splice_sites.txt"
exon = "--exon ${basenameGTF}.hisat2_exons.txt"

"""
if [ ! -e "${resultDir}/${basenameGenome}.8.ht2" ] ; then
    echo "${resultDir}/${basenameGenome}.8.ht2 Hisat2 index not found"
    rm -rf $tmpResultDir $resultDir && mkdir -p $tmpResultDir && cd $tmpResultDir 
    $extract_splice_sites
    $extract_exons
    hisat2-build ${hisat2_build_parameters} $ss $exon ${genome} ${basenameGenome}
    cd .. && mv $tmpResultDir $resultDir 
fi
"""

}

g114_8_genomeIndexPath_g116_11= g114_8_genomeIndexPath_g116_11.ifEmpty(file('hisat2index', type: 'any')) 

params.hisat2_index =  ""  //* @input


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 3
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 100
    $CPU  = 4
    $MEMORY = 32
    $QUEUE = "short"
} 
//* platform
//* autofill

process HISAT2_Module_Map_HISAT2 {

input:
 val mate from g_13_mate_g116_11
 set val(name), file(reads) from g_106_reads_g116_11
 val hisat2index from g114_8_genomeIndexPath_g116_11

output:
 set val(name), file("${newName}.bam")  into g116_11_mapped_reads_g116_13
 set val(name), file("${newName}.align_summary.txt")  into g116_11_outputFileTxt_g116_2
 set val(name), file("${newName}.flagstat.txt")  into g116_11_outputFileOut

when:
(params.run_HISAT2 && (params.run_HISAT2 == "yes")) || !params.run_HISAT2

script:
HISAT2_parameters = params.HISAT2_Module_Map_HISAT2.HISAT2_parameters
nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}

"""
$runGzip
if [ "${mate}" == "pair" ]; then
    hisat2 ${HISAT2_parameters} -x ${params.hisat2_index} -1 ${file1} -2 ${file2} -S ${newName}.sam &> ${newName}.align_summary.txt
else
    hisat2 ${HISAT2_parameters} -x ${params.hisat2_index} -U ${file1} -S ${newName}.sam &> ${newName}.align_summary.txt
fi
samtools view -bS ${newName}.sam > ${newName}.bam
samtools flagstat ${newName}.bam > ${newName}.flagstat.txt
"""

}


//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2000
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill

process HISAT2_Module_Merge_Bam {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_sorted.*bai$/) "sorted_bam_hisat2/$filename"
	else if (filename =~ /.*_sorted.*bam$/) "sorted_bam_hisat2/$filename"
}

input:
 set val(oldname), file(bamfiles) from g116_11_mapped_reads_g116_13.groupTuple()

output:
 set val(oldname), file("${oldname}.bam")  into g116_13_merged_bams
 set val(oldname), file("*_sorted*bai")  into g116_13_bam_index
 set val(oldname), file("*_sorted*bam")  into g116_13_sorted_bam_g123_92

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
    samtools merge !{oldname}.bam !{bamfiles.join(" ")} && samtools sort -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
else
    mv !{bamfiles.join(" ")} !{oldname}.bam 2>/dev/null || true
    samtools sort  -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
fi
'''
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill

process Single_Cell_Module_after_HISAT2_samtools_sort_index {

input:
 set val(name), file(bam) from g116_13_sorted_bam_g123_92

output:
 set val(name), file("bam/*.bam")  into g123_92_bam_file_g123_90
 set val(name), file("bam/*.bai")  into g123_92_bam_index_g123_90

script:
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = "samtools index " + bam 
} else {
    runSamtools = "samtools sort -o " + name +"_sorted " + bam +" && samtools index " + name + "_sorted.bam "
}
"""
$runSamtools
mkdir -p bam && mv *_sorted.ba* bam/.
"""
}

bowtie2_build_parameters = params.Check_and_Build_Module_Bowtie2_Index.bowtie2_build_parameters

process Check_and_Build_Module_Bowtie2_Index {

input:
 val genome from g114_15_genomePath_g114_6
 val gtf from g114_15_gtfPath_g114_6

output:
 val resultDir  into g114_6_genomeIndexPath_g114_18, g114_6_genomeIndexPath_g117_11

when:
(params.use_Bowtie2_Index == "yes") && (params.run_checkAndBuild == "yes")

script:

basedir  = genome.substring(0, genome.lastIndexOf('/'))
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
basename = genome.substring(genome.lastIndexOf('/')+1,genome.lastIndexOf('.'))
filename = genome.substring(genome.lastIndexOf('/')+1,genome.length())
newDirName = "Bowtie2Index"
resultDir = basedir.substring(0, basedir.lastIndexOf('/')) +"/"+ newDirName 
tmpResultDir = basedir.substring(0, basedir.lastIndexOf('/')) +"/_tmp_"+ newDirName 

"""
if [ ! -e "${resultDir}/${basename}.rev.1.bt2" ] ; then
    echo "${resultDir}/${basename}.rev.1.bt2 Bowtie2 index not found"
    rm -rf $tmpResultDir $resultDir && mkdir -p $tmpResultDir && cd $tmpResultDir
    ln -s ../main/${filename} ${filename}
    bowtie2-build ${bowtie2_build_parameters} ${filename} ${basename}
    cd .. && mv $tmpResultDir $resultDir 
fi
"""



}

g114_6_genomeIndexPath_g117_11= g114_6_genomeIndexPath_g117_11.ifEmpty(file('tophat2_index', type: 'any')) 

params.bowtie2_index =  ""  //* @input
params.gtf =  ""  //* @input


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 3
    $MEMORY = 24
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2500
    $CPU  = 4
    $MEMORY = 24
    $QUEUE = "long"
}
//* platform
//* autofill

process Tophat2_Module_Map_Tophat2 {

input:
 val mate from g_13_mate_g117_11
 set val(name), file(reads) from g_106_reads_g117_11
 val tophat2_index from g114_6_genomeIndexPath_g117_11

output:
 set val(name), file("${newName}.bam")  into g117_11_mapped_reads_g117_13
 set val(name), file("${newName}_unmapped.bam")  into g117_11_unmapped_reads
 set val(name), file("${newName}_align_summary.txt")  into g117_11_summary_g117_3

when:
(params.run_Tophat && (params.run_Tophat == "yes")) || !params.run_Tophat

script:
tophat2_parameters = params.Tophat2_Module_Map_Tophat2.tophat2_parameters
nameAll = reads.toString()
nameArray = nameAll.split(' ')

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file =  nameAll 
    runGzip = ''
}

"""
$runGzip
if [ "${mate}" == "pair" ]; then
    tophat2 ${tophat2_parameters}  --keep-tmp -G ${params.gtf} -o . ${params.bowtie2_index} $file
else
    tophat2 ${tophat2_parameters}  --keep-tmp -G ${params.gtf} -o . ${params.bowtie2_index} $file
fi

if [ -f unmapped.bam ]; then
    mv unmapped.bam ${newName}_unmapped.bam
else
    touch ${newName}_unmapped.bam
fi

mv accepted_hits.bam ${newName}.bam
mv align_summary.txt ${newName}_align_summary.txt
"""

}


//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2000
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill

process Tophat2_Module_Merge_Bam {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_sorted.*bai$/) "sorted_bam_tophat2/$filename"
	else if (filename =~ /.*_sorted.*bam$/) "sorted_bam_tophat2/$filename"
}

input:
 set val(oldname), file(bamfiles) from g117_11_mapped_reads_g117_13.groupTuple()

output:
 set val(oldname), file("${oldname}.bam")  into g117_13_merged_bams
 set val(oldname), file("*_sorted*bai")  into g117_13_bam_index
 set val(oldname), file("*_sorted*bam")  into g117_13_sorted_bam_g122_92

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
    samtools merge !{oldname}.bam !{bamfiles.join(" ")} && samtools sort -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
else
    mv !{bamfiles.join(" ")} !{oldname}.bam 2>/dev/null || true
    samtools sort  -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
fi
'''
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill

process Single_Cell_Module_after_Tophat2_samtools_sort_index {

input:
 set val(name), file(bam) from g117_13_sorted_bam_g122_92

output:
 set val(name), file("bam/*.bam")  into g122_92_bam_file_g122_90
 set val(name), file("bam/*.bai")  into g122_92_bam_index_g122_90

script:
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = "samtools index " + bam 
} else {
    runSamtools = "samtools sort -o " + name +"_sorted " + bam +" && samtools index " + name + "_sorted.bam "
}
"""
$runSamtools
mkdir -p bam && mv *_sorted.ba* bam/.
"""
}


process Check_and_Build_Module_STAR_Index_Check_Build {

input:
 val genome from g114_15_genomePath_g114_0
 val gtf from g114_15_gtfPath_g114_0

output:
 val resultDir  into g114_0_genomeIndexPath_g114_18, g114_0_genomeIndexPath_g115_16

when:
(params.use_STAR_Index == "yes") && (params.run_checkAndBuild == "yes")

script:
star_build_parameters = params.Check_and_Build_Module_STAR_Index_Check_Build.star_build_parameters
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
indexbasedir  = gtf_dir.substring(0, gtf_dir.lastIndexOf('/'))
genome_dir  = genome.substring(0, genome.lastIndexOf('/')) 
filename = genome.substring(genome.lastIndexOf('/')+1,genome.length())
newDirName = "STARIndex" 
resultDir = indexbasedir +"/"+ newDirName 
tmpResultDir = indexbasedir +"/_tmp_"+ newDirName
"""
if [ ! -e "${resultDir}/SA" ] ; then
    echo "STAR index not found"
    rm -rf $tmpResultDir ${resultDir} && mkdir -p $tmpResultDir && cd $tmpResultDir
    STAR --runMode genomeGenerate ${star_build_parameters} --genomeDir $tmpResultDir --genomeFastaFiles ${genome} --sjdbGTFfile ${gtf}
    cd .. && mv $tmpResultDir $resultDir && cd ${resultDir}
    ln -s ../../main/${filename} ${filename}
    ln -s ../../main/${filename}.fai ${filename}.fai
fi
"""





}

g114_0_genomeIndexPath_g115_16= g114_0_genomeIndexPath_g115_16.ifEmpty(file('star_index', type: 'any')) 

params.star_index =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 3
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 3
    $MEMORY = 32
    $QUEUE = "long"
} 
//* platform
//* autofill

process STAR_Module_Map_STAR {

input:
 val mate from g_13_mate_g115_16
 set val(name), file(reads) from g_106_reads_g115_16
 val star_index from g114_0_genomeIndexPath_g115_16

output:
 set val(name), file("${newName}Log.final.out")  into g115_16_outputFileOut_g115_18
 set val(name), file("${newName}.flagstat.txt")  into g115_16_outputFileTxt
 set val(name), file("${newName}Log.out")  into g115_16_logOut_g115_18
 set val(name), file("${newName}.bam")  into g115_16_mapped_reads_g115_14
 set val(name), file("${newName}SJ.out.tab")  into g115_16_outputFileTab_g115_18
 set val(name), file("${newName}Log.progress.out")  into g115_16_progressOut_g115_18
 set val(name), file("${newName}Aligned.toTranscriptome.out.bam") optional true  into g115_16_transcriptome_bam_g115_15

when:
(params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR

script:
params_STAR = params.STAR_Module_Map_STAR.params_STAR
nameAll = reads.toString()
nameArray = nameAll.split(' ')

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file =  nameAll 
    runGzip = ''
}

"""
$runGzip
STAR ${params_STAR}  --genomeDir ${params.star_index} --readFilesIn $file --outFileNamePrefix ${newName}
if [ ! -e "${newName}Aligned.toTranscriptome.out.bam" -a -e "${newName}Aligned.toTranscriptome.out.sam" ] ; then
    samtools view -S -b ${newName}Aligned.toTranscriptome.out.sam > ${newName}Aligned.toTranscriptome.out.bam
elif [ ! -e "${newName}Aligned.out.bam" -a -e "${newName}Aligned.out.sam" ] ; then
    samtools view -S -b ${newName}Aligned.out.sam > ${newName}Aligned.out.bam
fi
rm -rf *.sam
if [ -e "${newName}Aligned.sortedByCoord.out.bam" ] ; then
    mv ${newName}Aligned.sortedByCoord.out.bam ${newName}.bam
elif [ -e "${newName}Aligned.out.bam" ] ; then
    mv ${newName}Aligned.out.bam ${newName}.bam
fi

samtools flagstat ${newName}.bam > ${newName}.flagstat.txt
"""


}


//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2000
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill

process STAR_Module_Merge_Bam {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_sorted.*bai$/) "sorted_bam_star/$filename"
	else if (filename =~ /.*_sorted.*bam$/) "sorted_bam_star/$filename"
}

input:
 set val(oldname), file(bamfiles) from g115_16_mapped_reads_g115_14.groupTuple()

output:
 set val(oldname), file("${oldname}.bam")  into g115_14_merged_bams
 set val(oldname), file("*_sorted*bai")  into g115_14_bam_index
 set val(oldname), file("*_sorted*bam")  into g115_14_sorted_bam_g124_92

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
    samtools merge !{oldname}.bam !{bamfiles.join(" ")} && samtools sort -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
else
    mv !{bamfiles.join(" ")} !{oldname}.bam 2>/dev/null || true
    samtools sort  -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
fi
'''
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill

process Single_Cell_Module_after_STAR_samtools_sort_index {

input:
 set val(name), file(bam) from g115_14_sorted_bam_g124_92

output:
 set val(name), file("bam/*.bam")  into g124_92_bam_file_g124_90
 set val(name), file("bam/*.bai")  into g124_92_bam_index_g124_90

script:
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = "samtools index " + bam 
} else {
    runSamtools = "samtools sort -o " + name +"_sorted " + bam +" && samtools index " + name + "_sorted.bam "
}
"""
$runSamtools
mkdir -p bam && mv *_sorted.ba* bam/.
"""
}


//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2000
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill

process STAR_Module_merge_transcriptome_bam {

input:
 set val(oldname), file(bamfiles) from g115_16_transcriptome_bam_g115_15.groupTuple()

output:
 set val(oldname), file("${oldname}.bam")  into g115_15_merged_bams
 set val(oldname), file("*_sorted*bai")  into g115_15_bam_index
 set val(oldname), file("*_sorted*bam")  into g115_15_sorted_bam

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
    samtools merge !{oldname}.bam !{bamfiles.join(" ")} && samtools sort -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
else
    mv !{bamfiles.join(" ")} !{oldname}.bam 2>/dev/null || true
    samtools sort  -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
fi
'''
}

params.gtf =  ""  //* @input
params.genome =  ""  //* @input
params.commondb =  ""  //* @input
if (!(params.run_checkAndBuild == "yes" && params.run_Sequential_Mapping  == "yes")){
g114_15_commondb_path_g114_18.into{g114_18_commondb_path}
} else {


process Check_and_Build_Module_Check_Sequential_Mapping_Indexes {

input:
 val commondb from g114_15_commondb_path_g114_18
 val bowtieIndex from g114_13_genomeIndexPath_g114_18
 val bowtie2Index from g114_6_genomeIndexPath_g114_18
 val starIndex from g114_0_genomeIndexPath_g114_18

output:
 val commondb  into g114_18_commondb_path

when:
params.run_checkAndBuild == "yes" && params.run_Sequential_Mapping  == "yes"

script:
"""
"""
}
}



process STAR_Module_STAR_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.(out|tab)$/) "star/$filename"
}

input:
 set val(name), file(alignSum) from g115_16_outputFileOut_g115_18.groupTuple()
 set val(name), file(LogOut) from g115_16_logOut_g115_18.groupTuple()
 set val(name), file(progressOut) from g115_16_progressOut_g115_18.groupTuple()
 set val(name), file(TabOut) from g115_16_outputFileTab_g115_18.groupTuple()

output:
 file "*.tsv"  into g115_18_outputFile_g115_11
 file "*.{out,tab}"  into g115_18_logOut
 val "star_alignment_sum"  into g115_18_name_g115_11

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my %tsv;
my @headers = ();
my $name = "!{name}";

# merge output files 
`cat !{alignSum} >${name}_Merged_Log.final.out`;
`cat !{LogOut} >${name}_Merged_Log.out`;
`cat !{progressOut} >${name}_Merged_Log.progress.out`;
`cat !{TabOut} >${name}_Merged_SJ.out.tab`;

alteredAligned();

my @keys = keys %tsv;
my $summary = "$name"."_star_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
	my $values = join("\\t", @{ $tsv{$key} });
	`echo "$values" >> $summary`;
}


sub alteredAligned
{
	my @files = qw(!{alignSum});
	my $multimappedSum;
	my $alignedSum;
	my $inputCountSum;
	push(@headers, "Sample");
    push(@headers, "Total Reads");
	push(@headers, "Multimapped Reads Aligned (STAR)");
	push(@headers, "Unique Reads Aligned (STAR)");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($inputCount = `cat $file | grep 'Number of input reads' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($aligned = `cat $file | grep 'Uniquely mapped reads number' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($multimapped = `cat $file | grep 'Number of reads mapped to multiple loci' | awk '{sum+=\\$9} END {print sum}'`);
		$multimappedSum += int($multimapped);
        $alignedSum += int($aligned);
        $inputCountSum += int($inputCount);
	}
	$tsv{$name} = [$name, $inputCountSum];
	push(@{$tsv{$name}}, $multimappedSum);
	push(@{$tsv{$name}}, $alignedSum);
}
'''

}


process STAR_Module_merge_tsv_files_with_same_header {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.tsv$/) "star_alignment_summary/$filename"
}

input:
 file tsv from g115_18_outputFile_g115_11.collect()
 val outputFileName from g115_18_name_g115_11.collect()

output:
 file "${name}.tsv"  into g115_11_outputFileTSV

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}


process HISAT2_Module_HISAT2_Summary {

input:
 set val(name), file(alignSum) from g116_11_outputFileTxt_g116_2.groupTuple()

output:
 file "*.tsv"  into g116_2_outputFile_g116_10
 val "hisat2_alignment_sum"  into g116_2_name_g116_10

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my %tsv;
my @headers = ();
my $name = "!{name}";


alteredAligned();

my @keys = keys %tsv;
my $summary = "$name"."_hisat_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
	my $values = join("\\t", @{ $tsv{$key} });
	`echo "$values" >> $summary`;
}


sub alteredAligned
{
	my @files = qw(!{alignSum});
	my $multimappedSum;
	my $alignedSum;
	my $inputCountSum;
	push(@headers, "Sample");
    push(@headers, "Total Reads");
	push(@headers, "Multimapped Reads Aligned (HISAT2)");
	push(@headers, "Unique Reads Aligned (HISAT2)");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($inputCount = `cat $file | grep 'reads; of these:' | awk '{sum+=\\$1} END {print sum}'`);
		chomp($aligned = `cat $file | grep 'aligned.*exactly 1 time' | awk '{sum+=\\$1} END {print sum}'`);
		chomp($multimapped = `cat $file | grep 'aligned.*>1 times' | awk '{sum+=\\$1} END {print sum}'`);
		$multimappedSum += int($multimapped);
        $alignedSum += int($aligned);
        $inputCountSum += int($inputCount);
	}
	$tsv{$name} = [$name, $inputCountSum];
	push(@{$tsv{$name}}, $multimappedSum);
	push(@{$tsv{$name}}, $alignedSum);
}
'''

}


process HISAT2_Module_Merge_TSV_Files {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.tsv$/) "hisat2_alignment_summary/$filename"
}

input:
 file tsv from g116_2_outputFile_g116_10.collect()
 val outputFileName from g116_2_name_g116_10.collect()

output:
 file "${name}.tsv"  into g116_10_outputFileTSV

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill

process Tophat2_Module_Merge_Tophat_Summary {

input:
 set val(name), file(alignSum) from g117_11_summary_g117_3.groupTuple()
 val mate from g_13_mate_g117_3

output:
 set val(name), file("${name}_tophat_sum.tsv")  into g117_3_report_g117_9
 val "tophat2_alignment_sum"  into g117_3_name_g117_9

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my %tsv;
my @headers = ();
my $name = "!{name}";

alteredAligned();

my @keys = keys %tsv;
my $summary = "$name"."_tophat_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
	my $values = join("\\t", @{ $tsv{$key} });
	`echo "$values" >> $summary`;
}


sub alteredAligned
{
	my @files = qw(!{alignSum});
	my $multimappedSum;
	my $alignedSum;
	my $inputCountSum;
	push(@headers, "Sample");
    push(@headers, "Total Reads");
	push(@headers, "Multimapped Reads Aligned (Tophat2)");
	push(@headers, "Unique Reads Aligned (Tophat2)");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($aligned = `cat $file | grep 'Aligned pairs:' | awk '{sum=\\$3} END {print sum}'`);
		if ($aligned eq "") { # then it is single-end
		        chomp($inputCount = `cat $file | grep 'Input' | awk '{sum=\\$3} END {print sum}'`);
				chomp($aligned = `cat $file | grep 'Mapped' | awk '{sum=\\$3} END {print sum}'`);
				chomp($multimapped = `cat $file | grep 'multiple alignments' | awk '{sum+=\\$3} END {print sum}'`);
			}else{ # continue to pair end
			    chomp($inputCount = `cat $file | grep 'Input' | awk '{sum=\\$3} END {print sum}'`);
				chomp($multimapped = `cat $file | grep -A 1 'Aligned pairs:' | awk 'NR % 3 == 2 {sum+=\\$3} END {print sum}'`);
			}
        $multimappedSum += int($multimapped);
        $alignedSum += (int($aligned) - int($multimapped));
        $inputCountSum += int($inputCount);
        if ($alignedSum < 0){
            $alignedSum = 0;
        }
	}
	$tsv{$name} = [$name, $inputCountSum];
	push(@{$tsv{$name}}, $multimappedSum);
	push(@{$tsv{$name}}, $alignedSum);
}
'''

}


process Tophat2_Module_Merge_TSV_Files {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.tsv$/) "tophat_alignment_summary/$filename"
}

input:
 file tsv from g117_3_report_g117_9.collect()
 val outputFileName from g117_3_name_g117_9.collect()

output:
 file "${name}.tsv"  into g117_9_outputFileTSV

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}

params.countUniqueAlignedBarcodes_fromFile_filePath =  ""  //* @input

//* autofill
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 1
    $QUEUE = "long"
}
//*
if (!((params.run_Single_Cell_Module && (params.run_Single_Cell_Module == "yes")) || !params.run_Single_Cell_Module)){
g122_92_bam_file_g122_90.into{g122_90_sorted_bam_g122_91}
g122_92_bam_index_g122_90.into{g122_90_bam_index_g122_91}
g122_90_count_file_g122_91 = Channel.empty()
} else {


process Single_Cell_Module_after_Tophat2_Count_cells {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_count.txt$/) "cell_counts_after_tophat2/$filename"
}

input:
 set val(oldname), file(sorted_bams) from g122_92_bam_file_g122_90
 set val(oldname), file(bams_index) from g122_92_bam_index_g122_90
 val cutoff_reads_per_cell from g_111_cutoff_g122_90

output:
 set val(oldname), file("bam/*.bam")  into g122_90_sorted_bam_g122_91
 set val(oldname), file("bam/*.bam.bai")  into g122_90_bam_index_g122_91
 set val(oldname), file("*_count.txt")  into g122_90_count_file_g122_91

when:
(params.run_Single_Cell_Module && (params.run_Single_Cell_Module == "yes")) || !params.run_Single_Cell_Module

script:
"""
find  -name "*.bam" > filelist.txt
python ${params.countUniqueAlignedBarcodes_fromFile_filePath} -i filelist.txt -m ${cutoff_reads_per_cell} -o ${oldname}_count.txt
mkdir bam
mv $sorted_bams bam/.
mv $bams_index bam/.
"""
}
}


params.filter_lowCountBC_bam_print_py_filePath =  ""  //* @input
maxCellsForTmpFile = params.Single_Cell_Module_after_Tophat2_filter_lowCount.maxCellsForTmpFile

//* autofill
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 1
    $QUEUE = "long"
}
//*

process Single_Cell_Module_after_Tophat2_filter_lowCount {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}_filtered_.*.bam$/) "filtered_bam_after_tophat2/$filename"
}

input:
 set val(oldname), file(sorted_bams) from g122_90_sorted_bam_g122_91
 set val(name), file(count_file) from g122_90_count_file_g122_91
 set val(oldname), file(bam_index) from g122_90_bam_index_g122_91
 val cutoff_for_filter from g_111_cutoff_g122_91

output:
 set val(name), file("${name}_filtered_*.bam")  into g122_91_filtered_bam_g122_87

"""
python ${params.filter_lowCountBC_bam_print_py_filePath} -i ${sorted_bams} -b ${name}_count.txt -o ${name}_filtered.bam -n ${cutoff_for_filter} -c ${maxCellsForTmpFile}
"""
}

esat_parameters = params.Single_Cell_Module_after_Tophat2_ESAT.esat_parameters
params.ESAT_path =  ""  //* @input
params.gene_to_transcript_mapping_file =  ""  //* @input


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 40
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 40
    $QUEUE = "long"
}
//* platform
//* autofill

process Single_Cell_Module_after_Tophat2_ESAT {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.(txt|log)$/) "esat_after_tophat2/$filename"
	else if (filename =~ /.*umi.distributions.txt$/) "esat_after_tophat2/$filename"
}

input:
 set val(name), file(filtered_bam) from g122_91_filtered_bam_g122_87.transpose()

output:
 file "*.{txt,log}"  into g122_87_outputFileTxt
 set val(name), file("*umi.distributions.txt")  into g122_87_UMI_distributions_g122_88

script:
nameAll = filtered_bam.toString()
namePrefix = nameAll - ".bam"
"""    
find  -name "*.bam" | awk '{print "${namePrefix}\t"\$1 }' > ${namePrefix}_filelist.txt
java -Xmx40g -jar ${params.ESAT_path} -alignments ${namePrefix}_filelist.txt -out ${namePrefix}_esat.txt -geneMapping ${params.gene_to_transcript_mapping_file} ${esat_parameters}
mv scripture2.log ${namePrefix}_scripture2.log
"""
}

params.cleanLowEndUmis_path =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 30
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 700
    $CPU  = 1
    $MEMORY = 30
    $QUEUE = "long"
}
//* platform
//* autofill

process Single_Cell_Module_after_Tophat2_UMI_Trim {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_umiClean.txt$/) "UMI_count_final_after_tophat2/$filename"
}

input:
 set val(name), file(umi_dist) from g122_87_UMI_distributions_g122_88.groupTuple()

output:
 set val(name), file("*_umiClean.txt")  into g122_88_UMI_clean

"""
cat ${umi_dist} > ${name}_merged_umi.distributions.txt
python ${params.cleanLowEndUmis_path} \
-i ${name}_merged_umi.distributions.txt \
-o ${name}_umiClean.txt \
-n 2
"""
}

params.countUniqueAlignedBarcodes_fromFile_filePath =  ""  //* @input

//* autofill
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 1
    $QUEUE = "long"
}
//*
if (!((params.run_Single_Cell_Module && (params.run_Single_Cell_Module == "yes")) || !params.run_Single_Cell_Module)){
g123_92_bam_file_g123_90.into{g123_90_sorted_bam_g123_91}
g123_92_bam_index_g123_90.into{g123_90_bam_index_g123_91}
g123_90_count_file_g123_91 = Channel.empty()
} else {


process Single_Cell_Module_after_HISAT2_Count_cells {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_count.txt$/) "cell_counts_after_hisat2/$filename"
}

input:
 set val(oldname), file(sorted_bams) from g123_92_bam_file_g123_90
 set val(oldname), file(bams_index) from g123_92_bam_index_g123_90
 val cutoff_reads_per_cell from g_111_cutoff_g123_90

output:
 set val(oldname), file("bam/*.bam")  into g123_90_sorted_bam_g123_91
 set val(oldname), file("bam/*.bam.bai")  into g123_90_bam_index_g123_91
 set val(oldname), file("*_count.txt")  into g123_90_count_file_g123_91

when:
(params.run_Single_Cell_Module && (params.run_Single_Cell_Module == "yes")) || !params.run_Single_Cell_Module

script:
"""
find  -name "*.bam" > filelist.txt
python ${params.countUniqueAlignedBarcodes_fromFile_filePath} -i filelist.txt -m ${cutoff_reads_per_cell} -o ${oldname}_count.txt
mkdir bam
mv $sorted_bams bam/.
mv $bams_index bam/.
"""
}
}


params.filter_lowCountBC_bam_print_py_filePath =  ""  //* @input
maxCellsForTmpFile = params.Single_Cell_Module_after_HISAT2_filter_lowCount.maxCellsForTmpFile

//* autofill
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 1
    $QUEUE = "long"
}
//*

process Single_Cell_Module_after_HISAT2_filter_lowCount {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}_filtered_.*.bam$/) "filtered_bam_after_hisat2/$filename"
}

input:
 set val(oldname), file(sorted_bams) from g123_90_sorted_bam_g123_91
 set val(name), file(count_file) from g123_90_count_file_g123_91
 set val(oldname), file(bam_index) from g123_90_bam_index_g123_91
 val cutoff_for_filter from g_111_cutoff_g123_91

output:
 set val(name), file("${name}_filtered_*.bam")  into g123_91_filtered_bam_g123_87

"""
python ${params.filter_lowCountBC_bam_print_py_filePath} -i ${sorted_bams} -b ${name}_count.txt -o ${name}_filtered.bam -n ${cutoff_for_filter} -c ${maxCellsForTmpFile}
"""
}

esat_parameters = params.Single_Cell_Module_after_HISAT2_ESAT.esat_parameters
params.ESAT_path =  ""  //* @input
params.gene_to_transcript_mapping_file =  ""  //* @input


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 40
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 40
    $QUEUE = "long"
}
//* platform
//* autofill

process Single_Cell_Module_after_HISAT2_ESAT {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.(txt|log)$/) "esat_after_hisat2/$filename"
	else if (filename =~ /.*umi.distributions.txt$/) "esat_after_hisat2/$filename"
}

input:
 set val(name), file(filtered_bam) from g123_91_filtered_bam_g123_87.transpose()

output:
 file "*.{txt,log}"  into g123_87_outputFileTxt
 set val(name), file("*umi.distributions.txt")  into g123_87_UMI_distributions_g123_88

script:
nameAll = filtered_bam.toString()
namePrefix = nameAll - ".bam"
"""    
find  -name "*.bam" | awk '{print "${namePrefix}\t"\$1 }' > ${namePrefix}_filelist.txt
java -Xmx40g -jar ${params.ESAT_path} -alignments ${namePrefix}_filelist.txt -out ${namePrefix}_esat.txt -geneMapping ${params.gene_to_transcript_mapping_file} ${esat_parameters}
mv scripture2.log ${namePrefix}_scripture2.log
"""
}

params.cleanLowEndUmis_path =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 30
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 700
    $CPU  = 1
    $MEMORY = 30
    $QUEUE = "long"
}
//* platform
//* autofill

process Single_Cell_Module_after_HISAT2_UMI_Trim {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_umiClean.txt$/) "UMI_count_final_after_hisat2/$filename"
}

input:
 set val(name), file(umi_dist) from g123_87_UMI_distributions_g123_88.groupTuple()

output:
 set val(name), file("*_umiClean.txt")  into g123_88_UMI_clean

"""
cat ${umi_dist} > ${name}_merged_umi.distributions.txt
python ${params.cleanLowEndUmis_path} \
-i ${name}_merged_umi.distributions.txt \
-o ${name}_umiClean.txt \
-n 2
"""
}

params.countUniqueAlignedBarcodes_fromFile_filePath =  ""  //* @input

//* autofill
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 1
    $QUEUE = "long"
}
//*
if (!((params.run_Single_Cell_Module && (params.run_Single_Cell_Module == "yes")) || !params.run_Single_Cell_Module)){
g124_92_bam_file_g124_90.into{g124_90_sorted_bam_g124_91}
g124_92_bam_index_g124_90.into{g124_90_bam_index_g124_91}
g124_90_count_file_g124_91 = Channel.empty()
} else {


process Single_Cell_Module_after_STAR_Count_cells {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_count.txt$/) "cell_counts_after_star/$filename"
}

input:
 set val(oldname), file(sorted_bams) from g124_92_bam_file_g124_90
 set val(oldname), file(bams_index) from g124_92_bam_index_g124_90
 val cutoff_reads_per_cell from g_111_cutoff_g124_90

output:
 set val(oldname), file("bam/*.bam")  into g124_90_sorted_bam_g124_91
 set val(oldname), file("bam/*.bam.bai")  into g124_90_bam_index_g124_91
 set val(oldname), file("*_count.txt")  into g124_90_count_file_g124_91

when:
(params.run_Single_Cell_Module && (params.run_Single_Cell_Module == "yes")) || !params.run_Single_Cell_Module

script:
"""
find  -name "*.bam" > filelist.txt
python ${params.countUniqueAlignedBarcodes_fromFile_filePath} -i filelist.txt -m ${cutoff_reads_per_cell} -o ${oldname}_count.txt
mkdir bam
mv $sorted_bams bam/.
mv $bams_index bam/.
"""
}
}


params.filter_lowCountBC_bam_print_py_filePath =  ""  //* @input
maxCellsForTmpFile = params.Single_Cell_Module_after_STAR_filter_lowCount.maxCellsForTmpFile

//* autofill
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 1
    $QUEUE = "long"
}
//*

process Single_Cell_Module_after_STAR_filter_lowCount {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}_filtered_.*.bam$/) "filtered_bam_after_star/$filename"
}

input:
 set val(oldname), file(sorted_bams) from g124_90_sorted_bam_g124_91
 set val(name), file(count_file) from g124_90_count_file_g124_91
 set val(oldname), file(bam_index) from g124_90_bam_index_g124_91
 val cutoff_for_filter from g_111_cutoff_g124_91

output:
 set val(name), file("${name}_filtered_*.bam")  into g124_91_filtered_bam_g124_87

"""
python ${params.filter_lowCountBC_bam_print_py_filePath} -i ${sorted_bams} -b ${name}_count.txt -o ${name}_filtered.bam -n ${cutoff_for_filter} -c ${maxCellsForTmpFile}
"""
}

esat_parameters = params.Single_Cell_Module_after_STAR_ESAT.esat_parameters
params.ESAT_path =  ""  //* @input
params.gene_to_transcript_mapping_file =  ""  //* @input


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 40
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 40
    $QUEUE = "long"
}
//* platform
//* autofill

process Single_Cell_Module_after_STAR_ESAT {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.(txt|log)$/) "esat_after_star/$filename"
	else if (filename =~ /.*umi.distributions.txt$/) "esat_after_star/$filename"
}

input:
 set val(name), file(filtered_bam) from g124_91_filtered_bam_g124_87.transpose()

output:
 file "*.{txt,log}"  into g124_87_outputFileTxt
 set val(name), file("*umi.distributions.txt")  into g124_87_UMI_distributions_g124_88

script:
nameAll = filtered_bam.toString()
namePrefix = nameAll - ".bam"
"""    
find  -name "*.bam" | awk '{print "${namePrefix}\t"\$1 }' > ${namePrefix}_filelist.txt
java -Xmx40g -jar ${params.ESAT_path} -alignments ${namePrefix}_filelist.txt -out ${namePrefix}_esat.txt -geneMapping ${params.gene_to_transcript_mapping_file} ${esat_parameters}
mv scripture2.log ${namePrefix}_scripture2.log
"""
}

params.cleanLowEndUmis_path =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 30
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 700
    $CPU  = 1
    $MEMORY = 30
    $QUEUE = "long"
}
//* platform
//* autofill

process Single_Cell_Module_after_STAR_UMI_Trim {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_umiClean.txt$/) "UMI_count_final_after_star/$filename"
}

input:
 set val(name), file(umi_dist) from g124_87_UMI_distributions_g124_88.groupTuple()

output:
 set val(name), file("*_umiClean.txt")  into g124_88_UMI_clean

"""
cat ${umi_dist} > ${name}_merged_umi.distributions.txt
python ${params.cleanLowEndUmis_path} \
-i ${name}_merged_umi.distributions.txt \
-o ${name}_umiClean.txt \
-n 2
"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
