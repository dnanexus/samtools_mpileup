#!/usr/bin/env python

import dxpy
import subprocess, logging
import os
import sys
import re
import math
import operator

@dxpy.entry_point('main')
def main(**job_inputs):
    job_outputs = {}
    mappingsTable = dxpy.open_dxgtable(job_inputs['mappings']['$dnanexus_link'])
    mappingsTableId = mappingsTable.get_id()
    

    #This controls the degree of parallelism
    chunks = int(mappingsTable.describe()['length']/job_inputs['reads_per_job'])+1

    try:
        contigSetId = mappingsTable.get_details()['original_contigset']['$dnanexus_link']
        originalContigSet = mappingsTable.get_details()['original_contigset']
    except:
        raise Exception("The original reference genome must be attached as a detail")

    variants_schema = [
        {"name": "chr", "type": "string"}, 
        {"name": "lo", "type": "int32"},
        {"name": "hi", "type": "int32"},
        {"name": "ref", "type": "string"},
        {"name": "alt", "type": "string"},
        {"name": "qual", "type": "double"},
        {"name": "filter", "type": "string"},
        {"name": "ids", "type": "string"}
         ]

    #The information in these tags is elevated into specific columns, so additional columns for these tags will not be created
    elevatedTags = ['format_GT', 'format_DP', 'format_AD']
    
    #The info and format tags are extracted from the header printed by samtools
    #If additional code will add a tag to the output of the program, modify this header to include the tag.
    #TODO: Allow the table to be created by the first job that finishes to avoid this step.
    headerInfo = extractHeader("/tmp/header.txt", elevatedTags)
    description = {}
    samples = []


    indices = [dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')]
    formats = {}
    infos = {}
    filters = {}
    
    ##The following section creates the sample-specific table columns
    for k, v in headerInfo['tags']['info'].iteritems():
        variants_schema.append({"name": "info_"+k, "type":translateTagTypeToColumnType(v)})
        description[k] = {'name' : k, 'description' : v['description'], 'type' : v['type'], 'number' : v['number']}
    
    numSamples = 1
    #For each sample, write the sample-specific columns, at present only one sample is supported
    for i in range(numSamples):
      variants_schema.extend([
        {"name": "genotype_"+str(i), "type": "string"},
        {"name": "phasing_"+str(i), "type": "string"},
        {"name": "type_"+str(i), "type": "string"},
        {"name": "variation_qual_"+str(i), "type": "double"},
        {"name": "genotype_qual_"+str(i), "type": "double"},
        {"name": "coverage_"+str(i), "type": "string"},
        {"name": "total_coverage_"+str(i), "type": "int32"}
      ])
      indices.append(dxpy.DXGTable.lexicographic_index([["type_"+str(i), "ASC"]], 'type_'+str(i)))
      samples.append("Sample_0")
      for k, v in headerInfo['tags']['format'].iteritems():
        if "format_"+k not in elevatedTags:
          variants_schema.append({"name": "format_"+k+"_"+str(i), "type":translateTagTypeToColumnType(v)})

    #TODO: Add lexicographic indices when secondary indices are supported
    variants = dxpy.new_dxgtable(variants_schema, indices=[dxpy.DXGTable.genomic_range_index("chr", "lo", "hi", "gri")])
    tableId = variants.get_id()
    variants = dxpy.open_dxgtable(tableId)
    variants.add_types(["Variants", "gri"])
    
    details = {'samples':samples, 'original_contigset':job_inputs['reference'], 'original_mappings':job_inputs['mappings'], 'formats':headerInfo['tags']['format'], 'infos':headerInfo['tags']['info']}
    if headerInfo.get('filters') != {}:
      details['filters'] = headerInfo['filters']
    variants.set_details(details)

    if 'output name' in job_inputs:
        variants.rename(job_inputs['output name'])
    else:
        variants.rename(mappingsTable.describe()['name'] + " variant calls by Samtools mpileup")


    #Split the genome into evenly sized regions
    genomeRegions = splitGenomeLengthLargePieces(originalContigSet, chunks)

    #Generate the command line arguments needed to run samtools and bcftools
    samOptions = makeSamtoolsParameters(**job_inputs)
    bcfOptions = makeBcftoolsParameters(**job_inputs)


    #This 
    reduce_job_inputs = {}
    for i in range(len(commandList)):
        #print commandList[i]
        if len(commandList[i]) > 0:
            map_job_inputs = {
                'mappings_table_id':mappingsTableId,
                'original_contig_set': contigSetId,
                'interval': genomeRegions[i],
                'tableId': tableId,
                'compress_reference': job_inputs['compress_reference'],
                'compress_no_call': job_inputs['compress_no_call'],
                'infer_no_call' : job_inputs['infer_no_call'],
                'sam_options': samOptions,
                'bcf_options': bcfOptions,
                'part_number' : i
            }
            #Run a "map" job for each chunk, passing in the inputspec from above and looking for a function entry point given as "map" (@dxpy.entry_point('map'))
            map_job = dxpy.new_dxjob(map_job_inputs, "map")
            reduce_job_inputs["mapJob" + str(i) + "TableId"] = {'job': map_job.get_id(), 'field': 'ok'}

    reduce_job_inputs['tableId'] = tableId
    
    #Run a "reduce" job, which only begins once all of the map jobs singal they have completed by sending 'ok':True
    #The reduce job closes the table. This step is explicitly needed because table closing must wait till the completion of the map jobs
    #By giving the reduce job the map jobs as input, the reduce job will wait to start.
    reduce_job = dxpy.new_dxjob(reduce_job_inputs, "reduce")
    job_outputs = {'variants': {'job': reduce_job.get_id(), 'field': 'variants'}}
    
    return job_outputs

@dxpy.entry_point('map')
def mapPileup(**job_inputs):
    print "Downloading Reference Genome"
    subprocess.check_call("contigset2fasta %s ref.fa" % (job_inputs['original_contig_set']), shell=True)

    #The mappings-to-sam script takes a file with a list of regions in the form -L chrX:lo-hi, as produced by the genome splitting region.
    #This generates the file
    regionFile = open("regions.txt", 'w')
    regionFile.write(job_inputs['interval'])
    regionFile.close()

    print "Indexing Dictionary"
    subprocess.check_call("samtools faidx ref.fa", shell=True)

    #The sam-to-mappings script in dx-toolkit, the region_index_offset option 
    print "Converting Table to SAM"
    command = "dx-mappings-to-sam %s --output input.sam --region_index_offset -1 --region_file regions.txt" % (job_inputs['mappings_table_id'])
    print "Running: " + command
    subprocess.check_call(command, shell=True)
    if checkSamContainsRead("input.sam"):
        print "Converting to BAM"
        subprocess.check_call("samtools view -bS input.sam > input.bam", shell=True)
        print "Indexing"
        subprocess.check_call("samtools index input.bam", shell=True)

        variants = dxpy.open_dxgtable(job_inputs['tableId'])
        command = "samtools mpileup -uf ref.fa"
        #command += job_inputs['interval']
        bedFile = open("regions.bed", 'w')
        intervalMatch = re.findall("-L ([^:]*):(\d+)-(\d+)", job_inputs['interval'])
        if len(intervalMatch) > 0:
            for x in intervalMatch:
                bedFile.write(x[0]+"\t"+str(x[1])+"\t"+str(x[2])+"\n")
            bedFile.close()
            bedFile = open("regions.bed", 'r')
            bedFile.close()
            command += " -l regions.bed "
        command += job_inputs['sam_options']
        command += " input.bam | bcftools view "
        command += job_inputs['bcf_options']
        command += " - > output.vcf"
        
        print "Pileup Command: " + command
        subprocess.check_call(command ,shell=True)

        command = "dx_vcfToVariants2 --table_id %s --vcf_file output.vcf --region_file regions.txt" % (job_inputs['tableId'])
        if job_inputs['compress_reference']:
            command += " --compress_reference"
        if job_inputs['infer_no_call']:
            command += " --infer_no_call"
        if job_inputs['compress_no_call']:
            command += " --compress_no_call"

        print "Import variants command: " + command
        subprocess.check_call(command ,shell=True)
    job_outputs = {'ok':True}
    return job_outputs
        
@dxpy.entry_point('reduce')
def reducePileup(**job_inputs):
    t = dxpy.open_dxgtable(job_inputs['tableId'])
    print "Closing Table"
    t.close()
    job_outputs = {'variants': dxpy.dxlink(t.get_id())}
    return job_outputs
    
def makeSamtoolsParameters(**job_inputs):
    #Always output genotypes
    options = '-u'
    if job_inputs['disable_probabilistic_realignment'] == True:
        options += 'B'
    if job_inputs['extended_baq_computation'] == True:
        options += 'E'

    options += ' -C ' + str(job_inputs['excessive_mismatch_penalty'])
    options += ' -d ' + str(job_inputs['max_reads_per_position'])
    options += ' -q ' + str(job_inputs['minimum_mapping_quality'])
    options += ' -Q ' + str(job_inputs['minimum_base_quality'])

    return options

def makeBcftoolsParameters(**job_inputs):
    #Always output genotypes if possible
    options = '-g'

    if job_inputs['output_variants_only']:
        job_inputs['bayesian_inference'] = True
        options += 'v'
    if job_inputs['bayesian_inference']:
        options += 'c'
        job_inputs['perform_max_likelihood_inference'] = True
    if job_inputs['perform_max_likelihood_inference']:
        options += 'e'
    if job_inputs['retain_unlikely_alleles']:
        options += 'A'

    options += ' -i ' + str(job_inputs['indel_to_snp_ratio'])
    options += ' -p ' + str(job_inputs['variant_probability'])
    options += ' -t ' + str(job_inputs['scaled_mutation_rate'])

    return options

def extractHeader(vcfFileName, elevatedTags):
    result = {'columns': '', 'tags' : {'format' : {}, 'info' : {} }, 'filters' : {}}
    for line in open(vcfFileName):
        tag = re.findall("ID=(\w+),", line)
        if len(tag) > 0:
          tagType = ''
          if line.count("FORMAT") > 0:
            tagType = 'format'
          elif line.count("INFO") > 0:
            tagType = 'info'
          elif line.count("FILTER") > 0:
            result['filters'][re.findall("ID=(\w+),")[0]] = re.findall('Description="(.*)"')[0]
      
          typ = re.findall("Type=(\w+),", line)
          if tagType != '':
            number = re.findall("Number=(\w+)", line)
            description = re.findall('Description="(.*)"', line)
            if len(number) == 0:
              number = ['.']
            if len(description) == 0:
              description = ['']
            if "format_"+tag[0] not in elevatedTags:
                result['tags'][tagType][tag[0]] = {'type':typ[0], 'description' : description[0], 'number' : number[0]}
        if line[0] == "#" and line[1] != "#":
          result['columns'] = line.strip()
        if line == '' or line[0] != "#":
            break
    return result

def splitGenomeLengthLargePieces(contig_set, chunks):
    details = dxpy.DXRecord(contig_set).get_details()
    sizes = details['contigs']['sizes']
    names = details['contigs']['names']
    offsets = details['contigs']['offsets']
    
    for i in range(len(names)):
        print names[i]+":"+str(sizes[i])


    commandList = []
    for i in range(chunks):
        commandList.append('')
    position = 0
    chromosome = 0
    chunkSize = sum(sizes) / chunks
    currentChunk = 0
    currentLength = 0

    while chromosome < len(names):
        if position + (chunkSize - currentLength) >= sizes[chromosome]:
            commandList[currentChunk] += " -L %s:%d-%d" % (names[chromosome], position+1, sizes[chromosome])
            currentLength += sizes[chromosome] - position
            chromosome += 1
            position = 0
        else:
            commandList[currentChunk] += " -L %s:%d-%d" % (names[chromosome], position+1, position+(chunkSize-currentLength)+1)
            position += (chunkSize-currentLength) + 1
            if currentChunk < chunks-1:
                currentChunk += 1
            currentLength = 0
    return commandList

def checkSamContainsRead(samFileName):
    for line in open(samFileName, 'r'):
        if line[0] != "@":
            return True
    return False

def translateTagTypeToColumnType(tag):
  if tag['type'] == "Flag":
    return "boolean"
  if tag['number'] != '1':
    return 'string'
  if tag['type'] == "Integer":
    return 'int32'
  if tag['type'] == "Float":
    return "double"
  return "string"

dxpy.run()

