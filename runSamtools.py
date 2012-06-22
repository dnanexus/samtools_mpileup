#!/usr/bin/env python

import dxpy
import subprocess, logging
import os
import sys
import re
import math
import operator

def main():
    mappingsTable = dxpy.open_dxgtable(job['input']['mappings']['$dnanexus_link'])
    mappingsTableId = mappingsTable.get_id()
    try:
        contigSetId = mappingsTable.get_details()['original_contigset']['$dnanexus_link']
        originalContigSet = mappingsTable.get_details()['original_contigset']
    except:
        raise Exception("The original reference genome must be attached as a detail")

    #subprocess.check_call("contigset2fasta %s ref.fa" % (contigSetId), shell=True)
    #reference_sequence = dxpy.dxlink(dxpy.upload_local_file("ref.fa"))

    #print "Indexing Dictionary"
    #subprocess.check_call("samtools faidx ref.fa", shell=True)

    #referenceIndex = dxpy.dxlink(dxpy.upload_local_file("ref.fa.fai"))

    variants_schema = [{"name": "chr", "type": "string"},
                       {"name": "lo", "type": "int32"},
                       {"name": "hi", "type": "int32"},
                       {"name": "type", "type": "string"},     #change this type to uint once there is an abstraction method for enum
                       {"name": "ref", "type": "string"},
                       {"name": "alt", "type": "string"},
                       {"name": "qual", "type": "int32"},
                       {"name": "coverage", "type": "int32"},
                       {"name": "genotypeQuality", "type": "int32"}]
    if job['input']['store_full_vcf']:
        variants_schema.extend([{"name": "vcf_alt", "type": "string"}, {"name": "vcf_additional_data", "type": "string"}])

    simpleVar = dxpy.new_dxgtable(variants_schema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')])
    tableId = simpleVar.get_id()
    simpleVar = dxpy.open_dxgtable(tableId)
    simpleVar.set_details({'original_contigset':originalContigSet})
    simpleVar.add_types(["SimpleVar", "gri"])

    reduceInput = {}
    #commandList = splitGenomeLengthLargePieces(originalContigSet, job['input']['intervals_to_process'], job['input']['intervals_to_exclude'],  job['input']['minimum_chunk_size'], job['input']['maximum_chunks'])
    commandList = splitGenomeLengthLargePieces(originalContigSet, job['input']['maximum_chunks'])
    samOptions = makeSamtoolsParameters(job)
    bcfOptions = makeBcftoolsParameters(job)

    for i in range(len(commandList)):
        #print commandList[i]
        if len(commandList[i]) > 0:
            mapInput = {
                'mappings_table_id':mappingsTableId,
                'original_contig_set': contigSetId,
                'interval': commandList[i],
                'tableId': tableId,
                'compress_reference': job['input']['compress_reference'],
                'infer_no_call' : job['input']['infer_no_call'],
                'store_full_vcf' : job['input']['store_full_vcf'],
                'sam_options': samOptions,
                'bcf_options': bcfOptions,
                'part_number' : i
            }
            # Run a "map" job for each chunk
            mapJobId = dxpy.new_dxjob(fn_input=mapInput, fn_name="mapPileup").get_id()
            reduceInput["mapJob" + str(i) + "TableId"] = {'job': mapJobId, 'field': 'id'}

    reduceInput['tableId'] = tableId
    reduceJobId = dxpy.new_dxjob(fn_input=reduceInput, fn_name="reducePileup").get_id()
    #print "SimpleVar table" + json.dumps({'table_id':simpleVar.get_id()})
    job['output'] = {'simplevar': {'job': reduceJobId, 'field': 'simplevar'}}

def makeSamtoolsParameters(job):
    #Always output genotypes
    options = '-u'
    if job['input']['disable_probabilistic_realignment'] == True:
        options += 'B'
    if job['input']['extended_baq_computation'] == True:
        options += 'E'

    options += ' -C ' + str(job['input']['excessive_mismatch_penalty'])
    options += ' -d ' + str(job['input']['max_reads_per_position'])
    options += ' -q ' + str(job['input']['minimum_mapping_quality'])
    options += ' -Q ' + str(job['input']['minimum_base_quality'])

    return options

def makeBcftoolsParameters(job):
    #Always output genotypes if possible
    options = '-g'

    if job['input']['output_variants_only']:
        job['input']['bayesian_inference'] = True
        options += 'v'
    if job['input']['bayesian_inference']:
        options += 'c'
        job['input']['perform_max_likelihood_inference'] = True
    if job['input']['perform_max_likelihood_inference']:
        options += 'e'
    if job['input']['retain_unlikely_alleles']:
        options += 'A'

    options += ' -i ' + str(job['input']['indel_to_snp_ratio'])
    options += ' -p ' + str(job['input']['variant_probability'])
    options += ' -t ' + str(job['input']['scaled_mutation_rate'])

    return options

def mapPileup():
    print "Downloading Reference Genome"
    subprocess.check_call("contigset2fasta %s ref.fa" % (job['input']['original_contig_set']), shell=True)
    #subprocess.check_call("dx_writeReferenceIndex --contig_set %s --writeSamtoolsIndex ref.fa.fai" % (job['input']['original_contig_set']), shell=True)

    regionFile = open("regions.txt", 'w')
    regionFile.write(job['input']['interval'])
    regionFile.close()

    print "Indexing Dictionary"
    subprocess.check_call("samtools faidx ref.fa", shell=True)

    print "Converting Table to SAM"
    print "dx_mappingsTableToSam2 --table_id %s --output input.sam --region_index_offset -1 --region_file regions.txt" % (job['input']['mappings_table_id'])
    subprocess.check_call("dx_mappingsTableToSam2 --table_id %s --output input.sam --region_index_offset -1 --region_file regions.txt" % (job['input']['mappings_table_id']), shell=True)

    if checkSamContainsRead("input.sam"):
        print "Converting to BAM"
        subprocess.check_call("samtools view -bS input.sam > input.bam", shell=True)
        #print "Sorting"
        #subprocess.check_call("samtools sort input.bam input.sorted", shell=True)
        print "Indexing"
        subprocess.check_call("samtools index input.bam", shell=True)

        simpleVar = dxpy.open_dxgtable(job['input']['tableId'])

        #command = job['input']['command'] + job['input']['interval']
        #print command

        command = "samtools mpileup -uf ref.fa"
        #command += job['input']['interval']
        bedFile = open("regions.bed", 'w')
        intervalMatch = re.findall("(\w+):(\d+)-(\d+)", job['input']['interval'])
        if len(intervalMatch) > 0:
            for x in intervalMatch:
                bedFile.write(x[0]+"\t"+str(x[1])+"\t"+str(x[2])+"\n")
            bedFile.close()
            bedFile = open("regions.bed", 'r')
            bedFile.close()
            command += " -l regions.bed "
        command += job['input']['sam_options']
        command += " input.bam | bcftools view "
        command += job['input']['bcf_options']
        command += " - > output.vcf"
        subprocess.call(command, shell=True)

        command = "dx_vcfToSimplevar2 --table_id %s --vcf_file output.vcf --region_file regions.txt" % (job['input']['tableId'])
        if job['input']['compress_reference']:
            command += " --compress_reference"
        if job['input']['infer_no_call']:
            command += " --infer_no_call"
        if job['input']['store_full_vcf']:
            command += " --store_full_vcf"
        if job['input']['part_number'] == 0:
            command += " --extract_header"   

        subprocess.call(command ,shell=True)

def runTrivialTest(contig_set, command):
    details = dxpy.DXRecord(contig_set['$dnanexus_link']).get_details()
    sizes = details['contigs']['sizes']
    names = details['contigs']['names']
    offsets = details['contigs']['offsets']

    for i in range(len(names)):
        if sizes[i] > 0:
            chromosome = names[i]
            break
    command += ' -L ' + chromosome+':1-1'
    subprocess.call(command, shell=True)
    return extractHeader(open("output.vcf", 'r'))

def extractHeader(vcfFile):
    header = ''
    fileIter = vcfFile.__iter__()

    #Additional data will contain the extra format and info columns that are optional in VCF and may not be
    #   present in the VCF file. These are stored in an extended table
    additionalColumns = []
    while 1:
        try:
            input = fileIter.next()
            if input[0] == "#":
                header += input
                #extract additional column header data
                if(input[1] != "#"):
                    tabSplit = input.split("\t")
                    return header
        except StopIteration:
            break

def reducePileup():
    t = dxpy.open_dxgtable(job['input']['tableId'])
    print "Closing Table"
    t.close()
    job['output']['simplevar'] = dxpy.dxlink(t.get_id())

def splitGenomeLength(contig_set, includeInterval, excludeInterval, chunkSize, splits):
    details = dxpy.DXRecord(contig_set['$dnanexus_link']).get_details()
    sizes = details['contigs']['sizes']
    names = details['contigs']['names']
    offsets = details['contigs']['offsets']

    commandList = []
    position = 0
    chromosome = 0
    currentChunk = 0

    for i in range(splits):
        commandList.append(" "+excludeInterval)

    includeDictionary = {}
    includeMatch = re.findall("(\w+):(\d+)-(\d+)", includeInterval)
    for x in includeMatch:
        if includeDictionary.get(x[0]) == None:
            includeDictionary[x[0]] = []
            includeDictionary[x[0]].append([int(x[1]), int(x[2])])

    while chromosome < len(names):
        if position + chunkSize >= sizes[chromosome]:
            #print chromosome
            commandList[currentChunk] += checkIntervalRange(includeDictionary, names[chromosome], position+1, sizes[chromosome])
            chromosome += 1
            position = 0
        else:
            commandList[currentChunk] += checkIntervalRange(includeDictionary, names[chromosome], position+1, position+chunkSize)
            position += chunkSize
        currentChunk = (currentChunk+1)%splits

    return commandList

def checkIntervalRange(includeList, chromosome, lo, hi):
    included = False
    command = ''
    if len(includeList) == 0:
        return " -L %s:%d-%d" % (chromosome, lo, hi)
    if includeList.get(chromosome) != None:
        for x in includeList[chromosome]:
            min = lo
            max = hi
            if (lo >= x[0] and lo <= x[1]) or (hi <= x[1] and hi >= x[0]):
                if lo >= x[0] and lo <= x[1]:
                    min = lo
                elif lo <= x[0]:
                    min = x[0]
                if hi <= x[1] and hi >= x[0]:
                    max = hi
                elif hi >= x[1]:
                    max = x[1]
                command += " -L %s:%d-%d" % (chromosome, min, max)
    return command

def splitGenomeLengthLargePieces(contig_set, chunks):
    details = dxpy.DXRecord(contig_set).get_details()
    sizes = details['contigs']['sizes']
    names = details['contigs']['names']
    offsets = details['contigs']['offsets']

    commandList = []
    for i in range(chunks):
        commandList.append('')
    position = 0
    chromosome = 0
    chunkSize = sum(sizes)/chunks
    currentChunk = 0
    currentLength = 0

    while chromosome < len(names):
        if position + (chunkSize - currentLength) >= sizes[chromosome]:
            commandList[currentChunk] += checkIntervalRange({}, names[chromosome], position+1, sizes[chromosome])
            currentLength += sizes[chromosome] - position
            chromosome += 1
            position = 0
        else:
            commandList[currentChunk] += checkIntervalRange({}, names[chromosome], position+1, position+(chunkSize-currentLength)+1)
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
