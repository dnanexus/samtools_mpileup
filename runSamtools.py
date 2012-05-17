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
        ##NEED TO UPDATE TO NEW CONTIGSET SPELLING AFTER TESTING DONE
        #contigSetId = mappingsTable.get_details()['originalContigSet']['$dnanexus_link']
        #originalContigSet = mappingsTable.get_details()['originalContigSet']
        
        contigSetId = mappingsTable.get_details()['original_contigset']['$dnanexus_link']
        originalContigSet = mappingsTable.get_details()['original_contigset']
    except:
        raise Exception("The original reference genome must be attached as a detail")
    
    subprocess.check_call("contigset2fasta %s ref.fa" % (contigSetId), shell=True)
    reference_sequence = dxpy.dxlink(dxpy.upload_local_file("ref.fa"))

    print "Indexing Dictionary"
    subprocess.check_call("samtools faidx ref.fa", shell=True)
    
    referenceIndex = dxpy.dxlink(dxpy.upload_local_file("ref.fa.fai"))

    mappings_schema = [
            {"name": "chr", "type": "string"}, 
            {"name": "lo", "type": "int32"},
            {"name": "hi", "type": "int32"},
            {"name": "type", "type": "string"},     #change this type to uint once there is an abstraction method for enum
            {"name": "ref", "type": "string"},
            {"name": "alt", "type": "string"},
            {"name": "qual", "type": "int32"},
            {"name": "coverage", "type": "int32"},
            {"name": "genotypeQuality", "type": "int32"},    
        ]
    if job['input']['store_full_vcf']:
        mappings_schema.extend([{"name": "vcf_alt", "type": "string"}, {"name": "vcf_additional_data", "type": "string"}])
    
    simpleVar = dxpy.new_dxgtable(mappings_schema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')])
    tableId = simpleVar.get_id()
    simpleVar = dxpy.open_dxgtable(tableId)
    
    
    reduceInput = {}
    commandList = splitGenomeLength(originalContigSet, job['input']['intervals_to_process'], job['input']['intervals_to_exclude'],  job['input']['minimum_chunk_size'], job['input']['maximum_chunks'])
    
    samOptions = makeSamtoolsParameters(job)
    bcfOptions = makeBcftoolsParameters(job)
    
    
    
    for i in range(len(commandList)):
        print commandList[i]
        if len(commandList[i]) > 0:
            mapInput = {
                'mappings_table_id':mappingsTableId,
                'reference_sequence': reference_sequence,
                'reference_index': referenceIndex,
                'original_contig_set': contigSetId,
                'interval': commandList[i],
                'tableId': tableId,
                'compress_reference': job['input']['compress_reference'],
                'compress_no_call' : job['input']['compress_no_call'],
                'store_full_vcf' : job['input']['store_full_vcf'],
                'sam_options': samOptions,
                'bcf_options': bcfOptions, 
                'part_number' : i
            }
            # Run a "map" job for each chunk
            mapJobId = dxpy.new_dxjob(fn_input=mapInput, fn_name="mapPileup").get_id()
            reduceInput["mapJob" + str(i) + "TableId"] = {'job': mapJobId, 'field': 'id'}

    reduceInput['tableId'] = tableId
    reduceJobId = dxpy.new_dxjob(fn_input=reduceInput, fn_name="reduceGatk").get_id()
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
    
    print "Converting Table to SAM"
    print "dx_mappingsTableToSam --table_id %s --output input.sam --region_index_offset -1 %s" % (job['input']['mappings_table_id'], job['input']['interval'])
    subprocess.check_call("dx_mappingsTableToSam --table_id %s --output input.sam --region_index_offset -1 %s" % (job['input']['mappings_table_id'], job['input']['interval']), shell=True)
    print "Converting to BAM"
    subprocess.check_call("samtools view -bS input.sam > input.bam", shell=True)
    print "Sorting"
    subprocess.check_call("samtools sort input.bam input.sorted", shell=True)
    print "Indexing"
    subprocess.check_call("samtools index input.sorted.bam", shell=True)
    
    referenceFileName = dxpy.download_dxfile(job['input']['reference_sequence'], "ref.fa")
    indexFileName = dxpy.download_dxfile(job['input']['reference_index'], "ref.fa.fai")
    
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
            print x[0]+"\t"+str(x[1])+"\t"+str(x[2])+"\n"
        bedFile.close()
        command += " -l regions.bed"
    command += job['input']['sam_options']
    command += " input.sorted.bam | bcftools view "
    command += job['input']['bcf_options']
    command += " - > output.vcf"
    print command
    subprocess.call(command, shell=True)
    
    parseVcf(open("output.vcf", 'r'), simpleVar, job['input']['compress_reference'], job['input']['compress_no_call'], job['input']['store_full_vcf'])
    
    if job['input']['part_number'] == 0:
        header = extractHeader(open("output.vcf", 'r'))
        simpleVar.set_details({'header':header, 'original_contigset':dxpy.dxlink(job['input']['original_contig_set'])})


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
    
def reduceGatk():
    t = dxpy.open_dxgtable(job['input']['tableId'])
    t.close(block=True)
    print "Closing Table"
    job['output']['simplevar'] = dxpy.dxlink(t.get_id())
    
def parseVcf(vcfFile, simpleVar, compressNoCall, compressReference, storeFullVcf):

    #These prior variables are used for keeping track of contiguous reference/no-call
    #   in the event that compressReference or compressNoCall is True
    priorType = "None"
    priorPosition = -1

    fileIter = vcfFile.__iter__()
    count = 1

    #Additional data will contain the extra format and info columns that are optional in VCF and may not be
    #   present in the VCF file. These are stored in an extended table 
    additionalData = []
    
    while 1:
        try:
            input = fileIter.next()
            if count%100000 == 0:
                print "Processed count %i variants " % count
            count += 1
            
            if input[0] != "#":
                tabSplit = input.strip().split("\t")
                chr = tabSplit[0]
                lo = int(tabSplit[1])
                hi = lo + len(tabSplit[3])
                ref = tabSplit[3].replace(".","")
                
                #In VCF format, the ALT column holds possible candidate alleles. The actual call as to the
                #   variant and its zygosity is a combination of ALT and the genotype specified in the info field.
                #   We store all of the options (including ref) and calculated the actual calls later
                altOptions = [ref]
                altOptions.extend(tabSplit[4].split(","))
                qual = tabSplit[5]
                type = "Unknown"
                if qual == ".":
                    type = "No-call"
                else:
                    qual = int(float(tabSplit[5]))

                formatColumn = tabSplit[7]
                infoColumn = tabSplit[8]
                genotypeQuality = 0
                
                coverage = re.findall("DP=(\d+);", formatColumn)
                if(len(coverage) > 0):
                    coverage = int(coverage[0])
                else:
                    coverage = 0

                if altOptions == [ref, '.']:
                    if type == "No-call":
                        if compressNoCall == False:
                            entry = [chr, lo, hi, type, "", "", 0, 0, 0]
                            entry.append(tabSplit[4])
                            vcfSpecificData = ''
                            for x in tabSplit[7:]:
                                vcfSpecificData += x+"\t"
                            entry.append(tabSplit[7:].strip())
                    else:
                        type = "Ref"
                        if compressReference == False:
                            entry = [chr, lo, hi, type, "", "", 0, 0, 0]
                            entry.append(tabSplit[4])
                            vcfSpecificData = ''
                            for x in tabSplit[7:]:
                                vcfSpecificData += x+"\t"
                            entry.append(vcfSpecificData.strip())
                else:
                    #In VCF format, the prior character to a sequence change is given in some cases (Ins, Del)
                    #   we are removing this in our format, and so need to figure out which characters to filter   
                    overlap = findMatchingSequence(ref, altOptions, False)
                    ref = ref[overlap:]
                    for i in range(len(altOptions)):
                        altOptions[i] = altOptions[i][overlap:]
                    reverseOverlap = findMatchingSequence(ref, altOptions, True)

                    ref = ref[:len(ref)-reverseOverlap]
                    for i in range(len(altOptions)):
                        altOptions[i] = altOptions[i][:len(altOptions[i])-reverseOverlap]
                    
                    
                    
                    #Find all of the genotypes 
                    genotypePossibilities = {}
                    for x in tabSplit[9:]:
                        genotype = getInfoField("GT", infoColumn, x)
                        genotypeQuality = float(getInfoField("GQ", infoColumn, x))
                        if genotype != False and genotypeQuality != False:
                            if genotypePossibilities.get(genotype) == None:
                                genotypePossibilities[genotype] = float(genotypeQuality)
                            else:
                                genotypePossibilities[genotype] += float(genotypeQuality)
                        else:
                            genotypeQuality = 0
                    genotypePossibilities = sorted(genotypePossibilities.iteritems(), key=operator.itemgetter(1), reverse=True)
                    if len(genotypePossibilities) > 0:
                        genotype = genotypePossibilities[0][0]
                        genotypeQuality = genotypePossibilities[0][1]
                        if len(genotypePossibilities) > 1:
                            genotypeQuality -= genotypePossibilities[1][1]
                        alt = ""
                        if genotype == "0/0" or genotype == "0|0" or genotype == False:
                            if(len(genotypePossibilities) > 1):
                                genotype = genotypePossibilities[1][0]
                                genotypeQuality = 0
                                
                        genotypeSplit = re.split("[\|\/]", genotype)
                        for i in range(len(genotypeSplit)):
                            
                        #This is done to ensure the convention of placing the ref allele first
                        #   in practice, it seems that all VCFs already place the ref first
                            genotypeSplit[i] = int(genotypeSplit[i])
                        genotypeSplit.sort()
    
                        
                        alleles = []
                        for x in genotypeSplit:
                            if len(alt) > 0:
                                alt += "/"
                            alt += altOptions[x]
                            if x != 0:
                                alleles.append(altOptions[x])
                            if len(altOptions[x]) == 0:
                                alt += "-"
                    else:
                        alt = "?"
                        for x in altOptions:
                            alt += "/"+x
                        genotypeQuality = 0
                        
                    typeList = []                
                    #These rules determine how to characterize the type of change that has occurred
                    for x in alleles:
                        if len(x) == len(ref) and len(ref) == 1:
                            type = "SNP"
                        elif ref in x:
                            type = "Ins"
                        elif x in ref:
                            type = "Del"
                        else:
                            type = "Complex"
                        for x in typeList[1::]:
                            if typeList[0] != x:
                                type = "Mixed"

                    hi = lo+len(ref)-1
                    if len(ref) == 0:
                        ref = "-"
                    entry = [chr, lo+overlap-1, hi, type, ref.upper(), alt.upper(), qual, coverage, int(genotypeQuality)]
                    if storeFullVcf:
                        entry.append(tabSplit[4])
                        vcfSpecificData = ''
                        for x in tabSplit[7:]:
                            vcfSpecificData += x+"\t"
                        entry.append(vcfSpecificData.strip())
                    if type == "Ins" or type == "Complex":
                        print "Ins"
                        print entry
                        print input
                    #print entry
                    simpleVar.add_rows([entry])
                if compressReference:
                    if priorType == "Ref" and type != priorType:
                        entry = [chr, priorPosition, lo, type, "", "", 0, 0, 0]
                        if storeFullVcf:
                            entry.extend([".", ""])
                        simpleVar.add_rows([entry])                        
                if compressNoCall:
                    if priorType == "No-call" and type != priorType:
                        entry = [chr, priorPosition, lo, type, "", "", 0, 0, 0]
                        if storeFullVcf:
                            entry.extend([".",""])
                        simpleVar.add_rows([entry])
                if type != priorType:
                    priorType = type
                    priorPosition = lo-1
        except StopIteration:
            break

def findMatchingSequence(ref, altOptions, reverse):
    r = 1
    if reverse == True:
        r = -1
    position = 0
    minLength = len(ref)
    for x in altOptions:
        if len(x[::r]) < minLength:
            minLength = len(x[::r])
    for i in range(minLength):
        for x in altOptions:
            if ref[::r][i] != x[::r][i]:
                return i
    return minLength

def getInfoField(fieldName, infoColumn, infoContents):
    if infoColumn.count(fieldName) > 0:
        entrySplitColumn = infoColumn.split(":")
        position = -1
        for i in range(len(entrySplitColumn)):
            if entrySplitColumn[i] == fieldName:
                position = i
                entrySplitInfo = infoContents.split(":")
                if len(entrySplitInfo) == len(entrySplitColumn):
                    return entrySplitInfo[position]
    return False
    
def generateEmptyList(columns):
    result = []
    for i in range(columns):
        result.append('')
    return result
    
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
            print chromosome
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
        return " -l %s:%d-%d" % (chromosome, lo, hi)
    if includeList.get(chromosome) != None:
        for x in includeList[chromosome]:
            print "List"
            print x
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
                command += " -l %s:%d-%d" % (chromosome, min, max)
    return command

