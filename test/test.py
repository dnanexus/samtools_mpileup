#!/usr/bin/env python
import os, sys, unittest, json, subprocess

import dxpy, dxpy.program_builder
from dxpy.exceptions import *

src_dir = os.path.join(os.path.dirname(__file__), "..")
test_resources_dir = os.path.join(src_dir, "test", "resources")

from optparse import OptionParser




def makeInputsBwa():
    try:
        contigset_importer = dxpy.DXProgram(dxpy.find_data_objects(classname="applet", properties={"name": "fasta_contigset_importer"}).next()['id'])
        reads_importer = dxpy.DXProgram(dxpy.find_data_objects(classname="applet", properties={"name": "Letter Space FASTQ importer"}).next()['id'])
    except StopIteration:
        raise Exception("fasta_contigset_importer or Letter Space FASTQ importer not found, please upload them")
    
    genome_archive = dxpy.upload_local_file(os.path.join(test_resources_dir, "hg19_chr22.fa.xz"), wait_on_close=True)
    contigset_importer_input = {"name": "hg19_chr22", "sequence_file": dxpy.dxlink(genome_archive)}
    print "Running fasta_contigset_importer with", contigset_importer_input
    job = contigset_importer.run(contigset_importer_input)
    job.wait_on_done()
    contig_set = job.describe()["output"]["contig_set"]
    
    left_reads = dxpy.upload_local_file(os.path.join(test_resources_dir, "small_left.fq"), wait_on_close=True)
    right_reads = dxpy.upload_local_file(os.path.join(test_resources_dir, "small_right.fq"), wait_on_close=True)
    #left_reads = dxpy.upload_local_file(os.path.join(test_resources_dir, "SRR188205_1_1M.fastq.xz"), wait_on_close=True)
    #right_reads = dxpy.upload_local_file(os.path.join(test_resources_dir, "SRR188205_2_1M.fastq.xz"), wait_on_close=True)

    reads_importer_input = {"left_file": dxpy.dxlink(left_reads), "right_file": dxpy.dxlink(right_reads)}
    print "Running LetterSpaceFileObjectToReadsTable with", reads_importer_input
    job = reads_importer.run(reads_importer_input)
    job.wait_on_done()
    reads = job.describe()["output"]["reads"]
    
    return {"reads": [reads] * 3, "reference": contig_set}
    
class TestMyApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if mappingsId == False:
            try:
                dxpy.DXGTable('gtable-9yZkVKj0000PFpq67bpQ0003').describe()
    
                small_reads = {u'$dnanexus_link': u'gtable-9yZkVKj0000PFpq67bpQ0003'}
                big_reads = {u'$dnanexus_link': u'gtable-9yZYyj000002B13g3JQQ0015'}
    
                cls.base_input = {'indexed_reference': {u'$dnanexus_link': u'file-9yZ40880000PkPJJ3GFQ001q'}}
                cls.base_input['reads'] = [small_reads] * 3
                
                # cls.base_input['reads'] = [big_reads] * 3
    
    #            dxpy.DXGTable('zrecord-9yZ3Qj80000PkPJJ3GFQ000G').describe()
    #            cls.base_input = {'reads': [{u'$dnanexus_link': u'gtable-9yZYyj000002B13g3JQQ0015'},
    #                                        {u'$dnanexus_link': u'gtable-9yZYyj000002B13g3JQQ0015'},
    #                                        {u'$dnanexus_link': u'gtable-9yZYyj000002B13g3JQQ0015'}],
    #                              #'reference': {u'$dnanexus_link': u'record-9yZ3Qj80000PkPJJ3GFQ000G'},
    #                              'indexed_reference': {u'$dnanexus_link': u'file-9yZ40880000PkPJJ3GFQ001q'}}
    #            'reads': [{u'$dnanexus_link': u'gtable-9yZkVKj0000PFpq67bpQ0003'},
    #                                        {u'$dnanexus_link': u'gtable-9yZkVKj0000PFpq67bpQ0003'}],
                                  #'reference': {u'$dnanexus_link': u'record-9yZ3Qj80000PkPJJ3GFQ000G'},
                                  
            except:
                cls.base_input = makeInputsBwa()
            
           
            
        bundled_resources = dxpy.program_builder.upload_resources(src_dir)
        program_id = dxpy.program_builder.upload_program(src_dir, bundled_resources, overwrite=True)
        cls.program = dxpy.DXProgram(program_id)

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_full_pipeline(self):
        if mappingsId == False:
            input = self.base_input
            print "Running program with", input
            try:   
                bwa = dxpy.DXProgram(dxpy.find_data_objects(classname="applet", properties={"name": "BWA"}).next()['id'])
            except:
                print "BWA not found, please upload it"
            job = bwa.run(input)
            job.wait_on_done()
            print "Bwa output:"
            print json.dumps(job.describe()["output"])
            mappings = job.describe()["output"]["mappings"][0]
        else:
            mappings = {"$dnanexus_link":mappingsId}
            print mappings
        print {'mappings':mappings}
        #job = self.program.run({'mappings':mappings, 'reference':{"$dnanexus_link":"record-9ykz7KQ00006B3PXk1b00005"}})
        
        #job = self.program.run({'mappings':mappings, 'reference':{"$dnanexus_link":"record-9zPp07j000035P6yJ9kQ0006"}})
        
        job = self.program.run({'mappings':mappings, 'reference':{"$dnanexus_link":"record-9zV2FBQ0000293088JZ00005"}})
        
        job.wait_on_done()
        print "Mpileup output:"
        print json.dumps(job.describe()["output"])

if __name__ == '__main__':
    #If a mappings_id is provided, test suite will assume you want to resume to test at GATK with the output from bwa
    parser = OptionParser("Usage: % mappings_id")
    parser.add_option("--mappings_id", dest="mappings_id", default=False, help="gtable with the bwa output of a previous run")
    (opts, args) = parser.parse_args()
    mappingsId = opts.mappings_id
    del sys.argv[1:]
    unittest.main()
