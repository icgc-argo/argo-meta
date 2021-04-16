#!/usr/bin/env python3

import json
import os
import csv
import glob
from argparse import ArgumentParser
import sys
import subprocess
from collections import OrderedDict
import requests
import re
import tarfile
import time
from datetime import date


def prep_migrate_payloads(song_dump):
    date_str = date.today().strftime("%Y-%m-%d")
    output_dir = os.path.join(os.getcwd(), 'argo_song_migration', date_str)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    migrate_qc_metrics = open(os.path.join(output_dir, 'migrate_qc_metrics.jsonl'), 'w')
    migrate_aln = open(os.path.join(output_dir, 'migrate_sequencing_alignment.jsonl'), 'w')
    migrate_vc = open(os.path.join(output_dir, 'migrate_variant_calling.jsonl'), 'w')
    migrate_vc_supplement = open(os.path.join(output_dir, 'migrate_variant_calling_supplement.jsonl'), 'w')
    with open(song_dump, 'r') as fp:
        for fline in fp:
            analysis = json.loads(fline)
            analysis_state = analysis.pop('analysisState')
            if not analysis_state == 'PUBLISHED': continue
            analysisVersion = analysis['analysisType'].pop('version')
            for item in ['createdAt', 'createdAt', 'firstPublishedAt', 'publishedAt', 'analysisStateHistory', 'updatedAt']:
                if analysis.get(item):
                    analysis.pop(item)
            
            if analysis.get('info') and analysis['info'].get('origin') == 'ICGC-25K': continue
            analysis['info'] = {'origin': 'ICGC-25K'}
            change = True

            if analysis['analysisType']['name'] == 'qc_metrics':
                #change = False
                # files section
                for fl in analysis['files']:
                    if not fl.get('info') or not fl['info'].get('data_category') or not fl['info'].get('analysis_tools') or not fl['info'].get('description'): continue
                    if not fl['info']['data_category'] == "Quality Control Metrics": continue
                    if fl['dataType'] == 'Alignment QC':
                        fl['dataType'] = 'Aligned Reads QC'
                        fl['info']['data_subtype'] = ['Alignment Metrics']
                        change = True
                    elif fl['dataType'] == 'Cross Sample Contamination':
                        fl['dataType'] = 'Sample QC'
                        fl['info']['data_subtype'] = ['Cross Sample Contamination']
                        change = True
                    elif fl['dataType'] == 'Duplicates Metrics':
                        fl['dataType'] = 'Aligned Reads QC'
                        fl['info']['data_subtype'] = ['Duplicates Metrics']
                        change = True
                    elif fl['dataType'] == 'Genotyping Inferred Gender':
                        fl['dataType'] = 'Analysis QC'
                        fl['info']['data_subtype'] = ['Genotyping Stats']
                        change = True
                    elif fl['dataType'] in ['Mutect2 Callable Stats', 'Mutect2 Callabe Stats']:
                        fl['dataType'] = 'Analysis QC'
                        fl['info']['data_subtype'] = ['Variant Callable Stats']
                        change = True
                    elif fl['dataType'] == 'Mutect2 Filtering Stats':
                        fl['dataType'] = 'Analysis QC'
                        fl['info']['data_subtype'] = ['Variant Filtering Stats']
                        change = True
                    elif fl['dataType'] == 'OxoG Metrics':
                        fl['dataType'] = 'Aligned Reads QC'
                        fl['info']['data_subtype'] = ['OxoG Metrics']
                        change = True
                    elif fl['dataType'] == 'Ploidy and Purity Estimation':
                        fl['dataType'] = 'Analysis QC'
                        fl['info']['data_subtype'] = ['Ploidy', 'Tumour Purity']
                        change = True
                    elif fl['dataType'] == 'Read Group QC':
                        fl['dataType'] = 'Sequencing QC'
                        fl['info']['data_subtype'] = ['Read Group Metrics']
                        change = True
                    else:
                        pass

                    for index, tool in enumerate(fl['info']['analysis_tools']):
                        if tool in ['GATK-CalculateContamination', 'GATK-FilterMutectCalls', 'GATK-Mutect2']:
                            fl['info']['analysis_tools'][index] = tool.replace("-", ":")
                            change = True
                        elif tool in ['bas_stats']:
                            fl['info']['analysis_tools'][index] = 'Sanger:bam_stats'
                            change = True
                        elif tool in ['compareBamGenotypes']:
                            fl['info']['analysis_tools'][index] = 'Sanger:compareBamGenotypes'
                            change = True
                        elif tool in ['verifyBamHomChk']:
                            fl['info']['analysis_tools'][index] = 'Sanger:verifyBamHomChk'
                            change = True
                        else:
                            pass

                    if fl['info']['description'] == 'Alignment QC metrics generated by Sanger bas_stats.pl script':
                        fl['info']['description'] = 'Alignment QC metrics generated by Sanger bam_stats script'  
                        change = True
                
                if change:
                    migrate_qc_metrics.write(json.dumps(analysis)+"\n")
            
            elif analysis['analysisType']['name'] == 'variant_calling':
                #change = False
                for fl in analysis['files']:
                    if not fl.get('info') or not fl['info'].get('data_category') or not fl['info'].get('analysis_tools'): continue
                    if not fl['info']['data_category'] == "Simple Nucleotide Variation": continue
                    for index, tool in enumerate(fl['info']['analysis_tools']):
                      if tool == 'GATK-Mutect2':
                          fl['info']['analysis_tools'][index] = 'GATK:Mutect2'
                          change = True
                      else:
                          pass
                if change:
                    migrate_vc.write(json.dumps(analysis)+"\n")

            elif analysis['analysisType']['name'] == 'variant_calling_supplement':
                #change = False
                for fl in analysis['files']:
                    if not fl.get('info'): continue
                    if fl['fileName'].endswith('.timings-supplement.tgz'):
                        if not fl['info'].get('data_subtype'):                                
                            fl['dataType'] = 'Analysis QC'
                            fl['info']['data_category'] = 'Quality Control Metrics'
                            fl['info']['data_subtype'] = ['Runtime Stats']
                            fl['info']['analysis_tools'] = None
                            change = True
                    elif fl['fileName'].endswith('.pindel-supplement.tgz'):
                        if fl['dataType'] == 'Variant Calling Supplement':
                            fl['dataType'] = 'InDel Supplement'
                            change = True
                    elif fl['fileName'].endswith('.caveman-supplement.tgz'):
                        if fl['dataType'] == 'Variant Calling Supplement':
                            fl['dataType'] = 'SNV Supplement'
                            change = True
                    elif fl['fileName'].endswith('.brass-supplement.tgz'):
                        if fl['dataType'] == 'Variant Calling Supplement':
                            fl['dataType'] = 'SV Supplement'
                            change = True
                    elif fl['fileName'].endswith('.ascat-supplement.tgz'):
                        if fl['dataType'] == 'Variant Calling Supplement':
                            fl['dataType'] = 'CNV Supplement'
                            change = True
                    else:
                        pass
                if change:
                    migrate_vc_supplement.write(json.dumps(analysis)+"\n")

            elif analysis['analysisType']['name'] == 'sequencing_alignment':
                #change = False
                for fl in analysis['files']:
                    if not fl.get('info') or not fl['info'].get('data_category') or not fl['info'].get('analysis_tools'): continue
                    change = True

                if change:
                    migrate_aln.write(json.dumps(analysis)+"\n")
            else:
                pass

    migrate_qc_metrics.close()
    migrate_aln.close()
    migrate_vc.close()
    migrate_vc_supplement.close()

def main():
    parser = ArgumentParser()
    parser.add_argument("-d", "--dump_path", dest="dump_path", type=str, required=True, help="path to song dump jsonl file")
    args = parser.parse_args()
    
    # song dump analysis
    song_dump = args.dump_path

    # process the qc files and prep qc_metrics payload
    prep_migrate_payloads(song_dump)
    

if __name__ == "__main__":
    main()
