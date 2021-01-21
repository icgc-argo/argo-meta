#!/usr/bin/python
import os
import re
import json
import sys, getopt
import argparse
import csv

# Author: Hardeep Nahal-Bose hardeep.nahal@oicr.on.ca
# Script used to generate required data for GoogleDoc worksheet. This requires the intermediate song report generated after SONG payloads are imported (contains mapping information between legacy ICGC analysis ID and ARGO analysis IDs) and the summary file after processing (icgc_song_payloads/<ARGO program code>/<sequencing strategy>/summary_tables/<summary_table_file>).

## Example inputs:
## summary_file: https://raw.githubusercontent.com/icgc-argo/argo-meta/master/icgc_song_payloads/PACA-CA/WGS/PACA-CA_WGS_batch1_summary.tsv
## intermediate song mapping file: https://raw.githubusercontent.com/icgc-argo/argo-meta/master/icgc_song_payloads/PACA-CA/WGS/intermediate-song-rdpc-collab-report.PACA-CA.WGS.batch1.json
## output_file: worksheet.tsv
##
## Copy contents of worksheet.tsv into main GoogleDoc worksheet


parser = argparse.ArgumentParser()
parser.add_argument("-s", "-summary_file", dest="summary_file", required=True)
parser.add_argument("-i", "-intermediate_song_file", dest="intermediate_song_file", required=True)
parser.add_argument("-o", "-out_file", dest="out_file", required=True)
parser.add_argument("-p", "-project", dest="project", required=True)
args = parser.parse_args()

intermediateReport = args.intermediate_song_file

summaryFile = args.summary_file
outFile = open(args.out_file, "w")
outFile.write("submitter_donor_id\tsubmitter_sample_id\ttumour_normal_designation\tsequencing_strategy\tdonor_id\tsample_id\tstudy_id\tanalysis_id\n")
completed = {}



def get_analysisId_mapping():
   mapping = {}
   mapping_file = open("%s_analysis_id_mapping.tsv"%(args.project), "w")
   report = json.load(open(intermediateReport, "r"))
   for payload in report["success"]:
      legacyAnalysisIds = payload["legacyAnalysisIds"]
      for legacy_id in legacyAnalysisIds:
         mapping[legacy_id] = payload["targetAnalysisId"]

   for legacy_id in mapping:
      mapping_file.write("%s\t%s\n"%(legacy_id, mapping[legacy_id]))
   mapping_file.close()
   return mapping

mapping = get_analysisId_mapping()
reader = csv.DictReader(open(summaryFile), delimiter="\t")
for row in reader:
   icgc_analysis_id = row["bundle_id"]
   print(icgc_analysis_id)
   if icgc_analysis_id in mapping and icgc_analysis_id not in completed:
      print("Processing icgc_analysis_id=%s"%icgc_analysis_id)
      argo_analysis_id = mapping[icgc_analysis_id]
      dcc_donor_id = row["donor_id/donor_count"]
      submitter_donor_id = row["submitter_donor_id"]
      specimen_type = row["specimen_type"]
      dcc_sample_id = row["dcc_sample_id"]
      submitter_sample_id = row["submitter_sample_id"]
      sequencing_strategy = row["sequencing_strategy"]
      study_id = row["project_id/project_count"]
      if re.search("tumour", specimen_type, re.IGNORECASE):
         tumour_normal_designation = "tumour"
      if re.search("normal", specimen_type, re.IGNORECASE):
         tumour_normal_designation = "normal"

      outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(submitter_donor_id, submitter_sample_id, tumour_normal_designation, sequencing_strategy, dcc_donor_id, dcc_sample_id, study_id, argo_analysis_id))
   completed[icgc_analysis_id] = argo_analysis_id 
outFile.close()
