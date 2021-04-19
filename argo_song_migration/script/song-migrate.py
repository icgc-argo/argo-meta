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

def song_operation(endpoint, operation, headers, data=None):
    # wait a bit
    time.sleep(.5)
    try:
        if data is None:
            res = requests.put(endpoint, headers=headers)
        else:
            res = requests.put(endpoint, data=data, headers=headers)
        res.raise_for_status()
        
    except requests.exceptions.HTTPError as err:
        sys.exit("SONG %s failed, HTTPError: %s" % (operation, err))
    except requests.exceptions.RequestException as err:
        sys.exit("SONG %s failed, RequestException: %s" % (operation, err))

    if res.status_code != 200:
        sys.exit("SONG %s failed HTTP status code not 200: %s" % (operation, res))

def migrate_song_analysis(migrate_dump, song_url, token):

    headers = {
        "Authorization": "Bearer %s" % token,
        "Content-Type": "application/json",
        "Accept": "application/json"
    }

    with open(migrate_dump, 'r') as fp:
        for fline in fp:
            analysis = json.loads(fline)
            analysis.pop('analysisState', None)
            study = analysis.pop('studyId')
            analysis_id = analysis.pop('analysisId')
            files = analysis.pop('files')
            samples = analysis.pop('samples')
            
            print('Processing study: %s analysis: %s \n' % (study, analysis_id))
            # unpublish analysis
            endpoint = "%s/studies/%s/analysis/unpublish/%s" % (song_url, study, analysis_id)
            operation = 'analysis_unpublish'
            song_operation(endpoint, operation, headers)

            # update dynamic portion
            endpoint = "%s/studies/%s/analysis/%s" % (song_url, study, analysis_id)
            operation = 'analysis_update'
            data = json.dumps(analysis)
            song_operation(endpoint, operation, headers, data)

            # update file info
            for fl in files:
                object_id = fl.pop('objectId')
                endpoint = "%s/studies/%s/files/%s" % (song_url, study, object_id)
                operation = 'file_update'
                
                update_data = {
                    'dataType': fl['dataType'],
                    'info': {
                        'data_subtypes': fl['info'].get('data_subtypes', None),
                        'analysis_tools': fl['info'].get('analysis_tools', None),
                        'data_category': fl['info']['data_category'],
                        'description': fl['info'].get('description', None) 
                    }
                }

                data = json.dumps(update_data)
                song_operation(endpoint, operation, headers, data)

            # publish analysis
            endpoint = "%s/studies/%s/analysis/publish/%s" % (song_url, study, analysis_id)
            operation = 'analysis_publish'
            song_operation(endpoint, operation, headers)


def main():
    parser = ArgumentParser()
    parser.add_argument("-d", "--dump_path", dest="dump_path", type=str, required=True, help="path to migrated dump jsonl file")
    parser.add_argument("-m", "--song_url", dest="song_url", type=str, default="https://song.rdpc-qa.cancercollaboratory.org", help="SONG URL")
    parser.add_argument("-t", "--token", dest="token", type=str, required=True)
    args = parser.parse_args()

    # song dump analysis
    migrate_dump = args.dump_path

    # apply song migration
    migrate_song_analysis(migrate_dump, args.song_url, args.token)

    

if __name__ == "__main__":
    main()
