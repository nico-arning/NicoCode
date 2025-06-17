"""
Script to interact with the s3 bucket stored within the datapond.
This script is deployed by the bridge.nf nextflow pipleine. 
All neccessary functions are contained in the bridge_library.py file.
"""
from bridge_library import Datapond
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download files from BRIDGE datapond')
    parser.add_argument('-f','--folder', type=str, 
        help='folder to look into')
    parser.add_argument('-u','--upload', type=str, default = "",
        help='folder to upload')
    parser.add_argument('-d','--download', type=str, default = "",
        help='folder to download')                    
    datapond = Datapond()
    args = parser.parse_args()

    if args.upload != "":
        datapond.upload_files(
            args.upload, 
            args.folder)
    
    if args.download != "":
        datapond.download_files(
            args.download, 
            args.folder)

