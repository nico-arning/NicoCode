"""
Script to generate Powerpoint from the bridge.nf outputs.
This script is deployed by the bridge.nf nextflow pipleine. 
All neccessary functions are contained in the bridge_library.py file.
"""
from bridge_library import PowerPoint
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make a powerpoint presentation from the outputs of the BRIDGE pipeline')
    parser.add_argument(
        '-f',
        '--folder',
        type=str, 
        help="""Folder with tables stored by the DeSeq  analysis usually this is in 
        the project directory under <accession>/DeSeq/Tables"""
    )
    args = parser.parse_args()
    powerpointer = PowerPoint(args.folder)
    powerpointer.make_presentation()

