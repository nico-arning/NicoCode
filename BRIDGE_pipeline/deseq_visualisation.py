"""
Script to do some additional visualisation from DeSeq outputs.
This script is deployed by the bridge.nf nextflow pipleine. 
All neccessary functions are contained in the bridge_library.py file
"""
from bridge_library import DeSeqViz, DoseResponse
import argparse
import warnings
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Some extra python visualisation for deseq')
    parser.add_argument(
        '-f',
        '--folder',
        type=str, 
        help="""Folder with tables stored by the DeSeq  analysis usually this is in 
        the project directory under <accession>/DeSeq/Tables"""
    )
    parser.add_argument(
        '-c',
        '--cores',
        type=int, 
        default = 4,
        help="""how many cores to use"""
    )
    parser.add_argument(
        '-g',
        '--genome_dir',
        type=str, 
        default = "",
        help="""which genome directory to use for annotation"""
    )
                           
    args = parser.parse_args()

    visualiser = DeSeqViz(args.folder, args.genome_dir)
    #visualiser.well_pcas()
    visualiser.volcano_plot()
    visualiser.leiden_clustering()
    visualiser.read_decay()
    visualiser.replicate_correlations()

    doseresponser = DoseResponse(args.folder, args.genome_dir, args.cores)
    doseresponser.dose_response()


