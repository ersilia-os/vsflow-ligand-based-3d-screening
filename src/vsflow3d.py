import os
import tempfile
import argparse
import subprocess

root = os.path.dirname(os.path.abspath(__file__))

tmp_dir = tempfile.mkdtemp(prefix='vsflow3d-')


def main() -> None:
    args = parseArgs()


def parseArgs() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Generates conformers for a given list of SMILES strings.')

    parser.add_argument('-i',
                        dest='in_file',
                        required=True,
                        metavar='<file>',
                        help='Molecule input file in CSV format. A column named "smiles" is required.')



    return parser.parse_args()


if __name__ == '__main__':
    main()