import os
import tempfile
import argparse
import subprocess

root = os.path.dirname(os.path.abspath(__file__))

tmp_dir = tempfile.mkdtemp(prefix='vsflow3d-')


def main() -> None:
    args = parseArgs()

    query_file = args.query_file
    database_file = args.database_file
    output_folder = args.output_folder

    cmd = f'vsflow3d shape -q {query_file} -d {database_file} -o {output_folder} -t 1000'

    subprocess.Popen(cmd, shell=True).wait()

    results_file = os.path.join(output_folder, 'results.csv')
    
    with open(output_folder, 'r') as f:
        pass



def parseArgs() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Generates conformers for a given list of SMILES strings.')

    parser.add_argument('-q',
                        dest='query_file',
                        required=True,
                        metavar='<file>',
                        help='Query molecule in SDF format. Multiple conformers are permitted.')
    
    parser.add_argument('-d',
                        dest='database_file',
                        required=True,
                        metavar='<file>',
                        help='Database molecules in SDF format. Multiple conformers are permitted.')
    
    parser.add_argument('-o',
                        dest='output_folder',
                        required=True,
                        metavar='<folder>',
                        help='Output folder.')
    
    return parser.parse_args()


if __name__ == '__main__':
    main()