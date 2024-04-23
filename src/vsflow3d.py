import os
import tempfile
import argparse
import subprocess
import shutil
import csv
from rdkit import Chem


root = os.path.dirname(os.path.abspath(__file__))

tmp_dir = tempfile.mkdtemp(prefix='vsflow3d-')


def main() -> None:
    args = parseArgs()

    query_file = args.query_file
    database_file = args.database_file
    output_tmp = os.path.join(tmp_dir, "output")
    output_folder = args.output_folder
    top_hits = int(args.top)

    cmd = f'vsflow shape -i {query_file} -d {database_file} -o {output_tmp} -t {top_hits} --boost --nproc 12'

    subprocess.Popen(cmd, shell=True).wait()

    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    os.makedirs(output_folder)

    results_file = os.path.join(output_folder, 'results.csv')
    results_sdf = os.path.join(output_folder, 'results.sdf')
    
    shutil.copy(output_tmp+"_1.sdf", results_sdf)

    reader = Chem.SDMolSupplier(results_sdf)
    R = []
    for mol in reader:
        if mol is None:
            continue
        index = mol.GetProp("_InputIndex")
        smiles = Chem.MolToSmiles(mol)
        query_smiles = mol.GetProp("QuerySmiles")
        combo_score = mol.GetProp("Combo_Score")
        shape_sim = mol.GetProp("Shape_Similarity")
        fp_sim = mol.GetProp("3D_FP_Similarity")
        R += [[index, smiles, query_smiles, combo_score, shape_sim, fp_sim]]

    with open(results_file, "w") as f:
        writer = csv.writer(f)
        header = ["InputIndex", "Smiles", "QuerySmiles", "ComboScore", "ShapeSim", "Fp3dSim"]
        writer.writerow(header)
        for r in R:
            writer.writerow(r)
    

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
    
    parser.add_argument('-t',
                        dest='top',
                        required=False,
                        default=1000,
                        help='Top number of hits.')
    
    return parser.parse_args()


if __name__ == '__main__':
    main()