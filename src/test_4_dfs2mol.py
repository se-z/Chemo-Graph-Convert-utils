# -*- coding: utf-8 -*-

from utils import converter as cvt
import numpy as np
import sys
import argparse

CONFIG_PATH = "../config/atomMap.txt"
DB_PATH = "../demoDB/smiles"
WORK_DIR = "../out/work"

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="convert smiles format 2 adjancy matrix")
    parser.add_argument('id', help="data id for demo")

    args = parser.parse_args()
    data_id = args.id
    Converter = cvt.Converter(CONFIG_PATH)
    data_path = DB_PATH + "/" + data_id + ".txt"
    try:
        Converter.read(data_path)
        print("reader done")
    except Exception as e:
        print(e)
        sys.exit()

    A, v = Converter.mol2adjancy()

    id_label_dict = {int(k): int(v) for (k, v) in v}
    try:
        dfs_v, dfs_e = Converter.adjancy2dfs(A, id_label_dict)
    except Exception as e:
        print(e)
        sys.exit()

    # convert DFS2Mol
    incomplete_mol = Converter.dfs2Mol(dfs_v, dfs_e, WORK_DIR)
    print("---------------------------------------------------")
    print(incomplete_mol)
    # graphSVG(sfile + "/img-inverse.svg", inverse_mol)
