# -*- coding: utf-8 -*-

from utils import converter as cvt
import numpy as np
import sys
import argparse

CONFIG_PATH = "../config/atomMap.txt"
DB_PATH = "../demoDB/smiles"
OUT_PATH = "../out/adjacency"

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="convert smiles format 2 adjacency matrix")
    parser.add_argument('id', help="data_id ")

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

    A, v = Converter.mol2adjacency()
    # labeled adjacency matrix
    for i in range(len(A)):
        print(A[i, :])

    # [vertex id, vertex lavel ]
    print(v)
