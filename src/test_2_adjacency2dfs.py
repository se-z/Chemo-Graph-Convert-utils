# -*- coding: utf-8 -*-

from utils import converter as cvt
import numpy as np
import sys
import argparse

CONFIG_PATH = "../config/atomMap.txt"
DB_PATH = "../demoDB/smiles"
OUT_PATH = "../out/dfscode"

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="convert adjacency matrix to DFS Code")
    parser.add_argument('id', help="data_id")
    parser.add_argument(
        'out_path', default="../out/dfscode", help="output dir. default is ../out/dfscode ")

    args = parser.parse_args()
    data_id = args.id
    OUT_PATH = args.out_path

    Converter = cvt.Converter(CONFIG_PATH)
    data_path = DB_PATH + "/" + data_id + ".txt"
    try:
        Converter.read(data_path)
        print("reader done")
    except Exception as e:
        print(e)
        sys.exit()

    A, v = Converter.mol2adjacency()

    # convert adjacency to DFS Code
    id_label_dict = {int(k): int(v) for (k, v) in v}
    try:
        dfs_v, dfs_e = Converter.adjacency2dfs(A, id_label_dict)
    except Exception as e:
        print(e)
        sys.exit()

    print(dfs_v)
    print(dfs_e)

    with open(OUT_PATH + "/dfs.txt", mode="w") as f:
        for item in dfs_v:
            f.writelines("v {} {} \n".format(item[0], item[1]))
        for item in dfs_e:
            f.writelines("e {} {} {} \n".format(item[0], item[1], item[2]))
