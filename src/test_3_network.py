# -*- coding: utf-8 -*-

from utils import converter as cvt
import numpy as np
import sys
import argparse

CONFIG_PATH = "../config/atomMap.txt"
DB_PATH = "../demoDB/dfscode"
OUT_PATH = "../out/dfs2pic"

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="draw network picture")
    parser.add_argument('pic_num', default=5,
                        help="number of network pictures.")
    args = parser.parse_args()
    pic_mum = args.pic_mum

    Converter = cvt.Converter(CONFIG_PATH)
    Converter.graphDraw(OUT_PATH, DB_PATH + "/dfs.txt", pic_mum)
