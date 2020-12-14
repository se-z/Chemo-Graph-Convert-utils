# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# mpl.use('Agg')
import re
import os
import csv
import sys
from rdkit import rdBase, Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
import networkx as nx


class SMILESParseError_1(Exception):
    """
        this processing algorithm remove Asymmetric center expression
    """
    pass


class SMILESParseError_2(Exception):
    """
        this processing algorithm remove ion
    """
    pass


class SMILESParseError_3(Exception):
    """
        this processing algorithm remove > expression
    """
    pass


class SMILESParseError_4(Exception):
    """
        size, ion 
    """
    pass


class NeighborhoodRestrictionError(Exception):
    """
        gSpan's neighborhood Restriction Check
    """
    pass


class Converter:

    def __init__(self, config_file):
        self.ATOM2ID = {}
        self.ID2ATOM = {}
        with open(config_file, "r") as f:
            lines = f.readlines()

        for l in lines:
            atom = re.split(" +", l)[0]
            ID = (re.split(" +", l)[1]).rstrip("\n")
            self.ATOM2ID[atom] = ID
            self.ID2ATOM[ID] = atom

    def __isRemove(self):
        NUM_VER_EDGE_LINE = 3
        mol_list = self.MOL_str.split("\n")

        nums = re.split(" +", mol_list[NUM_VER_EDGE_LINE])
        v_num = nums[1]
        e_num = nums[2]

        # サイズが100を超えるものを削除しておく
        if int(v_num) >= 100 or int(e_num) >= 100:
            isRemove = True
            return isRemove

        # 境界値バグの原因になるので, 小さすぎるものも削除
        if int(v_num) == 1 or int(e_num) == 1:
            isRemove = True
            return isRemove

        VERTEX_FIRST_LINE = 4
        VERTEX_LAST_LINE = 4 + int(v_num)
        EDGE_FIRST_LINE = 4 + int(v_num)
        EDGE_LAST_LINE = EDGE_FIRST_LINE + int(e_num)

        vertex_list = mol_list[VERTEX_FIRST_LINE:VERTEX_LAST_LINE]
        edge_list = mol_list[EDGE_FIRST_LINE:EDGE_LAST_LINE]
        footer_list = mol_list[EDGE_LAST_LINE:]

        isRemove = False
        if not((int(v_num) == len(vertex_list)) and (int(e_num) == len(edge_list))):
            print("(vertex num, edge num, vertex list, edge list) = ({}, {}, {}, {})".format(
                v_num, e_num, len(vertex_list), len(edge_list)))
            isRemove = True
            return isRemove

        if len(footer_list) > 2:
            isRemove = True
            return True

        isIon = False
        for itr, item in enumerate(footer_list):
            marker = (re.split(" +", item))[1]
            if marker == "END":
                break

            if marker == "CHG":
                isION = True
                return isIon

        return isIon

    def __checkNeiborhoodRestriction(self, edge_set):

        p_e_1 = 0
        p_e_2 = 0  # a_k
        # for i in range(edge_set.shape[0]):
        for i in range(len(edge_set)):
            if i == 0:
                continue
            if i == 1:
                p_e_1 = edge_set[i][0]
                p_e_2 = edge_set[i][1]
                continue

            e_1 = edge_set[i][0]  # a_{k+1}
            e_2 = edge_set[i][1]

            if p_e_1 == p_e_2:
                print("[ERROR]: cyclic loop")
                return False

            if p_e_1 > p_e_2:
                # backward edge
                if e_1 < e_2:  # ak: backward edge, ak+1:forward edge
                    if not((e_1 <= p_e_1) and (e_2 == p_e_1 + 1)):
                        print("[ERROR]: {} 行目, 条件1".format(i + 1))
                        return False

                else:  # ak: backward edge, ak+1:backward edge
                    if not((p_e_1 == e_1) and (p_e_2 < e_2)):
                        print("[ERROR]: {} 行目, 条件2".format(i + 1))
                        return False

            else:
                # forward edge
                if e_1 < e_2:
                    if not((e_1 <= p_e_2) and (e_2 == p_e_2 + 1)):  # forward edge
                        print("[ERROR]: {} 行目, 条件3".format(i + 1))
                        return False
                else:
                    if not((e_1 == p_e_2) and (e_2 < p_e_1)):  # forward edge
                        print("[ERROR]: {} 行目, 条件4".format(i + 1))
                        return False

            p_e_1 = e_1
            p_e_2 = e_2

        # print("[operate sucessfully]")
        return True

    def read(self, file_path):
        smiles_str = ""
        with open(file_path) as f:
            smiles_str = f.read()

        print(smiles_str)

        """" 立体構造, イオン結合, 謎の記号を含むものをsmilesから削除 """
        if "@" in list(smiles_str):
            raise SMILESParseError_1(
                "[Error] this processing algorithm remove Asymmetric center expression")

        if "." in list(smiles_str):
            raise SMILESParseError_2(
                "[Error] this processing algorithm remove ion")

        if ">" in list(smiles_str):
            raise SMILESParseError_3(
                "[Error] this processing algorithm remove > expression")

        self.MOL_block = Chem.MolFromSmiles(smiles_str)
        self.MOL_str = Chem.MolToMolBlock(self.MOL_block)

        if self.__isRemove():
            raise SMILESParseError_4(
                "This molcule is size over, too small, or ion.")
        print(self.MOL_str)

    def mol2adjacency(self):
        NUM_VER_EDGE_LINE = 3
        mol_list = self.MOL_str.split("\n")
        nums = re.split(" +", mol_list[NUM_VER_EDGE_LINE])
        v_num = nums[1]
        e_num = nums[2]

        VERTEX_FIRST_LINE = 4
        VERTEX_LAST_LINE = 4 + int(v_num)
        EDGE_FIRST_LINE = 4 + int(v_num)
        EDGE_LAST_LINE = EDGE_FIRST_LINE + int(e_num)
        vertex_list = mol_list[VERTEX_FIRST_LINE:VERTEX_LAST_LINE]
        edge_list = mol_list[EDGE_FIRST_LINE:EDGE_LAST_LINE]
        footer_list = mol_list[EDGE_LAST_LINE:]

        vertex_arr = np.zeros((0, 2), dtype=np.int64)
        for itr, item in enumerate(vertex_list):
            v = (re.split(" +", item))[4]
            vertex_arr = np.append(
                vertex_arr, [[int(itr + 1), int(self.ATOM2ID[v])]], axis=0)

        edge_arr = np.zeros((0, 3))
        for itr, item in enumerate(edge_list):
            li = re.split(" +", item)
            v_from = int(li[1])
            v_to = int(li[2])
            connect_num = int(li[3])
            edge_arr = np.append(
                edge_arr, [[v_from, v_to, connect_num]], axis=0)

        col = len(vertex_arr)
        A = np.zeros((col, col), dtype=np.int64)
        for itr, item in enumerate(edge_arr):
            A[int(edge_arr[itr, 0]) - 1, int(edge_arr[itr, 1]) -
              1] = int(edge_arr[itr, 2])
            A[int(edge_arr[itr, 1]) - 1, int(edge_arr[itr, 0]) -
              1] = int(edge_arr[itr, 2])

        return A, vertex_arr

    def adjacency2dfs(self, old_A, id_labels):

        V_size = len(id_labels)
        DFS_VERTEX = []  # [頂点id, label]
        DFS_EDGE = []  # [頂点id, 頂点id, edge label]
        oID2nID = {}
        oID2oID_Edge_Stack = []  # old id from, old id to
        NEW_ID_V_CTR = int(0)

        for itr in range(1, V_size):
            if old_A[0, itr] == 0:
                continue
            oID2oID_Edge_Stack.append([1, itr + 1])

        oID2nID[1] = NEW_ID_V_CTR
        DFS_VERTEX.append([0, id_labels[1]])
        NEW_ID_V_CTR += 1

        while len(oID2oID_Edge_Stack) != 0:

            #  Forward Edgeをpop
            oEdge = oID2oID_Edge_Stack.pop()
            oEdgeFrom = oEdge[0]
            oEdgeTo = oEdge[1]
            if oEdgeTo in oID2nID:
                continue

            nFrom = oID2nID[oEdgeFrom]
            oID2nID[oEdgeTo] = NEW_ID_V_CTR
            nTo = NEW_ID_V_CTR
            DFS_VERTEX.append([nTo, id_labels[oEdgeTo]])  # ラベルを探したい
            NEW_ID_V_CTR += 1
            r, c = oEdgeFrom - 1, oEdgeTo - 1
            eLabel = old_A[r, c]
            DFS_EDGE.append([nFrom, nTo, int(eLabel)])

            #  toを主体に考えて, Edgeを探す
            nBack_ward_edge = []
            raw = oEdgeTo - 1
            for col in range(0, V_size):
                if col == raw:
                    continue
                if old_A[raw, col] == 0:
                    continue
                if col + 1 in oID2nID:
                    nBack_ward_edge.append(
                        [oID2nID[oEdgeTo], oID2nID[col + 1], int(old_A[raw, col])])
                    continue

                oID2oID_Edge_Stack.append([oEdgeTo, col + 1])

            # backward edge 追加
            nBack_ward_edge.sort(key=lambda x: x[1], reverse=True)
            while len(nBack_ward_edge) != 0:
                e = nBack_ward_edge.pop()
                nFrom = e[0]
                nTo = e[1]
                nlable = e[2]

                if [nFrom, nTo, nlable] in DFS_EDGE:
                    continue

                if [nTo, nFrom, nlable] in DFS_EDGE:
                    continue

                DFS_EDGE.append([nFrom, nTo, int(e[2])])

        if not(self.__checkNeiborhoodRestriction(DFS_EDGE)):
            raise NeighborhoodRestrictionError(
                "[Error] neighborhood restiction error")

        return DFS_VERTEX, DFS_EDGE

    def graphDraw(self, out_root_dir, file_path, pic_num):

        for idx in range(pic_num):
            vertexes = {}
            forward = []
            backward = []
            edges = {}
            graph = nx.Graph()
            for line in open(file_path):
                items = line.replace('\n', '').split(' ')
                if items[0] == 'v':
                    vertexes[items[1]] = self.ID2ATOM[items[2]]
                    # vertexes[items[1]] = items[2]
                    graph.add_node(items[1], label=self.ID2ATOM[items[2]])
                    # graph.add_node(items[1],label=items[2])

                elif items[0] == 'e':
                    graph.add_edge(items[1], items[2], label=items[3])
                    if int(items[1]) < int(items[2]):
                        forward.append((items[1], items[2]))
                    else:
                        backward.append((items[1], items[2]))
                    edges[(items[1], items[2])] = items[3]

            plt.figure()
            pos = nx.spring_layout(graph)
            nx.draw_networkx_nodes(graph, pos, node_size=200, node_color="w")
            nx.draw_networkx_edges(graph, pos, forward, width=1, style='-')
            nx.draw_networkx_edges(graph, pos, backward, width=1, style='-')
            nx.draw_networkx_edge_labels(graph, pos, edges, font_size=5)

            nx.draw_networkx_labels(graph, pos, vertexes,
                                    font_size=4, font_color="r")

            plt.xticks([])
            plt.yticks([])
            plt.savefig(out_root_dir + "/graph" + "-" + str(idx) + ".pdf")

            vertexes.clear()
            del forward[:]
            del backward[:]
            edges.clear()
            graph.clear()
            drawFlag = False

    def graphSVG(self, filedir, mol):
        tm = rdMolDraw2D.PrepareMolForDrawing(mol)
        view = rdMolDraw2D.MolDraw2DSVG(500, 500)
        view.DrawMolecule(tm)
        view.FinishDrawing()
        svg = view.GetDrawingText()
        with open(filedir, "w") as f:
            f.writelines(svg)

    def dfs2Mol(self, dfs_vertex, dfs_edge, work_path):

        molBody = "\n"
        molBody += "     RDKit          2D\n"
        molBody += "\n"
        v_num = str(len(dfs_vertex))
        e_num = str(len(dfs_edge))
        molBody += " " + v_num + " " + e_num + "  0  0  0  0  0  0  0  0999 V2000\n"

        for itr in range(len(dfs_vertex)):
            atom_info = self.ID2ATOM[str(dfs_vertex[itr][1])]

            if len(atom_info) == 1:
                molBody += "    0.0000    0.0000    0.0000 " + \
                    atom_info + "   0  0  0  0  0  0  0  0  0  0  0  0\n"
            else:
                molBody += "    0.0000    0.0000    0.0000 " + \
                    atom_info + "  0  0  0  0  0  0  0  0  0  0  0  0\n"

        for itr in range(len(dfs_edge)):
            if (dfs_edge[itr][0] + 1 < 10) and (dfs_edge[itr][1] + 1 < 10):
                molBody += '  ' + str(dfs_edge[itr][0] + 1) + '  ' + str(
                    dfs_edge[itr][1] + 1) + '  ' + str(dfs_edge[itr][2]) + "  0\n"

            elif (dfs_edge[itr][0] + 1 < 10) and (dfs_edge[itr][1] + 1 >= 10):
                molBody += '  ' + str(dfs_edge[itr][0] + 1) + ' ' + str(
                    dfs_edge[itr][1] + 1) + '  ' + str(dfs_edge[itr][2]) + "  0\n"

            elif (dfs_edge[itr][0] + 1 >= 10) and (dfs_edge[itr][1] + 1 < 10):
                molBody += ' ' + str(dfs_edge[itr][0] + 1) + '  ' + str(
                    dfs_edge[itr][1] + 1) + '  ' + str(dfs_edge[itr][2]) + "  0\n"

            else:
                molBody += ' ' + str(dfs_edge[itr][0] + 1) + ' ' + str(
                    dfs_edge[itr][1] + 1) + '  ' + str(dfs_edge[itr][2]) + "  0\n"

        molBody += "M  END\n"

        return Chem.MolToMolBlock(self.__incomp2complete(work_path, molBody))

    def dfs2Molbolck(self, dfs_vertex, dfs_edge, work_path):

        molBody = "\n"
        molBody += "     RDKit          2D\n"
        molBody += "\n"
        v_num = str(len(dfs_vertex))
        e_num = str(len(dfs_edge))
        molBody += " " + v_num + " " + e_num + "  0  0  0  0  0  0  0  0999 V2000\n"

        for itr in range(len(dfs_vertex)):
            atom_info = self.ID2ATOM[str(dfs_vertex[itr][1])]

            if len(atom_info) == 1:
                molBody += "    0.0000    0.0000    0.0000 " + \
                    atom_info + "   0  0  0  0  0  0  0  0  0  0  0  0\n"
            else:
                molBody += "    0.0000    0.0000    0.0000 " + \
                    atom_info + "  0  0  0  0  0  0  0  0  0  0  0  0\n"

        for itr in range(len(dfs_edge)):
            if (dfs_edge[itr][0] + 1 < 10) and (dfs_edge[itr][1] + 1 < 10):
                molBody += '  ' + str(dfs_edge[itr][0] + 1) + '  ' + str(
                    dfs_edge[itr][1] + 1) + '  ' + str(dfs_edge[itr][2]) + "  0\n"

            elif (dfs_edge[itr][0] + 1 < 10) and (dfs_edge[itr][1] + 1 >= 10):
                molBody += '  ' + str(dfs_edge[itr][0] + 1) + ' ' + str(
                    dfs_edge[itr][1] + 1) + '  ' + str(dfs_edge[itr][2]) + "  0\n"

            elif (dfs_edge[itr][0] + 1 >= 10) and (dfs_edge[itr][1] + 1 < 10):
                molBody += ' ' + str(dfs_edge[itr][0] + 1) + '  ' + str(
                    dfs_edge[itr][1] + 1) + '  ' + str(dfs_edge[itr][2]) + "  0\n"

            else:
                molBody += ' ' + str(dfs_edge[itr][0] + 1) + ' ' + str(
                    dfs_edge[itr][1] + 1) + '  ' + str(dfs_edge[itr][2]) + "  0\n"

        molBody += "M  END\n"

        return self.__incomp2complete(work_path, molBody)

    def __incomp2complete(self, path, incomplete_mol_str):
        with open(path + "/" + "incomplete-mol.txt", mode="w") as f:
            f.writelines(incomplete_mol_str)

        return Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromMolFile(path + "/" + "incomplete-mol.txt", strictParsing=False)))
