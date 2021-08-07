# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@file: pdb_parser.py
@time: 6/6/21 4:08 PM
@desc:
   fork from https://github.com/manasa711/Protein-Binding-Site-Visualization
   addition features to save new pdb
"""
import platform
import os
import numpy as np
from collections import defaultdict as ddict
# from utils.exception_message import ExceptionPassing

# from utils.log_output import Logger
# logger = Logger(log_path='//output/logs/select_residue.log', is_print=False)


def _continuity_check(res_list):
    """
    :param res_list: like ['2', '3', '3A', '4', '7', ...]
    :return: fix missing but remains special elements like 3A
            -> ['2', '3', '3A', '4', '5', '6', '7']
    """
    start = res_list[0]
    end = res_list[-1]

    if not start.isdigit():
        # start = re.findall(r'\d+', start)                  # op 1
        start = ''.join(e for e in start if e.isdigit())     # op 2 more pythonic
    if not end.isdigit():
        # end = re.findall(r'\d+', end)
        end = ''.join(e for e in end if e.isdigit())
    idx_list = range(int(start), int(end) + 1)  # [2, 3, 4, 5, 6, 7]
    point_res = 0
    point_idx = 0

    final_res = []
    while True:
        try:
            idx = str(idx_list[point_idx])
            res = res_list[point_res]
        except IndexError:
            break

        if res == idx:
            final_res.append(res)
            point_idx += 1
            point_res += 1
        elif not res.isdigit():  # or res not in idx_list
            final_res.append(res)
            point_res += 1
        else:
            final_res.append(idx)
            point_idx += 1
    return final_res


def _search_bind_sites(pdb_file, bind_radius, chain1, chain2):
    file = open(pdb_file, 'r')

    # creating lists with the coordinates of CA atoms from both the chains
    O_atoms_ch1 = []
    O_atoms_ch2 = []

    for line in file.readlines():
        if line.startswith('HETATM'):
            # elem = line.split()
            #
            # if elem[2] == 'CA':
            #     is_ch_id = elem[4]
            #     if is_ch_id.isalpha():
            #         chain_id = is_ch_id
            #         aa_idx = elem[5]
            #         x = elem[6]
            #         y = elem[7]
            #         z = elem[8]
            #     else:  # A1099A -> A 1099A
            #         chain_id = is_ch_id[0]
            #         aa_idx = is_ch_id[1:]
            #         x = elem[5]
            #         y = elem[6]
            #         z = elem[7]

            res_name = line[17:20].strip()
            chain_id = line[21].strip()
            res_seq = line[22:26].strip()

            x = line[30:38].strip()
            y = line[38:46].strip()
            z = line[46:54].strip()

            if chain_id == chain1 and res_name == 'HOH':
                O_atoms_ch1.append([res_seq, x, y, z])
                # print(line)
            elif chain_id == chain2 and res_name == 'HOH':
                O_atoms_ch2.append([res_seq, x, y, z])
                # print(line)

    file.close()
    # calculating Euclidean Distance between CA atoms of chain 1 and CA atoms of chain2

    bind_sites_1 = []  # list with interface atoms from chain 1
    bind_sites_2 = []  # list with interface atoms from chain 2

    for CA_i in O_atoms_ch1:
        for CA_j in O_atoms_ch2:
            x1 = float(CA_i[1])
            y1 = float(CA_i[2])
            z1 = float(CA_i[3])
            x2 = float(CA_j[1])
            y2 = float(CA_j[2])
            z2 = float(CA_j[3])
            e = ((x1 - x2) ** 2) + ((y1 - y2) ** 2) + ((z1 - z2) ** 2)
            euc = e ** 0.5
            if euc <= bind_radius:
                if CA_i[0] not in bind_sites_1:
                    bind_sites_1.append(CA_i[0])

    # We have to do this loop twice, to maintain the original order of amino acid
    for CA_j in O_atoms_ch2:
        for CA_i in O_atoms_ch1:
            x1 = float(CA_i[1])
            y1 = float(CA_i[2])
            z1 = float(CA_i[3])
            x2 = float(CA_j[1])
            y2 = float(CA_j[2])
            z2 = float(CA_j[3])
            e = ((x1 - x2) ** 2) + ((y1 - y2) ** 2) + ((z1 - z2) ** 2)
            euc = e ** 0.5
            if euc <= bind_radius:
                if CA_j[0] not in bind_sites_2:
                    bind_sites_2.append(CA_j[0])

    if len(bind_sites_1) == 0 or len(bind_sites_2) == 0:
        # raise ExceptionPassing('  can not find binding sites, maybe due to error chain selection: \n',
        #                        '   - bind_sites: ', bind_sites_1, bind_sites_2, pdb_file)
        pass
    # bind_sites_1 = _continuity_check(bind_sites_1)
    # bind_sites_2 = _continuity_check(bind_sites_2)
    return bind_sites_1, bind_sites_2


def save_bind_waters(pdb_file, save_dir, rec_chains, lig_chains, bind_radius):
    if platform.system() == 'Windows':
        pdb_id = pdb_file.split('\\')[-1].replace('.pdb', '')
    else:
        pdb_id = pdb_file.split('/')[-1].replace('.pdb', '')

    # get residues index of binding site
    final_res = ddict(set)
    for cr in rec_chains:
        for cl in lig_chains:
            try:
                bs1, bs2 = _search_bind_sites(pdb_file, bind_radius, cr, cl)
            except Exception as e:
                # logger.write(e.message, join_time=False)
                # logger.flush()
                break
            for res_seq in bs1:
                final_res[cr].add(res_seq)
            for res_seq in bs2:
                final_res[cl].add(res_seq)

    print('hydration water:')
    for k, v in final_res.items():
        print(k, v)

    # save pdb
    chain_ids = np.concatenate([rec_chains, lig_chains], axis=0)

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    new_file = open(save_dir + pdb_id + '_bind_site.pdb', 'a')

    for line in open(pdb_file, 'r'):
        if line.startswith('ATOM'):
            new_file.write(line)

        if line.startswith('HETATM'):
            chain_id = line[21].strip()
            res_seq = line[22:26].strip()
            if chain_id in chain_ids and res_seq in list(final_res[chain_id]):
                new_file.write(line)

    new_file.close()

    # logger.write(pdb_id, '_bind_sites: ', final_res.values(), join_time=False)
    # logger.flush()


if __name__ == "__main__":
    save_bind_waters(pdb_file='/media/xin/WinData/ACS/gmx/interaction/ding/7KGJ_2/renum.pdb',
                     save_dir='/media/xin/WinData/ACS/gmx/interaction/ding/7KGJ_2/',
                     bind_radius=5, rec_chains=['B'], lig_chains=['A'])
