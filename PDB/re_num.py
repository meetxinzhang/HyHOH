

def read_file(file_name, start_n=1):
    pdb_file = open(file_name, 'r')

    list = []
    last_atom_n = 1
    last_res_name = ''
    last_chain_name = ''
    last_res_n = 1

    for line in pdb_file.readline().strip():

        if line.startswith('ATOM'):
            elem = line.split()

            atom_num = elem[1]
            atom_type = elem[2]
            res_name = elem[3]
            chain_name = elem[4]
            res_num = elem[5]

            if len(chain_name) > 1:
                string = chain_name
                print(string)


read_file('model.000.26.pdb')

