import os
import re
import json
import sys
def grep_cp2k_basis_molopt(data):
    element = data['element']
    basis_type = data['basis_type']
    basis_file = os.path.join(data['cp2k_path'], 'data', 'BASIS_MOLOPT')
    pattern=fr'^\s*{element}\s*{basis_type}-MOLOPT-GTH.*'
    content=[]
    found=False
    with open(basis_file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if re.match(pattern, line):
            found=True
            content.append(line)
            continue
        if found:
            if not re.match(r'^\s*\d', line):
                break
            print("line: ", line)
            content.append(line)
    writefile=f"{element}-{basis_type}-MOLOPOT-GTH"
    with open(writefile, 'w') as f:
        f.writelines(content)
    return writefile

if __name__ == '__main__':
    file_path = sys.argv[1]
    with open(file_path, 'r') as file:
        data = json.load(file)
    writefile = grep_cp2k_basis_molopt(data)
    orb_submodule_pth=f"{os.path.dirname(os.path.abspath(__file__))}/deps/gaussian_orbital_for_ABACUS/cp2k2abacus_gaussian.py"
    os.system(f"python {orb_submodule_pth} {writefile} {writefile} {data['basis_rcut']} MOLOPT-GTH")