import re
electron_configurations = {
    'H': ['1s1'],
    'He': ['1s2'],
    'Li': ['He', '2s1'],
    'Be': ['He', '2s2'],
    'B': ['He', '2s2', '2p1'],
    'C': ['He', '2s2', '2p2'],
    'N': ['He', '2s2', '2p3'],
    'O': ['He', '2s2', '2p4'],
    'F': ['He', '2s2', '2p5'],
    'Ne': ['He', '2s2', '2p6'],
    'Na': ['Ne', '3s1'],
    'Mg': ['Ne', '3s2'],
    'Al': ['Ne', '3s2', '3p1'],
    'Si': ['Ne', '3s2', '3p2'],
    'P': ['Ne', '3s2', '3p3'],
    'S': ['Ne', '3s2', '3p4'],
    'Cl': ['Ne', '3s2', '3p5'],
    'Ar': ['Ne', '3s2', '3p6'],
    'K': ['Ar', '4s1'],
    'Ca': ['Ar', '4s2'],
    'Sc': ['Ar', '4s2', '3d1'],
    'Ti': ['Ar', '4s2', '3d2'],
    'V': ['Ar', '4s2', '3d3'],
    'Cr': ['Ar', '4s1', '3d5'],
    'Mn': ['Ar', '4s2', '3d5'],
    'Fe': ['Ar', '4s2', '3d6'],
    'Co': ['Ar', '4s2', '3d7'],
    'Ni': ['Ar', '4s2', '3d8'],
    'Cu': ['Ar', '4s1', '3d10'],
    'Zn': ['Ar', '4s2', '3d10'],
    'Ga': ['Ar', '4s2', '3d10', '4p1'],
    'Ge': ['Ar', '4s2', '3d10', '4p2'],
    'As': ['Ar', '4s2', '3d10', '4p3'],
    'Se': ['Ar', '4s2', '3d10', '4p4'],
    'Br': ['Ar', '4s2', '3d10', '4p5'],
    'Kr': ['Ar', '4s2', '3d10', '4p6'],
    'Rb': ['Kr', '5s1'],
    'Sr': ['Kr', '5s2'],
    'Y': ['Kr', '5s2', '4d1'],
    'Zr': ['Kr', '5s2', '4d2'],
    'Nb': ['Kr', '5s1', '4d4'],
    'Mo': ['Kr', '5s1', '4d5'],
    'Tc': ['Kr', '5s1', '4d5'],
    'Ru': ['Kr', '5s1', '4d7'],
    'Rh': ['Kr', '5s1', '4d8'],
    'Pd': ['Kr', '5s0', '4d10'],
    'Ag': ['Kr', '5s1', '4d10'],
    'Cd': ['Kr', '5s2', '4d10'],
    'In': ['Kr', '5s2', '4d10', '5p1'],
    'Sn': ['Kr', '5s2', '4d10', '5p2'],
    'Sb': ['Kr', '5s2', '4d10', '5p3'],
    'Te': ['Kr', '5s2', '4d10', '5p4'],
    'I': ['Kr', '5s2', '4d10', '5p5'],
    'Xe': ['Kr', '5s2', '4d10', '5p6'],
    'Cs': ['Xe', '6s1'],
    'Ba': ['Xe', '6s2'],
    'La': ['Xe', '6s2', '5d1'],
    'Ce': ['Xe', '6s2', '4f1', '5d1'],
    'Pr': ['Xe', '6s2', '4f3', '5d1'],
    'Nd': ['Xe', '6s2', '4f4', '5d1'],
    'Pm': ['Xe', '6s2', '4f5', '5d1'],
    'Sm': ['Xe', '6s2', '4f6', '5d1'],
    'Eu': ['Xe', '6s2', '4f7'],
    'Gd': ['Xe', '6s2', '4f7', '5d1'],
    'Tb': ['Xe', '6s2', '4f9', '5d1'],
    'Dy': ['Xe', '6s2', '4f10', '5d1'],
    'Ho': ['Xe', '6s2', '4f11', '5d1'],
    'Er': ['Xe', '6s2', '4f12', '5d1'],
    'Tm': ['Xe', '6s2', '4f13', '5d1'],
    'Yb': ['Xe', '6s2', '4f14'],
    'Lu': ['Xe', '6s2', '4f14', '5d1'],
    'Hf': ['Xe', '6s2', '4f14', '5d2'],
    'Ta': ['Xe', '6s2', '4f14', '5d3'],
    'W': ['Xe', '6s2', '4f14', '5d4'],
    'Re': ['Xe', '6s2', '4f14', '5d5'],
    'Os': ['Xe', '6s2', '4f14', '5d6'],
    'Ir': ['Xe', '6s2', '4f14', '5d7'],
    'Pt': ['Xe', '6s1', '4f14', '5d9'],
    'Au': ['Xe', '6s1', '4f14', '5d10'],
    'Hg': ['Xe', '6s2', '4f14', '5d10'],
    'Tl': ['Xe', '6s2', '4f14', '5d10', '6p1'],
    'Pb': ['Xe', '6s2', '4f14', '5d10', '6p2'],
    'Bi': ['Xe', '6s2', '4f14', '5d10', '6p3'],
    'Po': ['Xe', '6s2', '4f14', '5d10', '6p4'],
    'At': ['Xe', '6s2', '4f14', '5d10', '6p5'],
    'Rn': ['Xe', '6s2', '4f14', '5d10', '6p6'],
    'Fr': ['Rn', '7s1'],
    'Ra': ['Rn', '7s2'],
    # Add more elements as needed
}

def count_valence_num(lst):
    valence_num = 0
    for item in lst:
        if re.match(r'\d+', item):
            valence_num += int(re.search(r'\d+$', item).group())  #`\d+$` match the end number, and .group() extract the matched number
    return valence_num

def get_core_valence(element, nval):
    core_init = electron_configurations[element][0]
    vals_init =electron_configurations[element][1:]
    nval_init = count_valence_num(vals_init)
    if element == 'H' or element == 'He':
        return ["none"], electron_configurations[element]
    elif nval_init == nval:   #treat the core as another element to analyze
        return [core_init], vals_init
    elif nval_init < nval:
        core, val_inside = get_core_valence(core_init, nval-nval_init)
        return core, val_inside + vals_init
    elif nval_init > nval:   #try other possible valence numbers from the end of the list
        nval_try = 0
        core_left = electron_configurations[element]
        val=[]
        for str_try in reversed(vals_init):
            nval_try += count_valence_num([str_try])
            core_left.pop()
            val.insert(0, str_try)
            if(nval_try == nval):
                return core_left, val
    raise ValueError(f"Could not find core and valence for {element} with valence {nval}")