
import sys


class Atom:
    def __init__(self, vmd_index, chain, resname, resid, name, element):
        """
        Construct a new atom with the values, check the element (potentially infer from name) and set van der Waals
        radius based on it.

        Parameters
        ----------
        vmd_index: int
        chain: str
        resname: str
        resid: int
        name: str
        element: str
        """
        self.index = vmd_index
        self.chain = chain
        self.resid = resid
        self.resname = resname
        self.name = name
        self.element = element.upper()

        if element == "X":
            self.element = infer_element(resname, name)

        if self.element in vdw_radii:
            self.vdwradius = vdw_radii[self.element]
        else:
            sys.stderr.write("WARNING: Doesn't know van der Waals radius of "+self.element+". Defaulting to 1.7")
            self.vdwradius = 1.7

    def get_label(self):
        """
        Return a colon-separated representation of this atom: "chain:resname:resid:name:index" e.g. "A:ASP:114:CA:11205"

        Returns
        -------
        string
        """
        return "%s:%s:%d:%s:%d" % (self.chain, self.resname, self.resid, self.name, self.index)


def infer_element(resname, name):
    """
    Infer the element of an atom based on atom name. If no element can be inferred print a warning and return '?'.
    If there's doubt print a warning and resolve by returning the longer element in case the residue name is similar
    (e.g. "CLA/CLA" returns "CL") or the shortest element name otherwise (e.g. "CA" returns "C"). This is based on the
    observation that typically ions (that have frequent naming conflicts) have `name==resname`.

    Parameters
    ----------
    resname: string
        Residue name as specified in topology
    name: string
        Atom name as specified in topology

    Returns
    -------
    string
        Inferred element identity from the atom name
    """
    import re

    name_letters = re.search(r'[a-zA-Z]+', name).group(0)  # Only retain first word: "1HH2" -> "HH", "1H2S" -> "H"
    name_letters = name_letters.upper()
    if len(name_letters) == 0:
        sys.stderr.write("WARNING: Element can't be determined for atom '" + name + "'")
        return "?"

    # Single letter atom_name identifies element
    if len(name_letters) == 1:
        if name_letters in element_names:
            return name_letters
        else:
            sys.stderr.write("WARNING: Element can't be determined for atom '" + name + "'")
            return "?"

    elem_1ch = name_letters[0]
    elem_2ch = name_letters[0:2]
    if elem_1ch in element_names and elem_2ch in element_names:
        ret = elem_2ch if resname == name else elem_1ch
        # sys.stderr.write("WARNING: Ambiguous element of atom '" + name + "'. Assuming '" + ret + "'\n")
        return ret

    if elem_1ch in element_names:
        return elem_1ch

    if elem_2ch in element_names:
        return elem_2ch

    sys.stderr.write("WARNING: Element can't be determined for atom '" + name + "'")
    return "?"


# https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
# https://www.webelements.com/gadolinium/
vdw_radii = {
    'H': 1.20, 'HE': 1.40, 'LI': 1.82, 'BE': 1.53, 'B': 1.92, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
    'NE': 1.54, 'NA': 2.27, 'MG': 1.73, 'AL': 1.84, 'SI': 2.10, 'P': 1.80, 'S': 1.80, 'CL': 1.75, 'AR': 1.88,
    'K': 2.75, 'CA': 2.31, 'SC': 2.11, 'NI': 1.63, 'CU': 1.40, 'ZN': 1.39, 'GA': 1.87, 'GE': 2.11, 'AS': 1.85,
    'SE': 1.90, 'BR': 1.85, 'KR': 2.02, 'RB': 3.03, 'SR': 2.49, 'PD': 1.63, 'AG': 1.72, 'CD': 1.58, 'IN': 1.93,
    'SN': 2.17, 'SB': 2.06, 'TE': 2.06, 'I': 1.98, 'XE': 2.16, 'CS': 3.43, 'BA': 2.68, 'PT': 1.75, 'AU': 1.66,
    'HG': 1.55, 'TL': 1.96, 'PB': 2.02, 'BI': 2.07, 'PO': 1.97, 'AT': 2.02, 'RN': 2.20, 'FR': 3.48, 'RA': 2.83,
    'U': 1.86, 'GD': 2.79}


# https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
element_names = {
    'H': 'hydrogen', 'HE': 'helium', 'LI': 'lithium', 'BE': 'beryllium', 'B': 'boron',
    'C': 'carbon', 'N': 'nitrogen', 'O': 'oxygen', 'F': 'fluorine', 'NE': 'neon',
    'NA': 'sodium', 'MG': 'magnesium', 'AL': 'aluminium', 'SI': 'silicon', 'P': 'phosphorus',
    'S': 'sulfur', 'CL': 'chlorine', 'AR': 'argon', 'K': 'potassium', 'CA': 'calcium',
    'SC': 'scandium', 'TI': 'titanium', 'V': 'vanadium', 'CR': 'chromium', 'MN': 'manganese',
    'FE': 'iron', 'CO': 'cobalt', 'NI': 'nickel', 'CU': 'copper', 'ZN': 'zinc',
    'GA': 'gallium', 'GE': 'germanium', 'AS': 'arsenic', 'SE': 'selenium', 'BR': 'bromine',
    'KR': 'krypton', 'RB': 'rubidium', 'SR': 'strontium', 'Y': 'yttrium', 'ZR': 'zirconium',
    'NB': 'niobium', 'MO': 'molybdenum', 'TC': 'technetium', 'RU': 'ruthenium', 'RH': 'rhodium',
    'PD': 'palladium', 'AG': 'silver', 'CD': 'cadmium', 'IN': 'indium', 'SN': 'tin',
    'SB': 'antimony', 'TE': 'tellurium', 'I': 'iodine', 'XE': 'xenon', 'CS': 'caesium',
    'BA': 'barium', 'LA': 'lanthanum', 'CE': 'cerium', 'PR': 'praseodymium', 'ND': 'neodymium',
    'PM': 'promethium', 'SM': 'samarium', 'EU': 'europium', 'GD': 'gadolinium', 'TB': 'terbium',
    'DY': 'dysprosium', 'HO': 'holmium', 'ER': 'erbium', 'TM': 'thulium', 'YB': 'ytterbium',
    'LU': 'lutetium', 'HF': 'hafnium', 'TA': 'tantalum', 'W': 'tungsten', 'RE': 'rhenium',
    'OS': 'osmium', 'IR': 'iridium', 'PT': 'platinum', 'AU': 'gold', 'HG': 'mercury',
    'TL': 'thallium', 'PB': 'lead', 'BI': 'bismuth', 'PO': 'polonium', 'AT': 'astatine',
    'RN': 'radon', 'FR': 'francium', 'RA': 'radium', 'AC': 'actinium', 'TH': 'thorium',
    'PA': 'protactinium', 'U': 'uranium', 'NP': 'neptunium', 'PU': 'plutonium', 'AM': 'americium',
    'CM': 'curium', 'BK': 'berkelium', 'CF': 'californium', 'ES': 'einsteinium', 'FM': 'fermium',
    'MD': 'mendelevium', 'NO': 'nobelium', 'LR': 'lawrencium', 'RF': 'rutherfordium', 'DB': 'dubnium',
    'SG': 'seaborgium', 'BH': 'bohrium', 'HS': 'hassium', 'MT': 'meitnerium', 'DS': 'darmstadtium',
    'RG': 'roentgenium', 'CN': 'copernicium', 'NH': 'nihonium', 'FL': 'flerovium', 'MC': 'moscovium',
    'LV': 'livermorium', 'TS': 'tennessine', 'OG': 'oganesson'}
