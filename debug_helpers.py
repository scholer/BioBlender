import os
import json


def write_json_of_pdbtree(pdict):
    '''
    one might pass this function an instance of pdbIDmodelsDictionary
    write the dict as json to disc for debug
    '''
    this_root = os.path.dirname(__file__)
    destination_path = os.path.join(this_root, 'tmp', 'PDB_tree.json')
    m = json.dumps(pdict, sort_keys=True, indent=2)
    with open(destination_path, 'w') as PDB_tree:
        PDB_tree.writelines(m)
