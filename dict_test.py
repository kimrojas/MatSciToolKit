from pprint import pprint

# def ensure_key(dictionary, key, value):
    
#     def replace_key()

# def ensure_key(dictionary, key, value):
#     key_updated = False
    
#     for k, v in dictionary.items():        
#         if isinstance(v, dict):
#             v, key_updated = ensure_key(v, key, value)
        
#         if k == key and v != value:
#             v = value
#             key_updated = True
        
#         if key_updated:
#             break
#     return dictionary, key_updated
            
    



# def ensure_tprnfor(input_data):
#     updated = False
#     for key, value in input_data.items():
#         if isinstance(value, dict):
#             if "tprnfor" in value:
#                 if value["tprnfor"] is True:
#                     print("tprnfor is already set to True.")
#                 else:
#                     value["tprnfor"] = True
#                     print("tprnfor was set to True.")
#                     updated=True
#             else:
#                 value["tprnfor"] = True
#                 print("tprnfor added and set to True.")
                
#             # break  # Assuming only one nested dictionary needs to be checked
#     return input_data

# Test the function
# input_data_template = {
#     "control": {"tprnfor": True},
#     "pseudo_dir": "pseudo_dir",
#     "occupations": "smearing",
#     "smearing": "fermi-dirac",
#     "degauss": 0.02,
#     "ecutwfc": 20,
#     "ecutrho": 160,
# }

input_data_template = {
    "control": {"layer2": {"force":False}},
    "pseudo_dir": "pseudo_dir",
    "occupations": "smearing",
    "smearing": "fermi-dirac",
    "degauss": 0.02,
    "ecutwfc": 20,
    "ecutrho": 160,
}

pprint(input_data_template)
# updated_data = ensure_key(input_data_template, "force", True)

# pprint(updated_data)

# def update_force_to_true(d):
#     for key, value in d.items():
#         if key == "force":
#             d[key] = True
#         elif isinstance(value, dict):
#             update_force_to_true(value)

def update_force_to_true(d):
    # A flag to check if 'force' key is found
    found = {'exists': False}

    # Inner recursive function to update the value
    def update_recursive(sub_d):
        if "force" in sub_d:
            sub_d["force"] = True
            found['exists'] = True
        for key, value in sub_d.items():
            if isinstance(value, dict):
                update_recursive(value)
    
    # Update the dictionary recursively
    update_recursive(d)

    # If 'force' was not found, add it to the top-most level
    if not found['exists']:
        d["force"] = True
            
update_force_to_true(input_data_template)
pprint(input_data_template)