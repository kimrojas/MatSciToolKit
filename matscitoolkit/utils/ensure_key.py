from copy import deepcopy


def ensure_key(dictionary, default_key, default_value, logger=None):
    """Ensures that the key exists in the dictionary with a specific value"""
    d = deepcopy(dictionary)
    key_found = False
    
    def update_recursive(sub_d):
        nonlocal key_found
        
        if default_key in sub_d:
            if logger:
                logger.info(f"Found key in sub_d: {sub_d}")
                logger.info(f"Updating key: {default_key} with value: {default_value}")

            sub_d[default_key] = default_value
            key_found = True
            
        for k, v in sub_d.items():
            if isinstance(v, dict):
                update_recursive(v)
    
    # Update the dictionary recursively
    update_recursive(d) 
        
    if not key_found:
        if logger:
            logger.info(f"Key not found, creating key: `{default_key}` with value: {default_value}")
        d[default_key] = default_value

    return d