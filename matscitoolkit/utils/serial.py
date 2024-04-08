import json
import numpy as np


def convert_np_to_list(obj):
    """Convert numpy array to list."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def custom_deserializer(dct):
    """Custom deserializer for JSON."""
    if "__ndarray__" in dct:
        shape, dtype, data = dct["__ndarray__"]
        return np.array(data, dtype=dtype).reshape(shape)
    return dct


def deserialize(input_file):
    """Deserialize JSON file."""
    with open(input_file, "r") as f:
        return json.load(f, object_hook=custom_deserializer)


def custom_serializer(obj):
    """Custom serializer for JSON."""
    if isinstance(obj, np.ndarray):
        return {"__ndarray__": [obj.shape, obj.dtype.name, obj.flatten().tolist()]}
    raise TypeError(f"Type {type(obj)} not serializable")


# def serialize(self, force_arr, output_file):
def serialize(force_arr=None, dipole_arr=None, output_file=None):
    """Serialize data to JSON."""
    data = {}
    if force_arr is not None:
        data["forces"] = force_arr
    if dipole_arr is not None:
        data["dipole"] = dipole_arr
    
    if output_file is not None:
        with open(output_file, "w") as f:
            json.dump(data, f, default=custom_serializer)
    return json.dumps(data, default=custom_serializer)
