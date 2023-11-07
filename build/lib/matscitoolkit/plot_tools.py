def tick_combiner(tick_locations, tick_labels):
    """
    Combines tick labels if the tick labels are located in the same xcoord
    
    tick_locations = [0., 1.7321, 2.9568, 3.6639, 4.1639, 5.6639, 5.6639, 6.1639]
    tick_labels = ["G", "H", "N", "G", "P", "H", "P", "N"]
    to 
    [0.0, 1.7321, 2.9568, 3.6639, 4.1639, 5.6639, 6.1639],
    ['G', 'H', 'N', 'G', 'P', 'H | P', 'N']
    
    """
    
    
    processed_locations = []
    processed_labels = []

    for location, label in zip(tick_locations, tick_labels):
        if location in processed_locations:
            # merge the label with the existing one at that location
            idx = processed_locations.index(location)
            processed_labels[idx] += "|" + label
        else:
            processed_locations.append(location)
            processed_labels.append(label)

    return (processed_locations, processed_labels)