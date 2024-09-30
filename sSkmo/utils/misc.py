def intersection(lst1, lst2):
    '''returns intersection between two lists.'''
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def filter_dict(input_dict):
    return {k: v for k, v in input_dict.items() if isinstance(v, list) and len(v) > 1}