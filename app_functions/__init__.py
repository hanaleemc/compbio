import random
from tempfile import gettempdir, NamedTemporaryFile
import numpy as np
from flask import session
import biotite.sequence as seq
import biotite.sequence.io.fasta as fasta
import biotite.sequence.io.genbank as gb
import biotite.database.entrez as entrez
import matplotlib.pyplot as plt
import biotite.sequence.graphics as graphics
plt.switch_backend('Agg')
import os

__version__ = "1.0.1" # specify a version here in the backend

def random_id(length): # This function has two doctests inside it to give an example on how to make doctests
    """
    Creates a random configuration key for the session - for safety of session variables.

    >>> len(random_id(50)) == 50
    True

    >>> random_id('hello')
    Traceback (most recent call last):
        ...
    TypeError: The input must be a positive integer.

    """
    if type(length) != int or length < 1:
        raise TypeError('The input must be a positive integer.')

    choices = '0123456789abcdefghijklmnopqrstuvwxyz'

    id = ''
    for _ in range(length):
        id += random.choice(choices)
    return id

def make_feature_maps(gene):

    try:
        find_id = entrez.fetch( gene , gettempdir(), suffix="gb", db_name= "nuccore",ret_type= "gb" )
        read_file = gb.GenBankFile.read(find_id)
        file_annotation = gb.get_annotation(read_file)
    except:
        flash('The entered gene could not found. Please try again.', 'error')
        return None

    key_list=[]

    for feature in file_annotation:
        keys=feature.key
        key_list.append(keys)
        if feature.key == "source":
            # loc_range has exclusive stop
            loc = list(feature.locs)[0]
            loc_range = (loc.first, loc.last+1)
            Unique_key= np.unique(key_list)

    pwd = os.getcwd()

    Unique_key= np.unique(key_list)
    for j in range(len(Unique_key)):
        i = Unique_key[j]

        fig, ax = plt.subplots(figsize=(8.0, 2.0))
        graphics.plot_feature_map(ax,seq.Annotation( [feature for feature in file_annotation if feature.key == i]),
                                                multi_line=False, loc_range=loc_range,show_line_position=True)

        plt.title('This plot is for {} features'.format(i))
        plt.savefig(pwd + '/app/static/images/{}.png'.format(i), dpi=300)
        session['valid_gene'] = True

    return None
