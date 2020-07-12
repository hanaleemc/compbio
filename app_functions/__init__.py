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
import biotite.sequence as seq
import matplotlib.pyplot as plt
def RLFP (gene):
    #The cutting sequances. Notice that this is only two restriction enzymes
    #for forward and backward translation 
    gene=''.join(gene.split())#removing the spaces
    gene= seq.NucleotideSequence(gene)
    cut_TaqI= seq.NucleotideSequence("CGA")   
    cut_TaqI_rev= seq.NucleotideSequence("AGCT")
    cut_HpaII= seq.NucleotideSequence("CGG")
    cut_HpaII_rev= seq.NucleotideSequence("GGCT")
    #finding the indexes
    find_TaqI= list(seq.find_subsequence(gene, cut_TaqI)) 
    find_TaqI_rev= list(seq.find_subsequence(gene, cut_TaqI_rev))
    find_HpaII= list(seq.find_subsequence(gene, cut_HpaII))
    find_HpaII_rev= list(seq.find_subsequence(gene, cut_HpaII_rev))
    #lenghts of the cuts 
    #for Taqi
    passed_cut_TaqI=[] #list of the indexs that have the cut before it 
    Taqi_length=[] #final length
    for i in find_TaqI:
      if gene[i-1] == "T":
          passed_cut_TaqI.append(i)
    for i in passed_cut_TaqI:
      if i != passed_cut_TaqI[0]:
          if i != passed_cut_TaqI[-1]:
              Taqi_length.append(i- passed_cut_TaqI[passed_cut_TaqI.index(i)-1])
          else:
              Taqi_length.append(len(gene)-i)
      else: 
          Taqi_length.append(i)
    #------------------------------------------------------------------------------------------------
    #for Taqi reverse
    passed_cut_TaqI_rev=[] #list of the indexs that have the cut before it 
    TaqI_rev_length=[] #final length
    for i in find_TaqI_rev:
      if i != find_TaqI_rev[0]:
          if i != find_TaqI_rev[-1]:
              TaqI_rev_length.append(i- find_TaqI_rev[find_TaqI_rev.index(i)-1])
          else:
              TaqI_rev_length.append(len(gene)-i)
      else: 
          TaqI_rev_length.append(i)
    #------------------------------------------------------------------------------------------------
    #for HpaII 
    passed_cut_HpaII=[] #list of the indexs that have the cut before it 
    HpaII_length=[] #final length
    for i in find_HpaII:
      if gene[i-1] == "C":
          passed_cut_HpaII.append(i)
    for i in passed_cut_HpaII:
      if i != passed_cut_HpaII[0]:
          if i != passed_cut_HpaII[-1]:
              HpaII_length.append(i- passed_cut_HpaII[passed_cut_HpaII.index(i)-1])
          else:
              HpaII_length.append(len(gene)-i)
      else: 
          HpaII_length.append(i)
    #------------------------------------------------------------------------------------------------
    #for HpaII_rev
    passed_cut_HpaII_rev=[] #list of the indexs that hae the cut before it 
    HpaII_rev_length=[] #final length
    for i in find_HpaII_rev:
      if i != find_HpaII_rev[0]:
          if i != find_HpaII_rev[-1]:
              HpaII_rev_length.append(i- find_HpaII_rev[find_HpaII_rev.index(i)-1])
          else:
              HpaII_rev_length.append(len(gene)-i)
      else: 
          HpaII_rev_length.append(i)
    #---------------------------------------------------------------------------------------------
    #building the histograms
    pwd = os.getcwd()
    plt.hist(Taqi_length,  bins = 500, facecolor='blue', alpha=0.5)
    plt.ylabel("Number of strands")
    plt.xlabel('length of strands')
    plt.savefig(pwd + '/static/images/RLFP_Taqi_length.png', dpi=200)
    plt.hist(TaqI_rev_length, bins = 500, facecolor='blue', alpha=0.5)
    plt.savefig(pwd + '/static/images/RLFP_TaqI_rev_length.png', dpi=200)
    plt.hist(HpaII_length,  bins = 500, facecolor='blue', alpha=0.5)
    plt.savefig(pwd + '/static/images/RLFP_HpaII_length.png', dpi=200)
    plt.hist(HpaII_rev_length,  bins = 500, facecolor='blue', alpha=0.5)
    plt.savefig(pwd + '/static/images/HpaII_rev_length.png', dpi=200)
    session['valid_seq'] = True
