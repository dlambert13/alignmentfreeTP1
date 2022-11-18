
def kmer2str(val, k):
    """ Transform a kmer integer into its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for i in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)

def encode(nucleotide):
    """ Transform the string representation of a nucleotide into its
    integer representation
    :param str nucleotide: a string representation of a nucleotide
    :return 
    """
    char_dict = {'A': 0, 'C': 1, 'T': 2, 'G': 3}
    if nucleotide in char_dict.keys():
        return char_dict[nucleotide]
    else: # if nucleotide not in dictionary, consider it as equal to A
        return 0

def comp_encode(nucleotide):
    """ Transform the string representation of a nucleotide into the integer
    representation of its complementary nucleotide
    :param str nucleotide: a string representation of a nucleotide
    :return 
    """
    char_dict = {'A': 2, 'C': 3, 'T': 0, 'G': 1}
    if nucleotide in char_dict.keys():
        return char_dict[nucleotide]
    else: # if nucleotide not in dictionary, consider it as equal to A
        return 2

def stream_kmers(text, k):
    """ Transforms a text into a list of k-grams"""
    kmer_list = []
    rkmer_list = []
    mask = (1 << (k - 1) * 2) - 1
    rmask = mask
    kmer = 0
    rkmer = 0

    # encoding the first k-1 characters
    for i in range(k - 1):
        kmer = kmer << 2 # make room to the right to store the next item
        kmer += encode(text[i]) # encode then add the next item
        rkmer = rkmer << 2
        rkmer += comp_encode(text[len(text) - 1 - i])

    print("end of first loop: ", bin(kmer), bin(rkmer))
    
    # encoding 
    for nucleotide in text[(k - 1):]:
        kmer = kmer & mask
        kmer = kmer << 2
        kmer = kmer + encode(nucleotide)
        kmer_list.append(kmer)
        #r: mask operation not stricly necessary since we right bit shift
        rkmer = rkmer & rmask
        rkmer = rkmer >> 2
        rkmer = rkmer + comp_encode(nucleotide) << (k - 1) * 2
        rkmer_list.append(rkmer)

    return kmer_list, rkmer_list

if __name__ == '__main__':
    print(
        "A: ", encode("A"), ", ",
        "C: ", encode("C"), ", ",
        "T: ", encode("T"), ", ",
        "G: ", encode("G")
    )
    test_sequence = "ATTACGT"
    #test_sequence = "ATT"
    test_k = 3
    print("test sequence:", test_sequence)
    kmer_ints, rkmer_ints = stream_kmers(test_sequence, test_k)
    print(kmer_ints)
    for i in kmer_ints:
        print("kmer: ")
        print(i)
        print(kmer2str(i, test_k))
    for i in rkmer_ints:
        print("rkmer: ")
        print(i)
        print(kmer2str(i, test_k))
    #for r_i in rkmer_ints:
    #    print(kmer2str(r_i, test_k))
