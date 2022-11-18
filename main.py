from loading import load_directory
from kmers import stream_kmers, kmer2str


def similarity(A, inter, B):
    # --- To complete ---
    pass


def jaccard(A, inter, B):
    # --- To complete ---
    pass



if __name__ == "__main__":
    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    
    filenames = list(files.keys())
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            
            # --- Complete here ---

            A, inter, B = my_method(indexes[filenames[i]], files[filenames[j]], k)
            print(filenames[i], filenames[j], jaccard(A, inter, B), similarity(A, inter, B))
