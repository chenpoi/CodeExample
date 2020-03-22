import re
import subprocess
import time
import os
from multiprocessing.pool import ThreadPool as Pool
from functools import partial
import sys

###
# high_volume_MSA.py is to continue the work started by prefix_align.py. it read temporary files from prefix_align.py,
# which is separated fasta format files for each cluster defined by CD-HIT. Then clustalo command will get called and
# undergo Multiple Sequence Alignment(MSA). The final output is consensus sequences for each cluster. Multi-thread mode
# can be activated if a specific number has been taken as command line parameters


log_file = open("log.txt", "w")

# part 1: parameter interpreter
# ------------------------------------------------------------------------------------------------------------------

def parameter_interpreter():
    # check if multi-thread mode is called
    if len(sys.argv) > 2:
        return sys.argv[1], sys.argv[2]
    elif len(sys.argv) == 2:
        return sys.argv[1], 0
    else:
        return 0, 0

# part 2: clustal omega alignment
# ------------------------------------------------------------------------------------------------------------------


def clustalo_align(input_seq, output_file):
    # use subprocess module to do clustal omega
    p = subprocess.Popen(["/data/seed/Programs/clustal-omega-1.2.1/clustalo-1.2.0-Ubuntu-x86_64",
                         "--infile=-",
                         "--seqtype=DNA",
                         "--outfile=" + output_file,
                         "--outfmt=msf",
                         "--force",
                         "--use-kimura"],
                         stdin=subprocess.PIPE)
    p.stdin.write(input_seq.encode())
    p.stdin.close()
    p.communicate()
    # it's better that if you can find a way to print aln in standard out
    # /data/seed/Programs/clustal-omega-1.2.1/clustalo-1.2.0-Ubuntu-x86_64


# part 3: get consensus sequence
# ------------------------------------------------------------------------------------------------------------------

# generate consensus sequence from alignment file
def aln_file_process(aln_file, header):
    # get consensus sequence and then write it to consen.fasta
    os.system("cons -sequence "+aln_file+" -outseq "+aln_file+".fasta"+" -name '"+header+"'")
    os.system("cat "+aln_file+".fasta >> consen.fasta")
    # return consensus sequence


# helper function for merge_align_and_process()
def filename_to_cluster(filename):
    m = re.search("^(\S+)\.temp", filename)
    if m:
        return m.group(1)
    else:
        return None


def merge_align_and_process(i, f):
    temp_file = open(f[i], "r")
    cluster_id = filename_to_cluster(f[i])
    temp_content = temp_file.readlines()
    # check if fasta_map[i] contain more than one sequence
    if len(temp_content) > 2:
        log_file.write(cluster_id+"============="+f[i])
        temp_name = "temp2/"+cluster_id+"temp.msf"
        clustalo_align("".join(temp_content), temp_name)
        header = "Cluster"+cluster_id+", Number_of_elements: "+str(len(temp_content)/2)
        aln_file_process(temp_name, header)
    temp_file.close()


# part 5: main function
# ------------------------------------------------------------------------------------------------------------------

def main():
    start_time = time.time()

    print("test")

    # initialize consen.fasta
    exists = os.path.isfile('consen.fasta')
    if exists:
        os.system(" > consen.fasta")
    else:
        os.system("touch consen.fasta")

    multi_thread, initial = parameter_interpreter()

    f = []
    for (dirpath, dirnames, filenames) in os.walk(os.getcwd()):
        f.extend(filenames)
        break
    print(len(f))



    if multi_thread:
        pool_size = multi_thread
        pool = Pool(pool_size)
        pool_func = partial(merge_align_and_process, f=f)
        pool.map(pool_func, range(initial, len(f)))
    else:
        # single thread
        for i in range(initial, len(f)):
            merge_align_and_process(i, f)

    print("--- %s seconds ---" % (time.time() - start_time))
    # os.system("rm -r temp/")


if __name__ == '__main__':
    main()


