"""
Name: Au Yeung Tsz Lok
Date: 12/11/2023
Description: To output the combinations of kmers with two different method, and output the
             operation time with log file
"""

import time
import sys


def main():
    # Read the kmer values (a list of str or strs)
    km = sys.argv[2].split(",")
    km = [int(x) for x in km]

    # For non dictionary function
    sequence = load_fasta_file()
    answer, time_used_non_dict = count_kmers_non_dict(sequence, km)
    write2file_non_dict(answer, km)

    # For dictionary function
    sequence = load_fasta_file()
    answer, time_used_dict = count_kmers(sequence, km)
    write2file(answer, km)

    # Comparison if -c option is present and only have 1 requested kmer value
    if len(sys.argv) == 4 and len(km) == 1:
        with open("log_02836.txt", "a") as outfile:
            for x in range(len(km)):
                outfile.write(f"For {km[x]}-mers, time elapsed without the use of dictionary is \
{time_used_non_dict[x] :g} seconds.\n")
                outfile.write(f"For {km[x]}-mers, time elapsed with the use of dictionary is \
{time_used_dict[x] :g} seconds.\n")
                outfile.write(f"\n")
        outfile.close()
    elif len(sys.argv) == 4 and len(km) != 1:
        print("[ERROR] The -c option only applies to 1 requested k-value")


def load_fasta_file():
    """

    :return: Lists - sequence of each gene id(s) respectively according to order of the input file
    """
    file = sys.argv[1]
    with open(file, "r") as inputfile:
        # For storing the sequences from reading the file (multiple sequences if more than 1 id has read)
        seq = []
        # List of fasta file suffix
        fa = [".fasta", ".fa", ".fna"]
        # List of fastq file suffix
        fq = [".fastq", ".fq"]

        # Find the file whether is Fasta
        for suffix in fa:
            # Is fasta
            if file.endswith(suffix):
                # Skip the first line of the file, which is the name description
                lines = inputfile.readlines()[1:]
                # Empty list for storing the sequence
                sequence = []
                # Read sequence from every line of the file
                for line in lines:
                    # Skip the id/name line, which starts with ">"
                    if line.startswith(">"):
                        # Join the stored sequences to a single string and append to seq list
                        seq.append("".join(sequence))
                        # Clear the original sequences for storing new gene id sequences
                        sequence.clear()
                    else:
                        # Add the sequences for each line that is still within the same gene id
                        sequence.append(line.strip())
                # For adding the last gene id sequences, because there will be no new ">" containing lines after this
                seq.append(''.join(sequence))

        # Find the file whether is Fastq
        for suffixq in fq:
            # Is fastq
            if file.endswith(suffixq):
                # Skip the first line, add every sequence lines that appears every 4th row starting from the 1st seq
                lines = inputfile.readlines()[1::4]
                # Empty list for storing the sequence
                sequence = []
                # Assumption/limitation: The fastq file is/must be a single-end file
                # Concept: Suppose every gene id has two sequences in each fastq file\
                # so every 2 sequences can be joined as the same string for the gene id
                # a control value for counting 2 consecutive sequences
                i = 1
                # Read every sequence from the file
                for x in lines:
                    # add to the sequence list for every non-ending sequence\
                    # and the 1st sequence for every 2 consecutive sequences
                    if i != 2 and x != lines[-1]:
                        # Add the 1st sequence to the sequence list
                        sequence.append(x.strip())
                        # Add 1 to the control value indicating the following sequence is the 2nd sequence
                        i += 1
                    # If it is ending sequence, it maybe doesn't have a pair of sequences for the gene id\
                    # still need to add to sequence list
                    elif i != 2 and x == lines[-1]:
                        sequence.append(x.strip())
                        # Forming a string to append to the seq list
                        seq.append(''.join(sequence))
                    # For every 2nd sequence in every 2 consecutive sequences
                    elif i == 2:
                        # Add the sequence to the sequence list that contains the 1st sequence
                        sequence.append(x.strip())
                        # Join with the 1st sequence to become a string, then append to the seq list
                        seq.append(''.join(sequence))
                        # Clear the sequence list for storing the next gene id sequences
                        sequence.clear()
                        # Minus 1 to the control value indicating the following sequence is the 1st sequence\
                        # of the new gene id
                        i -= 1
    return seq


def count_kmers_non_dict(test, k):
    """

    :param test: str - the sequence(s) from user input file
    :param k: int - the value(s) of k
    :return: answer2: list - the occurrence of each kmer for each sequence in each k-value requested
             stime: list - the used time for running in each kmer value
    """
    # Empty list for storing kmer with occurrence information for each kmer value
    answer2 = []
    # Empty list for storing running time for each requested kmer value
    stime = []
    # Loop along each requested kmer value(s)
    for values in k:
        # Start timing
        start_time = time.time()
        # Empty list for storing occurrence of kmers for each gene id sequence
        ans = []
        # Loop the sequence of each gene id
        for seq in test:
            # Empty list for storing kmer with occurrence information respectively for each gene id
            answer = []
            # Empty list for storing kmers (included non-valid kmers)
            kmer = []
            # Transform all the bases in sequence to upper base
            seqs = seq.upper()
            # State the valid nucleotides
            valid_bases = ["A", "T", "C", "G"]
            # Empty list for storing the occurrence of each kmer
            value = []
            # Empty set for storing non-valid kmers
            non_valid = set()
            # Record each combination of kmers (includes duplicates and non-valid kmers)
            for base in range(len(seqs) - values + 1):
                kmer.append(seqs[base:base + values])
            # Sort the kmers list
            kmer.sort()
            # Record the non-valid kmers
            for nucleotides in kmer:
                for nucleo in nucleotides:
                    for x in range(len(nucleo)):
                        if nucleo[x] not in valid_bases:
                            non_valid.add(nucleotides)
            # Make a new list for removing non valid kmers
            kmer2 = [good for good in kmer if good not in non_valid]
            # Make a set to record valid kmers (remove duplicates)
            kmer3 = set(kmer2.copy())
            # Transform set to list, then sort it
            kmer3 = sorted(list(kmer3))
            # Record the occurrence of each kmer to the value list
            for item in kmer3:
                value.append(kmer2.count(item))
            # Add list of each kmer with occurrence information to answer list
            for x in range(len(kmer3)):
                answer.append([kmer3[x], value[x]])
            for info in answer:
                for x in ans:
                    for y in x:
                        if info[0] == y[0]:
                            y[1] += info[1]
            # Add list of every kmer with occurrence information to ans list
            ans.append(answer)
        # Stop the time
        end_time = time.time()
        # Count used time
        used_time = end_time - start_time
        # Add the occurrence of kmers information for every gene id to answer2 list
        answer2.append(ans[0])
        # Add the used time to stime list
        stime.append(int(used_time))
    return answer2, stime


def count_kmers(seqs, k):
    """

    :param seqs: str - the sequence(s) from user input file
    :param k: int - the value(s) of k
    :return: answer: list - the dictionaries of occurrence of each kmer for each sequence in each k-value requested
             time_used: list - the used time for running in each kmer value
    """
    # Empty list for storing kmer with occurrence information of each gene id(s) for each kmer value
    answer = []
    # Empty list for storing time used in running each kmer value
    time_used = []
    # Valid nucleotides
    valid = ["A", "T", "C", "G"]
    # Loop for each kmer value
    for ks in k:
        # Start timing
        start_time = time.time()
        # Empty list for storing non-valid kmers
        non_valid = []
        # Empty list for storing kmers with occurrence information
        ans = []
        # Loop for each gene id(s) sequence
        for seq in seqs:
            # Transform all bases in sequence to upper base
            seq = seq.upper()
            # Empty dictionary for storing the kmer and respective occurrence information
            dic = {}
            # Find possible kmer
            for base in range(len(seq) - ks + 1):
                # If the found kmer is not recorded, initialize it
                if seq[base:base + ks] not in dic:
                    dic[seq[base:base + ks]] = 1
                # If found, add 1 to the original occurrence value
                else:
                    dic[seq[base:base + ks]] += 1
            # Loop for the keys of the dictionary
            for key in dic.keys():
                # Loop for the nucleotides of each kmer, record the non_valid kmers to non_valid list
                for x in range(len(key)):
                    if key[x] not in valid:
                        non_valid.append(key)
            # Make a copy of dictionary for further modification
            dic2 = dic.copy()
            # Remove the non-valid kmers in dic2
            for key in dic.keys():
                if key in non_valid:
                    del dic2[key]
            # Sort the keys order in dic2
            keys = list(dic2.keys())
            keys.sort()
            dic2 = {i: dic2[i] for i in keys}
            # Add the dictionary of valid kmer occurrences of the gene id sequence to ans list
            if len(ans) == 0:
                ans.append(dic2)
            else:
                for key in dic2.keys():
                    if key in ans[0]:
                        (ans[0])[key] = (ans[0])[key] + dic2[key]
                    else:
                        pass

        # Sort the keys order in for the dictionary in ans
        keys = list(ans[0].keys())
        keys.sort()
        ans[0] = {i: (ans[0])[i] for i in keys}
        # Add the dictionary valid kmer occurrences of the gene id(s) sequence to answer list
        answer.append(ans)
        # Stop timing
        end_time = time.time()
        # Find the used time
        times = end_time - start_time
        # Add the used time to the time_used list
        time_used.append(int(times))  # Record the occurrence of kmer in different k value respectively
    return answer, time_used


def write2file_non_dict(seqs, k):
    """

    :param seqs: list - the occurrence of each kmer of gene id(s) for each kmer value
    :param k: list - the kmer value(s)
    :return: txt - testfile of the outputs of each kmer with frequencies respectively
    """
    i = 0
    for value in k:
        with open(f"testfile_{value}-mer-frequency.txt", "w") as outfile:
            for y in seqs[i]:
                outfile.write(f'{y[0]}:{y[1]}\n')
            outfile.write(f'\n')
        outfile.close()
        i += 1


def write2file(seqs, k):
    """

    :param seqs: list - the occurrence of each kmer of gene id(s) for each kmer value
    :param k: list - the kmer value(s)
    :return: txt - testfile of the outputs of each kmer with frequencies respectively
    """
    i = 0
    for value in k:
        with open(f"testfile_{value}-mer-frequency.txt", "w") as outfile:
            for y in seqs[i]:
                for key in y.keys():
                    outfile.write(f'{key}:{y[key]}\n')
            outfile.write(f'\n')
        outfile.close()
        i += 1


main()
