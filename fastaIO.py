"""This script contains several function for working with fasta formated sequence files. The main component is a fasta file iterator, FastaGeneralIterator, that is a modified version of the FastqGeneralIterator in the QualityIO.py script by Peter Cock that is is part of the Biopython distribution and is governed by its license. Please see the Biopython_License file that is included as part of this package for detail.

There is also code, FastaTitleStandardization, to standardize fasta record title lines output from third party software.

FastqGeneralIterator was modified to handle fasta files and FastaTitleStandardization was written by Brad Cavinder""" 


def FastaGeneralIterator(handle):
    """Iterate over Fasta records as string tuples (not as Biopython SeqRecord objects).
    """

    #Skip any text before the first record (e.g. blank lines, comments?)
    while True:
        line = handle.readline()
        if line == "" : return #Premature end of file, or just empty?
        if line[0] == ">":
            break

    while True:
        if line[0]!=">":
            raise ValueError("Records in Fasta files should start with '>' character")
        title_line = line[1:].rstrip()
        #Will now be at least one line of sequence data - string concatenation is used (if needed)
        #rather using than the "".join(...) trick just in case it is multiline:
        seq_string = handle.readline().rstrip()
        #There may now be more sequence lines, or the ">" for the next record:
        while True:
            line = handle.readline()
            if not line:
                break
            if len(line) == 0:
                continue
            if line[0] == ">":               
                break
            seq_string += line.rstrip() #removes trailing newlines
        #Assuming whitespace isn't allowed, any should be caught here:
        if " " in seq_string or "\t" in seq_string:
            raise ValueError("Whitespace is not allowed in the sequence.")
        seq_len = len(seq_string)

        #Return the record and then continue...
        yield (title_line, seq_string)
        if not line : return #StopIteration at end of file
    assert False, "Should not reach this line"
    
def FastaGeneralIterator2(handle):
    """Iterate over Fasta like records as string tuples. That is each entry starts with a ">" and header info and is followed by some type of text of one or more lines. This version does not remove newlines, keeping the input format in the output.
    """

    #Skip any text before the first record (e.g. blank lines, comments?)
    while True:
        line = handle.readline()
        if line == "" : return #Premature end of file, or just empty?
        if line[0] == ">":
            break

    while True:
        if line[0]!=">":
            raise ValueError("Records in Fasta files should start with '>' character")
        title_line = line
        #Will now be at least one line of sequence data - string concatenation is used (if needed)
        #rather using than the "".join(...) trick just in case it is multiline:
        seq_string = handle.readline()
        #There may now be more sequence lines, or the ">" for the next record:
        while True:
            line = handle.readline()
            if not line:
                break
            if line[0] == ">":               
                break
            seq_string += line

        #Return the record and then continue...
        yield (title_line, seq_string)
        if not line : return #StopIteration at end of file
    assert False, "Should not reach this line"
    
def FastaTitleStandardization(handle):
    """Use FastaGeneralIterator to iterate through records in a fasta file while standardizing the title of each record (third party programs are inconsitant to each other in title line naming), removing whitespace and special characters other than '-'. Also, ':' and '=' are changed to '-'. Title lines longer than 60 characters wil be shortened to 60.
    """
    import re
    pattern = re.compile('[^\w-]')
    for title, seq in FastaGeneralIterator(handle):
        #Remove whitespace and special characters except '-' from title and shorten it to no more than 100 characters
        title = title.replace(":", "-")
        title = title.replace("(", "-")
        title = title.replace(")", "-")
        title = title.replace("=", "-")
        title = pattern.sub("_", title)
        title = title.replace('___', '_')
        title = title.replace('__', '_')
        title = title.replace('_-_', '-')
        if len(title) >= 99:
            title = title[:99]
        title = title.replace(" ", "")
        yield (title, seq)

def FastqGeneralIterator2(handle):
    """Iterate over Fastq records as string tuples (not as SeqRecord objects).

    This code does not try to interpret the quality string numerically.  It
    just returns tuples of the title, sequence and quality as strings.  For
    the sequence and quality, any whitespace (such as new lines) is removed.

    Our SeqRecord based FASTQ iterators call this function internally, and then
    turn the strings into a SeqRecord objects, mapping the quality string into
    a list of numerical scores.  If you want to do a custom quality mapping,
    then you might consider calling this function directly.

    For parsing FASTQ files, the title string from the "@" line at the start
    of each record can optionally be omitted on the "+" lines.  If it is
    repeated, it must be identical.

    The sequence string and the quality string can optionally be split over
    multiple lines, although several sources discourage this.  In comparison,
    for the FASTA file format line breaks between 60 and 80 characters are
    the norm.

    WARNING - Because the "@" character can appear in the quality string,
    this can cause problems as this is also the marker for the start of
    a new sequence.  In fact, the "+" sign can also appear as well.  Some
    sources recommended having no line breaks in the  quality to avoid this,
    but even that is not enough, consider this example::

        @071113_EAS56_0053:1:1:998:236
        TTTCTTGCCCCCATAGACTGAGACCTTCCCTAAATA
        +071113_EAS56_0053:1:1:998:236
        IIIIIIIIIIIIIIIIIIIIIIIIIIIIICII+III
        @071113_EAS56_0053:1:1:182:712
        ACCCAGCTAATTTTTGTATTTTTGTTAGAGACAGTG
        +
        @IIIIIIIIIIIIIIICDIIIII<%<6&-*).(*%+
        @071113_EAS56_0053:1:1:153:10
        TGTTCTGAAGGAAGGTGTGCGTGCGTGTGTGTGTGT
        +
        IIIIIIIIIIIICIIGIIIII>IAIIIE65I=II:6
        @071113_EAS56_0053:1:3:990:501
        TGGGAGGTTTTATGTGGA
        AAGCAGCAATGTACAAGA
        +
        IIIIIII.IIIIII1@44
        @-7.%<&+/$/%4(++(%

    This is four PHRED encoded FASTQ entries originally from an NCBI source
    (given the read length of 36, these are probably Solexa Illumna reads where
    the quality has been mapped onto the PHRED values).

    This example has been edited to illustrate some of the nasty things allowed
    in the FASTQ format.  Firstly, on the "+" lines most but not all of the
    (redundant) identifiers are ommited.  In real files it is likely that all or
    none of these extra identifiers will be present.

    Secondly, while the first three sequences have been shown without line
    breaks, the last has been split over multiple lines.  In real files any line
    breaks are likely to be consistent.

    Thirdly, some of the quality string lines start with an "@" character.  For
    the second record this is unavoidable.  However for the fourth sequence this
    only happens because its quality string is split over two lines.  A naive
    parser could wrongly treat any line starting with an "@" as the beginning of
    a new sequence!  This code copes with this possible ambiguity by keeping
    track of the length of the sequence which gives the expected length of the
    quality string.

    Using this tricky example file as input, this short bit of code demonstrates
    what this parsing function would return:

    >>> handle = open("Quality/tricky.fastq", "rU")
    >>> for (title, sequence, quality) in FastqGeneralIterator(handle):
    ...     print title
    ...     print sequence, quality
    071113_EAS56_0053:1:1:998:236
    TTTCTTGCCCCCATAGACTGAGACCTTCCCTAAATA IIIIIIIIIIIIIIIIIIIIIIIIIIIIICII+III
    071113_EAS56_0053:1:1:182:712
    ACCCAGCTAATTTTTGTATTTTTGTTAGAGACAGTG @IIIIIIIIIIIIIIICDIIIII<%<6&-*).(*%+
    071113_EAS56_0053:1:1:153:10
    TGTTCTGAAGGAAGGTGTGCGTGCGTGTGTGTGTGT IIIIIIIIIIIICIIGIIIII>IAIIIE65I=II:6
    071113_EAS56_0053:1:3:990:501
    TGGGAGGTTTTATGTGGAAAGCAGCAATGTACAAGA IIIIIII.IIIIII1@44@-7.%<&+/$/%4(++(%
    >>> handle.close()

    Finally we note that some sources state that the quality string should
    start with "!" (which using the PHRED mapping means the first letter always
    has a quality score of zero).  This rather restrictive rule is not widely
    observed, so is therefore ignored here.  One plus point about this "!" rule
    is that (provided there are no line breaks in the quality sequence) it
    would prevent the above problem with the "@" character.
    """
    #We need to call handle.readline() at least four times per record,
    #so we'll save a property look up each time:
    handle_readline = handle.readline
    
    #Skip any text before the first record (e.g. blank lines, comments?)
    while True:
        line = handle_readline()
        if line == "" : return #Premature end of file, or just empty?
        if line[0] == "@":
            break

    while True:
        if line[0]!="@":
            raise ValueError("Records in Fastq files should start with '@' character")
        title_line = line[1:].rstrip()
        #Will now be at least one line of sequence data - in most FASTQ files
        #just one line! We therefore use string concatenation (if needed)
        #rather using than the "".join(...) trick just in case it is multiline:
        seq_string = handle_readline().rstrip()
        #There may now be more sequence lines, or the "+" quality marker line:
        while True:
            line = handle_readline()
            if not line:
                raise ValueError("End of file without quality information.")
            if line[0] == "+":
                #The title here is optional, but if present must match!
                second_title = line[1:].rstrip()
                #if second_title and second_title != title_line:
                #    raise ValueError("Sequence and quality captions differ.")
                break
            seq_string += line.rstrip() #removes trailing newlines
        #This is going to slow things down a little, but assuming
        #this isn't allowed we should try and catch it here:
        if " " in seq_string or "\t" in seq_string:
            raise ValueError("Whitespace is not allowed in the sequence.")
        seq_len = len(seq_string)

        #Will now be at least one line of quality data...
        quality_string = handle_readline().rstrip()
        #There may now be more quality data, or another sequence, or EOF
        while True:
            line = handle_readline()
            if not line : break #end of file
            if line[0] == "@":
                #This COULD be the start of a new sequence. However, it MAY just
                #be a line of quality data which starts with a "@" character.  We
                #should be able to check this by looking at the sequence length
                #and the amount of quality data found so far.
                if len(quality_string) >= seq_len:
                    #We expect it to be equal if this is the start of a new record.
                    #If the quality data is longer, we'll raise an error below.
                    break
                #Continue - its just some (more) quality data.
            quality_string += line.rstrip()
        
        #if seq_len != len(quality_string):
            #raise ValueError("Lengths of sequence and quality values differs "
                             #" for %s (%i and %i)." \
                             #% (title_line, seq_len, len(quality_string)))

        #Return the record and then continue...
        yield (title_line, seq_string, quality_string)
        if not line : return #StopIteration at end of file
    assert False, "Should not reach this line"

def median(in_list):
    #if list has even number of elements, take the avg of middle two
    #otherwise return middle elemenent
    in_list.sort()
    mid = len(in_list)/2
    if len(alist) % 2 == 0:  
        return (srtd[mid-1] + srtd[mid]) / 2.0
    else:
        return in_list[mid]            

def reverse_complement(seq):
    import string
    trans_table = string.maketrans("ATGCatgcNn", "TACGtacgNn")
    rev_seq = seq[::-1]
    rev_comp = rev_seq.translate(trans_table)
    return(rev_comp)
    
def complement_seq(seq):
    import string
    trans_table = string.maketrans("ATGCatgcNn", "TACGtacgNn")
    comp = seq.translate(trans_table)
    return(comp)
    
def individual_seq_len (in_path):
    in_handle = open(in_path, 'r') 
    out_path = in_path + ".lengths"
    out_handle = open(out_path, 'w')

    for title, seq in FastaGeneralIterator(in_handle):
        seq_len = len(seq)
        print>>out_handle, title, "\t", seq_len

    in_handle.close()
    out_handle.close()
    
def total_seq_len (in_path, write_mode, out_path, species_name):
    
    in_handle = open(in_path, "r")
    seq_len = 0
    for title, seq in FastaGeneralIterator(in_handle):
        seq_len +=  len(seq)
            
    in_handle.close()
    if write_mode == "append" or write_mode == "a" or write_mode == "A":
        out_handle = open(out_path, "a")
        print>>out_handle, "\t".join([species_name , str(seq_len)])
    else:
        out_handle = open(out_path, "w")
        print>>out_handle, "species\tgenome_length\n" + species_name + "\t" + str(seq_len)
        
    out_handle.close()

def sequence_retriever(contig, start, end, flank, genome_dict3):
    needed_left = 0
    needed_right = 0
    wanted_seq = ''
    add_left = ''
    add_right = ''
    left_coord = ''
    right_coord = ''
    if contig in genome_dict3:
        seq = genome_dict3[contig]
        contig_seq_len = len(seq)
        if int(flank) < int(start):
            left_coord = (int(start)-int(flank))
        else:
            needed_left = (int(flank) - int(start))
            left_coord = 0
        if (contig_seq_len - int(flank)) >= int(end):
            right_coord = (int(end+1) + int(flank))
        else:
            needed_right = int(end) - ((contig_seq_len - int(flank)))
            right_coord = contig_seq_len
        if needed_left > 0:
            add_left = "N" * needed_left
        if needed_right > 0:
            add_right = "N" * needed_right
        wanted_seq = add_left + seq[left_coord:right_coord] + add_right
    else:
        print "Didn't find " + contig + " in sequence dictionary."
    return wanted_seq
    
def sequence_retriever2(contig, start, end, left_flank, right_flank, genome_dict3):
    needed_left = 0
    needed_right = 0
    wanted_seq = ''
    add_left = ''
    add_right = ''
    left_coord = 0
    right_coord = 0
    if contig in genome_dict3:
        seq = genome_dict3[contig]
        contig_seq_len = len(seq)
        if int(left_flank) < int(start):
            left_coord = (int(start)-int(left_flank))
        else:
            needed_left = (int(left_flank) - int(start))
            left_coord = 0
        if (contig_seq_len - int(right_flank)) >= int(end):
            right_coord = (int(end)+1 + int(right_flank))
        else:
            needed_right = int(end) - ((contig_seq_len - int(right_flank)))
            right_coord = int(contig_seq_len)
        if needed_left > 0:
            add_left = "N" * needed_left
        if needed_right > 0:
            add_right = "N" * needed_right
        wanted_seq = add_left + seq[left_coord:right_coord] + add_right
    else:
        print "Didn't find " + contig + " in sequence dictionary."
    right_coord = right_coord - 1
    return wanted_seq, left_coord, right_coord

def SplitLongString(s, w):
    temp = ''.join(s[x:x+w] + '\n' for x in xrange(0, len(s), w))
    temp = temp.rstrip("\n")
    return temp
    
def shuffle_split(fpath, split_number):
    """Shuffle and split a fasta file into groups of ~350""" 
    
    import math
    import random
    
    copy_list = []
    copy_dict = {}
    group_list = []
    path_list = []
    
    in_handle = open(fpath, "r")
    for title, seq in FastaGeneralIterator(in_handle):
        title = title.strip("\n").strip()
        copy_list.append(title)
        copy_dict[title] = seq
    in_handle.close()

    copy_num = len(copy_list)
    groups = int(round(copy_num/float(split_number)))
    copies_to_group = int(math.ceil(float(copy_num)/groups))
    random.shuffle(copy_list)
    
    i = 0
    while i < groups:
        start = copies_to_group * i
        end = (start + copies_to_group)-1
        if start < copy_num:
            if end < copy_num:
                group_list.append(copy_list[start:end])
            else:
                group_list.append(copy_list[start:])
        i += 1
    c = 1
    for group in group_list:
        out_path = fpath + ".group" + str(c) + "_split"
        path_list.append(out_path)
        out_handle = open(out_path, "w")
        for title in group:
            print>>out_handle, ">" + title + "\n" + copy_dict[title]
        out_handle.close()
        c += 1
    return(path_list, copies_to_group)

class Vividict(dict):
    """Autovivification of dictionaries"""
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
