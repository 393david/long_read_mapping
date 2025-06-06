# This is used for not-standard SAM file, since the mapped position in reference is not 1-based indexed.
# Such as the mapAlign method.
# This reports all mapped records for each read, including the complement mapping
import argparse
import os, sys, cigar
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from aligntools import Cigar

if __name__ == "__main__":
    usage = '''  Translate SAM file to blast_like output file.
  
  Usage example:
  python ''' + sys.argv[0] + ' -i aligned.sam -r ref.fa -o aligned.txt'
    
    parser = argparse.ArgumentParser(description = usage, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input sam file [required]')
    parser.add_argument('-q', '--query', type = str, required = True, help = 'Query read file [required]')
    parser.add_argument('-r', '--ref', type = str, required = True, help = 'The reference fasta file [required]')
    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output similarity file [required]')
    parser.add_argument('-s', '--start', type = int, default = 1, help = 'Mapped start position in reference is 0-based or 1-based [1]')

    #parser.add_argument('-m', '--min', type = int, default = 0, help = 'Min length [0]')
    #parser.add_argument('-M', '--max', type = int, default = sys.maxsize, help = 'Max length [infinity]')
    args = parser.parse_args()
    based = args.start
    file = args.input
    outt = args.output
    queryy = args.query
    len_max = 0
    len_min = sys.maxsize
    ref_file = args.ref
    if based not in [0, 1]: #!= 0 and based != 1:
        print(' Start position can be only 0 or 1!')
        sys.exit(1)
    print()
    if not os.path.exists(file):
        print(' %s does not exist!' % file)
        sys.exit(1)
    if not os.path.exists(queryy):
        print(' %s does not exist!' % queryy)
        sys.exit(1)
    #print(' Reading reference file...')
    n_seq = 0
    n_bases = 0
    ref_seq = {}
    reads = {}
    reads_reverse = {}
    # reading the query reads
    for seq_record in SeqIO.parse(queryy, "fasta"):
        ll = len(seq_record)
        #print(seq_record.id)
        reads[seq_record.id] = str(seq_record.seq).upper()
        seq_class = Seq(str(seq_record.seq).upper())
        reads_reverse[seq_record.id] = str(seq_class.reverse_complement())
        n_seq += 1
        n_bases += ll
        sys.stdout.write(" Reading the %d-th query read, length: %d  \r" % (n_seq, ll))
        sys.stdout.flush()
    print('\n')
    # reading the reference file
    n_seq = 0
    n_bases = 0
    for seq_record in SeqIO.parse(ref_file, "fasta"):
        ll = len(seq_record)
        ref_seq[seq_record.id] = str(seq_record.seq).upper()
        n_seq += 1
        n_bases += ll
        if ll > len_max:
            len_max = ll
        if ll < len_min:
            len_min = ll
        sys.stdout.write(" Reading the %d-th reference, length: %d  \r" % (n_seq, ll))
        sys.stdout.flush()
    print('\n\n Reference file:')
    print(' Sequence number: %d' % n_seq)
    print(' Bases number:    %d' % n_bases)
    print(' Average length:  %f' % (n_bases/n_seq))
    print(' Max length:      %d' % len_max)
    print(' Min length:      %d\n' % len_min)
    
    # read the SAM file
    n_seq = 0
    average_match = 0
    average_insert = 0
    average_delete = 0
    average_mismatch = 0
    average_mapped_length = 0
    min_match = 1
    min_insert = 1
    min_delete = 1
    min_mismatch = 1
    min_mapped_length = sys.maxsize
    max_match = 0
    max_insert = 0
    max_delete = 0
    max_mismatch = 0
    max_mapped_length = 0
    #soft_left = 0
    #soft_right = 0
    aligned_header = []
    oo = open(outt, 'w')
    #f = pysam.AlignmentFile(file)
    f = open(file)
    ll = f.readline()
    while ll:
    #for rec in f:
        if ll.startswith('@'):
            ll = f.readline()
            continue
        soft = 0
        ll = ll.strip().split()
        query_name = ll[0]
        #if query_name != 'S1_996541':
        #    ll = f.readline()
        #    continue
        #print(query_name)
        strandd = 0
        if ll[5] == '*':
            ll = f.readline()
            continue
        if query_name not in aligned_header:
            aligned_header.append(query_name)
            nnn = int(ll[1])
            if nnn == 0:
                strandd = 'plus'
            elif nnn == 16:
                strandd = 'minus'
            elif nnn == 256:
                strandd = 'plus'
            elif nnn == 272:
                strandd = 'minus'
            elif nnn == 2048:
                strandd = 'plus'
            elif nnn == 2064:
                strandd = 'minus'
            else:
                print(' Query name: %s', query_name)
                print(' Can not recognize mapped flag: %d' % nnn)
                sys.exit(1)
            #if strandd == 'plus':
            #    reads[query_name] = ll[9]
            #    seq_class = Seq(ll[9])
            #    reverse = seq_class.reverse_complement()
            #   reads_reverse[query_name] = str(reverse)
            #elif strandd == 'minus':
            #    reads_reverse[query_name] = ll[9]
            #    seq_class = Seq(ll[9])
            #    pluss = seq_class.reverse_complement()
            #    reads[query_name] = str(pluss)
            if not (strandd == 'plus' or strandd == 'minus'):
                print( 'Query: ' % query_name)
                print(' Mapped strand can not be recognized!')
                sys.exit(1)
        # else:
        #     ll = f.readline()
        #     continue
        n_seq += 1
        middle = ''
        sys.stdout.write("Processing the %d-th mapped records, mapped unique reads: %d   \r" % (n_seq, len(aligned_header)))
        sys.stdout.flush()
        if strandd == 0:
            nnn = int(ll[1])
            if nnn == 0:
                strandd = 'plus'
            elif nnn == 16:
                strandd = 'minus'
            elif nnn == 256:
                strandd = 'plus'
            elif nnn == 272:
                strandd = 'minus'
            elif nnn == 2048:
                strandd = 'plus'
            elif nnn == 2064:
                strandd = 'minus'
            else:
                print(' Query name: %s', query_name)
                print(' Can not recognize mapped strand: %d' % nnn)
                sys.exit(1)
        if strandd == 'plus':
            query_is = reads[query_name]
        elif strandd == 'minus':
            query_is = reads_reverse[query_name]
        else:
            print(' Query name: %s', query_name)
            print(' Can not recognize mapped flag: %d' % nnn)
            sys.exit(1)
        #aligned_header.append(rec.query_name)
        len_query = len(query_is)
        ref_is = ref_seq[ll[2]]
        if based:
            ref_start = int(ll[3]) - 1 
        else:
            ref_start = int(ll[3])
        #ref_start = int(ll[3])#rec.reference_start
        cigar_class = Cigar.coerce(ll[5])#(rec.cigarstring)
        cc = cigar.Cigar(ll[5])#(rec.cigarstring)
        cigar_list = list(cc.items())
        soft_left = 0
        soft_right = 0
        if cigar_list[-1][1] == 'S':
            soft += cigar_list[-1][0]
            soft_right = cigar_list[-1][0]
        if cigar_list[0][1] == 'S':
            soft += cigar_list[0][0]
            soft_left = cigar_list[0][0]
        if cigar_list[0][1] == 'H':
            query_is = query_is[cigar_list[0][0] : ]
        len_ref_aligned   = cigar_class.ref_length
        len_query_aligned = cigar_class.query_length
        ref_msa, query_msa = cigar_class.to_msa(ref_is[ref_start : ref_start + len_ref_aligned], query_is)
        n_match  = 0
        n_delete = 0
        n_insert = 0
        n_mismatch = 0
        len_aligned = len(ref_msa)
        for i in range(len_aligned):
            if ref_msa[i] == query_msa[i]:
                n_match += 1
                middle += '|'
            elif ref_msa[i] == '-':
                n_insert += 1
                middle += ' '
            elif query_msa[i] == '-':
                n_delete += 1
                middle += ' '
            else:
                n_mismatch += 1
                middle += '*'
        n_insert = n_insert - soft      
        mapped_length = n_match + n_mismatch + n_insert  
        error_match = n_match/(len_aligned - soft)
        error_insert = n_insert/(len_aligned - soft)
        error_delete = n_delete/(len_aligned - soft)
        error_mismatch = n_mismatch/(len_aligned - soft)
        average_match += error_match
        average_insert += error_insert
        average_delete += error_delete
        average_mismatch += error_mismatch
        average_mapped_length += mapped_length
        if error_match > max_match:
            max_match = error_match
        if error_insert > max_insert:
            max_insert = error_insert
        if error_delete > max_delete:
            max_delete = error_delete
        if error_mismatch > max_mismatch:
            max_mismatch = error_mismatch
        if mapped_length > max_mapped_length:
            max_mapped_length = mapped_length
        if error_match < min_match:
            min_match = error_match
            #if error_match < 0.202:
            #    haha = 0
        if error_insert < min_insert:
            min_insert = error_insert
        if error_delete < min_delete:
            min_delete = error_delete
        if error_mismatch < min_mismatch:
            min_mismatch = error_mismatch
        if mapped_length < min_mapped_length:
            min_mapped_length = mapped_length
        oo.write(' Query:      ' + query_name + '\n')
        oo.write(' nMatch:     ' + str(n_match) + '(' + str(error_match)[0:8] + ')' + '\n')
        oo.write(' nInsert:    ' + str(n_insert) + '(' + str(error_insert)[0:8] + ')' + '\n')
        oo.write(' nDelete:    ' + str(n_delete) + '(' + str(error_delete)[0:8] + ')' + '\n')
        oo.write(' nMismatch:  ' + str(n_mismatch) + '(' + str(error_mismatch)[0:8] + ')' + '\n')
        oo.write(' Similarity: ' + str(error_match) + '\n')
        oo.write(' Target:     ' + ll[2] + '\n')
        oo.write(' Strand:     ' + strandd + '\n')
        if cigar_list[0][1] == 'H':
            oo.write(' Range:      ' + str(cigar_list[0][0]) + ' -> ' + str(cigar_list[0][0] + len_query_aligned - 1) + ' of ' + str(len_query) + ', target: ' + str(ref_start) + ' -> ' + str(ref_start + len_ref_aligned - 1) + ' of ' + str(len(ref_is)) + '\n\n')
        else:
            oo.write(' Range:      ' + str(soft_left) + ' -> ' + str(len_query - soft_right - 1) + ' of ' + str(len_query) + ', target: ' + str(ref_start) + ' -> ' + str(ref_start + len_ref_aligned - 1) + ' of ' + str(len(ref_is)) + '\n\n')
        per_line = 60
        pos_together_start = soft_left
        pos_together_end = 0
        pos_query_start = soft_left
        if cigar_list[0][1] == 'H':
            pos_query_start = cigar_list[0][0]
        pos_ref_start   = ref_start#rec.reference_start
        pos_query_end   = pos_query_start - 1
        pos_query_end_pre = pos_query_end
        pos_ref_end     = pos_ref_start - 1
        pos_ref_end_pre = pos_ref_end
        for i in range(soft_left, len(middle) - soft_right):
            if pos_together_start >= len(middle) - soft_right:
                oo.write('\n\n\n')
                break
            if pos_together_start + per_line <= len(middle) - soft_right:
                pos_together_end = pos_together_start + per_line
            else:
                pos_together_end = len(middle) - soft_right
            for jj in range(pos_together_start, pos_together_end):
                #oo.write(ref_msa[jj])
                if ref_msa[jj] != '-':
                    pos_ref_end += 1
            if pos_ref_end > pos_ref_end_pre:
                pos_ref_start = pos_ref_end_pre + 1
            else:
                pos_ref_start = pos_ref_end_pre
            oo.write('Target %11d %s %d\n' % (pos_ref_start, ref_msa[pos_together_start : pos_together_end], pos_ref_end))
            #oo.write(' %d\n' % pos_ref_end)
            # middle line
            oo.write('                   %s\n' % middle[pos_together_start : pos_together_end])
            #oo.write(' Query %11d ' % pos_query_start)
            for jj in range(pos_together_start, pos_together_end):
                #oo.write(query_msa[jj])
                if query_msa[jj] != '-':
                    pos_query_end += 1
            if pos_query_end > pos_query_end_pre:
                pos_query_start = pos_query_end_pre + 1
            else:
                pos_query_start = pos_query_end_pre
            oo.write(' Query %11d %s %d\n\n' %(pos_query_start, query_msa[pos_together_start : pos_together_end], pos_query_end))
            pos_query_end_pre = pos_query_end
            pos_ref_end_pre   = pos_ref_end
            #pos_ref_start = pos_ref_end + 1
            pos_together_start = pos_together_end
        ll = f.readline()

    f.close()
    oo.close()
    print('\n\n Original file:')
    print(' Aligned number:     %d' % n_seq)
    print(' Aligned Bases number:     %d' % n_bases)
    print(' Match average:   %f' % (average_match/n_seq))
    print('       max:       %f' % max_match)
    print('       min:       %f' % min_match)
    print(' Insert average:   %f' % (average_insert/n_seq))
    print('        max:       %f' % max_insert)
    print('        min:       %f' % min_insert)
    print(' Delete average:   %f' % (average_delete/n_seq))
    print('        max:       %f' % max_delete)
    print('        min:       %f' % min_delete)
    print(' Mismatch average: %f' % (average_mismatch/n_seq))
    print('          max:     %f' % max_mismatch)
    print('          min:     %f' % min_mismatch)
    print(' Mapped length average: %f' % (average_mapped_length/n_seq))
    print('               max:     %f' % max_mapped_length)
    print('               min:     %f' % min_mapped_length)

