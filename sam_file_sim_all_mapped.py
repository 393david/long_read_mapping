# Using SAM file to calculate similarity and different errors, does not use pysam to read SAM file.
# For standard sam file: the position in reference is 1-indexed.
# This script print the following information:
# match insert  delete  mismatch    length  mapped_based   matched_bases   score
# In the output file, following are outputed: 
# match insert  delete  mismatch    length  mapq
# This reports all mapped results for each read, including the complement mapping 
import argparse
import os, sys, cigar
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from aligntools import Cigar
if __name__ == "__main__":
    usage = '''  Compute all aligned similarity for a sam file.
  sys.maxsize
  Usage example:
  python ''' + sys.argv[0] + ' -s 1 -i aligned.sam -r ref.fa -o sim.txt'
    
    parser = argparse.ArgumentParser(description = usage, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input sam file [required]')
    parser.add_argument('-r', '--ref', type = str, required = True, help = 'The reference fasta file [required]')
    parser.add_argument('-q', '--query', type = str, required = True, help = 'The query read file [required]')
    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output similarity file [required]')
    parser.add_argument('-m', '--match', type = int, default = 1, help = 'match score [1]')
    parser.add_argument('-d', '--delete', type = int, default = -1, help = 'delete score [-1]')
    parser.add_argument('-I', '--insert', type = int, default = -1, help = 'insert score [-1]')
    parser.add_argument('-M', '--mismatch', type = int, default = -1, help = 'mismatch score [-1]')
    args = parser.parse_args()
    file = args.input
    outt = args.output
    query_file = args.query
    based = 1
    len_max = 0
    len_min = sys.maxsize
    ref_file = args.ref
    print()
    if not os.path.exists(file):
        print(' %s does not exist!' % file)
        sys.exit(1)
    if not os.path.exists(ref_file):
        print(' %s does not exist!' % ref_file)
        sys.exit(1)
    #print(' Reading reference file...')
    n_seq = 0
    n_bases = 0
    ref_seq = {}
    reads = {}
    reads_reverse = {}
    for seq_record in SeqIO.parse(ref_file, "fasta"):
        ll = len(seq_record)
        #print(seq_record.id)
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
    print(' Min length:      %d' % len_min)
    
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
    mapped_bases  = 0
    mapped_score  = 0
    matched_bases = 0
    #soft_right = 0
    aligned_header = []
    oo = open(outt, 'w')
    oo.write('match\tinsert\tdelete\tmismatch\tlength\tmapq\n')
    f = open(file)
    ll = f.readline()
    while ll:
        if ll.startswith('@'):
            ll = f.readline()
            continue
        ll = ll.strip().split()
        query_name = ll[0]
        mapq = ll[4]
        soft = 0
        if ll[5] == '*':
            ll = f.readline()
            continue
        strandd = 0
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
                print(' Can not recognize mapped flag: %d' % nnn)
                sys.exit(1)
            if strandd == 'plus':
                reads[query_name] = ll[9]
                seq_class = Seq(ll[9])
                reverse = seq_class.reverse_complement()
                reads_reverse[query_name] = str(reverse)
            elif strandd == 'minus':
                reads_reverse[query_name] = ll[9]
                seq_class = Seq(ll[9])
                pluss = seq_class.reverse_complement()
                reads[query_name] = str(pluss)
            else:
                print( 'Query: ' % query_name)
                print(' Mapped strand can not be recognized!')
                sys.exit(1)
        # else:
        #     ll = f.readline()
        #     continue
        n_seq += 1
        sys.stdout.write("Processing the %d-th mapped records, mapped unique reads: %d   \r" % (n_seq, len(aligned_header)))
        sys.stdout.flush()
        if len(aligned_header) < 22545:
            ll = f.readline()
            continue
        #aligned_header.append(rec.query_name)
        #print(query_name)
        #print(rec.reference_name)
        #if query_name == 'm140612_020550_42156_c100652082550000001823118110071460_s1_p0/7461/0_6631':
        #    bbbb = 0
        ref_is = ref_seq[ll[2]]
        if based:
            ref_start = int(ll[3]) - 1 
        else:
            ref_start = int(ll[3])
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
        cigar_class = Cigar.coerce(ll[5])
        cc = cigar.Cigar(ll[5])
        cigar_list = list(cc.items())
        if cigar_list[-1][1] == 'S':
            soft += cigar_list[-1][0]
        if cigar_list[0][1] == 'S':
            soft += cigar_list[0][0]
        if cigar_list[0][1] == 'H':
            query_is = query_is[cigar_list[0][0] : ]
        len_ref   = cigar_class.ref_length
        len_query = cigar_class.query_length
        ref_msa, query_msa = cigar_class.to_msa(ref_is[ref_start : ref_start + len_ref], query_is)
        n_match  = 0
        n_delete = 0
        n_insert = 0
        n_mismatch = 0
        len_aligned = len(ref_msa)
        for i in range(len_aligned):
            if ref_msa[i] == query_msa[i]:
                n_match += 1
            elif ref_msa[i] == '-':
                n_insert += 1
            elif query_msa[i] == '-':
                n_delete += 1
            else:
                n_mismatch += 1
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
        mapped_bases  += mapped_length
        matched_bases += n_match
        mapped_score  += n_match * args.match + n_insert * args.insert + n_delete * args.delete + n_mismatch * args.mismatch
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
            # if error_match < 0.62:
            #     haha = 0
        if error_insert < min_insert:
            min_insert = error_insert
        if error_delete < min_delete:
            min_delete = error_delete
        if error_mismatch < min_mismatch:
            min_mismatch = error_mismatch
        if mapped_length < min_mapped_length:
            min_mapped_length = mapped_length
        oo.write(str(error_match) + '\t')
        oo.write(str(error_insert) + '\t')
        oo.write(str(error_delete) + '\t')
        oo.write(str(error_mismatch) + '\t')
        oo.write(str(mapped_length) + '\t')
        oo.write(str(mapq) + '\n')
        ll = f.readline()
    f.close()
    oo.close()
    print('\n\n Original file:')
    print(' Mapped reads:     %d' % len(aligned_header))
    print(' Mapped records:   %d' % n_seq)
    print(' Mapped bases:     %d' % mapped_bases)
    print(' Matched bases:    %d' % matched_bases)
    print(' Mapped score:     %d' % mapped_score)
    print()
    print(' Match average:    %f' % (average_match/n_seq))
    print('       max:        %f' % max_match)
    print('       min:        %f' % min_match)
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

