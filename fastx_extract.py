import argparse
import os, sys
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

if __name__ == "__main__":
    usage = '''  Extract longer or shorter reads, write them to a fa file.
  
  Usage example:
  python ''' + sys.argv[0] + ' -m 500 -i reads.fa -o extracted.fa'
    
    parser = argparse.ArgumentParser(description = usage, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input fastx file [required]')
    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output fastx file [required]')
    parser.add_argument('-n', '--num', type = int, default = 0, help = 'Extract sequence number [0: all sequence that satified min and max length]')
    parser.add_argument('-m', '--min', type = int, default = 0, help = 'Min length [0]')
    parser.add_argument('-M', '--max', type = int, default = sys.maxsize, help = 'Max length [infinity]')
    args = parser.parse_args()
    file = args.input
    outt = args.output
    th_max = args.max
    th_min = args.min
    num2 = args.num
    print()
    if not os.path.exists(file):
        print('%s does not exist!' % file)
        sys.exit(1)
    name, extention = os.path.splitext(file)
    n_seq = 0
    n_seq2 = 0
    n_bases = 0
    n_bases2 = 0
    len_min = sys.maxsize
    len_max = 0
    len_min2 = sys.maxsize
    len_max2 = 0
    if extention in ['.fa', '.fasta']:
        oo = open(outt, 'w')
        for seq_record in SeqIO.parse(file, "fasta"):
            ll       = len(seq_record)
            n_bases += ll
            n_seq   += 1
            if ll > len_max:
                len_max = ll
            if ll < len_min:
                len_min = ll
            sys.stdout.write("Read number: %d   \r" % n_seq)
            sys.stdout.flush()
            if ll >= th_min and ll <= th_max:
                oo.write('>' + seq_record.id + '\n')
                oo.write(str(seq_record.seq) + '\n')
                n_seq2  += 1
                n_bases2 += ll
                if ll > len_max2:
                    len_max2 = ll
                if ll < len_min2:
                    len_min2 = ll
                if num2 > 0 and n_seq2 == num2:
                    break
        oo.close()
    elif extention in ['.fq', '.fastq']:
        oo = open(outt, 'w')
        in_handle = open(file)
        iterator = FastqGeneralIterator(in_handle)
        for triplet in iterator:
            (description, sequence, quality) = triplet
            ll = len(sequence)
            n_bases += ll
            n_seq   += 1
            if ll > len_max:
                len_max = ll
            if ll < len_min:
                len_min = ll
            sys.stdout.write("Reads number: %d   \r" % n_seq)
            sys.stdout.flush()
            if ll >= th_min and ll <= th_max:
                oo.write('>' + description + '\n')
                oo.write(sequence + '\n')
                n_seq2  += 1
                n_bases2 += ll
                if ll > len_max2:
                    len_max2 = ll
                if ll < len_min2:
                    len_min2 = ll
                if num2 > 0 and n_seq2 == num2:
                    break
        in_handle.close()
        oo.close()
    else:
        print('The file is not in fasta or fastq format!')
        sys.exit(1)
    print('\n\n Original file:')
    print(' Reads number:     %d' % n_seq)
    print(' Bases number:     %d' % n_bases)
    print(' Average length:   %f' % (n_bases/n_seq))
    print(' Max length:       %d' % len_max)
    print(' Min length:       %d' % len_min)
    print('\n\n After extraction:')
    print(' Reads number:     %d' % n_seq2)
    print(' Bases number:     %d' % n_bases2)
    print(' Average length:   %f' % (n_bases2/n_seq2))
    print(' Max length:       %d' % len_max2)
    print(' Min length:       %d' % len_min2)
