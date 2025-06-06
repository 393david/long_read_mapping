# Compute the accuracy for inversion
# 
import os, sys, argparse
from aligntools import Cigar

def strand_confirm(flag):
    # 1 plus and 0 reverse
    match flag:
        case 0 | 1 | 256 | 2048:
            return '+'
        case 16 | 272 | 2064:
            return '-'
        case _:
            print('SAM flag invalid: %d!' % flag)
            sys.exit(1)

if __name__ == "__main__":
    usage = ''' Compute the accuracy for simulated inversions.

  Usage example:
  python ''' + sys.argv[0] + ' -i mapped.sam -b inversion.groundtruth -o error_mapped.txt'

    parser = argparse.ArgumentParser(description = usage, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input the SAM file [required]')
    parser.add_argument('-b', '--groundtruth', type = str, required = True, help = 'The .groundtruth file [required]')
    parser.add_argument('-o', '--output', type = str, required = True, help = 'The output file not or wrong mapped[required]')
    #parser.add_argument('-o2', '--output2', type = str, required = True, help = 'The output file contained the output reads information[required]')
    parser.add_argument('-p', '--overlap', type = float, default = 0.1, help = 'The overlap ration that reads cover inversion [0.1]')
    args = parser.parse_args()
    p    = args.overlap
    sam  = args.input
    ground  = args.groundtruth
    outt = args.output
    #outt2= args.output2 
    if not os.path.exists(sam):
        print(' %s does not exist!' % sam)
        sys.exit(1)
    if not os.path.exists(ground):
        print(' %s does not exist!' % ground)
        sys.exit(1)
    if p < 0.0 or p > 1.0:
        print('p = %f, p should 0 < p <= 1' % p)
        sys.exit(1)
    seq_header      = []
    seq_inv_len     = {} # the inversion length in this read
    truth_pos_start = []
    truth_pos_end   = []
    truth_strand    = []
    truth_length    = []
    truth_set = []
    header_mapped  = {}
    n1 = 0
    ff = open(ground)
    line = ff.readline()
    line = ff.readline()
    line = ff.readline()
    while line:
        line = line.strip()
        aa = line.split()
        seq_header.append(aa[0])
        header_mapped[aa[0]] = 0
        truth_strand.append(aa[2])
        truth_pos_start.append(int(aa[5]))
        truth_pos_end.append(int(aa[6]))
        truth_length.append(int(aa[6]) - int(aa[5]) + 1)
        truth_set.append(set(range(int(aa[5]), int(aa[6]) + 1)))
        seq_inv_len[aa[0]] = int(aa[6]) - int(aa[5]) + 1
        n1 += 1
        #sys.stdout.write("Inversion in bed file: %d   \r" % n1)
        #sys.stdout.flush()
        line = ff.readline()
    ff.close()
    print()
    print('Sequences with inversion: %d' % n1)
    num2 = 0
    covered_inversion = []
    ff2 = open(sam)
    line = ff2.readline()
    ll = 0
    covered_num = 0
    mapped_seq_header = []
    while line:
        if line.startswith('@'):
            line = ff2.readline()
            continue
        ll += 1
        bb = line.strip().split()
        cigar_is = bb[5]
        if cigar_is == "*":
            line = ff2.readline()
            continue
        if bb[0] not in mapped_seq_header:
            mapped_seq_header.append(bb[0])
        if header_mapped[bb[0]] == 1:
            line = ff2.readline()
            continue
        aligned_start = int(bb[3])
        strandd = strand_confirm(int(bb[1]))
        indexx = seq_header.index(bb[0])
        if strandd == truth_strand[indexx]:
            line = ff2.readline()
            continue
        cigar = Cigar.coerce(cigar_is)
        aligned_end = aligned_start + cigar.ref_length - 1
        indexx = seq_header.index(bb[0])
        set_now = set(range(aligned_start, aligned_end + 1))
        set_true = set(range(truth_pos_start[indexx], truth_pos_end[indexx] + 1))
        intersection = set_now & set_true
        lenn = len(intersection)
        if lenn*1.0/len(set_true) >= p:
            covered_num += 1
            header_mapped[bb[0]] = 1
        line = ff2.readline()
        sys.stdout.write("Checking the %d -th mapped read, covered number: %d   \r" % (len(mapped_seq_header), covered_num))
        sys.stdout.flush()
    ff2.close()
    oo = open(outt, 'w')
    for key, value in header_mapped.items():
        if key in mapped_seq_header and value == 0:
            oo.write(key + '\t' + 'Wrong_mapped' + '\t' + str(seq_inv_len[key]) + '\n')
        if key not in mapped_seq_header:
            oo.write(key + '\t' + 'Not_mapped' + '\t' + str(seq_inv_len[key]) + '\n')
    oo.close()
    print('\n\n')
    print('Total reads:         %d' % n1)
    print('Mapped reads:        %d' % len(mapped_seq_header))
    print('Correctly mapped:    %d' % covered_num)
    print('Not or wrong mapped: %d' % (n1 - covered_num))
    print('Accuracy:            %f' % (1.0 * covered_num / n1))
    print('Precision:           %f' % (1.0 * covered_num / len(mapped_seq_header)))
    #print('Recordes in SAM1 covered by SAM2: %f' % (num1 * 1.0 / len(aligned_seq_header)))
    #print('Recordes in SAM2 covered by SAM1: %f' % (num2 * 1.0 / len(aligned_seq_header2)))
