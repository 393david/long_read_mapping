# report how best alignments from different methods cover each other for all records in SAM format by different methods
# all records aggrement.
import os, sys, argparse
from collections import Counter
from aligntools import Cigar

if __name__ == "__main__":
    usage = '''  Calculate the agreement between two SAM files for all records.

  Usage example:
  python ''' + sys.argv[0] + ' -i aligned1.sam -t aligned2.sam'

    parser = argparse.ArgumentParser(description = usage, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input the first sam file [required]')
    parser.add_argument('-t', '--target', type = str, required = True, help = 'The second sam file [required]')
    parser.add_argument('-p', '--overlap', type = float, default = 0.9, help = 'The overlap ration [0.9]')
    parser.add_argument('-q', '--mapq', type = int, default = 20, help = 'The mapping quality threshold [20]')
    args = parser.parse_args()
    p    = args.overlap
    sam1 = args.input
    sam2 = args.target
    mapq = args.mapq
    #if based not in [0, 1]: #!= 0 and based != 1:
    #    print(' Start position can be only 0 or 1!')
    #    sys.exit(1)
    #print()
    if not os.path.exists(sam1):
        print(' %s does not exist!' % sam1)
    if not os.path.exists(sam2):
        print(' %s does not exist!' % sam2)
    aligned_seq_header = []
    aligned_seq_positions = []
    aligned_seq_num = 0
    target_genome = []
    aligned_records_mapq = []
    ff = open(sam1)
    line = ff.readline()
    seq1 = ""
    seq2 = ""
    lengthIs = 0
    n1_records_than_mapq = 0
    n2_records_than_mapq = 0
    while line:
        line = line.strip()
        if line.startswith("@"):
            line = ff.readline()
            continue
        aa = line.split()
        aligned_is = aa[5]
        if aligned_is == "*":
            line = ff.readline()
            continue
        seq_header = aa[0]
        target_genome.append(aa[2])
        #startt = 0
        #if seq_header in aligned_seq_header: # just report the first aligned results
        #    line = ff.readline()
        #    continue
        aligned_seq_header.append(seq_header)
        aligned_records_mapq.append(int(aa[4]))
        if int(aa[4]) >= mapq:
            n1_records_than_mapq += 1
        position = int(aa[3]) - 1
        #seq_is = aa[9]
        aligned_seq_num += 1
        sys.stdout.write("Aligned records in SAM file 1: %d   \r" % aligned_seq_num)
        sys.stdout.flush()
        cigar_is = aa[5]
        cigar = Cigar.coerce(cigar_is)
        aligned_seq_positions.append(range(position, position + cigar.ref_length))
        line = ff.readline()
    ff.close()
    print()
    aligned_seq_num2 = 0
    aligned_seq_header2 = []
    aligned_seq_positions2 = []
    target_genome2 = []
    aligned_records_mapq2 = []
    ff2 = open(sam2)
    line = ff2.readline()
    while line:
        line = line.strip()
        if line.startswith("@"):
            line = ff2.readline()
            continue
        aa = line.split()
        aligned_is = aa[5]
        if aligned_is == "*":
            line = ff2.readline()
            continue
        seq_header = aa[0]
        target_genome2.append(aa[2])
        #startt = 0
        #if seq_header in aligned_seq_header2: # just report the first aligned results
        #    line = ff2.readline()
        #    continue
        aligned_seq_header2.append(seq_header)
        aligned_records_mapq2.append(int(aa[4]))
        if int(aa[4]) >= mapq:
            n2_records_than_mapq += 1
        position = int(aa[3]) - 1
        #seq_is = aa[9]
        aligned_seq_num2 += 1
        sys.stdout.write("Aligned records in SAM file 2: %d   \r" % aligned_seq_num2)
        sys.stdout.flush()
        cigar_is = aa[5]
        cigar = Cigar.coerce(cigar_is)
        aligned_seq_positions2.append(range(position, position + cigar.ref_length))
        line = ff2.readline()
    ff2.close()
    print()
    counter1 = Counter(aligned_seq_header)
    counter2 = Counter(aligned_seq_header2)
    num1 = 0
    num2 = 0
    num1_with_mapq = 0
    num2_with_mapq = 0
    iii = -1
    for ss in aligned_seq_header:
        iii += 1
        sys.stdout.write("Calculate the %d-th records in SAM file 1, total: %d   \r" % (iii + 1, aligned_seq_num))
        sys.stdout.flush()
        if ss in aligned_seq_header2:
            if counter2[ss] == 1:
                indexx = aligned_seq_header2.index(ss)
                if target_genome[iii] == target_genome2[indexx]:
                    pos1 = set(aligned_seq_positions[iii])
                    pos2 = set(aligned_seq_positions2[indexx])
                    intersection = pos1 & pos2
                    lenn = len(intersection)
                    if lenn*1.0/len(pos1) >= p:
                        num1 += 1
                        if aligned_records_mapq[iii] >= mapq and aligned_records_mapq2[indexx] >= mapq:
                            num1_with_mapq += 1
            else:
                startt = aligned_seq_header2.index(ss)
                flag1 = 0
                flag2 = 0
                for jj in range(startt, startt + counter2[ss]):
                    if target_genome[iii] == target_genome2[jj]:
                        pos1 = set(aligned_seq_positions[iii])
                        pos2 = set(aligned_seq_positions2[jj])
                        intersection = pos1 & pos2
                        lenn = len(intersection)
                        if lenn*1.0/len(pos1) >= p:
                            if flag1 == 0:
                                num1 += 1
                                flag1 = 1
                            if flag2 == 0:
                                if aligned_records_mapq[iii] >= mapq and aligned_records_mapq2[jj] >= mapq:
                                    num1_with_mapq += 1
                                    flag2 = 1
                            if flag1 ==1 and flag2 == 1:
                                break
    print()
    iii = -1
    for ss in aligned_seq_header2:
        iii += 1
        sys.stdout.write("Calculate the %d-th records in SAM file 2, total: %d   \r" % (iii + 1, aligned_seq_num2))
        sys.stdout.flush()
        if ss in aligned_seq_header:
            if counter1[ss] == 1:
                indexx = aligned_seq_header.index(ss)
                if target_genome[indexx] == target_genome2[iii]:
                    pos1 = set(aligned_seq_positions[indexx])
                    pos2 = set(aligned_seq_positions2[iii])
                    intersection = pos1 & pos2
                    lenn = len(intersection)
                    if lenn*1.0/len(pos2) >= p:
                        num2 += 1
                        if aligned_records_mapq[indexx] >= mapq and aligned_records_mapq2[iii] >= mapq:
                            num2_with_mapq += 1
            else:
                startt = aligned_seq_header.index(ss)
                flag1 = 0
                flag2 = 0
                for kk in range(startt, startt + counter1[ss]):
                    if target_genome[kk] == target_genome2[iii]:
                        pos1 = set(aligned_seq_positions[kk])
                        pos2 = set(aligned_seq_positions2[iii])
                        intersection = pos1 & pos2
                        lenn = len(intersection)
                        if lenn*1.0/len(pos2) >= p:
                            if flag1 == 0:
                                num2 += 1
                                flag1 = 1                            
                            if flag2 == 0:
                                if aligned_records_mapq[kk] >= mapq and aligned_records_mapq2[iii] >= mapq:
                                    num2_with_mapq += 1
                                    flag2 = 1
                            if flag1 ==1 and flag2 == 1:
                                break

    print('\n\n')
    print('Mapped recordes in SAM 1:         %d' % aligned_seq_num)
    print('Mapped recordes in SAM 2:         %d' % aligned_seq_num2)
    print('Recordes in SAM1 covered by SAM2: %f' % (num1 * 1.0 / len(aligned_seq_header)))
    print('Recordes in SAM2 covered by SAM1: %f' % (num2 * 1.0 / len(aligned_seq_header2)))
    print('\n\n')
    print('Mapping records larger than MAQ: %d'  % mapq)
    print('Mapped recordes in SAM 1:         %d' % n1_records_than_mapq)
    print('Mapped recordes in SAM 2:         %d' % n2_records_than_mapq)
    if n1_records_than_mapq > 0:
        print('Recordes in SAM1 covered by SAM2: %f' % (num1_with_mapq * 1.0 / n1_records_than_mapq))
    else:
        print('Recordes in SAM1 covered by SAM2: 0')
    if n2_records_than_mapq > 0:
        print('Recordes in SAM2 covered by SAM1: %f' % (num2_with_mapq * 1.0 / n2_records_than_mapq))
    else:
        print('Recordes in SAM1 covered by SAM2: 0')
