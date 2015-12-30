__author__ = 'Guoliang Lin'
Softwarename = 'sort_by_location_quik_version'
version = '2.0.1'
data = ""
bugfixs = '     version 2.0.1\n' \
          'Add GFF to filt out scafold'

import sys, getopt
import time


def trim(y):
    y = y.replace("[", '')
    y = y.replace(']', '')
    y = y.replace("',", '\t')
    y = y.replace("'", '')
    y = y.replace('\\n', '')
    y = y.replace(' ', '')
    # y=y.replace(',','')
    y = y.strip()
    y = y + '\n'
    return y


def trim_head(file, line=1):
    for i in range(0, line):
        file.readline()


def POSdetect(POS1, POS2):
    """
    :param POS1:
    :param POS2:
    :return:
    :rtype: int
    """
    if POS1.find(".") == -1:
        return int(POS1)
    else:
        assert isinstance(POS2, str)
        return int(POS2)


def GetBase(ALT1, ALT2):
    """


    :param ref:
    :param ALT1:
    :param ALT2:
    :return:
    :rtype: str
    """

    assert isinstance(ALT1, str)
    if ALT1.find('.') == -1:
        return ALT1
    else:
        return ALT2

def is_contain(x, segment):
    """

    :param x:
    :param segment:
    :return:
    :rtype: bool
    """
    assert isinstance(x, int)
    assert isinstance(segment, list)
    if int(segment[0]) < x and int(segment[1]) > x:
        return True
    else:
        return False


def typedet(listrefaalt,ref):
    return GetBase(listrefaalt[0], ref) + "->" + GetBase(listrefaalt[1], ref)

def f(a):
    if len(a)==2:
        return int(a[0])
    else:
        return POSdetect(a[0],a[1])


print('%s software version is %s' % (Softwarename, version))
print(bugfixs)
print('starts at :' + time.strftime('%Y-%m-%d %H:%M:%S'))
end=0
SegmentDict = {}
TypeDict = {}
snplist1=[]
Total = 0
opts, args = getopt.getopt(sys.argv[1:], 'i:g:s:h', ['inputfile=','GFF=', 'snp=', 'help'])
InputFileName = ''
for o, a in opts:
    if o in ['-i', '--inp,utfile']:
        InputFileName = a
    elif o in ['-s', '--snp']:
        snp = a
    elif o in ['-g','--GFF']:
        gff=a
    elif o in ['-h', '--help']:
        help = True
with open(InputFileName, 'r') as InputFile:
    with open(snp, 'r') as snpfile:
        with open(gff,'r') as gfffile:
            with open(InputFileName.split(".")[0] + ".diff-filter.out", 'w') as snpoutfile:
                with open(InputFileName.split(".")[0] + ".multivariation.out", 'w') as multiv:
                    with open("statistic.txt", 'w') as statistic:
                        statistic.write("Tpye\tNO.\tTotal\tRate\n")
                        snpoutfile.write("CHROM\tPOS1\tPOS2\tIN_FILE\tREF1\tREF2\tALT1\tALT2\tTYPE\n")
                        trim_head(InputFile)
                        for item in gfffile:
                            segmentlist=item.split('\t')
                            if not SegmentDict.has_key(segmentlist[0]):
                                SegmentDict[segmentlist[0]]=[]
                        for item in InputFile:
                            segmentlist = item.split("\t")
                            if SegmentDict.has_key(segmentlist[0]):
                                SegmentDict[segmentlist[0]].append(segmentlist[1:3])
                        for item in snpfile:
                            item=item.replace('\n', '')
                            snplist = item.split("\t")
                            if SegmentDict.has_key(snplist[0]):
                                SegmentDict[snplist[0]].append(snplist[1:])
                        for keys in SegmentDict.keys():
                            SegmentDict[keys].sort(key=f)
                            end=0
                            for snplist in SegmentDict[keys]:
                                if len(snplist)==2:
                                    end=int(snplist[1])
                                elif POSdetect(snplist[0],snplist[1])<end:
                                    tmplist=list(snplist[5])
                                    tmplist.reverse()
                                    if snplist[5]==snplist[6] or ''.join(tmplist)==snplist[6]:
                                        pass
                                    else:
                                        snplist.append(typedet(snplist[5:7],GetBase(snplist[3],snplist[4])))
                                        if snplist[-1] in TypeDict.keys():
                                            TypeDict[snplist[-1]] = TypeDict[snplist[-1]] + 1
                                        else:
                                            TypeDict[snplist[-1]] = 1
                                        if len(snplist[-1])==4:
                                            snpoutfile.write(keys+'\t'+trim(str(snplist)))
                                        else:
                                            multiv.write(keys+'\t'+trim(str(snplist)))
                                        Total = Total + 1
                        for key in TypeDict.keys():
                            statistic.write(key + '\t' + str(TypeDict[key]) + '\t' + str(Total) + '\t' + str(
                                TypeDict[key] * 100.0 / Total) + '\n')
print('ends at :' + time.strftime('%Y-%m-%d %H:%M:%S'))

