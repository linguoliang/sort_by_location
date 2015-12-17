__author__ = 'Guoliang Lin'
Softwarename = 'sort_by_location'
version = '1.0.1'
data = ""
bugfixs = ''

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



print('%s software version is %s' % (Softwarename, version))
print(bugfixs)
print('starts at :' + time.strftime('%Y-%m-%d %H:%M:%S'))

SegmentDict = {}
TypeDict = {}
snplist1=[]
Total = 0
opts, args = getopt.getopt(sys.argv[1:], 'i:s:h', ['inputfile=', 'snp=', 'help'])
InputFileName = ''
for o, a in opts:
    if o in ['-i', '--inp,utfile']:
        InputFileName = a
    elif o in ['-s', '--snp']:
        snp = a
    elif o in ['-h', '--help']:
        help = True
with open(InputFileName, 'r') as InputFile:
    with open(snp, 'r') as snpfile:
        with open(InputFileName.split(".")[0] + ".diff-filter.out", 'w') as snpoutfile:
            with open("statistic.txt", 'w') as statistic:
                statistic.write("Tpye\tNO.\tTotal\tRate\n")
                snpoutfile.write("CHROM\tPOS1\tPOS2\tIN_FILE\tREF1\tREF2\tALT1\tALT2\tTYPE\n")
                trim_head(InputFile)
                for item in InputFile:
                    segmentlist = item.split("\t")
                    if segmentlist[0] in SegmentDict.keys():
                        SegmentDict[segmentlist[0]].append(segmentlist[1:3])
                    else:
                        SegmentDict[segmentlist[0]] = []
                        SegmentDict[segmentlist[0]].append(segmentlist[1:3])
                for item in snpfile:
                    item=item.replace('\n', '')
                    snplist = item.split("\t")
                #     snplist1.append(snplist)
                # for snplist in snplist1:
                    if snplist[0] in SegmentDict.keys():
                        listseg = SegmentDict[snplist[0]]
                        for element in listseg:
                            if is_contain(POSdetect(snplist[1], snplist[2]), element):
                                tmplist=list(snplist[6])
                                tmplist.reverse()
                                if snplist[6]==snplist[7] or ''.join(tmplist)==snplist[7]:
                                    pass
                                else:
                                    snplist.append(typedet(snplist[6:8],GetBase(snplist[4],snplist[5])))
                                    if snplist[-1] in TypeDict.keys():
                                        TypeDict[snplist[-1]] = TypeDict[snplist[-1]] + 1
                                    else:
                                        TypeDict[snplist[-1]] = 1
                                    snpoutfile.write(trim(str(snplist)))
                                    Total = Total + 1
                                break;
                for key in TypeDict.keys():
                    statistic.write(key + '\t' + str(TypeDict[key]) + '\t' + str(Total) + '\t' + str(
                        TypeDict * 100.0 / Total) + '\n')
print('starts at :' + time.strftime('%Y-%m-%d %H:%M:%S'))
