__author__ = 'Guoliang Lin'
Softwarename = 'splitGFF'
version = '1.0.2'
data=""
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


def set_to_type(list_item, file,string,genestring,list_type):
    """

    :param listitem:listitem
    :param file:file intergenic or intron
    """
    list_item[2] = string
    list_item[3] = str(int(list_type[-2][1]) + 1)
    list_item[4] = str(int(list_type[-1][0]) - 1)
    list_item[-1] = genestring
    file.write(trim(str(list_item)))
def convert_to_type(list_item,list_type,file,string,genestring):
    """

    :param list_item: listitem
    :param list_type: ExonsList or GeneList
    :param file: file intergenic or intron
    """
    tmp=list_item[3:5]
    list_type.append(tmp)
    if len(list_type)>2:
        list_type.pop(0)
        set_to_type(list_item,file,string,genestring,list_type)
    elif len(list_type)==2:
        set_to_type(list_item,file,string,genestring,list_type)
ExonsList= []
GenesList=[]
scaffold=''
gene=''
genelist=[]
print('%s software version is %s' % (Softwarename, version))
print(bugfixs)
print('starts at :' + time.strftime('%Y-%m-%d %H:%M:%S'))

opts, args = getopt.getopt(sys.argv[1:], 'i:t:h', ['inputfile=','type=', 'help'])
InputFileName = ''
for o, a in opts:
    if o in ['-i', '--inputfile']:
        InputFileName = a
    elif o in ['-t','--type']:
        etype=a
    elif o in ['-h', '--help']:
        help = True
with open(InputFileName, 'r') as InputFile:
    with open(InputFileName+'-'+etype,'w') as cds:
        with open(InputFileName+'-intron','w') as intron:
            with open(InputFileName+'-utr','w') as utr:
                with open(InputFileName+'-intergenic','w') as intergenic:
                    for item in InputFile:
                        listitem=item.split();
                        if len(listitem)!=0:
                            if scaffold!=listitem[0]:
                                scaffold=listitem[0]
                                GenesList=[]
                                ExonsList=[]
                            if listitem[2]==etype:
                                cds.write(trim(str(listitem)))
                            elif listitem[2].find('UTR')!=-1:
                                utr.write(trim(str(listitem)))
                            elif listitem[2]=='gene':
                                tmpname=gene
                                gene=listitem[-1]
                                convert_to_type(listitem,GenesList,intergenic,'intergenic',tmpname+';'+gene)
                                ExonsList=[]
                            if listitem[2]=='exon':
                                convert_to_type(listitem,ExonsList,intron,'intron',gene)
print('starts at :' + time.strftime('%Y-%m-%d %H:%M:%S'))