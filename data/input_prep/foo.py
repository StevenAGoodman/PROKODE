import re
import pandas as pd

def annotation():
    with open('genomic.gbff','r') as file:
        out = 'geneid\tstart\tend\tsynonyms\n'
        get_synonym = False
        spec = False
        prev_lin = 0
        for line in file:
            if spec:
                stri += ' ' + line.strip()[:-1]
                out += f"{stri.split('; ')}\n"
                spec = False
            if get_synonym:
                search2 = re.search('/gene_synonym=".+"', line)
                search3 = re.search('/gene_synonym="[^"]+\n', line)
                if search2 != None:
                    get_synonym = False
                    out += f'{search2.group()[15:-1].split('; ')}\n'
                elif search2 == None and search3 != None:
                    spec = True
                    stri = search3.group()[15:-1]


            search = re.search('/gene=".+"', line)
            if search != None:
                geneid = search.group()[7:-1]
                start = re.search('\s(complement\()?\d+\.\.\d+', prev_lin)
                if start == None:
                    continue
                start = start.group()
                end = start[start.find('..')+2:]
                start = start[12:start.find('..')] if start.find('nt(') >= 0 else start[1:start.find('..')]
                prev_lin = line
                out += f'{geneid}\t{start}\t{end}\t'
                get_synonym = True
                continue
            else:
                prev_lin = line
                continue


    out = out.split('\n')
    res = []
    foo_prev = 'na'
    for lin in out:
        foo = lin[:lin.find('\t')]
        if foo == foo_prev:
            foo_prev = foo
            continue
        else:
            res.append(lin + '\n')
            foo_prev = foo
            continue

    open('../../src/inputs/annotation.tsv', 'w').writelines(res)

def operons():
    annotation_df = pd.read_csv('../../operon_src/inputs/annotation.tsv', delimiter='\t')
    out = 'operonid\tstart\tgeneids\n'
    started = False
    genes = []
    with open("operons.txt", 'r') as operonsfile:
        for line in operonsfile:
            search1 = re.search('.+\t', line)
            if search1 == None:
                if started:
                    print(genes)
                    print(genes[0])
                    try:
                        start = annotation_df.loc[annotation_df['geneid']==genes[0], 'start'].values[0]
                        print(start)
                        out += f'{operonid}\t{start}\t{genes}\n'
                        genes = []
                    except:
                        genes = [] 
                operonid = re.search('\d+\n', line).group().replace('\n', '')
                print(operonid)
                started = True
                continue
            elif search1 != None and started:
                genes.append(search1.group()[1:search1.group()[1:].find('\t')+1].replace('\t',''))


    open('../../operon_src/inputs/operons.tsv', 'w').writelines(out)

operons()
