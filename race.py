#!/usr/bin/env python

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


def complement(N):
    NN = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return(NN[N])

def self_dimers(seq):
    rev = seq[::-1]
    start = 0
    end = 4
    while end != len(seq):
        rev1 = ''
        for i in rev[start:end]:
            rev1 += complement(i)
        if seq.find(rev1) == -1:
            start += 1
            end += 1
        else:
            start1 = seq.find(rev1)
            end1 = start1 + 4
            alseq = start * '-' + seq
            alref = start1 * '-' + rev
            if len(alseq) > len(alref):
                alref = alref + '-'*(len(alseq) - len(alref))
                scrstr = start*'-' + len(rev1) * '*' 
                scrstr = scrstr + '-'*(len(alseq)-len(scrstr))
                return(1)
            elif len(alseq) < len(alref):
                alseq = alseq + '-'*(len(alref) - len(alseq))
                scrstr = start1 *'-' + len(rev1) * '*' 
                scrstr = scrstr + '-'*(len(alseq)-len(scrstr))
                return(1)
            else:
                scrstr = start * '-' + '-' * start + len(rev1) * '*'
                scrstr = scrstr + '-'*(len(alseq)-len(scrstr))
                return(1)
    if end == len(seq):
        return(0)

def hairpin(seq):
    count = -4
    while count >= -((len(seq)-3)):
        end_seq = seq[count:]
        end_seq = end_seq[::-1]
        start_seq = seq[:-len(end_seq)]
        if len(end_seq) >= len(start_seq):
            end_seq, start_seq = start_seq, end_seq
        alig = SWalig(SWvlmtrx(end_seq, start_seq), end_seq, start_seq)
        if alig[2].find('****') >= 1 or alig[2].find('**-***') >= 1 or alig[2].find('***-**') >= 1:
            return(0)
        count -= 1
    return(alig)

def GC(pr):
    return((pr.count('G') + pr.count('C'))/len(pr)*100)

def GCend(pr):
    if pr[-5:].count('G') + pr[-5:].count('C') <= 3:
        return(1)
    else:
        return(0)
#Возвращает вес комплементарных нуклеотидов для матрицы
def Vlpoint(refN, seqN):
    global NN
    NN = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    m = 10
    if refN == NN[seqN]:
        return m
    else:
        return 0

def SWvlmtrx(REF, seq, p=-10):
    # Функция возвращает матрицу весов и таблицу пути в виде кортежа
    if len(REF) < len(seq):
        REF, seq = seq, REF
    nref = len(REF) + 1
    nseq = len(seq) + 1
    kref = 1
    k2 = 0
    hseq = 1
    c = {}
    d = {}
    for i in range(nref):
        d[i] = dict.fromkeys([j for j in range(nseq)], 0)
    for j in range(len(REF)):
        c[j] = dict.fromkeys([i for i in range(len(seq))], 0)
    for j in seq:
        d[0][hseq] = hseq * 0
        hseq += 1
    hseq = 0
    for i in REF:
        for j in seq:
            d[kref][0] = kref * 0
            match = d[k2][hseq] + Vlpoint(i, j)

            c2 = d[k2 + 1][(hseq)] + p

            c3 = d[k2][hseq + 1] + p

            if k2+1 < len(c):
                if match > c2 and match > c3:
                    c[k2+1][hseq+1] = k2, hseq
                elif c2 > match and c2 > c3:
                    c[k2+1][hseq+1] = k2 + 1, hseq
                elif c3 > match and c3 > c2:
                    c[k2+1][hseq+1] = k2, hseq + 1
                else:
                    if match == c2:
                        c[k2+1][hseq+1] = k2, hseq, k2 + 1, hseq
                    elif match == c3:
                        c[k2+1][hseq+1] = k2, hseq, k2, hseq + 1
                    elif c2 == c3:
                        c[k2+1][hseq+1] = k2 + 1, hseq, k2, hseq + 1
                    else:
                        c[k2+1][hseq+1] = k2, hseq, k2 + 1, hseq, k2, hseq + 1
            else:
                break
            d[kref][hseq + 1] = max(match, c2, c3)
            hseq += 1
        k2 += 1
        kref += 1
        hseq = 0
    return d, c

def SWalig(matrx, REF, seq):
    if len(REF) < len(seq):
        REF, seq = seq, REF
    #Функция возвращает выравненные последовательности
    alref = ''
    alseq = ''
    scrstr = ''
    d = matrx[0]
    c = matrx[1]
    corr_R = 0
    corr_s = 0
    #Поиск начальной точки сбора выравниваия в словаре c на выходе получаем координаты отсчёта
    max_into_d = 0
    counter = 10
    for R in d:
        if R == 0:
            continue
        else:
            for s in d[R]:
                if s == 0:
                    continue
                else:
                    if (d[R][s]) > max_into_d:
                        max_into_d = d[R][s]
                        coord_R = R-1
                        coord_s = s-1
                    elif (d[R][s]) == 0 and counter == 0:
                        coord_R = len(REF)-1
                        coord_s = len(seq)-1 
                        break                 
                    else:
                        counter -= 1
                        if counter == -1:
                            break
                        else:
                            continue
    alref = REF[coord_R] + alref
    alseq = seq[coord_s] + alseq
    if Vlpoint(REF[coord_R], seq[coord_s]) > 0:
        scrstr = '*' + scrstr
    else:
        scrstr = '-' + scrstr
    
    #Цикл собирает выравненные участки последовательностей
    while coord_R != 0 and coord_s != 0:
      if coord_s == 0 and coord_R == 0:
            alref = REF[coord_R] + alref
            alseq = seq[coord_s] + alseq
            if Vlpoint(REF[coord_R], seq[coord_s]) > 0:
                scrstr = '*' + scrstr
            else:                    
                scrstr = '-' + scrstr
            coord_R = 0
            coord_s = 0
      else:
        if len(c[coord_R][coord_s]) == 2:
            if coord_R-1 == c[coord_R][coord_s][0] and coord_s-1 == c[coord_R][coord_s][1]:
                coord_R, coord_s = c[coord_R][coord_s][0], c[coord_R][coord_s][1]
                alref = REF[coord_R] + alref
                alseq = seq[coord_s] + alseq

                if Vlpoint(REF[coord_R], seq[coord_s]) > 0:
                    scrstr = '*' + scrstr
                else:
                    scrstr = '-' + scrstr
            elif coord_R == c[coord_R][coord_s][0] and coord_s-1 == c[coord_R][coord_s][1]:
                coord_R, coord_s = c[coord_R][coord_s][0], c[coord_R][coord_s][1]
                alref = '-' + alref
                corr_R += 1
                alseq = seq[coord_s] + alseq
                scrstr = '-' + scrstr
            elif  coord_R-1 == c[coord_R][coord_s][0] and coord_s == c[coord_R][coord_s][1]:
                coord_R, coord_s = c[coord_R][coord_s][0], c[coord_R][coord_s][1]
                alref = REF[coord_R] + alref
                alseq = '-' + alseq
                corr_s += 1
                scrstr = '-' + scrstr
        elif len(c[coord_R][coord_s]) == 4: 
            if  coord_R-1 == c[coord_R][coord_s][0] and coord_s-1 == c[coord_R][coord_s][1] or coord_R-1 == c[coord_R][coord_s][2] and coord_s-1 == c[coord_R][coord_s][3]:
                if coord_R-1 == c[coord_R][coord_s][0]:
                    coord_R, coord_s = c[coord_R][coord_s][0], c[coord_R][coord_s][1]
                elif coord_R-1 == c[coord_R][coord_s][2]:
                    coord_R, coord_s = c[coord_R][coord_s][2], c[coord_R][coord_s][3]               
                alref = REF[coord_R] + alref
                alseq = seq[coord_s] + alseq

                if Vlpoint(REF[coord_R], seq[coord_s]) > 0:
                    scrstr = '*' + scrstr
                else:
                    scrstr = '-' + scrstr
            else:
                if coord_R == c[coord_R][coord_s][0] and coord_s-1 == c[coord_R][coord_s][1] or coord_R == c[coord_R][coord_s][2] and coord_s-1 == c[coord_R][coord_s][3]:
                    if coord_R == c[coord_R][coord_s][0]:
                        coord_R, coord_s = c[coord_R][coord_s][0], c[coord_R][coord_s][1]
                    elif coord_R == c[coord_R][coord_s][2]:
                        coord_R, coord_s = c[coord_R][coord_s][2], c[coord_R][coord_s][3]                    
                    alref = '-' + alref
                    corr_R += 1
                    alseq = seq[coord_s] + alseq
                    scrstr = '-' + scrstr
                elif  coord_R-1 == c[coord_R][coord_s][0] and coord_s == c[coord_R][coord_s][1] or coord_R-1 == c[coord_R][coord_s][2] and coord_s == c[coord_R][coord_s][3]:
                    if coord_R == c[coord_R][coord_s][0]:
                        coord_R, coord_s = c[coord_R][coord_s][0], c[coord_R][coord_s][1]
                    elif coord_R == c[coord_R][coord_s][2]:
                        coord_R, coord_s = c[coord_R][coord_s][2], c[coord_R][coord_s][3]                    
                    alref = REF[coord_R] + alref
                    alseq = '-' + alseq
                    corr_s += 1
                    scrstr = '-' + scrstr
    #Последний кусок который добавляет учаски начала и конца последовательности если таковые не выравненны
    if len(alref)-corr_R < len(REF):
        alseq_const = len(alseq)
        alref_const = len(alref)
        alref = REF[:coord_R] + alref + REF[alref_const+coord_R-corr_R:]
        alseq = seq[:coord_s] + alseq + seq[alseq_const+coord_s-corr_s:]
        scrstr = '-'*max(len(REF[:coord_R]), len(seq[:coord_s])) + scrstr
        if coord_R != 0:
            if len(alseq)-corr_s < len(seq):
                alseq = seq[:coord_s] + alseq + seq[len(alseq)+coord_s:]
            else:
                alseq = '-'*len(REF[:coord_R]) + alseq
                alseq = alseq + '-'* (len(alref)-len(alseq))
            scrstr = scrstr + '-'*len(REF[alref_const+coord_R:])
        else:
            alref = '-'*len(seq[:coord_s]) +alref
            alseq_const = len(alseq)
            scrstr_const = len(alseq)-alseq_const
            alseq =  alseq + '-'*(len(alref)-len(alseq))
            alref = alref + '-' * (len(alseq)-len(alref))
            scrstr = scrstr + '-' * len(REF[alref_const-corr_R:])
    return(alref, alseq, scrstr)

def nested(seq):
    shrt = 18
    lng = 25
    pr = ''
    maxa = 0    
    ok = 0
    while ok != 1:
        if shrt == lng+1:
            shrt = 18
            seq = seq[1:]
        pr = seq[:shrt] 
        if GC(pr) >= 60 or GC(pr) <= 40:
            shrt += 1
            continue
        else:
            if GCend(pr) == 0:
                shrt += 1 
            else:
                if mt.Tm_Wallace(pr) > 60 or mt.Tm_Wallace(pr) < 56:
                    shrt += 1
                else:
                    pr_revers = pr[::-1]
                    complalig = SWalig(SWvlmtrx(pr, pr_revers), pr, pr_revers)
                    if complalig[2].find("****") >= 0  or self_dimers(pr) == 1 or complalig[2].find("***-**") >= 0 or complalig[2].find("**-***") >= 0:
                        shrt += 1
                    else: 
                        hpin = hairpin(pr)
                        if hpin == 0:
                            shrt += 1
                        else:
                            return pr

def compare(pr1, pr2):
    self_dimers = SWalig(SWvlmtrx(pr1, pr2[::-1]), pr1, pr2[::-1])
    if self_dimers[2].find('****') == -1 and self_dimers[2].find('***-**') == -1 and self_dimers[2].find('**-***') == -1:
        return(1)
    else:
        return(0)

def primers(end3, end5):
    list_prF = []
    list_prR = []
    re_end3 = end3
    re_end5 = end5
    n = 10
    maxa = 0
    while maxa != n:
        try:
            prF = nested(re_end3)
            prR = nested(re_end5)
            list_prF.append(prF)
            list_prR.append(prR)
            re_end3 = re_end3[re_end3.find(prF)+len(prF):]
            re_end5 = re_end5[re_end5.find(prR)+len(prR):]
            maxa += 1
        except:
            n -= 1
            

    for i in list_prR:
        for j in list_prF:
            if abs(mt.Tm_Wallace(i) - mt.Tm_Wallace(j)) > 3:
               continue
            else:
                if compare(i, j) == 0:

                    continue
                else:
                    prF = j 
                    prR = i
                    if len(end3[end3.find(prF)+len(prF):])*4 < len(end5[end5.find(prR)+len(prR):]):
                        continue
                    elif len(end3[end3.find(prF)+len(prF):]) > 4*len(end5[end5.find(prR)+len(prR):]):
                        continue
                    else:
                        end3 = end3[end3.find(prF)+len(prF)-3:]
                        end5 = end5[end5.find(prR)+len(prR)-3:]
                        return(prF, prR, end3, end5)

def rcr_vis(prF1, prF2, prF3, prR1, prR2, prR3, end3, end5):
    prR1 = prR1[::-1]
    prR2 = prR2[::-1]
    prR3 = prR3[::-1]
    prF1_vis = ' ' * (len(end3 [:end3.find(prF1)]) - 3) + "5'-" + prF1 + "-3' --->" 
    prF2_vis = ' ' * (len(end3 [:end3.find(prF2)]) - 3) + "5'-" + prF2 + "-3' --->"  
    prF3_vis = ' ' * (len(end3 [:end3.find(prF3)]) - 3) + "5'-" + prF3 + "-3' --->"  
    prR1_vis = ' ' * (len(end5 [:end5.find(prR1)]) - 8) + "<--- 3'-" + prR1 + "-5'" 
    prR2_vis = ' ' * (len(end5 [:end5.find(prR2)]) - 8) + "<--- 3'-" + prR2 + "-5'"
    prR3_vis = ' ' * (len(end5 [:end5.find(prR3)]) - 8) + "<--- 3'-" + prR3 + "-5'"
    return(str('' + prF1_vis +  '\n' + prF2_vis + '\n' + prF3_vis + '\n' + end3 + '\n' + end5 + '\n' + prR3_vis + '\n' + prR2_vis + '\n' + prR1_vis) + '\n')

end3 = 'tggctgaccAtAgGaCgAgAGCttCCTGGtgaaGtgtgTTtCTtgaAATCatCACcaCCATGGACAgCAaaGGtTcGTCgcagaaAGGGTCCcGCCTGCTCCTGCTGCTGGTggTGTCAAATCTaCTCTTGTGCCAgGGtgTggTcTcCaCccCCgTCTGTCCcaAtGGGccTGgCaaCTGCCaGGtaTCCCTTcgAGACCTGTTTGaCCggGCagtCaTggTGTCCcaCtaCATccATgACctCTCCtcggAAATGTTCAacGAaTTTGATaAACgGTATGCCCAGGGcAAAGgGTtCaTTAcCAtgGcCctCaACAgcTGCCAtACCtccTCCCTtCcTacCCCtGAAGAtAaagAaCAAgcCCaACaGAcccACAtgAAGtCCTTATgAgcTtGATtCTtggGTTgCTgcgCTCCTGGaAtgacCCTCTgTATCAcCTAGTCACcGAggTgCGGgGTATGAAAGgAGcCcCAgATgCTATCCTATCgAGgGCcAtAGAgAtTGaGgAAgAAaacAAAcgACTTCtgGaAggCATgGAGAtGataTTtgGCCAGGTTATTccTggAGCCAAaGAgacTGagCCctaCccTgtgTGGTcaGGACTCCcgTCCCTgCaaaCtAaggATGAAGATgcaCGTtATTCTGCtTTTTATAaCCTGcTCcaCTGCCTGCGCAGgGATtCaaGcAAgaTTGACAcTTACcttAAGcTCcTGaatTGCaGAaTcATCTacaAcAAcAAcTGcTAAgcCCACATccATcCtATCCAttTCTGAGAtGGtTCtTAATGATCcATtCCcTgGCAaacttctctgagctttatAGCTTTgTAATGCATGCTtgGcTctAATGGGTtTCaTCTTAAATAAAaACAgAcTCTGTAGcGATGTCAAAAtct'



n = int((len(end3))/3)*2-50

end3 = end3.upper()
end3 = Seq(end3, IUPAC.unambiguous_dna)
end5 = end3.reverse_complement()

pcr_end3 = end3
pcr_end5 = end5[::-1]

end5 = end5[:n]
end3 = end3[:n]

primers_and_seq = primers(end3, end5)
prF1 = primers_and_seq[0]
prR1 = primers_and_seq[1]
end3 = primers_and_seq[2]
end5 = primers_and_seq[3]


print(' >Прямой праймер 1: ', '\n', prF1, ': Темпертура отжига', '%0.2f' % mt.Tm_Wallace(prF1), '\n', ">Обратный праймер 1: ", '\n',prR1, ': Темпертура отжига', '%0.2f' % mt.Tm_Wallace(prR1))


primers_and_seq1 = primers(end3, end5)
prF2 = primers_and_seq1[0]
prR2 = primers_and_seq1[1]
end3 = primers_and_seq1[2]
end5 = primers_and_seq1[3]

print(' >Прямой праймер 2: ', '\n', prF2, ': Темпертура отжига', '%0.2f' % mt.Tm_Wallace(prF2), '\n', ">Обратный праймер 2: ", '\n',prR2, ': Темпертура отжига', '%0.2f' % mt.Tm_Wallace(prR2))

primers_and_seq2 = primers(end3, end5)
prF3 = primers_and_seq2[0]
prR3 = primers_and_seq2[1]
end3 = primers_and_seq2[2]
end5 = primers_and_seq2[3]

print(' >Прямой праймер 3: ', '\n', prF3, ': Темпертура отжига', '%0.2f' % mt.Tm_Wallace(prF3), '\n', ">Обратный праймер 3: ",'\n',prR3, ': Темпертура отжига', '%0.2f' % mt.Tm_Wallace(prR3))

with open('file.txt', 'a') as file:
    file.write(rcr_vis(prF1, prF2, prF3, prR1, prR2, prR3, pcr_end3, pcr_end5))
