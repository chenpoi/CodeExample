import time
import os
import re
import pandas as pd
import numpy as np
import math
#读取蛋白对应dna结合位点的ChIP-seq频数表
sequence = {}
filepath = '/Users/mengqichen/maybe_useful/zheda/data/'
foldernames = os.listdir(path=filepath)[1:]
for i in foldernames:
      filenames = os.listdir(path=filepath+i)   # string adding will concatenate into a longer one
      for j in filenames:
            fo = open(filepath+i+'/'+j,'r')
            data = fo.readlines()
            sequence[re.sub('.jaspar','',j)] = pd.DataFrame([[int(s) for s in re.split(' +',data[i])[2:-1]] for i in range(1,5)],index=('a','c','g','t'))

#annotation
annotationpath = '/Users/mengqichen/maybe_useful/zheda/annotation/'
annotation = pd.read_csv(annotationpath+'JASPAR-Browse Collection core.csv', encoding='latin-1')
#必须要加上这encoding才能够正常读取这个csv文件。
vertebratesID = sequence.keys()
annotation = annotation[pd.Series([i in vertebratesID for i in annotation.ID])]
annotation = pd.concat([annotation,pd.read_csv(annotationpath+'JASPAR-Browse Collection pbm.csv', encoding='latin-1'),pd.read_csv(annotationpath+'JASPAR-Browse Collection pbm_hlh.csv', encoding='latin-1'),pd.read_csv(annotationpath+'JASPAR-Browse Collection pbm_homeo.csv',encoding='latin-1')],axis=0)
jasparpath = '/Users/mengqichen/maybe_useful/zheda/jaspar.txt'
fo = open(jasparpath,'r')
annotation['jaspar'] = fo.readlines()


#构建pwm
#function for working out the position weight matrix value

def pwm(freq, bg=0.25):
    #using the formulae above
    total = freq.sum(axis=0)
    p = (freq + (np.sqrt(total) * 1/4)) / (total + (4 * (np.sqrt(total) * 1/4)))
    return pd.DataFrame(np.log2(p/bg),index = ('a','c','g','t'))

pwmsequence = {}
pwmcolmax = {}
pwmsum = {}
bestseq = {}

for i in sequence.keys():
    pwmsequence[i] = pwm(sequence[i])
    pwmcolmax[i] = np.max(pwmsequence[i], axis=0)
    pwmsum[i] = np.sum(pwmcolmax[i])
    bestseq[i] = [np.argmax(pwmsequence[i][j]) for j in range(0, pwmsequence[i].shape[1])]

pwmlength={}
for i in range(1,23):
    pwmlength[i] = [j for j in pwmsequence.keys() if pwmsequence[j].shape[1] == i]
#和R不同，这里保存的是motif名字而不是编号

#检验一个序列的pwm分数



def seqsearch(seq):
    # timing = {}
    # timing['secondtime'] = list()
    # timing['thirdtime'] = list()
    starttime = time.time()  #计时用
    seq = re.sub('[^A-Za-z]','',seq).lower()
    seq = re.sub('u','t',seq)    #两步预处理
    x = list(seq)     #讲string转化为character的list
    outputvalue = []   #保存结果
    last = -100   #last-check系统用于快速筛选，这一步是last的初始化
    for i in range(6,min(len(x),23)): #最短的识别序列是6 bases.最长的识别序列是23 bases
        for j in range(0,(len(x)-i)):
            y = x[j:(j+i)]
            if (last < -2) & (j != 0):
                last = last + 1  #如果last过小，说明这部分区域相当不符合motif，可以向后跳过
                continue
            last = -100    #last重新初始化
            scorelevel = pd.Series()  #scorelevel用于存储相关motif的pwm计算结果
            length = len(y)
            for k in pwmlength[length]:
                #secondtime = time.time()
                check = sum([x == y for x, y in zip(y, bestseq[k])]) / length  #check意为sequence中有多少位是和bestseq一致的
                #timing['secondtime'] = timing['secondtime'] + [time.time() - secondtime]
                #thirdtime = time.time()
                if check <= 0.5:
                    last = max(last, length * ( check - 0.5 )) #阈值2为0.5,意为如果有一半的位点不是最优的话，那就直接跳过。
                    #timing['thirdtime'] = timing['thirdtime'] + [time.time() - thirdtime]
                    continue
                #timing['thirdtime'] = timing['thirdtime'] + [time.time() - thirdtime]
                scorelevel[k] = sum([pwmsequence[k].loc[y[l],l] for l in range(0,length)])/pwmsum[k]
                # python中对于index的操作很麻烦。。。有没有简单的方法————————————化成pd.series可以方便的加上index
            #scorelevel = pd.Series([sum(scorelist[l])/pwmsum[l] for l in scorelist.keys()],index=scorelist.keys())
            passcore = scorelevel[scorelevel > -0.25/4913*(i-5)**3+0.9]   #阈值1为-0.25/4913*(i-5)^3+0.9，超过这个阈值就可以输出了
            if len(passcore) > 0:
                outputvalue = outputvalue + [[l, list(annotation.loc[annotation.ID == l,'Name'].values),list(annotation.loc[annotation.ID == l,'jaspar'].values), passcore[l],''.join(y)] for l in passcore.keys()]
                #cat(sprintf('the TF and sequence score are %s\t%s\nthe partial sequence is %s, from %s base to %s base\n',names(passcore),passcore,paste(y, collapse = ""), j, j+i))
    print(time.time()-starttime)
    if len(outputvalue) == 0:
        return('No TFs Meeting Requirement!')
    else:
        return(outputvalue)