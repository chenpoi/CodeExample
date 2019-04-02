#---------------------------------------------------------------------------------------------------
#download protein data from jaspar and uniprot database

import re
import urllib.error
import urllib.request

geneid = ['MA0002.2', 'MA0003.3', 'MA0004.1', 'MA0006.1', 'MA0007.3', 'MA0009.2', 'MA0014.3']
geneentry = []


def genesearch(gene):
    response = urllib.request.urlopen('http://www.uniprot.org/uniprot/?query={0}&sort=score'.format(gene))
    html = str(response.read())
    response.close()
    htmlsplit = html.split('_MOUSE')
    returnvalue = [htmlsplit[i][-24:-16] for i in range(0,len(htmlsplit)-1) if '><strong>{0}</strong>'.format(gene) in htmlsplit[i+1][350:850]]
    print(str(returnvalue[0]))
    return str(returnvalue[0])

def check(gene):
    response = urllib.request.urlopen('http://www.uniprot.org/uniprot/?query={0}&sort=score'.format(gene))
    # print(response.info())
    html = str(response.read())
    # do something
    response.close()  # best practice to close the file
    # file = open('urlq.txt','x')
    # file.write(str(html))
    # file.close()
    htmlsplit = html.split('_MOUSE')
    print(htmlsplit[1][350:850])
    print(htmlsplit[0][-24:-16])
    print(len(htmlsplit))
    print([htmlsplit[i][-24:-16] for i in range(0, 1) if '<strong>{0}</strong>'.format(gene) in htmlsplit[i + 1][350:850]])
    returnvalue = [htmlsplit[i][-24:-16] for i in range(0, len(htmlsplit) - 1) if '<span><strong>{0}</strong>'.format(gene) in htmlsplit[i + 1][350:850]]
    return str(returnvalue[0])

def jaspargenesearch(ID):
    response = urllib.request.urlopen('http://jaspar.genereg.net/matrix/{0}'.format(ID))
    html = str(response.read())
    response.close()
    # file = open('urlq.txt','x')
    # file.write(str(html))
    # file.close()
    htmlsplit = html.split('http://www.uniprot.org/uniprot/')
    if len(htmlsplit) > 1:
        returnvalue = htmlsplit[1][0:10]
    else:
        returnvalue = 'NONE'
    print(str(returnvalue))
    return str(returnvalue)

def main():
    for i in geneid :
        print(i)
        geneentry.append(re.sub('[^0-9A-Z]+','',jaspargenesearch(i)))
    file = open('jaspar.txt','x')
    for item in geneentry:
        file.write("%s\n" % item)
    file.close()

    for i in geneid :
        geneentry.append(re.sub('[^0-9A-Z]+','',genesearch(i)))
        print(i)
    file = open('urlmouse.txt','x')
    for item in geneentry:
      file.write("%s\n" % item)
    file.close()

if __name__ == '__main__':
    main()