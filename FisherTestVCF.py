import scipy
import getopt
import sys
from scipy.stats import fisher_exact
def as_fisher(line='',pheno1=63,pheno2=215,refPos=3,altPos=4,colRem=5,colPheno=9,callrate=0.9,sep='\t',showGT=True):
    'Fisher exact test in association study such GWAS'
    line_s=line.strip().replace('|','/')
    info=line_s.split(sep)
    gt_l=info[:colRem]
    rg=range(colPheno,colPheno+pheno1+pheno2)
    refs=info[refPos]
    alts=info[altPos]
    for index in rg:
        gt_i=info[index].split(':')[0]
        if showGT:
            gt_i=gt2str(gt_i,refs,alts)
        gt_l.append(gt_i)
    missRateCase=''.join(gt_l[colRem:colRem+pheno1]).count('./.')/float(pheno1)
    missRateCtrl=''.join(gt_l[colRem+pheno1:]).count('./.')/float(pheno2)
    if (1-missRateCase)<callrate:
        return('Too low calling rate!')
    if not showGT:
        refs='0'
        alts='1'
    pt1Ref=''.join(gt_l[colRem:colRem+pheno1]).count(refs)
    pt1Alt=''.join(gt_l[colRem:colRem+pheno1]).count(alts)
    pt2Ref=''.join(gt_l[colRem+pheno1:]).count(refs)
    pt2Alt=''.join(gt_l[colRem+pheno1:]).count(alts)
    cont_mat=[[pt1Alt,pt2Alt],[pt1Ref,pt2Ref]]
    gt_l.extend([str(pt1Alt),str(pt2Alt),str(pt1Ref),str(pt2Ref),str(1-missRateCase),str(1-missRateCtrl),str(fisher_exact(cont_mat)[1])])
    return(gt_l)
def gt2str(gt,ref='',alt=''):
    if gt=='0/0':
        gt_s='%s%s'%(ref,ref)
    elif gt=='1/1':
        gt_s='%s%s'%(alt,alt)
    elif gt=='./.':
        gt_s='./.'
    else:
        gt_s='%s%s'%(ref,alt)
    return(gt_s)
def usage():
    pass
def recParams(params):
    try:
        opts,args=getopt.getopt(params,'i:p:P:r:a:R:s:f:d:D:',['vcf=','pheno1=','pheno2=','refPos=','altPos=','colRem=','colPheno=','callrate=','sep=','showGT='])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)
    if not opts:
        usage()
        sys.exit()
    for opt, value in opts:
        if opt in ('-i','--vcf'):
            vcf=value
        elif opt in ('-p','--pheno1'):
            pheno1=int(value)
        elif opt in ('-P','--pheno2'):
            pheno2=int(value)
        elif opt in ('-r','--refPos'):
            refPos=int(value)
        elif opt in ('-a','--altPos'):
            altPos=int(value)
        elif opt in ('-R','--colRem'):
            colRem=int(value)
        elif opt in ('-s','--colPheno'):
            colPheno=int(value)
        elif opt in ('-f','--callrate'):
            callrate=float(value)
        elif opt in ('-d','--sep'):
            sep=value
        elif opt in ('-D','--showGT'):
            showGT=eval(value)
        else:
            raise getopt.GetoptError('Invalid options')
    arg=locals()
    #for key in tuple(locals().keys())[5::-1]:
    #    arg[key]=vars()[key]
    arg.pop('args')#remove the odd arg item
    arg.pop('params')
    arg.pop('value')
    arg.pop('opts')
    arg.pop('opt')
    return(arg)
if __name__=='__main__':
    params=recParams(sys.argv[1:])
    print(params)
    if not set(['vcf','pheno1','pheno2']).issubset(list(params.keys())):
        raise Exception('VCF files and number of phenotypes must bu supplied')
    vcf=open(params['vcf'],'r').readlines()
    params.pop('vcf')
    for line in vcf:
        if line.startswith('#'):
            continue
        res=as_fisher(line,**params)
        #gtl='\t'.join(res[0])
        print('\t'.join(res))
  #  vcf=open('example','r').readlines()

#    for line in vcf:
 #       res=as_fisher(line,63,215,showGT=True,callrate=.1)
  #      print(res[1])


