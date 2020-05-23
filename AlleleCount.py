""" This python script is applied with counting the number of allele
    counts in a standard vcf
    Now only biallelic locus is supported"""
import getopt
import sys
import os
def countAllele(line = '', samGenotypeStartsCol = 9, showGT = False, delim = '\t'):
    """ Count the allele number on each locus """
    if line.strip().startswith('#'):
        return ['1','1'] 
    line_splited = line.strip().split(delim)
    genotype_info_concatenate = ''.join(line_splited[samGenotypeStartsCol:])
    refCount = genotype_info_concatenate.count('0')
    altCount = genotype_info_concatenate.count('1')
    returnValue = line_splited[:5]
    returnValue.extend([str(altCount), str(refCount)])
    return returnValue
def recParams(params):
    """ Receive the command line arguments """
    try:
        opts, args = getopt.getopt(params, 'i:', ['vcf='])
    except getopt.GetoptError as e:
        print(e)
        sys.exit(2)
    if not opts :
        usage()
        sys.exit(0)
    for opt, value in opts :
        if opt in ('-i', '--vcf'):
            if not os.access(value, os.R_OK):
                raise getopt.GetoptError('File not readable')
            else:
                return value
        else:
            raise getopt.GetoptError('Invalid options')
def usage():
    """ The usage of this script """
if __name__ == '__main__' :
    vcf = recParams(sys.argv[1:])
    with open(vcf, 'r') as f :
        for line in f:
            print('\t'.join(countAllele(line)))
