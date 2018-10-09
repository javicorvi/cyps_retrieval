import sys
import argparse
import ConfigParser
import urllib, urllib2
import os
import logging
import re
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

parser=argparse.ArgumentParser()
parser.add_argument('-p', help='Path Parameters')
args=parser.parse_args()
parameters={}
if __name__ == '__main__':
    import cyps_retrieval
    parameters = cyps_retrieval.ReadParameters(args)     
    cyps_retrieval.Main(parameters)


ROMAN = [
    (1000, "M"),
    ( 900, "CM"),
    ( 500, "D"),
    ( 400, "CD"),
    ( 100, "C"),
    (  90, "XC"),
    (  50, "L"),
    (  40, "XL"),
    (  10, "X"),
    (   9, "IX"),
    (   5, "V"),
    (   4, "IV"),
    (   1, "I"),
]


def Main(parameters):
    output=parameters['output']
    outputDict=parameters['outputDict']
    uniprot_search_query= parameters['uniprot_search_query']
    if not os.path.exists(output):
        os.makedirs(output)
    outputFileUniprot = output+"/cyps_uniprot_search.txt"
    #download_cyps(uniprot_search_query, outputFileUniprot)
    outputFileFilter = outputDict
    filter_cyps(outputFileUniprot, outputFileFilter)
    
    
def ReadParameters(args):
    if(args.p!=None):
        Config = ConfigParser.ConfigParser()
        Config.read(args.p)
        parameters['output']=Config.get('MAIN', 'output')
        parameters['outputDict']=Config.get('MAIN', 'outputDict')
        parameters['uniprot_search_query']=Config.get('MAIN', 'uniprot_search_query')
    else:
        logging.error("Please send the correct parameters config.properties --help ")
        sys.exit(1)
    return parameters   

def filter_cyps(unputUniprotFile, outputFile):
    logging.info("Filter and Generate Variations  " )
    with open(unputUniprotFile,'r') as result: 
        with open(outputFile,'w') as new_result: 
            for line in result:
                data = line.split('\t')
                if(data[4].upper().startswith("CYP")):
                    generate_terms_for_gene_name(data[4], data, new_result)
                    generate_terms_for_entry_name(data[2], data, new_result)
                    #new_result.write(line)
                    #new_result.flush()
    logging.info(" Process end" )
#TODO pregu
def generate_terms_for_entry_name(entryName, data, new_result):
    new_result.write(data[0]+'\t'+data[1]+'\t null \t'+entryName+'\n')     
    new_result.flush()

def generate_terms_for_gene_name(geneNames, data, new_result):
    geneNames = geneNames.replace(" ",",")
    geneNames = geneNames.replace(";",",")
    geneNames = geneNames.replace("/",",")
    geneNames_data = geneNames.split(",")
    #[os.path.join(input_file, f) for f in os.listdir(input_file) if (os.path.isfile(os.path.join(input_file, f)) & f.endswith('.xml.txt') & (os.path.basename(f) not in ids_list))]
    #geneNames_data = [(x.replace("-","")).upper() for x in geneNames_data if x!=""]
    geneNames_data = [(x.replace("-","")) for x in geneNames_data if x!=""]
    geneNames_data = set(geneNames_data)
    roots = ['CYP','Cyp','cyp','Cytochrome P450','CYTOCHROME P450','cytochrome P450','Cytochrome P-450','CYTOCHROME P-450','cytochrome P-450','Cytochrome-P450','CYTOCHROME-P450','cytochrome-P450','P450','p450','P-450','p-450','P450 (CYP)','P450 (cyp)','P450 (Cyp)','Cytochrome','CYTOCHROME','cytochrome']
    #subroot = ['P450','p450','P-450','p-450']
    for root in roots:
        for gn in geneNames_data:
            new_result.write('\tnull\tnull\t'+root+'\n')
            if(gn.startswith("CYP") or gn.startswith("Cyp")):
                generate_variations(gn, root, data, new_result)
                new_result.flush()   
        
def generate_variations(gn, root, data, new_result):
    try:
        pattern = re.compile(r'([a-zA-Z]+)([0-9]+)([a-zA-Z]+)([0-9]+)')
        divi = [""," ","-","/","*"]
        variants = []
        for (cyp, number, fam, sub ) in re.findall(pattern, gn):
            for d in divi:
                for d2 in divi:
                    for d3 in divi:
                        variants.append(''.join((root, d, number, d2,fam, d3, sub )))
                        #variants.append(''.join((root, d, number, d2,fam.lower(), d3, sub )))#no se puede diferentes case pertenecen a otras species.
                        variants.append(''.join((root, d, int_to_roman(int(number)), d2,fam, d3, int_to_roman(int(sub)) )))
                        variants.append(''.join((root, d, int_to_roman(int(number)), d2,fam, d3, sub ))) 
                        variants.append(''.join((root, d, number, d2,fam, d3, int_to_roman(int(sub)) )))
        for vari in variants:
            #new_result.write(data[0]+'\t'+data[1]+'\t'+gn+'\t'+vari+'\n')
            new_result.write(data[0]+'\t'+data[1]+'\t'+gn+'\t'+vari+'\n')     
        new_result.flush()
    except Exception as inst:
        print "Error generating variations for CYP : " + gn 
    
def download_cyps(uniprot_search_query, outputFile):
    logging.info("Downloading CYPs Query : " + uniprot_search_query )
    url = 'https://www.uniprot.org/uniprot/'
    params = urllib.urlencode({'query':uniprot_search_query,'format':'tab','force':'true','cols':'id,entry name,reviewed,protein names,genes,organism,length','sort':'score','compress':'no'})
    request = urllib2.Request(url, params)
    response = urllib2.urlopen(request)
    response_cyps = response.read()
    with open(outputFile,'w') as result: 
        result.write(response_cyps)
        result.flush()
    logging.info("Download End ")  


def int_to_roman(number):
    result = []
    for (arabic, roman) in ROMAN:
        (factor, number) = divmod(number, arabic)
        result.append(roman * factor)
        if number == 0:
            break
    return "".join(result)
