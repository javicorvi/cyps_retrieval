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

def Main(parameters):
    output=parameters['output']
    uniprot_search_query= parameters['uniprot_search_query']
    if not os.path.exists(output):
        os.makedirs(output)
    outputFileUniprot = output+"/cyps_uniprot_search.txt"
    #download_cyps(uniprot_search_query, outputFileUniprot)
    outputFileFilter = output+"/cyps_uniprot_search_filter.txt"
    filter_cyps(outputFileUniprot, outputFileFilter)
    
    
def ReadParameters(args):
    if(args.p!=None):
        Config = ConfigParser.ConfigParser()
        Config.read(args.p)
        parameters['output']=Config.get('MAIN', 'output')
        parameters['uniprot_search_query']=Config.get('MAIN', 'uniprot_search_query')
    else:
        logging.error("Please send the correct parameters config.properties --help ")
        sys.exit(1)
    return parameters   

def filter_cyps(unputUniprotFile, outputFile):
    with open(unputUniprotFile,'r') as result: 
        with open(outputFile,'w') as new_result: 
            for line in result:
                data = line.split('\t')
                if(data[4].upper().startswith("CYP")):
                    generate_terms_for_gene_name(data[4], data, new_result)
                    #new_result.write(line)
                    #new_result.flush()

def generate_terms_for_gene_name(geneNames, data, new_result):
    geneNames = geneNames.replace(" ",",")
    geneNames = geneNames.replace(";",",")
    geneNames = geneNames.replace("/",",")
    geneNames_data = geneNames.split(",")
    #[os.path.join(input_file, f) for f in os.listdir(input_file) if (os.path.isfile(os.path.join(input_file, f)) & f.endswith('.xml.txt') & (os.path.basename(f) not in ids_list))]
    geneNames_data = [(x.replace("-","")).upper() for x in geneNames_data if x!=""]
    geneNames_data = set(geneNames_data)
    for gn in geneNames_data:
        #new_result.write(data[0]+'\t'+data[1]+'\t'+data[2]+'\t'+data[3]+'\t'+gn+'\t'+data[5]+'\t'+data[6])
        #new_result.write(data[0]+'\t'+data[1]+'\t'+gn+'\n')
        if(gn.startswith("CYP")):
            generate_variations(gn, data, new_result)
            gn=gn.lower()
            generate_variations(gn, data, new_result)
            new_result.flush()
            gn = gn[:0]+'C'+gn[1:]
            generate_variations(gn, data, new_result)
        new_result.flush()   
        
        
def generate_variations(gn, data, new_result):
    try:
        pattern = re.compile(r'([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)')
        divi = [""," ","-","/","*"]
        variants = []
        for (cyp, number, fam, sub ) in re.findall(pattern, gn):
            for d in divi:
                for d2 in divi:
                    for d3 in divi:
                        variants.append(''.join((cyp, d, number, d2,fam, d3, sub ))) 
                        
        #print variants   
        for vari in variants:
            new_result.write(data[0]+'\t'+data[1]+'\t'+vari+'\n')     
        new_result.flush()
    except Exception as inst:
        print gn     
    '''
    s = re.search(r"\d+(\D+\d+\D+)?", gn)
    d = s.group(0)
    divi = [" ","-","/","*"]
    for d in divi:
        with_1 = gn[:3] + d + gn[3:]
        with_2 = gn[:4] + d + gn[4:]
        with_3 = gn[:5] + d + gn[5:]
        new_result.write(data[0]+'\t'+data[1]+'\t'+with_1+'\n')        
        new_result.write(data[0]+'\t'+data[1]+'\t'+with_2+'\n')   
        new_result.write(data[0]+'\t'+data[1]+'\t'+with_3+'\n') 
        for d2 in divi:
            with_11 = gn[:3] + d + gn[3:4] + d2 + gn[4:] 
            with_12 = gn[:3] + d2 + gn[3:4] + d + gn[4:] 
            new_result.write(data[0]+'\t'+data[1]+'\t'+with_11+'\n')  
            new_result.write(data[0]+'\t'+data[1]+'\t'+with_12+'\n')
            for d3 in divi:
                with_11 = gn[:3] + d + gn[3:4] + d2 + gn[4:5] + d2 + gn[4:]
                with_11 = gn[:3] + d + gn[3:4] + d2 + gn[4:]
                new_result.write(data[0]+'\t'+data[1]+'\t'+with_11+'\n')        
                new_result.write(data[0]+'\t'+data[1]+'\t'+with_12+'\n')  
    new_result.flush()   
    print ""'''       
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



