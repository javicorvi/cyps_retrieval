import sys
import argparse
import ConfigParser
import urllib, urllib2
import os
import logging
import re
import xml.etree.ElementTree as ET

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

roots = ['CYP','Cyp','cyp','Cytochrome P450','CYTOCHROME P450','cytochrome P450','Cytochrome P-450','CYTOCHROME P-450','cytochrome P-450','Cytochrome-P450','CYTOCHROME-P450','cytochrome-P450','P450','p450','P-450','p-450','P450 (CYP)','P450 (cyp)','P450 (Cyp)','Cytochrome','CYTOCHROME','cytochrome']

def Main(parameters):
    output=parameters['output']
    outputDict=parameters['outputDict']
    uniprot_search_query= parameters['uniprot_search_query']
    if not os.path.exists(output):
        os.makedirs(output)
    #outputFileUniprot = output+"/cyps_uniprot_search.txt"
    outputFileUniprot = output+"/cyps_uniprot_search.xml"
    download_cyps(uniprot_search_query, outputFileUniprot)
    outputFileFilter = outputDict
    #filter_cyps_tab(outputFileUniprot, outputFileFilter)
    filter_cyps_xml(outputFileUniprot, outputFileFilter)
    
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

def filter_cyps_tab(unputUniprotFile, outputFile):
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




def filter_cyps_xml(unputUniprotFile, outputFile):
    logging.info("Filter and Generate Variations  " )
    name_space = "{http://uniprot.org/uniprot}"
    docXml = ET.parse(unputUniprotFile)
    root = docXml.getroot()
    with open(outputFile,'w') as new_result: 
        new_result.write('uniprot_entry_name\torganism\ttype\tkeyword\tkeyword_variant\n')
        for root_cyps in roots:
            new_result.write('\t\t'+ 'generic_keyword'+'\t'+root_cyps+'\n')
        new_result.flush()
        for entry in root.findall(name_space+"entry"):
            try:
                
                name=entry.find(name_space+"name").text
                if(name=="CP21A_HUMAN"):
                    print ""
                #First read the gene and verify if starts with CYP
                gene=entry.find(name_space+"gene")
                if(gene is not None):
                    
                   
                    
                    geneNames = gene.findall(name_space+"name")
                    genes_alter = []
                    for geneName in geneNames:
                        if(geneName.attrib.get("type")=="primary"):
                            gene_primary = geneName.text
                        else:
                            genes_alter.append(geneName.text)
                    
                    if(gene_primary is not None and gene_primary.upper().startswith("CYP")):
                        #Name
                        name=entry.find(name_space+"name").text
                        #read organism information
                        organism=entry.find(name_space+"organism")
                        organismNames = organism.findall(name_space+"name")
                        for organismName in organismNames:
                            organismName.text
                        dbReference = organism.find(name_space+"dbReference")
                        db_organism = dbReference.attrib.get("type")
                        id_organism = dbReference.attrib.get("id")
                        
                        new_result.write(name+'\t'+ db_organism+':'+id_organism +'\t'+ 'entryName'+'\t'+name+'\n')
                        
                        accessions=entry.findall(name_space+"accession")
                        accessions_keys = []
                        #Accession Data
                        for access in accessions:
                            accessions_keys.append(access.text)
                            new_result.write(name+'\t'+ db_organism+':'+id_organism +'\t'+ 'accession'+'\t'+access.text+'\n')
                        
                        #Protein
                        protein = entry.find(name_space+"protein")
                        #Recomended Name
                        recomendedName = protein.find(name_space+"recommendedName")
                        if(recomendedName is not None):
                            #Full Name
                            fullname=recomendedName.find(name_space+"fullName").text
                            new_result.write(name+'\t'+ db_organism+':'+id_organism +'\t'+ 'recomendedName'+'\t'+fullname+'\n')
                            #EC Number
                            ecNumber=recomendedName.find(name_space+"ecNumber")
                            if(ecNumber is not None):
                                ecNumber_str = ecNumber.text
                        alternativeNames = protein.findall(name_space+"alternativeName")
                        fullNameAlterNames = []
                        for alternativeName in alternativeNames:
                            fullNameAlterNames.append(alternativeName.find(name_space+"fullName").text)
                            new_result.write(name+'\t'+ db_organism+':'+id_organism +'\t'+ 'alternativeName'+'\t'+alternativeName.find(name_space+"fullName").text+'\n')
                        new_result.flush()
                        new_result.write(name+'\t'+ db_organism+':'+id_organism +'\t'+ 'gene_primary'+'\t'+gene_primary+'\n')
                        #genesprocessed is for not repeat variations in the same family
                        genesToProcess = []
                        generate_terms_for_gene_name_singular(gene_primary,  name,'gene_primary_variant',db_organism+':'+id_organism , new_result, genesToProcess)
                        for genes_alt in genes_alter:
                            new_result.write(name+'\t'+ db_organism+':'+id_organism +'\t'+ 'gene_synonym'+'\t'+genes_alt+'\n')
                            generate_terms_for_gene_name_singular(genes_alt,  name, 'gene_synonym_variant',db_organism+':'+id_organism , new_result, genesToProcess)
                else:
                    print "Entry without gene data : " +  entry.find(name_space+"name").text          
            except Exception as inst:
                print "error reading " + entry.find(name_space+"name").text
                print str(inst)
    logging.info("filter_cyps_xml Process end " )

    
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

def generate_terms_for_gene_name_singular(geneName, name, type, organism, new_result,genesProcessed):
    #subroot = ['P450','p450','P-450','p-450']
    #The gene can contain - and spaces inside the name.
    #delete repeated families
    geneName = geneName.upper()
    geneName = geneName.replace("-","")
    if not any(geneName == s for s in genesProcessed): 
        for root in roots:
            generate_variations(geneName, root, name, type, organism, new_result)
            new_result.flush()
        genesProcessed.append(geneName)    
    else:
        print "The family is already processed: " +  geneName   
'''
def generate_variations(gn, root, data, new_result):
    try:
        pattern = re.compile(r'([a-zA-Z]+)([0-9]+)([a-zA-Z]+)([0-9]+)')
        divi = [""," ","-"]
        #divi = [""," ","-","/","*"]
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
'''  
def generate_variations(geneName, root, name, type, organism, new_result):
    try:
        pattern = re.compile(r'([a-zA-Z]+)([0-9]+)([a-zA-Z]+)([0-9]+)')
        divi = [""," ","-"]
        #divi = [""," ","-","/","*"]
        variants = []
        for (cyp, number, fam, sub ) in re.findall(pattern, geneName):
            for d in divi:
                for d2 in divi:
                    for d3 in divi:
                        variants.append(''.join((root, d, number, d2,fam, d3, sub )))
                        
                        #variants.append(''.join((root, d, number, d2,fam.lower(), d3, sub )))#no se puede diferentes case pertenecen a otras species.
                        
                        #variants.append(''.join((root, d, int_to_roman(int(number)), d2,fam, d3, int_to_roman(int(sub)) )))
                        #variants.append(''.join((root, d, int_to_roman(int(number)), d2,fam, d3, sub ))) 
                        #variants.append(''.join((root, d, number, d2,fam, d3, int_to_roman(int(sub)) )))
        for vari in variants:
            #new_result.write(data[0]+'\t'+data[1]+'\t'+gn+'\t'+vari+'\n')
            new_result.write(name+'\t'+organism+'\t'+type+'\t'+geneName+'\t'+vari+'\n')     
        new_result.flush()
    except Exception as inst:
        print "Error generating variations for CYP : " + geneName 
    
def download_cyps(uniprot_search_query, outputFile):
    logging.info("Downloading CYPs Query : " + uniprot_search_query )
    url = 'https://www.uniprot.org/uniprot/'
    #params = urllib.urlencode({'query':uniprot_search_query,'format':'tab','force':'true','cols':'id,entry name,reviewed,protein names,genes(PREFERRED),genes(ALTERNATIVE),genes(OLN),genes(ORF),organism,organism-id','sort':'score','compress':'no'})
    params = urllib.urlencode({'query':uniprot_search_query,'format':'xml','offset':1,'cols':'id,entry name,reviewed,protein names,genes,organism,length','sort':'score','compress':'no'})
    
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
