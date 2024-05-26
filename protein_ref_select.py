#! /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/02.apps/miniconda3/envs/scical/bin/python3
# -*- coding: utf-8 -*-
# #Author: guanxiangyu
# #Date:   2022-06-22
# #Last Modified by:   guanxiangyu
# #Last Modified time: 2022-06-27

import pandas as pd
import os, sys, gzip
from Bio import SeqIO
import re

def FilteUniProDB(uniprot_json_file, csv_file):
    UniProdb = pd.read_json(uniprot_json_file, typ='series')
    read_in = pd.read_csv(csv_file, index_col=2, usecols=range(1,6))
    read_in = read_in.loc[read_in['clade_name'].str.startswith('s'),:]
    UniProdb = UniProdb.loc[UniProdb.index.intersection(read_in.index)]
    UniProdb = pd.concat([UniProdb, read_in['clade_name']], axis=1, join='inner')
    UniProdb.columns = pd.Index(['Path','species'])
    UniProdb['species'] = UniProdb['species'].str.slice(start=3).str.replace("_", " ")
    UniProdb['Proteome_ID'] = UniProdb['Path'].apply(lambda x: os.path.basename(x))
    return(UniProdb)

def GetUniProDBHandle(Series):
    taxID = Series.name
    Path = Series.Path
    Species = Series.species
    ProteomeID = Series.Proteome_ID
    Pro_fasta = os.path.join(Path, "{}_{}.fasta.gz".format(ProteomeID, taxID))
    Pro_addi_fasta = os.path.join(Path, "{}_{}_additional.fasta.gz".format(ProteomeID, taxID))
    if os.path.isfile(Pro_addi_fasta):
        return(Pro_fasta, Pro_addi_fasta)
    else:
        return(Pro_fasta, 0)

def ReadFasta(filename):
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'rt')
    else:
        f = open(filename)
    return(f)

def ReadLine(record):
    '''
    [UniProtKB]
    >db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
    Where:
    db is 'sp' for UniProtKB/Swiss-Prot and 'tr' for UniProtKB/TrEMBL.
    UniqueIdentifier is the primary accession number of the UniProtKB entry.
    EntryName is the entry name of the UniProtKB entry.
    ProteinName is the recommended name of the UniProtKB entry as annotated in the RecName field. For UniProtKB/TrEMBL entries without a RecName field, the SubName field is used. In case of multiple SubNames, the first one is used. The 'precursor' attribute is excluded, 'Fragment' is included with the name if applicable.
    OrganismName is the scientific name of the organism of the UniProtKB entry.
    OrganismIdentifier is the unique identifier of the source organism, assigned by the NCBI.
    GeneName is the first gene name of the UniProtKB entry. If there is no gene name, OrderedLocusName or ORFname, the GN field is not listed.
    ProteinExistence is the numerical value describing the evidence for the existence of the protein.
    SequenceVersion is the version number of the sequence.
    Examples:
    >sp|Q8I6R7|ACN2_ACAGO Acanthoscurrin-2 (Fragment) OS=Acanthoscurria gomesiana OX=115339 GN=acantho2 PE=1 SV=1
    >sp|P27748|ACOX_CUPNH Acetoin catabolism protein X OS=Cupriavidus necator (strain ATCC 17699 / H16 / DSM 428 / Stanier 337) OX=381666 GN=acoX PE=4 SV=2
    >sp|P04224|HA22_MOUSE H-2 class II histocompatibility antigen, E-K alpha chain OS=Mus musculus OX=10090 PE=1 SV=1
    >tr|Q3SA23|Q3SA23_9HIV1 Protein Nef (Fragment) OS=Human immunodeficiency virus 1  OX=11676 GN=nef PE=3 SV=1
    >tr|Q8N2H2|Q8N2H2_HUMAN cDNA FLJ90785 fis, clone THYRO1001457, moderately similar to H.sapiens protein kinase C mu OS=Homo sapiens OX=9606 PE=2 SV=1
    Alternative isoforms (this only applies to UniProtKB/Swiss-Prot):
    >sp|IsoID|EntryName Isoform IsoformName of ProteinName OS=OrganismName OX=OrganismIdentifier[ GN=GeneName]
    Where:
    IsoID is the isoform identifier as assigned in the ALTERNATIVE PRODUCTS section of the UniProtKB entry.
    IsoformName is the isoform name as annotated in the ALTERNATIVE PRODUCTS Name field of the UniProtKB entry.
    ProteinExistence and SequenceVersion do not apply to alternative isoforms (ProteinExistence is dependent on the number of cDNA sequences, which is not known for individual isoforms).
    Example:
    >sp|Q4R572-2|1433B_MACFA Isoform Short of 14-3-3 protein beta/alpha OS=Macaca fascicularis OX=9541 GN=YWHAB

    SeqRecord(seq=Seq('MEEMYNVSSNPHVRDKMTTSRIMQLVVIALLPATLFGIWNFGFHALLVVLVTVI...KKA'), 
            id='tr|A0A173W9M9|A0A173W9M9_9FIRM', 
            name='tr|A0A173W9M9|A0A173W9M9_9FIRM',
            description='tr|A0A173W9M9|A0A173W9M9_9FIRM Ion-translocating oxidoreductasecomplex subunit D OS=Blautia obeum OX=40520 GN=nqrB PE=3 SV=1',
            dbxrefs=[])
    '''
    name = record.name
    # db, UniqueIdentifier, EntryName = name.split("|")
    seq = str(record.seq)
    desc = record.description
    desc = desc.replace("{} ".format(name), "", 1)
    # print(desc)
    desc = list(re.search(r'(.+)\sOS=([\w\s]+)\sOX=(\w+)\s(GN=(.+)\s)?PE=(\d)\sSV=(\d)', desc).groups())
    if len(desc) == 7:
        desc.pop(3)
        ProteinName, OrganismName, OrganismIdentifier, GeneName, ProteinExistence, SequenceVersion = desc
    return(name, seq, ProteinName, OrganismName, OrganismIdentifier, GeneName, ProteinExistence, SequenceVersion)
    # return(name, seq, desc)

class ProDB(object):
    'fasta merge'
    def __init__(self, seq, Id=None, ProteinName=None, OrganismName=None, OrganismIdentifier=None, GeneName=None, ProteinExistence=None, SequenceVersion=None):
        self.seq = seq
        self.Id = []
        self.ProteinName = []
        self.OrganismName = []
        self.OrganismIdentifier = []
        self.GeneName = []
        self.ProteinExistence = []
        self.SequenceVersion = []
        if Id is not None:
            self.add_Id(Id)
            self.add_ProName(ProteinName)
            self.add_OrgName(OrganismName)
            self.add_OrgaId(OrganismIdentifier)
            self.add_GenName(GeneName)
            self.add_ProExis(ProteinExistence)
            self.add_SeqVer(SequenceVersion)
    def add_Id(self, Id):
        self.Id.append(Id)
    def add_ProName(self, ProteinName):
        self.ProteinName.append(ProteinName)
    def add_OrgName(self, OrganismName):
        self.OrganismName.append(OrganismName)
    def add_OrgaId(self,OrganismIdentifier):
        self.OrganismIdentifier.append(OrganismIdentifier)
    def add_GenName(self, GeneName):
        if GeneName is not None:
            self.GeneName.append(GeneName)
        else:
            self.GeneName.append(" ")
    def add_ProExis(self, ProteinExistence):
        self.ProteinExistence.append(ProteinExistence)
    def add_SeqVer(self, SequenceVersion):
        self.SequenceVersion.append(SequenceVersion)

def DictAdd(record, database):
    name, seq, ProteinName, OrganismName, OrganismIdentifier, GeneName, ProteinExistence, SequenceVersion = ReadLine(record)
    my_key = seq
    if my_key not in database:
        database[my_key] = ProDB(seq, name, ProteinName, OrganismName, OrganismIdentifier, GeneName, ProteinExistence, SequenceVersion)
    else:
        # database[my_key].add_entry(name)
        # database[my_key].add_desc(desc)
        database[my_key].add_Id(name)
        database[my_key].add_ProName(ProteinName)
        database[my_key].add_OrgName(OrganismName)
        database[my_key].add_OrgaId(OrganismIdentifier)
        database[my_key].add_GenName(GeneName)
        database[my_key].add_ProExis(ProteinExistence)
        database[my_key].add_SeqVer(SequenceVersion)

def DictOutPut(database, OutPutPath, sample):
    filename = "{}.microbe.prodb.fasta".format(sample)
    f = open(os.path.join(OutPutPath, filename), 'w')
    for key in sorted(database.keys()):
        # f.write('>{}|{}| {}\n'.format(key, ";".join(database[key].EntryName), ";".join(database[key].desc)))
        f.write('>{Id}||{ProteinName}||OS={OrganismName}||OX={OrganismIdentifier}||GN={GeneName}||PE={ProteinExistence}||SV={SequenceVersion}\n'.format(
            Id=";".join(database[key].Id), ProteinName=";".join(database[key].ProteinName), OrganismName=";".join(database[key].OrganismName), OrganismIdentifier=";".join(database[key].OrganismIdentifier), 
            GeneName=";".join(database[key].GeneName), ProteinExistence=";".join(database[key].ProteinExistence), SequenceVersion=";".join(database[key].SequenceVersion)))
        f.write('{}\n'.format(database[key].seq))
    f.close()
        
def main():
    uniprot_json_file = "/hwfssz5/ST_SUPERCELLS/P20Z10200N0067/guanxiangyu/source_data/uniprot/previous_releases/release_2021_04/knowledgebase/reference_proteomes/taxid.json"
    csv_file =sys.argv[1]
    OutPutPath = sys.argv[2]
    if not os.path.exists(OutPutPath):
        os.makedirs(OutPutPath, 755)
    sample = os.path.basename(csv_file).split(".")[0]
    UniProdb = FilteUniProDB(uniprot_json_file, csv_file)
    UniProdb.to_csv(os.path.join(OutPutPath, "{}.selected.uniprotdb.csv".format(sample)), index_label="Tax_id")
    database = dict()
    for i in range(UniProdb.shape[0]):
        Pro_fasta, Pro_addi_fasta = GetUniProDBHandle(UniProdb.iloc[i,])
        fasta_file = ReadFasta(Pro_fasta)
        for record in SeqIO.parse(fasta_file,'fasta'):
            DictAdd(record, database)
        if Pro_addi_fasta:
            add_fasta_file = ReadFasta(Pro_addi_fasta)
            for record in SeqIO.parse(add_fasta_file, 'fasta'):
                DictAdd(record, database)
    DictOutPut(database, OutPutPath, sample)

if __name__ == '__main__':
    main()