import sys
from pathlib import Path
import requests
from hashlib import md5
import tarfile
import warnings

def download_url(url : str, target_dir : Path) -> Path:
    '''
        downloads a file from the web.
        Args: 
            url : path to the file to be downloaded
            target_dir : local directory to save the file into
        Returns:
            Path object pointing to the newly downloaded file
    '''
    # based on: https://stackoverflow.com/questions/16694907/download-large-file-in-python-with-requests
    local_file_path = Path(target_dir) / url.split('/')[-1]

    with requests.get(url, stream=True) as req:
        req.raise_for_status() #raises HTTPError if there is an error.
        with open(local_file_path, 'wb') as local_file:
            for chunk in req.iter_content(chunk_size=12288): #chunk size has not been optimized
                local_file.write(chunk)
    return local_file_path

def unpack(filename, path='.'):
    """
        uncompress a tar.gz archive
        from: https://github.com/HadrienG/taxadb
    """

    with tarfile.open(filename, "r:gz") as archive:
        archive.extractall(path,)
        archive.close()

class NCBITaxonomy:
    def __init__(self, local_path, overwrite=True):
        local_path = Path(local_path)
        

        # check if local copy of taxdump.tar.gz exists and has the same md5 as the most recent one from ncbi
        # partly derived from https://github.com/etetoolkit/ete/blob/1f587a315f3c61140e3bdbe697e3e86eda6d2eca/ete3/ncbi_taxonomy/ncbiquery.py#L667

        redownload = True
        local_taxdump_tar_gz_path = local_path / "taxdump.tar.gz"
        if (local_taxdump_tar_gz_path).exists():
            if overwrite == False:
                redownload = False
            else:
                taxdmp_md5 = md5(open(local_taxdump_tar_gz_path, "rb").read()).hexdigest()
                try:
                    md5_path = download_url(
                        "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5", local_path)
                    with open(md5_path, "r") as md5_file:
                        md5_check = md5_file.readline().split()[0]
                    if taxdmp_md5 == md5_check:
                        print("Local taxdump md5 matches server", file=sys.stderr)
                        redownload = False
                    else:
                        print(f"Local taxdump md5 does not match server {taxdmp_md5}, {md5_check}", file=sys.stderr)
                except Exception as e:
                    warnings.warn(f"Could not check md5 of taxdump.tar.gz: {e}")
                    redownload = False
        
        # download and extract taxdump.tar.gz from ncbi if the one from ncbi is different. 
        if redownload:
            print("Downloading taxdump.tar.gz", file=sys.stderr)
            local_taxdump_tar_gz_path = download_url("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", local_path)
            unpack(local_taxdump_tar_gz_path, local_path)

        # load the taxdmp file into memory
        # print("Loading nodes.dmp", file=sys.stderr)
        print("Loading NCBI taxonomy database", file=sys.stderr)
        self._nodes_tab = {int(x['tax_id']) : x for x in self.parse_taxdmp(local_path / 'nodes.dmp', 'nodes')}
        self._merged_tab = {int(x["old_tax_id"]) : x for x in self.parse_taxdmp(local_path / 'merged.dmp', 'merged')}
        self._deleted_tab = {int(x["tax_id"]) : x for x in self.parse_taxdmp(local_path / 'delnodes.dmp', 'delnodes')}
        self._names = self._names_dict(local_path / 'names.dmp')
        print("Loaded NCBI taxonomy database", file=sys.stderr)

    @classmethod
    def _names_dict(cls, path):
        '''
            returns a dict of all names for each taxid
        '''
        out_dict = dict()
        for line in cls.parse_taxdmp(path, 'names'):
            taxid = int(line['tax_id'])
            if taxid not in out_dict:
                out_dict[taxid] = {'names': [], 'unique_names': [], 'name_classes': [], 'scientific_name': None} 
            out_dict[taxid]["names"].append(line['name_txt'])
            out_dict[taxid]["unique_names"].append(line['unique name'])
            out_dict[taxid]["name_classes"].append(line['name class'])
            if line['name class'] == 'scientific name':
                out_dict[taxid]["scientific_name"] = line['name_txt']
        
        return out_dict
    

    @staticmethod
    def parse_taxdmp(path, filetype):
        headers = {
        "nodes": ["tax_id","parent_tax_id","rank","embl_code",'division_id', 'inherited_div_flag','genetic_code_id', 'inherited_GC_flag','mitochondrial_genetic_code_id','inherited_MGC_flag','GenBank_hidden_flag','hidden_subtree_root_flag','comments'],
        "names": ["tax_id","name_txt","unique name","name class"],
        "gencode": ['genetic code id', 'abbreviation', 'name', 'cde', 'starts'],
        "delnodes": ['tax_id'],
        "merged": ["old_tax_id","new_tax_id"],
        "citations": ["cit_id","cit_key","pubmed_id","medline_id",'url','text', 'tax_id_list'],
        "division" : ["division_id","division_cde","pubmed_id","division_name",'comments']
        }

        if filetype not in headers:
            raise ValueError(f"invalid filedtype {filetype}")

        with open(path, 'r') as infile:
            for line in infile:
                line = line.rstrip("\n|\t")
                parts = line.split("\t|\t")
                yield dict(zip(headers[filetype],parts))

    
    @staticmethod
    def _lineage(nodes_tab, merged_tab, deleted_tab, tax_id):
        '''
            returns all parent nodes of tax_id (starting with tax_id itself)
        '''
        out_list = list()

        if tax_id in merged_tab:
            tax_id = int(merged_tab[tax_id]["new_tax_id"])

        if tax_id in deleted_tab: #TODO: Not sure what the best way to handle this is, we may want to do some kind of custom handling for deleted taxids to ensure desired behavior
            warnings.warn(f"Attempted to calculate lineage for a deleted taxid: {tax_id}, lineage search prematurely stopped, which might confuse lineage exclusions.")
            return out_list
        
        out_list.append(tax_id)

        while nodes_tab[tax_id]['parent_tax_id'] != "1":
            # print(tax_id)
            parent_taxid = int(nodes_tab[tax_id]['parent_tax_id'])
            if parent_taxid in merged_tab:
                parent_taxid = int(merged_tab[parent_taxid]["new_tax_id"])
            
            out_list.append(parent_taxid) 
            tax_id = parent_taxid

        out_list.append(1)
        return out_list

    def normalized_taxid(self, tax_id):
        '''
            returns the normalized taxid of tax_id
        '''
        if tax_id in self._merged_tab:
            return int(self._merged_tab[tax_id]["new_tax_id"])
        elif tax_id in self._deleted_tab: #TODO: Not sure what the best way to handle this is, we may want to do some kind of custom handling for deleted taxids to ensure desired behavior
            return 32644 # this is the normalized taxid for "unidentified sequences"
        else:
            return tax_id

    def lineage(self, tax_id):
        '''
            returns all parent nodes of tax_id
        '''
        return self._lineage(self._nodes_tab, self._merged_tab, self._deleted_tab, tax_id)

    def name(self, tax_id):
        '''
            returns the scientific name of tax_id
        '''
        if tax_id not in self._names:
            return None
        elif self._names[tax_id]['scientific_name'] is not None:
            return self._names[tax_id]['scientific_name']
        else:
            return self._names[tax_id]['names'][0]

    def rank(self, tax_id):
        '''
            returns the rank of tax_id
        '''
        return self._nodes_tab[tax_id]['rank']

    # def named_lineage(self, tax_id):
    #     '''
    #         returns all parent nodes of tax_id with their names
    #     '''
    #     return {self.rank(tax_id): {"taxid": int(tax_id), "name": self.name(tax_id)} for tax_id in self.lineage(tax_id)}


