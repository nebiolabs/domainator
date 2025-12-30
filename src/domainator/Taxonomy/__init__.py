import sys
from pathlib import Path
import requests
from hashlib import md5
import tarfile
import warnings
import io
import os
import fcntl
from contextlib import contextmanager
from typing import Iterator, Union, TextIO

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
    target_dir = Path(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    local_file_path = target_dir / url.split('/')[-1]
    tmp_path = local_file_path.with_suffix(local_file_path.suffix + ".part")

    with requests.get(url, stream=True, timeout=60) as req:
        req.raise_for_status()  # raises HTTPError if there is an error.
        with open(tmp_path, 'wb') as local_file:
            for chunk in req.iter_content(chunk_size=12288):  # chunk size has not been optimized
                if chunk:
                    local_file.write(chunk)

    os.replace(tmp_path, local_file_path)
    return local_file_path


@contextmanager
def _locked_file(lock_path: Path, lock_flag: int):
    """Context manager that acquires an advisory lock on a lockfile.

    Uses fcntl.flock (Linux/Unix). LOCK_SH allows parallel readers; LOCK_EX is exclusive.
    """
    lock_path = Path(lock_path)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "a+") as lock_handle:
        try:
            fcntl.flock(lock_handle.fileno(), lock_flag | fcntl.LOCK_NB)
        except BlockingIOError:
            lock_type = "exclusive" if (lock_flag & fcntl.LOCK_EX) else "shared"
            print(f"Waiting for {lock_type} lock: {lock_path}", file=sys.stderr)
            fcntl.flock(lock_handle.fileno(), lock_flag)
        try:
            yield lock_handle
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def _md5_path(path: Path) -> str:
    hasher = md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _fetch_remote_md5(url: str) -> str:
    # NCBI .md5 file contains: "<md5sum>  <filename>"
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    return resp.text.strip().split()[0]

class NCBITaxonomy:
    def __init__(self, local_path, overwrite=True):
        local_path = Path(local_path)
        local_path.mkdir(parents=True, exist_ok=True)
        

        # check if local copy of taxdump.tar.gz exists and has the same md5 as the most recent one from ncbi
        # partly derived from https://github.com/etetoolkit/ete/blob/1f587a315f3c61140e3bdbe697e3e86eda6d2eca/ete3/ncbi_taxonomy/ncbiquery.py#L667

        local_taxdump_tar_gz_path = local_path / "taxdump.tar.gz"
        lock_path = local_path / "taxdump.tar.gz.lock"

        # First take a shared lock to allow concurrent readers. Only escalate to exclusive when
        # a download/overwrite is required.
        with _locked_file(lock_path, fcntl.LOCK_SH):
            needs_download = False
            if not local_taxdump_tar_gz_path.exists():
                needs_download = True
            elif overwrite is False:
                needs_download = False
            else:
                try:
                    taxdmp_md5 = _md5_path(local_taxdump_tar_gz_path)
                    md5_check = _fetch_remote_md5(
                        "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5"
                    )
                    if taxdmp_md5 == md5_check:
                        print("Local taxdump md5 matches server", file=sys.stderr)
                        needs_download = False
                    else:
                        print(
                            f"Local taxdump md5 does not match server {taxdmp_md5}, {md5_check}",
                            file=sys.stderr,
                        )
                        needs_download = True
                except Exception as e:
                    warnings.warn(f"Could not check md5 of taxdump.tar.gz: {e}")
                    needs_download = False

        if needs_download:
            with _locked_file(lock_path, fcntl.LOCK_EX):
                # Re-check under exclusive lock in case another process already updated it.
                still_needs_download = False
                if not local_taxdump_tar_gz_path.exists():
                    still_needs_download = True
                elif overwrite is False:
                    still_needs_download = False
                else:
                    try:
                        taxdmp_md5 = _md5_path(local_taxdump_tar_gz_path)
                        md5_check = _fetch_remote_md5(
                            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5"
                        )
                        still_needs_download = taxdmp_md5 != md5_check
                    except Exception:
                        still_needs_download = False

                if still_needs_download:
                    print("Downloading taxdump.tar.gz", file=sys.stderr)
                    download_url(
                        "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", local_path
                    )

        # Load the taxonomy tables into memory. Prefer reading directly from taxdump.tar.gz
        # (no unpacking). Hold a shared lock during reads so another process cannot overwrite
        # the tarball mid-parse.
        print("Loading NCBI taxonomy database", file=sys.stderr)
        with _locked_file(lock_path, fcntl.LOCK_SH):
            loaded = False
            if local_taxdump_tar_gz_path.exists():
                try:
                    with tarfile.open(local_taxdump_tar_gz_path, "r:gz") as archive:
                        nodes_raw = archive.extractfile("nodes.dmp")
                        merged_raw = archive.extractfile("merged.dmp")
                        delnodes_raw = archive.extractfile("delnodes.dmp")
                        names_raw = archive.extractfile("names.dmp")
                        if any(x is None for x in (nodes_raw, merged_raw, delnodes_raw, names_raw)):
                            raise ValueError("taxdump.tar.gz is missing required .dmp members")

                        nodes_fh = io.TextIOWrapper(nodes_raw, encoding="utf-8", errors="replace")
                        merged_fh = io.TextIOWrapper(merged_raw, encoding="utf-8", errors="replace")
                        delnodes_fh = io.TextIOWrapper(delnodes_raw, encoding="utf-8", errors="replace")
                        names_fh = io.TextIOWrapper(names_raw, encoding="utf-8", errors="replace")
                        try:
                            self._nodes_tab = {int(x['tax_id']) : x for x in self.parse_taxdmp(nodes_fh, 'nodes')}
                            self._merged_tab = {int(x["old_tax_id"]) : x for x in self.parse_taxdmp(merged_fh, 'merged')}
                            self._deleted_tab = {int(x["tax_id"]) : x for x in self.parse_taxdmp(delnodes_fh, 'delnodes')}
                            self._names = self._names_dict(names_fh)
                            loaded = True
                        finally:
                            nodes_fh.close()
                            merged_fh.close()
                            delnodes_fh.close()
                            names_fh.close()
                except (tarfile.ReadError, ValueError) as e:
                    warnings.warn(
                        f"Could not read required taxonomy files from {local_taxdump_tar_gz_path}: {e}. "
                        "Falling back to extracted .dmp files if available."
                    )

            if not loaded:
                # Backwards-compatible fallback for directories that already contain extracted .dmp files.
                nodes_path = local_path / 'nodes.dmp'
                merged_path = local_path / 'merged.dmp'
                delnodes_path = local_path / 'delnodes.dmp'
                names_path = local_path / 'names.dmp'
                if all(p.exists() for p in (nodes_path, merged_path, delnodes_path, names_path)):
                    self._nodes_tab = {int(x['tax_id']) : x for x in self.parse_taxdmp(nodes_path, 'nodes')}
                    self._merged_tab = {int(x["old_tax_id"]) : x for x in self.parse_taxdmp(merged_path, 'merged')}
                    self._deleted_tab = {int(x["tax_id"]) : x for x in self.parse_taxdmp(delnodes_path, 'delnodes')}
                    self._names = self._names_dict(names_path)
                    loaded = True

            if not loaded:
                # If we get here, neither the tarball nor extracted dmps were readable.
                raise FileNotFoundError(
                    f"NCBI taxonomy database not found or unreadable in {local_path}. "
                    "Expected taxdump.tar.gz (with nodes.dmp/names.dmp/merged.dmp/delnodes.dmp) "
                    "or extracted .dmp files."
                )
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

        def _iter_lines(handle: TextIO) -> Iterator[str]:
            for raw_line in handle:
                yield raw_line

        if isinstance(path, (str, Path)):
            with open(path, 'r') as infile:
                for line in _iter_lines(infile):
                    line = line.rstrip("\n|\t")
                    parts = line.split("\t|\t")
                    yield dict(zip(headers[filetype],parts))
        else:
            # file-like object (e.g., streamed from tarfile.extractfile)
            for line in _iter_lines(path):
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


