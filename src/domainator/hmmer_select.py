"""Find profiles in an hmm file using text-based searches

Search through an hmm file to find profiles with names or descriptions matching a specified regex or string, or containing a specified string.

All specifications are treated as "OR" with relation to each other. So any profile that matches any of the specifications will be returned.

"""
#subset an hmm file based on name, description or other properties. use pyhmmer.
from jsonargparse import ArgumentParser, ActionConfigFile
import pyhmmer
import re
import sys
from domainator import __version__, RawAndDefaultsFormatter
from domainator.utils import pyhmmer_decode
from typing import Optional, List

def hmmer_select(hmm_path=None, query_regex:Optional[List[str]]=None, query_exact:Optional[List[str]]=None, query_contains:Optional[List[str]]=None, hmm_iterator=None, search_name=True, search_description=True, search_accession=True, case_sensitive=False):
    """
        hmm_path: path to an hmmer file
        query_regex: regular expression to try to find in the hmmer profiles
        query_exact: exact string to try to find in the hmmer profiles
        query_contains: string to try to find in the hmmer profiles
        hmm_iterator: instead of an hmm file path, you can supply an iterator over pyhmmer HMM objects
        search_name: If True then search in the NAME field of the hmm profiles
        search_description: If True then search in the DESC field of the hmm profiles
        search_accession: If True then search in the ACC field of the hmm profiles
        case_sensitive: If True then the search is case sensitive
    """

    # match_types = 0
    # if query_regex is not None:
    #     match_types += 1
    # if query_exact is not None:
    #     match_types += 1
    # if query_contains is not None:
    #     match_types += 1
    # if match_types != 1:
    #     raise RuntimeError("Exactly one of query_regex or query_exact or query_contains must be specified.")

    if not case_sensitive:
        if query_exact is not None:
            for i in range(len(query_exact)):
                query_exact[i] = query_exact[i].lower()
        if query_contains is not None:
            for i in range(len(query_contains)):
                query_contains[i] = query_contains[i].lower()

    if query_regex is not None:
        regexes = []
        for i in range(len(query_regex)):
            if case_sensitive:
                regexes.append(re.compile(query_regex[i]))
            else:
                regexes.append(re.compile(query_regex[i], re.IGNORECASE))

    inputs_specified = 0 
    if hmm_path is not None:
        inputs_specified += 1
    if hmm_iterator is not None:
        inputs_specified += 1
    if inputs_specified != 1:
        raise RuntimeError("Exactly one of input_path or hmm_iterator must be specified.")

    if ((not search_name) and (not search_description) and (not search_accession)):
        raise RuntimeError("at least one of search_name or search_description must be True")
    
    def inner(hmm_iterator):
        for model in hmm_iterator:
            model_name = pyhmmer_decode(model.name)
            model_description = model.description
            model_accession = model.accession
  
            if model_description == None:
                model_description = ""
            else:
                model_description = pyhmmer_decode(model_description)
            if model_accession == None:
                model_accession = ""
            else:
                model_accession = pyhmmer_decode(model_accession)

            if not case_sensitive:
                model_name = model_name.lower()
                model_description = model_description.lower()
                model_accession = model_accession.lower()


            found = False
            if query_regex is not None:
                for regex in regexes:
                    if (search_name and regex.search(model_name)):
                        found = True
                    elif (search_description and regex.search(model_description)):
                        found = True
                    elif (search_accession and regex.search(model_accession)):
                        found = True
            if query_exact is not None:
                for exact in query_exact:
                    if (search_name and model_name == exact):
                        found = True
                    elif (search_description and model_description == exact):
                        found = True
                    elif (search_accession and model_accession == exact):
                        found = True
            if query_contains is not None:
                for contains in query_contains:
                    if (search_name and contains in model_name):
                        found = True
                    elif (search_description and contains in model_description):
                        found = True
                    elif (search_accession and contains in model_accession):
                        found = True

            if found:
                yield  model
    
    if hmm_iterator is None:
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_iterator:
            yield from inner(hmm_iterator)
    else:
        yield from inner(hmm_iterator)
    


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument("-i", "--input", default=None, required=False,
                          nargs="+", type=str,
                          help="names of input hmm files. If not supplied, reads from stdin.")

    parser.add_argument("-o", "--output", default=None, required=False,  type=str,
                        help="Hmm output file name. If not supplied writes to stdout.")

    parser.add_argument("--regex", default=None, required=False, type=str, nargs="+",
                        help="What regex(es) to search.")
    parser.add_argument("--exact", default=None, required=False, type=str, nargs="+",
                            help="What exact string(s) to search. Target must match exactly with no additions at either end.")
    parser.add_argument("--contains", default=None, required=False, type=str, nargs="+",
                            help="What string(s) to search. Target must contain this as a substring.")
    
    parser.add_argument("--field", default="all", required=False, type=str.lower, choices={"name", "desc", "acc", "all"},
                        help="whether to search the NAME field, the DESC field, the ACC field, or all three.")
    
    parser.add_argument("--case_sensitive", default=False, required=False, action="store_true",
                        help="If set, then the search is case sensitive. Otherwise it is case insensitive.")

    parser.add_argument("--config", action=ActionConfigFile)

    params = parser.parse_args(argv)


    search_name = True
    search_description = True
    search_accession = True

    
    if params.field == "name":
        search_name = True
        search_description = False
        search_accession = False

    if params.field == "desc":
        search_name = False
        search_description = True
        search_accession = False
    
    if params.field == "acc":
        search_name = False
        search_description = False
        search_accession = True
    
    if params.input is None:
        input_files = [sys.stdin.buffer]
    else:
        input_files = params.input

    if params.output is None:
        output_handle = sys.stdout.buffer
    else:
        output_handle = open(params.output, "wb")


    found_count = 0
    for input_file in input_files:
        with pyhmmer.plan7.HMMFile(input_file) as hmm_iterator:
            for found_profile in hmmer_select(query_regex=params.regex, query_exact=params.exact, query_contains=params.contains, hmm_iterator=hmm_iterator, search_name=search_name, search_description=search_description, search_accession=search_accession, case_sensitive=params.case_sensitive):
                found_profile.write(output_handle)
                found_count += 1

    print(f"Found {found_count} matching profiles.", file=sys.stderr)        
    if params.output is not None:
        output_handle.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == "__main__":
    main(sys.argv[1:])