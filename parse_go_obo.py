#!/usr/bin/env python3

"""
> parse_go_obo.py <

Python script is a helper function that can be called by other
scripts to parse and return GO annotation data .

Script requires a _file object_ (not filename).
"""

# --- example entry ---
# [Term]
# id: GO:0008156
# name: negative regulation of DNA replication
# namespace: biological_process
# def: "Any process that stops, prevents, or reduces the frequency, rate or
#      extent of DNA replication." [GOC:go_curators]
# subset: gosubset_prok
# synonym: "DNA replication inhibitor" RELATED []
# synonym: "down regulation of DNA replication" EXACT []
# synonym: "down-regulation of DNA replication" EXACT []
# synonym: "downregulation of DNA replication" EXACT []
# synonym: "inhibition of DNA replication" NARROW []
# is_a: GO:0006275 ! regulation of DNA replication
# is_a: GO:0051053 ! negative regulation of DNA metabolic process
# is_a: GO:2000113 ! negative regulation of cellular macromolecule biosynthetic
#                    process
# relationship: negatively_regulates GO:0006260 ! DNA replication

# imports from standard library
import re

# hardcode filename - it's not often that this gets changed anyway *shrug*
go_term_file = open('/lithium/data_repo/goa/go-basic.obo')

go_entries = go_term_file.read().split('[Term]')
go_entries.pop(0)

# construct a dictionary to facilitate conversion of alt_id --> go_id
alt_to_go_id = {}

# construct a dictionary to store GO data
#  go_data[go_id] = {'name': go_name, 'namespace': go_namespace, ...}
go_data = {}
for e in go_entries:
    go_id = re.search('id: (GO:\d+)', e).group(1)

    if go_id not in go_data:
        go_data[go_id] = {}

    go_data[go_id]['name'] = re.search('name: (.*?)\n', e).group(1)
    go_data[go_id]['namespace'] = re.search('namespace: (.*?)\n', e).group(1)
    go_data[go_id]['alt_id'] = re.findall('alt_id: (GO:\d+)', e)
    go_data[go_id]['is_a'] = re.findall('is_a: (GO:\d+)', e)

    if go_data[go_id]['alt_id']:
        for a in go_data[go_id]['alt_id']:
            alt_to_go_id[a] = go_id

def translate_alt_id(alt_id):
    """
    This function attempts to translate the alt_id into its proper GO term;
    else it'd just return the original value back.
    """
    return alt_to_go_id[alt_id] if alt_id in alt_to_go_id else alt_id

def get_property(query_go_id, prop='name'):
    """
    This function returns the property of interest of the GO term, and raises
    an error if an invalid GO term is used, or if an invalid property is
    queried for. TRANSLATE ALT_IDs BEFORE USE!
    """
    if query_go_id in alt_to_go_id:
        raise KeyError(query_go_id + ' is an alt_id of ' + \
                       alt_to_go_id[query_go_id]+ '. Translate it!')

    if query_go_id in go_data:
        if prop in go_data[query_go_id]:
            return go_data[query_go_id][prop]
        else:
            raise KeyError(prop + ' is not a valid property!')
    else:
        raise KeyError(query_go_id + ' does not exist in the GO database!')

def get_all_parent_terms(query_go_id):
    """
    This function returns all parent terms of the query GO term (including
    itself!), and checks all parent terms (whether they are deprecated) before
    returning it in the form of a list.
    """
    if query_go_id in alt_to_go_id:
        raise KeyError(query_go_id + ' is an alt_id of ' + \
                       alt_to_go_id[query_go_id]+ '. Translate it!')

    parent_go_terms = [query_go_id]
    unresolved_terms = [query_go_id]

    while unresolved_terms:
        u = unresolved_terms.pop(0)
        if go_data[u]['is_a']:
            parent_go_terms += go_data[u]['is_a']
            unresolved_terms += go_data[u]['is_a']

            # debug:
            # print (u, go_data[u]['namespace'], go_data[u]['is_a'])

    parent_go_terms = sorted(list(set(parent_go_terms)))

    # QC, will slow down performance, but can be disabled:
    for p in parent_go_terms:
        if p in alt_to_go_id:
            raise KeyError(p + ' is an alt_id of ' + \
                           alt_to_go_id[p]+ ". This shouldn't happen!")

    return parent_go_terms

def bin_go_terms_by_namespace(list_go_terms):
    """
    This function bins the GO terms based on their namespaces. i.e.
    [GO_term_1, GO_term_2, ...] --> {'molecular_function': [...],
                                     'biological_process': [...],
                                     'cellular_component': [...]}
    """

    binned_terms = {'biological_process': [], 'cellular_component': [],
                    'molecular_function': []}
    for g in list_go_terms:
        binned_terms[get_property(g, 'namespace')].append(g)

    for b in binned_terms:
        binned_terms[b] = sorted(list(set(binned_terms[b])))

    return binned_terms


if __name__ == '__main__':
    print (get_property('GO:0036105', 'namespace')) # 'molecular_function'
    print (get_all_parent_terms('GO:0000003'))      # ['GO:0000003', 'GO:0008150']
    print (get_all_parent_terms('GO:0000113'))      # ['GO:0000109', 'GO:0000113', 'GO:0005575', 'GO:0032991', 'GO:0043234', 'GO:0044422', 'GO:0044424', 'GO:0044428', 'GO:0044446', 'GO:0044464', 'GO:1990391']
    print (get_all_parent_terms('GO:0000001'))
    print (bin_go_terms_by_namespace(get_all_parent_terms('GO:0000003') +
                                     get_all_parent_terms('GO:0000113') +
                                     get_all_parent_terms('GO:0000001')))