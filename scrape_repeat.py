import physcraper
import dendropy
import pickle
import sys
import os
import json
from peyotl.nexson_syntax import (
    extract_tree,
    get_subtree_otus,
    extract_otu_nexson,
    PhyloSchema,
)
import dendropy


study_id = "pg_873"
tree_id = "tree1679"
configfi = "full_run.config"

conf = physcraper.ConfigObj(configfi)

#dataset = physcraper.get_dataset_from_treebase(study_id,
#                                phylesystem_loc='api')

#aln = dataset.char_matrices[0]

#if len(aln) == 42:
#  pass
#else:
#  aln = dataset.char_matrices[1]


#aln.write(path="before.aln", schema="nexus")

aln = dendropy.DnaCharacterMatrix.get(file=open("before.aln"), schema="nexus")

nexson = physcraper.get_nexson(study_id, 'api')
newick = extract_tree(nexson,
                          tree_id,
                          PhyloSchema('newick',
                                      output_nexml2json='1.2.1',
                                      content="tree",
                                      tip_label="ot:originalLabel"))
tre = dendropy.Tree.get(data=newick,
                   schema="newick",
                   preserve_underscores=True,
                   taxon_namespace=aln.taxon_namespace)

tre.write(path="before.tre", schema="nexus")

data_obj = physcraper.generate_ATT_from_phylesystem(aln=aln,
                                         workdir='treebase',
                                         config_obj=conf,
                                         study_id=study_id,
                                         tree_id=tree_id)

data_obj.write_files()
json.dump(data_obj.otu_dict, open('treebase/otu_dict.json', 'wb'))

sys.stdout.write("{} taxa in alignement and tree\n".format(len(data_obj.aln)))

ids = physcraper.IdDicts(conf, workdir='treebase')

scraper = physcraper.PhyscraperScrape(data_obj, ids)

scraper.run_blast_wrapper()
scraper.read_blast_wrapper()
scraper.remove_identical_seqs()
scraper.generate_streamed_alignment()

sys.stdout.write("round 1 complete\n")
sys.stdout.write("{} taxa in alignment and tree\n".format(len(data_obj.aln)))

sys.stdout.write("round2\n")

scraper.run_blast_wrapper()
scraper.read_blast_wrapper()
scraper.remove_identical_seqs()
scraper.generate_streamed_alignment()
