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


study_id = "pg_1711"
tree_id = "tree3450"
configfi = "local.config"
workdir = "run_3"

conf = physcraper.ConfigObj(configfi)
conf.blast_loc = "local"


nexson = physcraper.opentree_helpers.get_nexson(study_id, 'api')
newick = extract_tree(nexson,
                          tree_id,
                          PhyloSchema('newick',
                                      output_nexml2json='1.2.1',
                                      content="tree",
                                      tip_label="ot:originalLabel"))
tre = dendropy.Tree.get(data=newick,
                   schema="newick",
                   preserve_underscores=True)



#dataset = physcraper.opentree_helpers.get_dataset_from_treebase(study_id,
#                                phylesystem_loc='api')

#aln = dataset.char_matrices[0]

##order of data matrices is arbitratry!!!
#if len(aln) == len(tre.taxon_namespace):
#  pass
#else:
#  aln = dataset.char_matrices[1]


#aln.write(path="{}{}.aln".format(study_id, tree_id), schema="nexus")

aln = dendropy.DnaCharacterMatrix.get(file=open("{}{}.aln".format(study_id, tree_id)), schema="nexus", taxon_namespace=tre.taxon_namespace)


tre.write(path="{}{}.tre".format(study_id, tree_id), schema="nexus")

data_obj = physcraper.generate_ATT_from_phylesystem(aln=aln,
                                         workdir=workdir,
                                         config_obj=conf,
                                         study_id=study_id,
                                         tree_id=tree_id)

data_obj.write_files()
json.dump(data_obj.otu_dict, open('{}/otu_dict.json'.format(workdir), 'wb'))

sys.stdout.write("{} taxa in alignement and tree\n".format(len(data_obj.aln)))

ids = physcraper.IdDicts(conf, workdir=workdir)

scraper = physcraper.PhyscraperScrape(data_obj, ids)
scraper.threshold=25
scraper.mrca_ncbi =40996
#scraper.read_blast_wrapper()
scraper.est_full_tree()

'''scraper.run_blast_wrapper()
scraper.read_blast_wrapper()
scraper.remove_identical_seqs()
scraper.write_all_unaligned(filename="combo.fas")

json.dump(data_obj.otu_dict, open('treebase/otu_dict2.json', 'wb'))


scraper.generate_streamed_alignment()

sys.stdout.write("round 1 complete\n")
sys.stdout.write("{} taxa in alignment and tree\n".format(len(data_obj.aln)))
#
'''
'''
sys.stdout.write("round2\n")

scraper.run_blast_wrapper()
scraper.read_blast_wrapper()
scraper.remove_identical_seqs()
scraper.remove_identical_seqs()
#scraper.generate_streamed_alignment()
'''