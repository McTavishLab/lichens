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

configfi = "aws.config"


study_id = "pg_873"
tree_id = "tree1679"
workdir ="local_shoch"


conf = physcraper.ConfigObj(configfi)

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



dataset = physcraper.opentree_helpers.get_dataset_from_treebase(study_id,
                                phylesystem_loc='api')

aln = dataset.char_matrices[0]

##order of data matrices is arbitratry!!!
if len(aln) == len(tre.taxon_namespace):
  pass
else:
  aln = dataset.char_matrices[1]


aln.write(path="{}{}.aln".format(study_id, tree_id), schema="nexus")

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