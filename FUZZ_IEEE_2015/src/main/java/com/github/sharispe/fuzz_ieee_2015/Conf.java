/*
 *  Copyright or © or Copr. Ecole des Mines d'Alès (2012-2014) 
 *  
 *  This software is a computer program whose purpose is to provide 
 *  several functionalities for the processing of semantic data 
 *  sources such as ontologies or text corpora.
 *  
 *  This software is governed by the CeCILL  license under French law and
 *  abiding by the rules of distribution of free software.  You can  use, 
 *  modify and/ or redistribute the software under the terms of the CeCILL
 *  license as circulated by CEA, CNRS and INRIA at the following URL
 *  "http://www.cecill.info". 
 * 
 *  As a counterpart to the access to the source code and  rights to copy,
 *  modify and redistribute granted by the license, users are provided only
 *  with a limited warranty  and the software's author,  the holder of the
 *  economic rights,  and the successive licensors  have only  limited
 *  liability. 

 *  In this respect, the user's attention is drawn to the risks associated
 *  with loading,  using,  modifying and/or developing or reproducing the
 *  software by the user in light of its specific status of free software,
 *  that may mean  that it is complicated to manipulate,  and  that  also
 *  therefore means  that it is reserved for developers  and  experienced
 *  professionals having in-depth computer knowledge. Users are therefore
 *  encouraged to load and test the software's suitability as regards their
 *  requirements in conditions enabling the security of their systems and/or 
 *  data to be ensured and,  more generally, to use and operate it in the 
 *  same conditions as regards security. 
 * 
 *  The fact that you are presently reading this means that you have had
 *  knowledge of the CeCILL license and that you accept its terms.
 */
package com.github.sharispe.fuzz_ieee_2015;

/**
 *
 * @author Sébastien Harispe <sebastien.harispe@gmail.com>
 */
public class Conf {

    // Windows users must change the file separator 
    final static String root = System.getProperty("user.dir");
    final static String dataloc = root+"/data/WordNet-3.1/dict/";
    final static String mappingWordNet = dataloc+"index.noun.mapping";
    final static String mappingMC_Wordnet = root+"/data/benchmarks/semantic_similarity/Miller_Charles_index.txt";
    final static String benchmarkMC = root+"/data/benchmarks/semantic_similarity/Miller_Charles.csv";
    final static String benchmarkRG = root+"/data/benchmarks/semantic_similarity/ConceptSim/RG_sense_idWordnet.csv";
    final static String file_concept_representation = root+"/data/xp/ic_belief_theory/concept_representations.csv";
    final static String annotations = root+"/data/xp/ic_belief_theory/annot/wnid-3.1-uniq";
    final static String benchmarkSimlexNoun999 = root+"/data/benchmarks/semantic_similarity/SimLex-999/SimLex-999-nouns.txt";
    
    final static String prefixWordNet = "http://graph/wordnet/";
    

}
