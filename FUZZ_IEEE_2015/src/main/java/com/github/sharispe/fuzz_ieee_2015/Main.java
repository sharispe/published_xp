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

import java.util.Map;
import java.util.Set;
import org.openrdf.model.URI;
import slib.graph.model.graph.G;
import slib.graph.model.impl.repo.URIFactoryMemory;
import slib.graph.model.repo.URIFactory;

/**
 *
 * @author Sébastien Harispe <sebastien.harispe@gmail.com>
 */
public class Main {

    public static void main(String[] args) throws Exception {
        
        // Please refer to the classes XP_UTILS_* to reproduce the results

        boolean buildConceptRepresentation = false;
        boolean useOnlyNouns = true;
        URIFactory factory = URIFactoryMemory.getSingleton();
        G wordNet = DataLoaderUtils.loadWordNet(factory, Conf.dataloc, useOnlyNouns);
        DataLoaderUtils.loadAnnotations(factory, wordNet, Conf.annotations);
        
        Map<URI, String> labels = DataLoaderUtils.loadURILabels(factory, wordNet, Conf.dataloc, useOnlyNouns);

        FuzzIEEE_XPBrain computer = new FuzzIEEE_XPBrain(factory, wordNet, labels, Conf.file_concept_representation, buildConceptRepresentation, 1);

        Set<URI> leaves = computer.getLeaves();
        
        for (URI u : wordNet.getV()) {
            
            boolean isLeaf = leaves.contains(u);
            
            System.out.println(labels.get(u) 
                    + "\tleaf: "+isLeaf
                    + "\t" + u);
            
            System.out.println("\t[" + computer.getBelief(u) + ";" + computer.getPlausibility(u) + "]  "
                    + "eq:"+(computer.getBelief(u).equals(computer.getPlausibility(u)))
                    + "\tQ: " + computer.getQ_quantity(u)
                    + "\tPign: " + computer.getPignisticProbabilities(u));
            
            System.out.println("\tIC: [" + computer.getIC_Plausibility(u) + ";" + computer.getIC_Belief(u) + "]  "
                    + "eq:"+(computer.getIC_Belief(u).equals(computer.getIC_Plausibility(u)))
                    + "\tIC_Q: " + computer.getIC_Q_quantity(u)
                    + "\tIC_Pign: " + computer.getIC_PignisticProbabilities(u));
        }
    }
}
