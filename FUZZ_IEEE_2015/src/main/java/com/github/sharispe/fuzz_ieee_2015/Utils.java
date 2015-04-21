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

import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.openrdf.model.URI;
import org.openrdf.model.vocabulary.RDFS;
import slib.graph.model.graph.G;
import slib.graph.model.graph.elements.E;
import slib.graph.model.graph.utils.Direction;
import slib.graph.model.impl.graph.elements.Edge;
import slib.graph.model.repo.URIFactory;

/**
 *
 * @author Sébastien Harispe <sebastien.harispe@gmail.com>
 */
public class Utils {

    public static Set<URI> addRequiredLeaves(URIFactory factory, G wordnet, Map<URI, String> labels) {

        Set newLeaves = new HashSet();
        // For each node that is not a leaf
        // and that do not have more than one children, 
        // we add a topic children
        Set newEdges = new HashSet();
        for (URI u : wordnet.getV()) {

            Set<E> childrens = wordnet.getE(RDFS.SUBCLASSOF, u, Direction.IN);
            Set<E> fathers = wordnet.getE(RDFS.SUBCLASSOF, u, Direction.OUT);
//            System.out.println(u+"\t"+nbChildren);
            if (childrens.size() == 1 && fathers.size() == 1) {

                String label = labels.get(u);
//                System.out.println("Adding Topic node to " + u + "\t" + label);
//                URI uriChild = childrens.iterator().next().getSource();
//                URI uriFather = fathers.iterator().next().getTarget();
//                System.out.println("\tFather: " + labels.get(uriFather));
//                System.out.println("\t *    : " + label);
//                System.out.println("\tChild : " + labels.get(uriChild));

                URI topicURI = factory.getURI(u.toString() + "__TOPIC");
                newLeaves.add(topicURI);
                newEdges.add(new Edge(topicURI, RDFS.SUBCLASSOF, u));
                labels.put(topicURI, labels.get(u) + " [TOPIC]");
            }
        }
        System.out.println("New edges: " + newEdges.size());
        wordnet.addE(newEdges);
        return newLeaves;
        
    }

    static double[] extractColumn(int i, double[][] scoreMatrix) {
        
        double[] col = new double[scoreMatrix[i].length];
        for (int j = 0; j < scoreMatrix[i].length; j++) {
            col[j] = scoreMatrix[i][j];
        }
        return col;
    }

}
