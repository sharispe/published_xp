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

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.openrdf.model.URI;
import org.openrdf.model.vocabulary.OWL;
import org.openrdf.model.vocabulary.RDF;
import org.openrdf.model.vocabulary.RDFS;
import slib.graph.algo.utils.GAction;
import slib.graph.algo.utils.GActionType;
import slib.graph.algo.utils.GraphActionExecutor;
import slib.graph.algo.validator.dag.ValidatorDAG;
import slib.graph.io.conf.GDataConf;
import slib.graph.io.loader.wordnet.GraphLoader_Wordnet;
import slib.graph.io.util.GFormat;
import slib.graph.model.graph.G;
import slib.graph.model.graph.elements.E;
import slib.graph.model.graph.utils.Direction;
import slib.graph.model.impl.graph.elements.Edge;
import slib.graph.model.impl.graph.memory.GraphMemory;
import slib.graph.model.repo.URIFactory;
import slib.indexer.wordnet.IndexerWordNetBasic;
import slib.utils.ex.SLIB_Ex_Critic;
import slib.utils.ex.SLIB_Exception;

/**
 *
 * @author Sébastien Harispe <sebastien.harispe@gmail.com>
 */
public class DataLoaderUtils {

    public static G loadWordNet(URIFactory factory, String dataloc, boolean useOnlyNouns) throws SLIB_Exception {

        URI guri = factory.getURI(Conf.prefixWordNet);
        G wordnet = new GraphMemory(guri);

        // We load the data into the graph
        GraphLoader_Wordnet loader = new GraphLoader_Wordnet();

        GDataConf dataNoun = new GDataConf(GFormat.WORDNET_DATA, dataloc + "data.noun");
        loader.populate(dataNoun, wordnet);

        if (!useOnlyNouns) {

            GDataConf dataVerb = new GDataConf(GFormat.WORDNET_DATA, dataloc + "data.verb");
            GDataConf dataAdj = new GDataConf(GFormat.WORDNET_DATA, dataloc + "data.adj");
            GDataConf dataAdv = new GDataConf(GFormat.WORDNET_DATA, dataloc + "data.adv");

            loader.populate(dataVerb, wordnet);
            loader.populate(dataAdj, wordnet);
            loader.populate(dataAdv, wordnet);

            Set< E> nEdges = new HashSet();
            for (URI u : wordnet.getV()) {

                if (wordnet.getE(u, Direction.BOTH).isEmpty()) {
                    System.out.println("Add relationship " + u + "\t" + OWL.THING);
                    nEdges.add(new Edge(u, RDFS.SUBCLASSOF, OWL.THING));
                }
            }
            System.out.println("New Edges: " + nEdges.size());
            wordnet.addE(nEdges);

            // We root the graph which has been loaded (this is optional but may be required to compare synset which do not share common ancestors).
            GAction addRoot = new GAction(GActionType.REROOTING);
            GraphActionExecutor.applyAction(addRoot, wordnet);
        }

        // This is optional. It just shows which are the synsets which are not subsumed
        ValidatorDAG validatorDAG = new ValidatorDAG();
        Set<URI> roots = validatorDAG.getTaxonomicRoots(wordnet);
        System.out.println("Roots: " + roots);

        return wordnet;
    }

    /**
     * Load WordNet mapping that will be able to link a URI to a specific label
     *
     * @param factory
     * @param wordnet
     * @param datadir
     * @return
     * @throws SLIB_Ex_Critic
     */
    public static Map<URI, String> loadURILabels(URIFactory factory, G wordnet, String datadir, boolean useOnlyNouns) throws SLIB_Ex_Critic {

        // We create an index to map the nouns to the vertices of the graph
        // We only build an index for the nouns in this example
        String data_noun = datadir + "/index.noun";
        String data_verb = datadir + "/index.verb";
        String data_adj = datadir + "/index.adj";
        String data_adv = datadir + "/index.adv";

        Set<String> indexes = new HashSet();
        indexes.add(data_noun);

        if (!useOnlyNouns) {
            indexes.add(data_verb);
            indexes.add(data_adj);
            indexes.add(data_adv);
        }
        Map<URI, String> mapURIToLabel = new HashMap();
        for (String s : indexes) {

            IndexerWordNetBasic indexWordnetNoun = new IndexerWordNetBasic(factory, wordnet, s);

            // uncomment if you want to show the index, i.e. nouns and associated URIs (identifiers)
            for (Map.Entry<String, Set<URI>> entry : indexWordnetNoun.getIndex().entrySet()) {
                for (URI u : entry.getValue()) {
                    mapURIToLabel.put(u, entry.getKey());
                }
            }
        }
        System.out.println("Mapping loaded");

        return mapURIToLabel;
    }

    static Map<String, Set<URI>> loadLabelToURIs(G g, URIFactory factory) throws Exception {

        Map<String, Set<URI>> mapping = new HashMap();

        BufferedReader br = new BufferedReader(new FileReader(Conf.mappingWordNet));
        try {
            String line = br.readLine();

            String prefix = Conf.prefixWordNet;

            while (line != null) {

                String[] data = line.split("\t");
                String label = data[0];
                String[] ids = data[1].split(";");

                Set<URI> uris = new HashSet();

                for (String id : ids) {
                    URI u = factory.getURI(prefix + id);
                    if (!g.containsVertex(u)) {
                        System.out.println("Warning building the mapping, cannot found: " + u);
                    }
                    uris.add(u);
                }

                System.out.println(line);
                System.out.println(label + "\t" + uris.size());
                mapping.put(label, uris);

                if (uris.isEmpty()) {
                    throw new SLIB_Exception("Error in building the mapping, cannot found any URI associated to " + label);
                }

                line = br.readLine();
            }
        } finally {
            br.close();
        }
        return mapping;
    }

    static Set<BenchmarkEntry> loadBenchmarkEntries_MC(String benchmark) throws Exception {
        int posW1 = 0;
        int posW2 = 1;
        int posScore = 2;
        boolean header = true;
        return loadBenchmarkEntries(benchmark, posW1, posW2, posScore, header);
    }

    static Set<BenchmarkEntry> loadBenchmarkEntries(String benchmark, int posW1, int posW2, int posScore, boolean header) throws Exception {

        System.out.println("Loading benchmark entries");

        BufferedReader br = new BufferedReader(new FileReader(benchmark));

        Set<BenchmarkEntry> entries = new HashSet();
        try {
            String line = br.readLine();
            if (header) {
                line = br.readLine(); // skip header
            }
            while (line != null) {

                String[] data = line.split("\t");
                String label_A = data[posW1];
                String label_B = data[posW2];
                double score = Double.parseDouble(data[posScore]);

//                System.out.println(label_A+"\t"+label_B+"\t"+score);
                BenchmarkEntry e = new BenchmarkEntry(label_A, label_B, score);
                entries.add(e);

                line = br.readLine();
            }
        } finally {
            br.close();
        }

        System.out.println(entries.size() + " entries loaded");
        return entries;
    }

    static void loadAnnotations(URIFactory factory, G wordNet, String annotations) throws Exception {

        System.out.println("Loading Annotations");

        BufferedReader br = new BufferedReader(new FileReader(annotations));

        String prefix = Conf.prefixWordNet;
        int count = 0;
        int skipped = 0;
        try {
            String line = br.readLine();

            while (line != null) {

                String[] data = line.split(";");

                if (data.length != 2) {
                    continue;
                }

                String id_Instance = data[0];
                String id_Concept = data[1];

                URI i = factory.getURI("http:/i/" + id_Instance);
                URI c = factory.getURI(prefix + id_Concept);

                if (wordNet.containsVertex(c)) {
                    E e = new Edge(i, RDF.TYPE, c);
//                    System.out.println("Annot: "+e);
                    wordNet.addE(e);
                    count++;
                } else {
                    skipped++;
                }

                line = br.readLine();
            }
        } finally {
            br.close();
        }

        System.out.println(count + " annotations loaded  (skipped:" + skipped + ")");
    }

}
